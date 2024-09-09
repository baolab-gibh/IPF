#!/usr/bin/env python3
# File: create_bigwig_file.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 21, 2024
# Updated: May 21, 2024

import sys, logging, time
from argparse import ArgumentParser
from multiprocessing import Pool

import pysam
import pyBigWig


CHROM_ORDER = dict(zip([f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y', 'M']], list(range(1, 26))))


def get_cli_args():
    """Parse command line arguments"""
    parser = ArgumentParser(description="Create bigwig file from bam file")

    parser.add_argument('-b', '--bamfiles', required=True, nargs='+', help='Input bam file')
    parser.add_argument('-f', '--incl-flags', default=4095, type=int, help='Read flags to be included. Default: %(default)s')
    parser.add_argument('-F', '--excl-flags', default=0, type=int, help='Read flags to be excluded. Default: %(default)s')
    # parser.add_argument('-B', '--bin-size', default=50, help="Bin size (in bases) for bigwig file. Default: %(default)s")
    # parser.add_argument('-n', '--normalize', default=None, action='store_true', help="Normalize bigwig/bedgraph file. Default: %(default)s")
    # parser.add_argument('-d', '--dynamic-chunks', default=False, action='store_true', help="Determine chunks dynamically. Default: %(default)s")
    parser.add_argument('-t', '--strand-spec', default='*', choices=['*', '+', '-'], help="Strand specificity. Default: %(default)s")
    parser.add_argument('-p', '--which-end', default='five_prime', choices=['five_prime', 'three_prime'], help="Which end to count the read. Default: %(default)s")
    parser.add_argument('-c', '--chunks', default=None, nargs='+', help="Genomic regions, available format: chrom:start-end, chrom:start, chrom. Default: %(default)s")
    parser.add_argument('-C', '--chunks-file', default=None, help="Multiple genomic regions in a file. One region per line, available region format: chrom:start-end, chrom:start, chrom. Default: %(default)s")
    parser.add_argument('-n', '--n-cpus', default=1, type=int, help="Number of CPUs. Default: %(default)s")
    parser.add_argument('-D', '--deflation-cpus', default=1, type=int, help="Number of CPUs to decompress the bam file. Default: %(default)s")
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help="Verbose mode. Default: %(default)s")
    parser.add_argument('-O', '--out-format', default='ucsc_bw', choices=['ucsc_bw', 'dreg_bw', 'ucsc_bg'], help="Output format, default: %(default)s")
    parser.add_argument('-o', '--out-prefix', default='output', help='Output prefix. Default: %(default)s')

    return parser.parse_args()


def generate_dym_chunks(coverage_list):
    """Generate dynamic chunks."""
    pass


def check_chunk(chunk_str):
    """Check if the chunk is valid.

    Check a given region string is valid or not.
    Available region format: chrom:start-end, chrom:start, chrom.
    The chrom could be a mixture of digits and number.
    The start and end (if available) should be digits, meanwhile the start should be greated than 0 and less than end.
    """
    if chunk_str is None: return False

    chunk_str = chunk_str.strip()
    is_valid = chunk_str.count(':') <= 1 and chunk_str.count('-') <= 1
    if ':' in chunk_str and '-' in chunk_str:
        chrom, start_end = chunk_str.split(':')
        start, end = start_end.split('-')
        is_valid = is_valid and chrom.isalnum() and start.isdigit() and end.isdigit() and int(start) <= int(end) and int(start) > 0
    elif ':' in chunk_str:
        chrom, start = chunk_str.split(':')
        is_valid = is_valid and start.isdigit() and chrom.isalnum() and int(start) > 0
    else:
        is_valid = is_valid and chunk_str.isalnum()

    return is_valid


def parse_chunks(chunks=None, chunks_file=None):
    """Obtain the chunks files.

    Collecting chunks from -c/--chunks and -C/--chunks-file options.
    Both of them will be used if they are specified at the same time.
    """
    chunk_list = []
    if chunks is not None:
        for per_chunk in chunks:
            if check_chunk(per_chunk): # Obtain chunks from -c/--chunks option
                chunk_list.append(per_chunk) 

    if chunks_file is not None: # Obtain chunks from -C/--chunks-file option
        with open(chunks_file, 'r') as f:
            chunk_list.extend([line.strip() for line in f.readlines() if check_chunk(line)])

    if len(chunk_list) == 0: # Generate chunks dynamically
        raise ValueError("Please specify the chunks with -c/--chunks or -C/--chunks-file.")

    return chunk_list


def scatter(args):
    """Get genome coverage"""
    bam_file, per_chunk, options, *_ = args

    coverage = {}
    bam_records = pysam.AlignmentFile(bam_file, 'rb', threads=options.deflation_cpus)
    _, chrom, start, end = bam_records.parse_region(region=per_chunk)
    chrom = bam_records.get_reference_name(chrom)
    for per_read in bam_records.fetch(chrom, start + 1, end, multiple_iterators=True):
        tb_incl = (per_read.flag & (options.excl_flags | 1796) == 0) and (per_read.flag & options.incl_flags != 0)
        if not tb_incl: continue

        blocks = per_read.get_blocks()
        if blocks:
            rd_start, rd_end = blocks[0]
        else:
            continue
        strand = '+' if per_read.is_forward else '-'
        if options.out_format.lower() == "dreg_bw":
            if options.which_end == 'five_prime' and per_read.is_forward:
                key = (chrom, rd_start, rd_start + 1, strand)
            elif options.which_end == 'five_prime' or per_read.is_forward:
                key = (chrom, rd_end - 2, rd_end - 1, strand)
            else:
                key = (chrom, rd_start, rd_start + 1, strand)

        elif options.out_format.lower() in ["ucsc_bw", "ucsc_bg"]:
            if options.which_end == 'five_prime':
                key = (chrom, rd_start, rd_start + 1, '*')
            else:
                key = (chrom, rd_end - 2, rd_end - 1, "*")
        else:
            raise ValueError("Unsupported output format: {}".format(options.out_format))
            
        val = 1 if per_read.is_forward else -1
        if key in coverage:
            coverage[key] += val
        else:
            coverage[key] = val
    bam_records.close()

    return coverage


def collect_basic_info(bamfiles):
    """Generate the chromosome sizes"""
    chrom_sizes = {}
    for per_bam in bamfiles:
        bam = pysam.AlignmentFile(per_bam, 'rb')
        for chrom, length in zip(bam.references, bam.lengths):
            if chrom in chrom_sizes:
                chrom_sizes[chrom] = max(chrom_sizes[chrom], length)
            else:
                chrom_sizes[chrom] = length
        bam.close()

    return sorted([(k, v) for k, v in chrom_sizes.items()], key=lambda x: CHROM_ORDER[x[0]])


def collect(coverage_list):
    """Collect the coverage list."""
    cov_list = []
    for per_cov_dict in coverage_list:
        for (chrom, start, end, strand), count in per_cov_dict.items():
            cov_list.append([chrom, start, end, strand, count])
    return sorted(cov_list, key=lambda x: (CHROM_ORDER[x[0]], x[1]))


def write_bigwig(cov_list, header, out_prefix, outfmt, strand_spec, bin_size=200, normalize=False):
    """Write bigwig file."""
    strand_corr = 1 if strand_spec in ['+', '*'] else -1
    if outfmt.lower() == "dreg_bw":
        plus_cov, minus_cov = {'chrom':[], 'start':[], 'end':[], 'value':[]}, {'chrom':[], 'start':[], 'end':[], 'value':[]}
        with pyBigWig.open(f"{out_prefix}.plus_strand.bw", 'wb') as ps_bw, pyBigWig.open(f"{out_prefix}.minus_strand.bw", 'wb') as ms_bw:
            ps_bw.addHeader(header)
            ms_bw.addHeader(header)
            for chrom, start, end, strand, value in cov_list:
                if strand == '+': # Write the plus strand
                    plus_cov['chrom'].append(chrom)
                    plus_cov['start'].append(start)
                    plus_cov['end'].append(end)
                    plus_cov['value'].append(float(value) * strand_corr)
                elif strand == '-': # Write the minus strand
                    minus_cov['chrom'].append(chrom)
                    minus_cov['start'].append(start)
                    minus_cov['end'].append(end)
                    minus_cov['value'].append(float(value) * strand_corr)
                else:
                    continue

            ps_bw.addEntries(plus_cov['chrom'], plus_cov['start'], ends=plus_cov['end'], values=plus_cov['value'])
            ms_bw.addEntries(minus_cov['chrom'], minus_cov['start'], ends=minus_cov['end'], values=minus_cov['value'])
    elif outfmt.lower() in ["ucsc_bw", "ucsc_bg"]:
        print("Not implemented yet")


def main():
    """Main function"""
    options = get_cli_args()
    bamfiles = options.bamfiles
    chunk_list = parse_chunks(options.chunks, options.chunks_file)

    log_level = logging.INFO if options.verbose else logging.ERROR
    logging.basicConfig(level=log_level, format="[%(asctime)s] %(levelname)s: %(message)s")
    logman = logging.getLogger()

    inputs = [[bam_file, per_chunk, options] for bam_file in bamfiles for per_chunk in chunk_list]
    logman.info(f"Obtaining coverages from {len(bamfiles)} bam files for {len(chunk_list)} chunks.")
    pool = Pool(options.n_cpus)
    results = pool.map_async(scatter, inputs)
    pool.close()
    pool.join()

    try:
        coverage_results = results.get()
        coverages = collect(coverage_results)
    except Exception as e:
        coverages = None
        print(f"Failed to get the coverages. Error: {e}")
        sys.exit(1)
    logman.info(f"Obtained {len(coverages)} coverages.")

    time.sleep(1)
    logman.info(f"Writing bigwig files into disk in {options.out_format} format. Check {options.out_prefix}.* files.")
    if coverages:
        header = collect_basic_info(bamfiles)
        write_bigwig(coverages, header, options.out_prefix, options.out_format, options.strand_spec)
    logman.info("Done!")


if __name__ == "__main__":
    main()
