#!/usr/bin/env pyton3
# File: create_chromhmm_binfile.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 03, 2024
# Updated:

import sys
from pathlib import Path
from argparse import ArgumentParser

import numpy as np
import polars as pls


def chunknize(start: int, end: int, step: int = 500000):
    starts = list(range(start, end + 1, step))
    ends = [x - 1 for x in starts[1:]] + [end]
    return zip(starts, ends)


def create_indice_tab(step: int = 500000, which_chrom: str = "all") -> dict:
    CHROM_LENGTH = {
         "chr1": 248956422,  "chr2": 242193529,  "chr3": 198295559,  "chr4": 190214555,  "chr5": 181538259,  "chr6": 170805979,
         "chr7": 159345973,  "chr8": 145138636,  "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16":  90338345, "chr17":  83257441, "chr18":  80373285,
        "chr19":  58617616, "chr20":  64444167, "chr21":  46709983, "chr22":  50818468,  "chrX": 156040895,  "chrY":  57227415,
    }

    indice_dict = {}
    for k, v in CHROM_LENGTH.items():
        if k in which_chrom or which_chrom == "all":
            per_chunk = {(k, s): 0 for s, _ in chunknize(1, v, step)}
            indice_dict.update(per_chunk)

    return indice_dict


def get_chunk_stats(start, end):
    return [x[0] for x in chunknize(start, end, 200)]


def subpartition_intervals(bed_tab, bin_size=200, min_overlap=10, drop_dups=True, sort=True):
    itvl_tab = bed_tab.with_columns(
        pls.when(pls.col("starts") % bin_size >= min_overlap)
           .then(pls.col("starts").floordiv(bin_size))
           .otherwise(pls.col("starts").floordiv(bin_size).add(1))
           .mul(bin_size)
           .add(1)
           .alias("starts_binned"),
        pls.when(pls.col("ends") % bin_size >= min_overlap)
           .then(pls.col("ends").floordiv(bin_size).add(1))
           .otherwise(pls.col("ends").floordiv(bin_size))
           .mul(bin_size)
           .alias("ends_binned"),
    ).select("seqnames", "starts_binned", "ends_binned")

    signal_bins = []
    for seqname, starts_binned, ends_binned in itvl_tab.iter_rows():
        if ends_binned - starts_binned > 200:
            chunks = get_chunk_stats(starts_binned, ends_binned)
            signal_bins.extend([(seqname, per_bin_start) for per_bin_start in chunks])
        else:
            signal_bins.append((seqname, starts_binned))

    if drop_dups:
        signal_bins = list(set(signal_bins))

    if sort:
        signal_bins = sorted(signal_bins)

    return signal_bins


def get_cliopt():
    parser = ArgumentParser()
    parser.add_argument("in_file", help="The input .tsv file. Required")
    parser.add_argument("-t", "--cell-type", required=True, help="The cell type. Required")
    parser.add_argument("-i", "--sample-id", required=True, help="Sample ID. Required")
    parser.add_argument("-b", "--bin-size", default=200, type=int, help="The bin size used to partition the genome. Default: %(default)s")
    parser.add_argument("-O", "--min-overlap", default=10, type=int, help="The minimum overlap between genome-wide bin and functional bin. Default: %(default)s")
    parser.add_argument("-c", "--which-seqnames", default="all", nargs="+", help="The chromosome names. Default: %(default)s")
    parser.add_argument("-o", "--out_dir", default="ChromHMM_bins", help="The output file. Default: %(default)s")
    return parser.parse_args()


def create_chromhmm_binfile(in_file, bin_size=200, min_overlap=10):
    bed_tab = pls.read_csv(in_file, has_header=False, new_columns=["seqnames", "starts", "ends", "score", "p_value"], separator="\t")

    ts_bins = subpartition_intervals(bed_tab, bin_size=bin_size, min_overlap=min_overlap, drop_dups=True, sort=True)
    gw_bins = create_indice_tab(bin_size)

    for per_sig_bin in ts_bins: gw_bins[per_sig_bin] = 1
    signal_list = [[seqname, start, signal] for (seqname, start), signal in gw_bins.items()]
    gw_sig_tab = pls.DataFrame(signal_list).rename(
        {"column_0": "seqnames", "column_1": "starts", "column_2": "signal"}
    ).with_columns(
        pls.col("starts").add(bin_size - 1).alias("ends")
    ).select("seqnames", "starts", "ends", "signal")

    return gw_sig_tab


def write_results(gw_sig_tab: pls.DataFrame, out_dir, sample_id=None, cell_type=None):
    all_seqnames = gw_sig_tab.select("seqnames").unique().to_numpy().squeeze().tolist()
    for per_sn in all_seqnames:
        out_file = out_dir / f"{cell_type}/{per_sn}/{sample_id}.{per_sn}.chromhmm_binary.txt"
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True, exist_ok=True)

        with open(out_file, "w") as out:
            out.write(f"genome\t{per_sn}\n")
            out.write(f"{cell_type}_Enhancer\n")
            signal = "\n".join(gw_sig_tab.filter(pls.col("seqnames") == per_sn).with_columns(pls.col("signal").cast(str)).select("signal").to_numpy().squeeze().tolist())
            out.write(signal+"\n")


def main():
    args = get_cliopt()
    in_file = args.in_file
    cell_type = args.cell_type
    bin_size = args.bin_size
    min_overlap = args.min_overlap
    out_dir = Path(args.out_dir)
    sample_id = args.sample_id

    gw_sig_tab = create_chromhmm_binfile(in_file, bin_size=bin_size, min_overlap=min_overlap)
    write_results(gw_sig_tab, out_dir, sample_id=sample_id, cell_type=cell_type)


if __name__ == "__main__":
    main()
