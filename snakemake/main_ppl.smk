#!/usr/bin/env snakemake
# File: main_ppl.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 23, 2024
# Updated: Jun 13, 2024

# TODO:
# 1. Add md5sum check for the uploaded files.


import sys, csv, pathlib, yaml, benedict
from snakemake.utils import min_version

# Ensure the snakemake knows how to handle moduliazations.
min_version('6.0')


IGNORE_ME = [
  # Ignored due to reads number less than 5M. Check the 
  "SRR13076596", "SRR13076597", "SRR13076598", "SRR13076599", "SRR13076600", "SRR13076601", "SRR13076602", "SRR13076603",
  "SRR13076604", "SRR13076605", "SRR13076606", "SRR13076607", "SRR13076608", "SRR13076609", "SRR13076610", "SRR13076611",
  "SRR13076612", "SRR13076613",
]


# Configuration using YAML. Please pass the custom configuration file by `--config myconfig=PATH/TO/CONFIG/FILE`
class AsPath:
  def __init__(self, *values):
    self.values = values
  def parse(self):
    return pathlib.Path('/'.join([str(x) for x in self.values])).expanduser()


def as_path(loader, node):
  '''A custom YAML constructor for Path'''
  if isinstance(node, yaml.ScalarNode):
    values = loader.construct_scalar(node)
    return AsPath(values).parse()
  elif isinstance(node, yaml.SequenceNode):
    values = loader.construct_sequence(node)
    return AsPath(*values).parse()
  raise NotImplementedError('Unsupported node type: {}'.format(type(node)))


# Custom YAML constructor
yaml.add_constructor('!AsPath', as_path)


# Utility function to collect FASTQ files
def prepare_samples(sample_table, out_dir, create_link=False):
  '''Collect FASTQ files from a sample-to-path file.

  The sample_table should be plain text file contain 4 or 3 columns.
  The 1st and 2ed columns are the run ID and sample ID, respectively.
  The 3rd column indicates which sequencing technology was used. Options are 'pro-seq' and 'gro-seq', case-insensitive.
  The 4rd and 5th columns are the path to the R1 and R2 FASTQ file, respectively, if the sample is paired-end.
  Otherwise, the 3rd column is the path to the FASTQ file for single-end samples.
  '''
  out_dir = pathlib.Path(out_dir)
  if not out_dir.exists(): out_dir.mkdir(parents=True)

  _sample_info_dict = {}
  with open(sample_table) as fh:
    for per_run in fh:
      group_id, sample_id, seq_tech, strand_spec, src_read_r1, *flex_col = per_run.strip().split(',')

      if group_id in IGNORE_ME or sample_id in IGNORE_ME: continue
      if sample_id in _sample_info_dict: raise KeyError(f'Duplicate sample ID: {sample_id}. Ensure they are unique (2ed column).')

      if len(flex_col) == 0:
        is_pe = False
        src_read_r2 = None
      elif len(flex_col) == 1:
        is_pe = flex_col[0] != ''
        src_read_r2 = flex_col[0]
      else:
        raise ValueError('Invalid # of columns!!! Expected 4 or 3, but found {}.'.format(3 + len(flex_col)))

      if create_link:
        if is_pe:
          dst_read_r1 = out_dir / '.'.join([sample_id, group_id, 'R1.fastq.gz'])
          if not pathlib.Path(dst_read_r1).is_symlink(): os.symlink(src_read_r1, dst_read_r1)

          dst_read_r2 = out_dir / '.'.join([sample_id, group_id, 'R2.fastq.gz']) 
          if not pathlib.Path(dst_read_r2).is_symlink(): os.symlink(src_read_r2, dst_read_r2)
        else:
          dst_read_r1 = out_dir / '.'.join([sample_id, group_id, 'fastq.gz'])
          if not pathlib.Path(dst_read_r1).is_symlink(): os.symlink(src_read_r1, dst_read_r1)

      which_prime = 'three_prime' if seq_tech.lower() in ['pro-seq', 'proseq'] else 'five_prime'
      _sample_info_dict[sample_id] = {
        'group_id': group_id, 'paired_end': is_pe, 'which_prime': which_prime, 'strand_spec': strand_spec
      }

  return _sample_info_dict


# Configurations
bconfig = None
with open(config['myconfigfile'], 'r') as f:
  config = yaml.load(f, Loader=yaml.FullLoader)
  bconfig = benedict.BeneDict(config, keypath_separator = '/')

chunk_token = bconfig.get('chunk_token')

output_dir = bconfig.get('output_dir')
alignment_dir = output_dir / 'read_alignment'
variant_dir = output_dir / 'variant_calling'
quant_dir = output_dir / 'quantify_expression'
tre_dir = output_dir / 'tre_identification'
upload_dir = output_dir / 'uploading'

sample_table = bconfig.get('sample_table')
sample_info_dict = prepare_samples(sample_table, alignment_dir / 'fastqs')

remote_dir = bconfig.get("remote_dir")
remote_user = bconfig.get("remote_user")
remote_ipaddr = bconfig.get("remote_ipaddr")

# Singularity image
singularity_image = bconfig.get('singularity_image')
assert singularity_image.exists()

python_scripts = bconfig.get('python_scripts')
assert python_scripts.exists()

# Wildcard constraints
wildcard_constraints:
  sample_id = r'\w+', group_id = r'\w+'

localrules:
  all,
  main_s01_download_single_end, main_s01_download_paired_end,
  it_s02_bigwig_is_ready, main_s03_upload_tres, main_s04_upload_tres_all,
  cv_s02_haplotype_caller_ready, main_s03_upload_variants, main_s04_upload_variants_all,
  qt_s02_count_reads_ready, main_s03_upload_read_counts, main_s04_upload_readcounts_all,
  main_s03_upload_alignments, main_s04_upload_alignments_all,


#
## Rules
#
# Pipeline to control read alignment
module align_reads:
  snakefile: 'align_reads_ppl.smk'
  config: (bconfig, sample_info_dict)


module quantify_expression:
  snakefile: 'quantify_expression_ppl.smk'
  config: (bconfig, sample_info_dict)


# Pipeline to control variant calling
module call_variants:
  snakefile: 'call_variants_ppl.smk'
  config: (bconfig, sample_info_dict)


# Pipeline to control transcriptional regulatory element identification.
module identify_tre:
  snakefile: 'identify_tre_ppl.smk'
  config: (bconfig, sample_info_dict)


use rule * from align_reads as ar_*
use rule * from call_variants as cv_*
use rule * from identify_tre as it_*
use rule * from quantify_expression as qt_*


# Final outputs
all_groups = list(sample_info_dict.keys())
align_reads_done = upload_dir / ('align_reads.' + chunk_token + '.done')
call_variants_done = upload_dir / ('call_variant.' + chunk_token + '.done')
identify_tres_done = upload_dir / ('identify_tre.' + chunk_token + '.done')
read_counts_done = upload_dir / ('quantify_expression.' + chunk_token + '.done')
final_output = [align_reads_done, call_variants_done, read_counts_done, identify_tres_done]

rule all:
  input: final_output
  default_target: True


rule main_s01_download_paired_end: # Download paired-end samples
  input:
    sample_table
  output:
    read_r1 = alignment_dir / 'fastqs/{group_id}.{sample_id}.R1.fastq.gz',
    read_r2 = alignment_dir / 'fastqs/{group_id}.{sample_id}.R2.fastq.gz',
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    IFS=','
    key_file=${{HOME}}/.ssh/baolab_id_rsa 
    while read -r gid sid _ _ r1 r2; do
      rsync -avzhL -e "ssh -i ${{key_file}}" {remote_user}@{remote_ipaddr}:${{r1}} {output.read_r1}
      rsync -avzhL -e "ssh -i ${{key_file}}" {remote_user}@{remote_ipaddr}:${{r2}} {output.read_r2}
    done < <(grep -w {wildcards.sample_id} {input})
    '''


rule main_s01_download_single_end: # Download single-end samples
  input:
    sample_table
  output:
    read_r1 = alignment_dir / 'fastqs/{group_id}.{sample_id}.fastq.gz'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    IFS=','
    key_file=${{HOME}}/.ssh/baolab_id_rsa 
    while read -r gid sid _ _ r1 r2; do
      rsync -avzhL -e "ssh -i ${{key_file}}" {remote_user}@{remote_ipaddr}:${{r1}} {output.read_r1}
    done < <(grep -w {wildcards.sample_id} {input})
    '''


rule main_s02_compress_alignments:
  input:
    rules.it_s02_create_bigwig_done.output.bw_ready,
    rules.cv_s02_haplotype_caller_ready.output,
    rules.qt_s02_count_reads_ready.output,
    fasta_file = bconfig.reference_genome,
    bam_file = alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.bqsr.bam',
    bai_file = alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.bqsr.bam.bai'
  output:
    cram_file = alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.bqsr.cram',
    cram_md5_file = alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.bqsr.cram.md5'
  resources:
    time = '1:00:00', mem = '4G', cpus_per_task = 8
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} samtools view \
      -C -@ {resources.cpus_per_task} -T {input.fasta_file} -o {output.cram_file} {input.bam_file}
    md5sum {output.cram_file} > {output.cram_md5_file}
    rm -f {input.bam_file} && touch {input.bam_file}
    rm -f {input.bai_file} && touch {input.bai_file}
    '''


rule main_s03_upload_alignments:
  input:
    alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.bqsr.cram',
    alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.bqsr.cram.md5',
    alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.recal.table',
    alignment_dir / 'alignment/{group_id}/{group_id}.mkdup.sncr.recal.table.md5',
  output:
    upload_dir / 'alignment/{group_id}.upload_alignment.done'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent,
    remote_dir = lambda wc: remote_dir / 'read_alignment/alignment' / wc.group_id,
  shell:
    '''
    mkdir -p {params.output_dir}

    key_file=${{HOME}}/.ssh/baolab_id_rsa 
    for per_file in {input}; do
      rsync -avzhL -e "ssh -i ${{key_file}}" ${{per_file}} {remote_user}@{remote_ipaddr}:{params.remote_dir}/
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output}
    '''


rule main_s04_upload_alignments_all:
  input: expand(upload_dir / 'alignment/{group_id}.upload_alignment.done', group_id=all_groups),
  output: align_reads_done
  shell: '''touch {output}'''


rule main_s03_upload_tres:
  input:
    bw_pl_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_pl.bw',
    bw_mn_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_mn.bw',
    dv_tre_file = tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_divergent_peaks.bed',
    bd_tre_file = tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_bidirectional_peaks.bed',
    ud_tre_file = tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_unidirectional_peaks.bed',
  output:
    upload_dir / 'tres/{group_id}.upload_tre.done'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent,
    remote_dir = lambda wc: remote_dir / 'tre_identification/tre_results' / wc.group_id,
  shell:
    '''
    mkdir -p {params.output_dir}

    key_file=${{HOME}}/.ssh/baolab_id_rsa 
    for per_file in {input}; do
      rsync -avzhL -e "ssh -i ${{key_file}}" ${{per_file}} {remote_user}@{remote_ipaddr}:{params.remote_dir}/
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output}
    '''


rule main_s04_upload_tres_all:
  input: expand(upload_dir / 'tres/{group_id}.upload_tre.done', group_id=all_groups),
  output: identify_tres_done
  shell: '''touch {output}'''


rule main_s03_upload_variants:
  input:
    variant_dir / 'variant' / '{group_id}/{group_id}.g.vcf.gz',
    variant_dir / 'variant' / '{group_id}/{group_id}.g.vcf.gz.md5',
    variant_dir / 'variant' / '{group_id}/{group_id}.g.vcf.gz.tbi',
    variant_dir / 'variant' / '{group_id}/{group_id}.g.vcf.gz.tbi.md5',
  output:
    upload_dir / 'variants/{group_id}.upload_variants.done',
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent,
    remote_dir = lambda wc: remote_dir / 'variant_calling/variant' / wc.group_id,
  shell:
    '''
    mkdir -p {params.output_dir}

    key_file=${{HOME}}/.ssh/baolab_id_rsa 
    for per_file in {input}; do
      rsync -avzhL -e "ssh -i ${{key_file}}" ${{per_file}} {remote_user}@{remote_ipaddr}:{params.remote_dir}/
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output}
    '''


rule main_s04_upload_variants_all:
  input: expand(upload_dir / 'variants/{group_id}.upload_variants.done', group_id=all_groups),
  output: call_variants_done
  shell: '''touch {output}'''


rule main_s03_upload_read_counts:
  input:
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.md5',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.exonic_regions.tsv.gz',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.intronic_regions.tsv.gz',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.intergenic_enhancers.tsv.gz',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.transcript_regions.tsv.gz',
  output:
    upload_dir / 'quantify_expression/{group_id}.upload_count_tables.done',
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent,
    remote_dir = lambda wc: remote_dir / 'quantify_expression/featureCounts' / wc.group_id,
  shell:
    '''
    mkdir -p {params.output_dir}

    key_file=${{HOME}}/.ssh/baolab_id_rsa 
    for per_file in {input}; do
      rsync -avzhL -e "ssh -i ${{key_file}}" ${{per_file}} {remote_user}@{remote_ipaddr}:{params.remote_dir}/
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output}
    '''


rule main_s04_upload_readcounts_all:
  input: expand(upload_dir / 'quantify_expression/{group_id}.upload_count_tables.done', group_id=all_groups),
  output: read_counts_done
  shell: '''touch {output}'''


# rule main_s03_upload_tres:
#   input:
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.infp.bed.gz',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.infp.bed.gz.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.infp.bed.gz.tbi',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.infp.bed.gz.tbi.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.raw.bed.gz',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.raw.bed.gz.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.raw.bed.gz.tbi',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.raw.bed.gz.tbi.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.full.bed.gz',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.full.bed.gz.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.full.bed.gz.tbi',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.full.bed.gz.tbi.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.prob.bed.gz',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.prob.bed.gz.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.prob.bed.gz.tbi',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.prob.bed.gz.tbi.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.score.bed.gz',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.score.bed.gz.md5',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.score.bed.gz.tbi',
#     tre_dir / 'tre_results/{group_id}/{group_id}.dREG.peak.score.bed.gz.tbi.md5',
#   output:
#     upload_dir / 'tres/{group_id}.upload_tre.done'
#   params:
#     output_dir = lambda wildcards, output: Path(output[0]).parent,
#     remote_dir = lambda wc: remote_dir / 'tre_identification/tre_results' / wc.group_id,
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     key_file=${{HOME}}/.ssh/baolab_id_rsa 
#     for per_file in {input}; do
#       rsync -avzhL -e "ssh -i ${{key_file}}" ${{per_file}} {remote_user}@{remote_ipaddr}:{params.remote_dir}/
#       rm -f ${{per_file}} && touch ${{per_file}}
#     done
#     touch {output}
#     '''
