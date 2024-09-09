#!/usr/bin/env snakemake
# File: identify_tre_ppl.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 23, 2024
# Updated: Jun 13, 2024


import sys, csv, pathlib, yaml, benedict
from snakemake.utils import min_version

# Ensure the snakemake knows how to handle moduliazations.
min_version('6.0')


def determine_time_by_chrom(wc):
  return '4:59:59'

def which_gpu(wc):
  chrom = wc.per_chunk.replace('chr', '').lower()
  if chrom in ['x', 'y']:
    return '1'
  elif int(chrom) % 3 == 1:
    return '1'
  elif int(chrom) % 3 == 2:
    return '2'
  return '3'


#
## Resources
#
# Project directory
bconfig, sample_info_dict = config
project_dir = bconfig.get('project_dir')

reference_dir = project_dir / 'outputs/references'
reference_genome = reference_dir / 'genome/GRCh38.primary_assembly.genome.only_chroms.fa'
reference_genome_dict = reference_dir / 'genome/GRCh38.primary_assembly.genome.only_chroms.dict'
reference_genome_index = reference_dir / 'genome/GRCh38.primary_assembly.genome.only_chroms.fa.fai'
reference_genome_chrom_length = bconfig.get('reference_genome_chrom_length')

dreg_pretrained_model = bconfig.get('dreg_pretrained_model')

# I/O controlling
input_dir = bconfig.get('output_dir') / 'read_alignment'
alignment_dir = input_dir / 'alignment'

output_dir = bconfig.get('output_dir')
tre_dir = output_dir / 'tre_identification'


# Singularity image
singularity_image = bconfig.get('singularity_image')
assert singularity_image.exists()

python_scripts = bconfig.get('python_scripts')
assert python_scripts.exists()


# Outputs
# all_chunks = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY'] # + ['chrM']
all_groups = list(sample_info_dict.keys())
all_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']
all_chunks = ["chr" + "_".join(all_chroms[s::6]) for s in range(6)]
final_output = expand(tre_dir / 'tre_results/{group_id}/{group_id}.dREG.{per_ext}.bed.gz', group_id=all_groups, per_ext=['infp', 'peak.full', 'peak.score'])


# Wildcard constraints
wildcard_constraints:
  per_run = r'\w+', group_id = r'\w+', per_chunk = r'\w+'


#
## Rules
#
rule all:
  input: final_output


rule s01_create_bigwig: # Create BigWig files, which will be used for visualization
  input:
    bam_file = alignment_dir / 'alignment/{sample_id}/{sample_id}.mkdup.sncr.bqsr.bam',
  output:
    bw_pl_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_pl.bw',
    bw_mn_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_mn.bw',
  resources:
    time = '10:00:00', mem = '16G'
  params:
    output_prefix = lambda wc, output: Path(output[0]).parent / (wc.sample_id + '.pints'),
    exp_type = lambda wc: 'PROseq' if seq_tech_dict[wc.sample_id].lower() in ['pro-seq', 'proseq'] else 'GROseq',
    input_dir = lambda wildcards, output, input: Path(input[0]).parent
  shell:
    '''
    singularity exec {singularity_image} pints_visualizer --exp-type {params.exp_type} --bam {input.bam_file} --output-prefix {params.output_prefix}
    rm -f {input.bam_file} {params.input_dir}/*.npy && touch {input.bam_file}
    touch {output.bw_ready}
    '''


rule s02_create_bigwig_done:
  input:
    bw_pl_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_pl.bw',
    bw_mn_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_mn.bw',
  output:
    tre_dir / 'read_coord/{sample_id}/{sample_id}.create_bigwig.done',
  shell:
    ''' touch {output} '''


rule s02_call_tres: # Call TREs using PINTs from BigWig files (check rule s03_create_bigwig first)
  input:
    bw_pl_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_pl.bw',
    bw_mn_file = tre_dir / 'read_coord/{sample_id}/{sample_id}.pints_mn.bw',
  output:
    dv_tre_file = tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_divergent_peaks.bed',
    bd_tre_file = tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_bidirectional_peaks.bed',
    ud_tre_file = tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_unidirectional_peaks.bed',
  priority: 200
  resources:
    cpus_per_task = 10, mem = '32G', time = '16:00:00'
  params:
    out_dir = lambda wc, output: Path(output[0]).parent,
    out_pre = lambda wc, output: wc.sample_id + '.pints',
    exp_type = lambda wc: 'PROseq' if seq_tech_dict[wc.sample_id].lower() in ['pro-seq', 'proseq'] else 'GROseq'
  shell:
    '''
    mkdir -p {params.out_dir}
    singularity exec {singularity_image} pints_caller \
        --bw-pl {input.bw_pl_file} \
        --bw-mn {input.bw_mn_file} \
        --exp-type {params.exp_type} \
        --thread {resources.cpus_per_task} \
        --save-to {params.out_dir} \
        --file-prefix {params.out_pre} \
        --dont-check-updates \
        --disable-small
     rm -fr {params.out_dir}/*.csi
    '''


rule s03_call_tre_done:
  input:
    expand(tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_divergent_peaks.bed', sample_id=all_samples),
    expand(tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_bidirectional_peaks.bed', sample_id=all_samples),
    expand(tre_dir / 'tre_results/{sample_id}/{sample_id}.pints_1_unidirectional_peaks.bed', sample_id=all_samples),
  output:
    tre_dir / 'tre_results/{group_id}/{group_id}.call_tres.done',
  shell:
    ''' touch {output} '''
