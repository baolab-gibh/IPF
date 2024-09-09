#!/usr/bin/env snakemake
# File: quantify_expression.ppl.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 26, 2024
# Updated: Jun 26, 2024


import sys, csv, pathlib, yaml, benedict
from snakemake.utils import min_version

# Ensure the snakemake knows how to handle moduliazations.
min_version('6.0')


#
## Resources
#
# Project directory
bconfig, sample_info_dict = config
project_dir = bconfig.get('project_dir')

reference_dir = project_dir / 'outputs/references'
exonic_regions = reference_dir / 'genomic_feature/gencode.v46.basic.annotation.exon_and_utr.gtf.gz'
intronic_regions = reference_dir / 'genomic_feature/gencode.v46.basic.annotation.intron_regions.gtf.gz'
enhancer_regions = reference_dir / 'genomic_feature/homo_sapiens.GRCh38.Regulatory_Build.intergenic_enhancers.20231016.gff.gz'
transcript_regions = reference_dir / 'genomic_feature/gencode.v46.basic.annotation.transcript_regions.gtf.gz'


# I/O controlling
input_dir = bconfig.get('output_dir') / 'read_alignment'
alignment_dir = input_dir / 'alignment'

output_dir = bconfig.get('output_dir')
quant_dir = output_dir / 'quantify_expression'


# Singularity image
singularity_image = bconfig.get('singularity_image')
assert singularity_image.exists()


# Outputs
all_groups = list(sample_info_dict.keys())
final_output = expand(quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.ready', group_id=all_groups)


# Wildcard constraints
wildcard_constraints:
  group_id = r'\w+', per_chunk = r'\w+'


#
## Rules
#
rule all:
  input: final_output


# Count reads using featureCounts
rule s01_count_reads:
  input:
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.bam.bai',
    exons = exonic_regions,
    introns = intronic_regions,
    enhancers = enhancer_regions,
    transcripts = transcript_regions,
  output:
    count_table_md5 = quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.md5',
    exon_counttab_file = quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.exonic_regions.tsv.gz',
    intron_counttab_file = quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.intronic_regions.tsv.gz',
    enhancers_counttab_file = quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.intergenic_enhancers.tsv.gz',
    transcript_counttab_file = quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.transcript_regions.tsv.gz',
  resources:
    mem = '8G', time = '2:00:00', cpus_per_task = 2
  params:
    output_dir = lambda _, output: Path(output[0]).parent,
    exon_counttab_file = lambda wc: quant_dir / 'featureCounts' / wc.group_id / (wc.group_id + '.read_count.exonic_regions.tsv'),
    intron_counttab_file = lambda wc: quant_dir / 'featureCounts' / wc.group_id / (wc.group_id + '.read_count.intronic_regions.tsv'),
    enhancers_counttab_file = lambda wc: quant_dir / 'featureCounts' / wc.group_id / (wc.group_id + '.read_count.intergenic_enhancers.tsv'),
    transcript_counttab_file = lambda wc: quant_dir / 'featureCounts' / wc.group_id / (wc.group_id + '.read_count.transcript_regions.tsv'),
  shell:
    '''
    mkdir -p {params.output_dir}

    # Count reads mapped to exons
    singularity exec {singularity_image} featureCounts \
        -p --ignoreDup -t exon -g exon_id --extraAttributes gene_name,transcript_id,gene_type -T {resources.cpus_per_task} \
        -a {input.exons} -o {params.exon_counttab_file} {input.bam_file}
    gzip {params.exon_counttab_file}

    # Count reads mapped to introns
    singularity exec {singularity_image} featureCounts \
        -p --ignoreDup -t intron -g intron_id --extraAttributes gene_name,transcript_id,gene_type -T {resources.cpus_per_task} \
        -a {input.introns} -o {params.intron_counttab_file} {input.bam_file}
    gzip {params.intron_counttab_file}

    # Count reads mapped to enhancers
    singularity exec {singularity_image} featureCounts \
        -p --ignoreDup -t enhancer -g ID -F GFF -T {resources.cpus_per_task} \
        -a {input.enhancers} -o {params.enhancers_counttab_file} {input.bam_file}
    gzip {params.enhancers_counttab_file}

    # Count reads mapped to transcripts
    singularity exec {singularity_image} featureCounts \
        -p --ignoreDup -t transcript -g transcript_id --extraAttributes gene_name,transcript_id,gene_type -T {resources.cpus_per_task} \
        -a {input.transcripts} -o {params.transcript_counttab_file} {input.bam_file}
    gzip {params.transcript_counttab_file}

    md5sum {output.exon_counttab_file} {output.intron_counttab_file} {output.enhancers_counttab_file} {output.transcript_counttab_file} > {output.count_table_md5}
    '''


rule s02_count_reads_ready:
  input: 
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.md5',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.exonic_regions.tsv.gz',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.intronic_regions.tsv.gz',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.intergenic_enhancers.tsv.gz',
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.transcript_regions.tsv.gz',
  output:
    quant_dir / 'featureCounts/{group_id}/{group_id}.read_count.ready'
  shell:
    '''touch {output}'''
