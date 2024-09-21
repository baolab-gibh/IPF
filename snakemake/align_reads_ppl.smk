#!/usr/bin/env snakemake
# File: align_reads_ppl.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 23, 2024
# Updated: Jun 13, 2024


import sys, csv, pathlib, yaml, benedict
from snakemake.utils import min_version

# Ensure the snakemake knows how to handle moduliazations.
min_version("6.0")


#
## Resources
#
bconfig, sample_info_dict = config
project_dir = bconfig.get('project_dir')

# Container
singularity_image = bconfig.get('singularity_image')
assert singularity_image.exists()

# Output controlling
output_dir = bconfig.get('output_dir')
fastq_dir = output_dir / 'read_alignment/fastqs'
quality_control_dir = output_dir / 'read_alignment/quality_control'
alignment_dir = output_dir / 'read_alignment/alignment'

# Other parameters
zero_down_bams = bconfig.get('zero_down_bams')


# Final outputs
all_groups = list(sample_info_dict.keys())
final_output = expand(alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.{per_ext}', group_id=all_groups, per_ext=['bam', 'bam.bai'])


# Wildcard constraints
wildcard_constraints:
  sample_id = r'\w+', group_id = r'\w+'


#
## Rules
#
rule all:
  input: final_output


rule s01_build_genome_index:
  input:
    fasta_file = bconfig.reference_genome, gff_file = bconfig.genomic_feature,
  output:
    genome_dir = directory(bconfig.get("reference_genome_star_index"))
  resources:
    cpus_per_task = 4, mem = '64G'
  shell:
    '''
    singularity exec {singularity_image} STAR \
        --runMode genomeGenerate \
        --runThreadN {resources.cpus_per_task} \
        --genomeFastaFiles {input.fasta_file} \
        --genomeDir {output.genome_dir} \
        --sjdbGTFfile {input.gff_file} \
        --sjdbOverhang 99
    '''


rule s01_preproc_fastq:
  input:
    lambda wc:
      [fastq_dir / '{group_id}.{sample_id}.R1.fastq.gz', fastq_dir / '{group_id}.{sample_id}.R2.fastq.gz'] \
      if sample_info_dict[wc.sample_id]['paired_end'] else [fastq_dir / '{group_id}.{sample_id}.fastq.gz']
  output:
    clean_fastq_list = quality_control_dir / '{group_id}/{group_id}.{sample_id}.clean_fastq.txt',
    html_report = quality_control_dir / '{group_id}/{group_id}.{sample_id}.html',
    json_report = quality_control_dir / '{group_id}/{group_id}.{sample_id}.json'
  resources:
    cpus_per_task = 2, mem = '4G', time = '2:00:00'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    read -r raw_read_r1 raw_read_r2 <<<'{input}'
    if [[ -n ${{raw_read_r1}} && -n ${{raw_read_r2}} ]]; then
      clean_read_r1=${{raw_read_r1/.R1.fastq.gz/.clean.R1.fastq.gz}}
      clean_read_r2=${{raw_read_r2/.R2.fastq.gz/.clean.R2.fastq.gz}}
      singularity exec {singularity_image} fastp \
          -i ${{raw_read_r1}} \
          -o ${{clean_read_r1}} \
          -I ${{raw_read_r2}} \
          -O ${{clean_read_r2}} \
          -h {output.html_report} \
          -j {output.json_report} \
          -w {resources.cpus_per_task}
      echo "${{clean_read_r1}} ${{clean_read_r2}}" > {output.clean_fastq_list}
    elif [[ -n ${{raw_read_r1}} ]]; then
      clean_read_r1=${{raw_read_r1/.fastq.gz/.clean.fastq.gz}}
      singularity exec {singularity_image} fastp \
          -i ${{raw_read_r1}} \
          -o ${{clean_read_r1}} \
          -h {output.html_report} \
          -j {output.json_report} \
          -w {resources.cpus_per_task}
      echo ${{clean_read_r1}} > {output.clean_fastq_list}
    else
      echo "Error: Both R1 and R2 reads are missing."
      exit 1
    fi

    for per_file in {input}; do
      real_path=$(readlink -f ${{per_file}})
      rm -f ${{real_path}} && touch ${{real_path}}
    done
    touch {output.clean_fastq_list} {output.html_report} {output.json_report}
    '''


rule s02_align_reads:
  input:
    clean_fastq_list = quality_control_dir / '{group_id}/{group_id}.{sample_id}.clean_fastq.txt',
    genome_dir = bconfig.reference_genome_star_index,
  output:
    bam_file = alignment_dir / '{group_id}/{group_id}.{sample_id}.Aligned.sortedByCoord.out.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.{sample_id}.Aligned.sortedByCoord.out.bam.bai',
    log_out = alignment_dir / '{group_id}/{group_id}.{sample_id}.Log.out',
    log_final_out = alignment_dir / '{group_id}/{group_id}.{sample_id}.Log.final.out',
    log_progress_out = alignment_dir / '{group_id}/{group_id}.{sample_id}.Log.progress.out',
    log_sj_out_tab = alignment_dir / '{group_id}/{group_id}.{sample_id}.SJ.out.tab',
  resources:
    cpus_per_task = 4, mem = '64G', time = '5:59:00'
  params:
    out_prefix = lambda wc: alignment_dir / wc.group_id / (wc.group_id + '.' + wc.sample_id + '.'),
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    read -r clean_r1 clean_r2 < {input.clean_fastq_list}
    singularity exec {singularity_image} STAR \
        --runMode alignReads \
        --genomeDir {input.genome_dir} \
        --twopassMode Basic \
        --readFilesIn ${{clean_r1}} ${{clean_r2}} \
        --readFilesCommand gunzip -c \
        --runThreadN {resources.cpus_per_task} \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN {resources.cpus_per_task} \
        --outSAMattrRGline ID:{wildcards.sample_id} SM:{wildcards.group_id} PL:illumina \
        --outSAMattributes NH HI AS nM NM \
        --outSAMstrandField intronMotif \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outFilterMatchNmin 3 \
        --outBAMcompression 6 \
        --outFileNamePrefix {params.out_prefix}
    singularity exec {singularity_image} samtools index -@ {resources.cpus_per_task} {output.bam_file}

    for per_file in ${{clean_r1}} ${{clean_r2}}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output.bai_file} {output.bam_file} {output.log_out} {output.log_final_out} {output.log_progress_out} {output.log_sj_out_tab}
    '''


rule s03_merge_bam_alignment:
  input:
    lambda wc: expand(
      alignment_dir / wc.group_id / (wc.group_id + '.{sample_id}.Aligned.sortedByCoord.out.bam'),
      sample_id=[x for x in sample_info_dict.keys() if sample_info_dict[x]['group_id'] == wc.group_id]
    )
  output:
    bam_file = alignment_dir / '{group_id}/{group_id}.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.bam.bai'
  resources:
    time = '59:59', mem = '8G', cpus_per_task = 4
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    n_bams=$(echo {input} | wc -w)
    if [[ ${{n_bams}} -eq 1 ]]; then
      singularity exec {singularity_image} samtools sort --no-PG -@ {resources.cpus_per_task} -o {output.bam_file} {input}
    else
      singularity exec {singularity_image} samtools merge --no-PG -c --no-PG -@ {resources.cpus_per_task} {input} \
        | singularity exec {singularity_image} samtools sort --no-PG -@ 2 -o {output.bam_file}
    fi
    singularity exec {singularity_image} samtools index -@ {resources.cpus_per_task} {output.bam_file}

    for per_file in {input}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output.bai_file} {output.bam_file}
    '''


rule s04_mark_duplicates:
  input:
    bam_file = alignment_dir / '{group_id}/{group_id}.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.bam.bai'
  output:
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.bam.bai',
    metrics_file = alignment_dir / '{group_id}/{group_id}.mkdup.metrics.txt'
  resources:
    time = '24:00:00', mem = '8G', cpus_per_task = 8
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} samtools sort -n -@ 2 -O SAM {input.bam_file} \
      | singularity exec {singularity_image} samtools fixmate -@ 2 -cm -O SAM - - \
      | singularity exec {singularity_image} samtools sort -@ 2 -O SAM \
      | singularity exec {singularity_image} samtools markdup -c --no-PG -@ 2 -f {output.metrics_file} - {output.bam_file}
    singularity exec {singularity_image} samtools index -@ {resources.cpus_per_task} {output.bam_file}

    for per_file in {input.bam_file} {input.bai_file}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output.bai_file} {output.bam_file}
    '''


rule s05_split_n_cigar_reads:
  input:
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.bam.bai',
    reference_genome = bconfig.reference_genome,
    reference_genome_dict = bconfig.reference_genome_dict,
    reference_genome_index = bconfig.reference_genome_index,
  output:
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bam.bai'
  resources:
    cpus_per_task = 1, time = '24:00:00'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} gatk SplitNCigarReads \
        -OBI false \
        -R {input.reference_genome} \
        -I {input.bam_file} \
        -O {output.bam_file}
    singularity exec {singularity_image} samtools index -@ {resources.cpus_per_task} {output.bam_file}

    for per_file in {input.bam_file} {input.bai_file}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output.bai_file} {output.bam_file}
    '''


rule s06_base_recalibrator:
  input:
    reference_genome = bconfig.reference_genome,
    reference_variants_vcf = bconfig.reference_variants_vcf,
    reference_variants_vcf_tbi = bconfig.reference_variants_vcf_tbi,
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bam.bai',
  output:
    recal_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.recal.table',
    recal_md5_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.recal.table.md5',
  resources:
    time = '3:00:00', mem = '8G'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} gatk BaseRecalibrator \
        -R {input.reference_genome} \
        -I {input.bam_file} \
        -O {output.recal_file} \
        --known-sites {input.reference_variants_vcf}
    md5sum {output.recal_file} > {output.recal_md5_file}
    '''


rule s07_apply_bqsr:
  input:
    reference_genome = bconfig.reference_genome,
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bam.bai',
    recal_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.recal.table'
  output:
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.bam.bai'
  resources:
    time = '6:00:00', mem = '4G', cpus_per_task = 2
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} gatk ApplyBQSR \
        -OBI false \
        -R {input.reference_genome} \
        -I {input.bam_file} \
        --use-original-qualities \
        --bqsr-recal-file {input.recal_file} \
        -O {output.bam_file}
    singularity exec {singularity_image} samtools index -@ {resources.cpus_per_task} {output.bam_file}

    for per_file in {input.bam_file} {input.bai_file}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output.bai_file} {output.bam_file}
    '''
