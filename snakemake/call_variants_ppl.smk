#!/usr/bin/env snakemake
# File: call_variants_ppl.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 23, 2024
# Updated: Jun 13, 2024


import sys, csv, pathlib, yaml, benedict
from snakemake.utils import min_version

# Ensure the snakemake knows how to handle moduliazations.
min_version('6.0')


# Utility functions
def time_per_chunk(wc):
  chunk = wc.per_chunk
  if chunk in ['chr' + str(i) for i in range(1, 4)]:
    return '5:00:00'
  elif chunk in ['chr' + str(i) for i in range(4, 8)]:
    return '4:00:00'
  elif chunk in ['chrY', 'Y']:
    return '30:00'
  return '3:00:00'


#
## Resources
#
bconfig, sample_info_dict = config
project_dir = bconfig.get('project_dir')
chunk_token = bconfig.get('chunk_token')

# Container
singularity_image = bconfig.get('singularity_image')
assert singularity_image.exists()

# I/O controlling
input_dir = bconfig.get('output_dir') / 'read_alignment'
alignment_dir = input_dir / 'alignment'

output_dir = bconfig.get('output_dir') / 'variant_calling'
variant_dir = output_dir / 'variant'
# annotation_dir = output_dir / 'annotation'
# plink_dir = output_dir / 'plink'


# Final outputs
all_groups = list(sample_info_dict.keys())
all_chunks = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
final_output = [
    variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz'),
    variant_dir / (chunk_token + '.gatk_ppl.snps.clean.variant_calling_detail_metrics'),
    variant_dir / (chunk_token + '.gatk_ppl.snps.clean.variant_calling_summary_metrics'),
]

final_output = variant_dir / '{group_id}/{group_id}.haplotype_caller_ready'


# Wildcard constraints
wildcard_constraints:
  sample_id = r'\w+', group_id = r'\w+'


#
## Rules
#
rule all:
  input: final_output


# Not included in the current pipeline, due to missing R/gplot in the current version of singularity image (May 20, 2024)
rule s01_analyze_covariates:
  input:
    recal_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.recal.table'
  output:
    covar_plot_file = variant_dir / '{group_id}/{group_id}.covariates_report.pdf'
  resources:
    time = '10:00', mem = '1G'
  shell:
    '''
    singularity exec {singularity_image} gatk AnalyzeCovariates -bqsr {input.recal_file} -plots {output.covar_plot_file}
    '''


rule s01_haplotype_caller_perchunk:
  input:
    bam_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.bam',
    bai_file = alignment_dir / '{group_id}/{group_id}.mkdup.sncr.bqsr.bam.bai',
    reference_genome = bconfig.reference_genome,
    reference_variants_vcf = bconfig.reference_variants_vcf,
    reference_variants_vcf_tbi = bconfig.reference_variants_vcf_tbi,
  output:
    vcf_file = variant_dir / '{group_id}/per_chunk/{group_id}.{per_chunk}.g.vcf.gz',
    tbi_file = variant_dir / '{group_id}/per_chunk/{group_id}.{per_chunk}.g.vcf.gz.tbi',
    md5_file = variant_dir / '{group_id}/per_chunk/{group_id}.{per_chunk}.g.vcf.gz.md5',
  resources:
    time = time_per_chunk, mem = '8G', cpus_per_task = 1
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} gatk HaplotypeCaller \
        -R {input.reference_genome} \
        -I {input.bam_file} \
        -L {wildcards.per_chunk} \
        --emit-ref-confidence GVCF \
        --dont-use-soft-clipped-bases \
        --dbsnp {input.reference_variants_vcf} \
        -O {output.vcf_file}
    md5sum {output.vcf_file} > {output.md5_file}
    '''


rule s02_haplotype_caller_ready:
  input:
    lambda wc: expand(variant_dir / wc.group_id / 'per_chunk' / (wc.group_id + '.{per_chunk}.g.vcf.gz'), per_chunk=all_chunks),
  output:
    variant_dir / '{group_id}/{group_id}.haplotype_caller_ready'
  shell:
    '''touch {output}'''


rule s02_concat_gvcfs:
  input:
    lambda wc: expand(variant_dir / wc.group_id / 'per_chunk' / (wc.group_id + '.{per_chunk}.g.vcf.gz'), per_chunk=all_chunks),
  output:
    vcf_file = variant_dir / '{group_id}/{group_id}.g.vcf.gz',
    vcf_md5_file = variant_dir / '{group_id}/{group_id}.g.vcf.gz.md5',
    tbi_file = variant_dir / '{group_id}/{group_id}.g.vcf.gz.tbi',
    tbi_md5_file = variant_dir / '{group_id}/{group_id}.g.vcf.gz.tbi.md5',
  resources:
    time = '30:00', mem = '4G'
  params:
    output_dir = lambda wildcards, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.output_dir}

    singularity exec {singularity_image} bcftools concat --no-version -o {output.vcf_file} {input}
    singularity exec {singularity_image} bcftools index -ft {output.vcf_file}

    for per_file in {input}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done

    md5sum {output.vcf_file} > {output.vcf_md5_file}
    md5sum {output.tbi_file} > {output.tbi_md5_file}
    touch {output.vcf_file} {output.vcf_md5_file} {output.tbi_file} {output.tbi_md5_file}
    '''


# ---------- Previous version to estimate genotypes per chunk or chromosome. -------------------------------------------
# rule s03_combine_gvcfs_perchunk:
#   input:
#     vcf_file = lambda wc: expand(variant_dir / '{group_id}' / 'per_chunk' / ('{group_id}.' + wc.per_chunk + '.g.vcf.gz'), group_id=all_groups),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.g.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.g.vcf.gz.tbi'),
#   resources:
#     time = '5:00:00', mem = '4G'
#   params:
#     input_string = lambda wc, input: ' -V '.join(input.vcf_file),
#     output_dir = lambda wildcards, output: Path(output[0]).parent
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     singularity exec {singularity_image} gatk CombineGVCFs \
#         -R {input.reference_genome} \
#         -V {params.input_string} \
#         -O {output.vcf_file}
#     singularity exec {singularity_image} bcftools index -ft {output.vcf_file}
# 
#     for per_file in {input.vcf_file}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     done
#     touch {output.vcf_file} {output.tbi_file}
#     '''
# 
# 
# rule s04_genotype_gvcfs_perchunk:
#   input:
#     vcf_file = variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.g.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.g.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.vcf.gz.tbi'),
#   resources:
#     time = '12:00:00', mem = '4G'
#   params:
#     output_dir = lambda wildcards, output: Path(output[0]).parent
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     singularity exec {singularity_image} gatk GenotypeGVCFs \
#         -R {input.reference_genome} \
#         -V {input.vcf_file} \
#         -O {output.vcf_file}
#     singularity exec {singularity_image} bcftools index -ft {output.vcf_file}
# 
#     for per_file in {input.vcf_file} {input.tbi_file}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     fi
#     touch {output.vcf_file} {output.tbi_file}
#     '''
# 
# 
# rule s05_concat_gvcfs:
#   input:
#     vcf_file = expand(variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.vcf.gz'), per_chunk=all_chunks),
#     tbi_file = expand(variant_dir / (chunk_token + '{per_chunk}.gatk_ppl.vcf.gz.tbi'), per_chunk=all_chunks),
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz.tbi'),
#   resources:
#     time = '30:00', mem = '4G'
#   params:
#     output_dir = lambda wildcards, output: Path(output[0]).parent
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     singularity exec {singularity_image} bcftools concat --no-version -o {output.vcf_file} {input}
#     singularity exec {singularity_image} bcftools index -ft {output.vcf_file}
# 
#     for per_file in {input}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     done
#     touch {output.tbi_file} {output.vcf_file}
#     '''
# 
# 
# rule s06_select_variants_snp:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz.tbi'),
#   resources:
#     time = '1:30:00', mem = '2G'
#   shell:
#     '''
#     singularity exec {singularity_image} gatk SelectVariants \
#         -V {input.vcf_file} \
#         --select-type-to-include SNP \
#         -O {output.vcf_file}
#     '''
# 
# 
# rule s07_variant_filtration_snp:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz.tbi'),
#   resources:
#     time = '1:30:00', mem = '2G'
#   shell:
#     '''
#     singularity exec {singularity_image} gatk VariantFiltration \
#         -R {input.reference_genome} \
#         -V {input.vcf_file} \
#         --window 35 \
#         --cluster 3 \
#         --filter-expression "QD < 2.0" --filter-name "QD2" \
#         --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#         --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#         --filter-expression "FS > 60.0" --filter-name "FS60" \
#         --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#         --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#         --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#         -O {output.vcf_file}
# 
#     for per_file in {input.vcf_file} {input.tbi_file}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     fi
#     touch {output.vcf_file} {output.tbi_file}
#     '''
# 
# 
# rule s08_collect_variant_call_metrics_snp:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz.tbi'),
#     reference_variants_vcf = bconfig.reference_variants_vcf,
#     reference_variants_vcf_tbi = bconfig.reference_variants_vcf_tbi,
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#   output:
#     detail_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.variant_calling_detail_metrics'),
#     summary_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.variant_calling_summary_metrics'),
#   params:
#     out_prefix = variant_dir / (chunk_token + '.gatk_ppl.snps.clean')
#   resources:
#     time = '30:00', mem = '2G'
#   shell:
#     '''
#     singularity exec {singularity_image} gatk CollectVariantCallingMetrics \
#         -I {input.vcf_file} \
#         --DBSNP {input.reference_variants_vcf} \
#         --SEQUENCE_DICTIONARY {input.reference_genome_dict} \
#         -O {params.out_prefix}
#     '''
# 
# 
# rule s06_select_variants_indel: # Indels were not included in the current pipeline (May 20, 2024)
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz.tbi'),
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz.tbi'),
#   shell:
#     '''
#     singularity exec {singularity_image} gatk SelectVariants \
#         -V {input.vcf_file} \
#         --select-type-to-include INDEL \
#         -O {output.vcf_file}
#     '''
# 
# 
# rule s07_variant_filtration_indel:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz.tbi'),
#   shell:
#     '''
#     singularity exec {singularity_image} gatk VariantFiltration \
#         -R {input.reference_genome} \
#         -V {input.vcf_file} \
#         --window 35 \
#         --cluster 3 \
#         --filter-expression "QD < 2.0" --filter-name "QD2" \
#         --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#         --filter-expression "FS > 200.0" --filter-name "FS200" \
#         --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
#         -O {output.vcf_file}
#     '''
# 
# 
# rule s08_collect_variant_call_metrics_indel:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz.tbi'),
#     reference_variants_vcf = bconfig.reference_variants_vcf,
#     reference_variants_vcf_tbi = bconfig.reference_variants_vcf_tbi,
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#   output:
#     detail_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_detail_metrics'),
#     summary_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_summary_metrics'),
#   params:
#     detail_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_detail_metrics'),
#     summary_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_summary_metrics'),
#   shell:
#     '''
#     singularity exec {singularity_image} gatk CollectVariantCallingMetrics \
#         -I {input.vcf_file} \
#         --DBSNP {input.reference_variants_vcf} \
#         --SEQUENCE_DICTIONARY {input.reference_genome_dict} \
#         -O {params.out_prefix}
#     '''



# --------- Previouse version of calling variants for each sample. -----------------------------------------------------
# rule s02_concat_gvcfs:
#   input:
#     lambda wc: expand(variant_dir / wc.group_id / 'per_chunk' / (wc.group_id + '.{per_chunk}.g.vcf.gz'), per_chunk=all_chunks),
#   output:
#     vcf_file = variant_dir / '{group_id}/{group_id}.g.vcf.gz',
#     tbi_file = variant_dir / '{group_id}/{group_id}.g.vcf.gz.tbi',
#   resources:
#     time = '30:00', mem = '4G'
#   params:
#     output_dir = lambda wildcards, output: Path(output[0]).parent
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     singularity exec {singularity_image} bcftools concat --no-version -o {output.vcf_file} {input}
#     singularity exec {singularity_image} bcftools index -ft {output.vcf_file}
# 
#     for per_file in {input}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     done
#     touch {output.tbi_file} {output.vcf_file}
#     '''
# 
# 
# rule s03_combine_gvcfs:
#   input:
#     vcf_file = expand(variant_dir / '{group_id}/{group_id}.g.vcf.gz', group_id=all_groups),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.g.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.g.vcf.gz.tbi'),
#   resources:
#     time = '24:00:00', mem = '4G'
#   params:
#     input_string = lambda wc, input: ' -V '.join(input.vcf_file),
#     output_dir = lambda wildcards, output: Path(output[0]).parent
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     singularity exec {singularity_image} gatk CombineGVCFs \
#         -R {input.reference_genome} \
#         -V {params.input_string} \
#         -O {output.vcf_file}
#     singularity exec {singularity_image} bcftools index -ft {output.vcf_file}
# 
#     for per_file in {input.vcf_file}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     done
#     touch {output.vcf_file} {output.tbi_file}
#     '''
# 
# 
# rule s04_genotype_gvcfs:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.g.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.g.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz.tbi'),
#   resources:
#     time = '12:00:00', mem = '4G'
#   params:
#     output_dir = lambda wildcards, output: Path(output[0]).parent
#   shell:
#     '''
#     mkdir -p {params.output_dir}
# 
#     singularity exec {singularity_image} gatk GenotypeGVCFs \
#         -R {input.reference_genome} \
#         -V {input.vcf_file} \
#         -O {output.vcf_file}
#     singularity exec {singularity_image} bcftools index -ft {output.vcf_file}
# 
#     for per_file in {input.vcf_file} {input.tbi_file}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     fi
#     touch {output.vcf_file} {output.tbi_file}
#     '''
# 
# 
# rule s05_select_variants_snp:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz.tbi'),
#   resources:
#     time = '1:30:00', mem = '2G'
#   shell:
#     '''
#     singularity exec {singularity_image} gatk SelectVariants \
#         -V {input.vcf_file} \
#         --select-type-to-include SNP \
#         -O {output.vcf_file}
#     '''
# 
# 
# rule s06_variant_filtration_snp:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz.tbi'),
#   resources:
#     time = '1:30:00', mem = '2G'
#   shell:
#     '''
#     singularity exec {singularity_image} gatk VariantFiltration \
#         -R {input.reference_genome} \
#         -V {input.vcf_file} \
#         --window 35 \
#         --cluster 3 \
#         --filter-expression "QD < 2.0" --filter-name "QD2" \
#         --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#         --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#         --filter-expression "FS > 60.0" --filter-name "FS60" \
#         --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#         --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#         --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#         -O {output.vcf_file}
# 
#     for per_file in {input.vcf_file} {input.tbi_file}; do
#       rm -f ${{per_file}} && touch ${{per_file}}
#     fi
#     touch {output.vcf_file} {output.tbi_file}
#     '''
# 
# 
# rule s07_collect_variant_call_metrics_snp:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.vcf.gz.tbi'),
#     reference_variants_vcf = bconfig.reference_variants_vcf,
#     reference_variants_vcf_tbi = bconfig.reference_variants_vcf_tbi,
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#   output:
#     detail_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.variant_calling_detail_metrics'),
#     summary_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.snps.clean.variant_calling_summary_metrics'),
#   params:
#     out_prefix = variant_dir / (chunk_token + '.gatk_ppl.snps.clean')
#   resources:
#     time = '30:00', mem = '2G'
#   shell:
#     '''
#     singularity exec {singularity_image} gatk CollectVariantCallingMetrics \
#         -I {input.vcf_file} \
#         --DBSNP {input.reference_variants_vcf} \
#         --SEQUENCE_DICTIONARY {input.reference_genome_dict} \
#         -O {params.out_prefix}
#     '''
# 
# 
# rule s05_select_variants_indel: # Indels were not included in the current pipeline (May 20, 2024)
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.vcf.gz.tbi'),
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz.tbi'),
#   shell:
#     '''
#     singularity exec {singularity_image} gatk SelectVariants \
#         -V {input.vcf_file} \
#         --select-type-to-include INDEL \
#         -O {output.vcf_file}
#     '''
# 
# 
# rule s06_variant_filtration_indel:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.vcf.gz.tbi'),
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#     reference_genome_index = bconfig.reference_genome_index,
#   output:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz.tbi'),
#   shell:
#     '''
#     singularity exec {singularity_image} gatk VariantFiltration \
#         -R {input.reference_genome} \
#         -V {input.vcf_file} \
#         --window 35 \
#         --cluster 3 \
#         --filter-expression "QD < 2.0" --filter-name "QD2" \
#         --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#         --filter-expression "FS > 200.0" --filter-name "FS200" \
#         --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
#         -O {output.vcf_file}
#     '''
# 
# 
# rule s07_collect_variant_call_metrics_indel:
#   input:
#     vcf_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz'),
#     tbi_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.vcf.gz.tbi'),
#     reference_variants_vcf = bconfig.reference_variants_vcf,
#     reference_variants_vcf_tbi = bconfig.reference_variants_vcf_tbi,
#     reference_genome = bconfig.reference_genome,
#     reference_genome_dict = bconfig.reference_genome_dict,
#   output:
#     detail_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_detail_metrics'),
#     summary_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_summary_metrics'),
#   params:
#     detail_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_detail_metrics'),
#     summary_metrics_file = variant_dir / (chunk_token + '.gatk_ppl.indels.clean.variant_calling_summary_metrics'),
#   shell:
#     '''
#     singularity exec {singularity_image} gatk CollectVariantCallingMetrics \
#         -I {input.vcf_file} \
#         --DBSNP {input.reference_variants_vcf} \
#         --SEQUENCE_DICTIONARY {input.reference_genome_dict} \
#         -O {params.out_prefix}
#     '''
