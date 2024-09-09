#!/bin/bash
# File: workflow.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: May 11, 2024
# Updated: Jun 03, 2024


project_dir=~/Documents/projects/wp_ipf

if [[ ${USER} =~ 'zhzhang_gibh' ]]; then
  conda deactivate
  source ~/tools/miniconda3/bin/activate
  mamba activate sm732 # Snakemake v7.32.4
elif [[ ${USER} =~ 'zzhang' && -d ${project_dir}/scripts/.env ]]; then
  source ${project_dir}/scripts/.env/bin/activate
fi

cd ${project_dir}/temps || return 1

#
## Main pipeline
#
n_jobs=20
# per_run=identify_tre_pints
prof_dir=${project_dir}/scripts/snakemake/configs
log_dir=${project_dir}/outputs/analysis/logs
for per_run in main_ppl identify_tre_pints; do
  if [[ ${per_run} == "main_ppl" ]]; then
    # Modulaizations
    work_dir=${project_dir}/outputs/analysis
    smk_file=${project_dir}/scripts/snakemake/main_ppl.smk
    cfg_file=${project_dir}/scripts/snakemake/main_ppl.yaml
  elif [[ ${per_run} == "identify_tre_pints" ]]; then
    # Identify TREs by PINTS
    work_dir=${project_dir}/outputs/analysis/tre_identification
    smk_file=${project_dir}/scripts/snakemake/identify_tre_pints.smk
    cfg_file=${project_dir}/scripts/snakemake/identify_tre_pints.yaml
  else
    echo "Error: Unknown module ${per_run}."
  fi

  # Job graph, not able to show the DAG of modules
  snakemake -d ${work_dir} -s ${smk_file} -j ${n_jobs} -C myconfigfile=${cfg_file} --profile ${prof_dir} --rulegraph --forceall | dot -Tpdf >| ${log_dir}/${per_run}.job_graph.pdf

  # Touch existing files
  snakemake -d ${work_dir} -s ${smk_file} -j ${n_jobs} -C myconfigfile=${cfg_file} --profile ${prof_dir} --touch 2>/dev/null

  # Dry-run
  snakemake -d ${work_dir} -s ${smk_file} -j ${n_jobs} -C myconfigfile=${cfg_file} --profile ${prof_dir} --dry-run >| ${log_dir}/${per_run}.dry_run.txt
  awk '/^Job stats:$/ {count++} count >= 2 {print}' ${log_dir}/${per_run}.dry_run.txt
  # cat ${log_dir}/${per_run}.dry_run.txt | tail -n 100
  # less -X ${log_dir}/${per_run}.dry_run.txt

  # grep -oE 'SRR[0-9]+.mkdup.sncr.bqsr.cram' ${log_dir}/${per_run}.dry_run.txt | cut -f1 -d. | sort -u | xargs -I% rm -r %/%.mkdup.sncr.bqsr.cram
  # Unlock current .snakemake folder
  snakemake -d ${work_dir} -s ${smk_file} -j ${n_jobs} -C myconfigfile=${cfg_file} --profile ${prof_dir} --unlock

 # Run
  if [[ ${USER} == 'zhzhang_gibh' ]]; then
    snakemake -d ${work_dir} -s ${smk_file} -j ${n_jobs} -C myconfigfile=${cfg_file} --profile ${prof_dir} 2>&1 | grep -v 'sacct: error:'
  elif [[ ${USER} == 'zzhang' ]]; then
    snakemake -k -d ${work_dir} -s ${smk_file} -j ${n_jobs} -C myconfigfile=${cfg_file} | grep -v 'sacct: error:'
  fi
done


#
## Prepare genotypes for cell line
#
# The genotyping array were download from https://www.encodeproject.org
# Download the IDAT files
work_dir=~/Documents/projects/resources/Encode/genotyping_array.cell_line
bpm_file=~/Documents/projects/resources/Illumina/InfiniumOmni5_4/1.2/HumanOmni5Exome-4v1-2_A.bpm
egt_file=~/Documents/projects/resources/Illumina/InfiniumOmni5_4/1.2/HumanOmni5Exome-4v1-2_A_ClusterFile.egt

cd "${work_dir}" || echo 'Fail to change directory'
xargs -L1 curl -O -J -L < ${work_dir}/Richard_Myers_etal.genotyping_array.cell_lines.txt

while read -r sid channel gid; do
  channel=$(sed 's/idat green channel/Grn/g; s/idat red channel/Red/g' <<<"${channel}")
  [[ ! -d ${work_dir}/idat/"${gid}" ]] && echo mkdir -p ${work_dir}/idat/"${gid}"
  echo mv "${sid}".idat ${work_dir}/idat/"${gid}/${gid}_${sid}_${channel}".idat
done < <(cut -f1,5,7 -d$'\t' ${work_dir}/metadata.tsv | grep -v File)

for per_sample in *; do
  n_idats=$(ls ${per_sample}/*.idat | wc -l)
  if [[ $n_idats -eq 2 ]]; then
    barcode=$(bcftools plugin gtc2vcf -i -g ${per_sample} 2>/dev/null | cut -d $'\t' -f 6 | grep -E '[0-9]+' | sort -u)
    position=$(bcftools plugin gtc2vcf -i -g ${per_sample} 2>/dev/null | cut -d $'\t' -f 8 | grep -E '[A-Z0-9]+' | sort -u)
    mv ${per_sample}/$(ls ${per_sample} | grep Grn.idat) ${per_sample}/${barcode}_${position}_Grn.idat
    mv ${per_sample}/$(ls ${per_sample} | grep Red.idat) ${per_sample}/${barcode}_${position}_Red.idat
    iaap-cli gencall ${bpm_file} ${egt_file} ../ped/${per_sample} -p -f ${per_sample}
  else
    continue
  fi
done


#
## Copy samples to baolab PC
#
if [[ ! -d ${project_dir}/inputs/fastqs ]]; then echo Making dir ... && mkdir -p ${project_dir}/inputs/fastqs; fi
cd "${project_dir}/inputs/fastqs"
IFS=","; while read -r gid x x x rr1 rr2; do
  if [[ -d ${gid} ]]; then
    echo existing "${gid}", skipping
  else
    mkdir -p "${gid}" && ln -s "${rr1}" "${rr2}" "${gid}" && echo "linked ${gid}"
  fi
done < <(grep -v -e read_r1 ${project_dir}/outputs/overview/sample_information/sample_info.human.paired_end.txt | head -n 800)


ls --color=never ${project_dir}/inputs/fastqs/ \
  | xargs -P8 -I% bash -c 'echo %,$(zgrep -c ^@ %/*fastq.gz),$(zgrep -c ^@ %/*fastq.2.gz),$(du -h %|cut -f1)' \
  >| ${project_dir}/inputs/fastqs/sample_size.csv


# Split the samples in to small batches
mkdir -p ${project_dir}/outputs/overview/sample_information/chunks
cd ${project_dir}/outputs/overview/sample_information/chunks
rm -f x*.txt
grep -v group_id ${project_dir}/outputs/overview/sample_information/sample_info.human.paired_end.txt \
  | split -l 25 --additional-suffix=.txt --numeric-suffixes=1
for x in x*.txt; do mv ${x} ${x/x/chunk_}; done


# Collect failed jobs
for per_log_dir in ${project_dir}/outputs/analysis/logs/*; do
  if [[ -d ${per_log_dir} ]]; then
    grep -iw error ${per_log_dir}/*.out | cut -f1 -d: | sort -u
  fi
done | xargs -I% rm -f %

watch squeue -O NAME:40,TimeLimit:30,TimeUsed \| tail -n 40
cat /public/home/zhzhang_gibh/Documents/projects/wp_ipf/outputs/overview/sample_information/sample_info.human.paired_end.cram_path.txt

key_file=${HOME}/.ssh/baolab_id_rsa
IFS=','; read -r _ _ _ src_fp < <(grep --color=none SRR17735640 /public/home/zhzhang_gibh/Documents/projects/wp_ipf/outputs/overview/sample_information/sample_info.human.paired_end.cram_path.txt)
rsync -avzhL -e "ssh -i ${key_file}" zzhang@192.168.143.54:${src_fp} /public/home/zhzhang_gibh/Documents/projects/wp_ipf/outputs/analysis/read_alignment/alignment/SRR17735640/SRR17735640.mkdup.sncr.bqsr.cram
