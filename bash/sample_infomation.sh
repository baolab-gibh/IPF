#!/bin/bash
# File: sample_infomation.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Aug 15, 2024
# Updated:

project_dir=/home/zzhang/Documents/projects/wp_ipf

if [[ ${USER} =~ 'zhzhang_gibh' ]]; then
  conda deactivate
  source ~/tools/miniconda3/bin/activate
  mamba activate sm732 # Snakemake v7.32.4
elif [[ ${USER} =~ 'zzhang' && -d ${project_dir}/scripts/.env_py_3.10 ]]; then
  source ${project_dir}/scripts/.env_py_3.10/bin/activate
fi


while read -r start end; do
  pysradb metadata --saveto "sample_bam_info_batch_${start}-${end}.csv" $(awk -F, -v start="$start" -v end="$end" 'start < NR && NR <= end {print $1}' sample_bam_info.txt)
done < <(paste <(seq 0 200 "$(grep -c . sample_bam_info.txt)") <(seq 200 200 $(( "$(grep -c . sample_bam_info.txt)" + 200 ))))
