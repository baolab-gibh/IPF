#!/bin/bash
# File: misc.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 16, 2024
# Updated:

project_dir=${HOME}/Documents/projects/wp_ipf

if [[ -e ${project_dir}/scripts/.env ]]; then
  source ${project_dir}/scripts/.env/bin/activate
fi

# Covert CRAM into BAM
samtools view -T ${project_dir}/outputs/references/genome/GRCh38.primary_assembly.genome.fa \
  -bo ${project_dir}/temps/test.bam ${project_dir}/temps/SRR25074724.mkdup.sncr.bqsr.cram


pints_visualizer \
  --bam ${project_dir}/temps/test.bam \
  --output-prefix ${project_dir}/temps/pints \
  --exp-type PROseq

pints_caller \
  --bw-mn ${project_dir}/temps/pints_mn.bw \
  --bw-pl ${project_dir}/temps/pints_pl.bw \
  --save-to ${project_dir}/temps \
  --file-prefix pints \
  --disable-small \
  --exp-type PROseq \
  --thread 4
