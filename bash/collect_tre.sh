#!/bin/bash
# File: collect_tre.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 03, 2024
# Updated: Aug 03, 2024

# Steps to collect transcriptional regulatory elements (TREs) estimated using dREG.
# ChromHMM were adopted to merge the TREs in a region.
# ChromHMM is a tool learns chromatin-state signatures using a multivariate hidden Markov model (HMM) that explicitly 
# models the combinatorial presence or absence of each mark.
# More detailes are availabe at [Ernst and Kellis](http://dx.doi.org/10.1038/nprot.2017.124) or the official website 
# [ChromHMM](https://compbio.mit.edu/ChromHMM)

project_dir=/home/zzhang/Documents/projects/wp_ipf
chromhmm_dir=/home/zzhang/tools/ChromHMM/

if [[ ${USER} =~ 'zhzhang_gibh' ]]; then
  conda deactivate
  source ~/tools/miniconda3/bin/activate
  mamba activate sm732 # Snakemake v7.32.4
elif [[ ${USER} =~ 'zzhang' && -d ${project_dir}/scripts/.env_py_3.10 ]]; then
  source ${project_dir}/scripts/.env_py_3.10/bin/activate
fi

# Create binarized_bed. The dREG BEDs contains a score indicating the TRE potential
# sample_id=SRR11785045
create_binned_bed() {
  local info_file="$1"; local bin_size=100
  local project_dir=/home/zzhang/Documents/projects/wp_ipf
  local chromhmm_dir=/home/zzhang/tools/ChromHMM
  IFS=","; while read -r sample_id _ cell_type _; do
    local out_dir="${project_dir}/outputs/analysis/tre_identification/harmonization/binarized_bed"
    local in_file="${project_dir}/outputs/analysis/tre_identification/tre_results/${sample_id}/${sample_id}.dREG.peak.full.bed.gz"
    if [[ ! -e "${in_file}" ]]; then continue; fi
    if [[ -e ${out_dir}/chr1/${sample_id}.chr1.chromhmm_binary.txt ]]; then continue; fi
    python "${project_dir}/scripts/py3/create_chromhmm_binfile.py" \
      -i "${sample_id}" -b "${bin_size}" -t "${cell_type}" -o "${out_dir}" "${in_file}"

    # while [[ 1 -le 2 ]]; do
    #   n_jobs=$(pgrep python | wc -l)
    #   echo "${n_jobs} running ..."
    #   if [[ ${n_jobs} -lt 2 ]]; then break; else echo "${n_jobs} running ..."; fi
    # done
  done < "${info_file}"
}

create_binned_bed "${project_dir}/outputs/overview/sample_information/sample_info.human.paired_end.sample_cellline_celltype_tissue.txt"

# Learn chromatin states. NB: this step needs to be justified
bin_size=100
for per_cell in Acute_myeloid_leukemia Burkitt_lymphoma Embryonic_stem_cells Erythroid_progenitor_cells Histiocytic_cells Rhabdoid_tumor Breast_adenocarcinoma Chronic_myeloid_leukemia Epithelial_cells Fibroblast Monocyte Stromal_cells Breast_carcinoma Colorectal_carcinoma Erythroid Fibrosarcoma Neuroblastoma; do
  for per_chr in chr{{1..22},X,Y}; do
    echo "[I]: Working on ${per_cell}/${per_chr}"
    mkdir -p ${project_dir}/outputs/analysis/tre_identification/harmonization/chromatin_states/${per_cell}/${per_chr}
    java -Xmx16G -jar ${chromhmm_dir}/ChromHMM.jar LearnModel \
      -b ${bin_size} \
      -nobrowser -noautoopen -noenrich -noimage \
      ${project_dir}/outputs/analysis/tre_identification/harmonization/binarized_bed/${per_cell}/${per_chr} \
      ${project_dir}/outputs/analysis/tre_identification/harmonization/chromatin_states/${per_cell}/${per_chr} \
      2 hg38
  done
done


# Estimate enrichment of identified TRE in IPF/COPD GWAS loci
#QTLtools fdensity --help
for per_gwas in IPF COPD; do
  if [[ ${per_gwas} == 'IPF' ]]; then
    gwas_sumstats=${project_dir}/inputs/public_gwas/GCST90399721.h.tsv.gz
  elif [[ ${per_gwas} == 'COPD' ]]; then
    gwas_sumstats=${project_dir}/inputs/public_gwas/GCST90399695.h.tsv.gz
  else
    continue
  fi

  selected_gwas_loci=${project_dir}/temps/${per_gwas}.significant_vars.bed
  if [[ ! -e ${selected_gwas_loci} ]]; then
    zcat ${gwas_sumstats} | awk -F$'\t' 'NR>1 && $8<5e-8 { OFS="\t"; if ($10~/rs/) {rsid=$9} else {rsid=$1"_"$2"_"$3"_"$4}; print $1, $2-1, $2, rsid, "NA", "-"}' | sed 's/^/chr/g' > ${selected_gwas_loci}
  fi

  for per_anno in Breast_carcinoma Colorectal_carcinoma Erythroid Fibrosarcoma Neuroblastoma Acute_myeloid_leukemia Burkitt_lymphoma Embryonic_stem_cells Erythroid_progenitor_cells Histiocytic_cells Rhabdoid_tumor Breast_adenocarcinoma Chronic_myeloid_leukemia Epithelial_cells Fibroblast Monocyte Stromal_cells; do
    func_annots=${project_dir}/outputs/analysis/tre_identification/harmonization/chromatin_states/${per_anno}/${per_anno}.segments.bed

    selected_func_annots=${project_dir}/temps/${per_anno}.bed
    if [[ ! -e ${selected_func_annots} ]]; then cut -f1-3 ${func_annots} > ${selected_func_annots}; fi

    # QTLtools fdensity \
    #   --qtl ${selected_gwas_loci} \
    #   --bed ${selected_func_annots} \
    #   --out ${project_dir}/outputs/analysis/functional_annotation/TREs/${per_anno}_${per_gwas}.erichment.txt

    bedtools intersect \
      -a ${selected_gwas_loci} \
      -b ${selected_func_annots} \
      -wa -wb \
      > ${project_dir}/outputs/analysis/functional_annotation/TREs/${per_anno}_${per_gwas}.overlaps.txt
  done
done


# Misc
for per_chr in chr{{2..22},X,Y}; do
  cd ${per_chr}
  for x in *.txt; do
    celltype=$(grep -m1 "Enhancer" "$x" | rev | cut -f2- -d_ | rev)
    mkdir -p "${celltype}" && mv "$x" "${celltype}"
  done
  cd ..
done

for per_chr in chr{{1..22},X,Y}; do
  for per_celltype in Acute_myeloid_leukemia Chronic_myeloid_leukemia Epithelial_cells Fibroblast Monocyte Stromal_cells Breast_adenocarcinoma Colorectal_carcinoma Erythroid Fibrosarcoma Neuroblastoma Unknown Burkitt_lymphoma Embryonic_stem_cells Erythroid_progenitor_cells Histiocytic_cells Rhabdoid_tumor; do
    # mkdir -p ${per_celltype}/${per_chr}
    mv ${per_chr}/${per_celltype}/*.txt ${per_celltype}/${per_chr}/
  done
done
