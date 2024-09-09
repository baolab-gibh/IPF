#!/usr/bin/env bash
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: May 08, 2024

#
## Functional enrichment of IPF GWAS loci
#

project_dir=~/Documents/projects/wp_ipf
source ${project_dir}/scripts/.env/bin/activate

singularity_image=~/Documents/git/singularity/base_tools.sif

func_annall=${project_dir}/outputs/analysis/functional_annotation/encodeCcreCombined.bed

# 1. Convert functional annotations from BigBed to BED
bigBedToBed ${project_dir}/inputs/functional_annotation/encodeCcreCombined.bb ${func_annall}

for x in CTCF K4m3 enhD enhP prom; do
  awk -v ANNOTYPE=${x} '$13 == ANNOTYPE {print}' ${func_annall} >| ${func_annall/.bed/.${x}.bed}
done

# sum([388814, 25537, 667599, 141830, 34803])


# 2. Clump GWAS summary statistics
plink --bfile ${project_dir}/inputs/public_gwas/gwas_1KG_phase3_plink \
  --clump ${project_dir}/outputs/analysis/public_gwas/GCST90399721.5e-8.bed \
  --out


# 2. Reformant the GWAS summary statistics
gwas_id=GCST90399695
gwas_id=GCST90399721

gwas_sumstat=${project_dir}/outputs/analysis/public_gwas/${gwas_id}.5e-8.bed
awk -F $'\t' -v GWAS_ID=${gwas_id} -f- <<'EOF' \
  <(zcat ${project_dir}/inputs/public_gwas/${gwas_id}.tsv.gz) \
  | sort -k1,1g -k2,2g \
  >| ${gwas_sumstat}
BEGIN {
  gwas_dict["GCST90399721"] = "IPF.Zhou_etal.36777996"
  gwas_dict["GCST90399695"] = "COPD.Zhou_etal.36777996"
}
NR > 1 && $8 < 5e-8 {
  OFS = "\t"; if ($10 == "") { RSID = $1":"$2":"$3":"$4 } else { RSID = $10 }
  print $1,$2 - 1,$2,RSID,gwas_dict[GWAS_ID],"*"
}
EOF

ipf_func_enrich=${project_dir}/outputs/analysis/functional_enrichment/ipf_func_enrich.${gwas_id}.all.txt
# QTLtools fdensity --help
QTLtools fdensity \
  --qtl "${gwas_sumstat:?Missing GWAS summary statistics!}" \
  --bed "${func_annall:?Missing functional annotation!}" \
  --out "${ipf_func_enrich:-ipf_func_enrich.txt}"


for per_anno in CTCF K4m3 enhD enhP prom; do
  ipf_func_enrich=${project_dir}/outputs/analysis/functional_enrichment/ipf_func_enrich.${gwas_id}.${per_anno}.txt
  func_annper=${func_annall/.bed/.${per_anno}.bed}
  # QTLtools fdensity --help
  QTLtools fdensity \
    --qtl "${gwas_sumstat:?Missing GWAS summary statistics!}" \
    --bed "${func_annper:?Missing functional annotation!}" \
    --out "${ipf_func_enrich:-ipf_func_enrich.txt}"
done

# Intersection between GWAS snp and functional annotation
func_per_snp=${project_dir}/outputs/analysis/functional_enrichment/snp_annotation.${gwas_id}.txt
bedtools intersect -a ${gwas_sumstat} -b ${func_annall} -wa -wb >| ${func_per_snp}

# Estimation of heritability
gwas_id=GCST90399721
gwas_trait=IPF
for gwas_id in GCST90399721 GCST90399695; do
  if [[ ${gwas_id} == GCST90399721 ]]; then
    gwas_trait=IPF
  elif [[ ${gwas_id} == GCST90399695 ]]; then
    gwas_trait=COPD
  fi

  awk -f- <<'EOFx' <(zcat ${project_dir}/inputs/public_gwas/${gwas_id}.tsv.gz)  \
    >| ${project_dir}/outputs/analysis/public_gwas/${gwas_trait}.${gwas_id}.tsv
  NR == 1 {print $10,$3,$4,$8,$5,$17,$11; next}
  NR >= 2 && length($3) == 1 && length($4) == 1 {
    OFS = "\t"
    if ($10 == "") { rs_id = $1":"$2":"$3":"$4 } else { rs_id = $10 }
    print rs_id,$3,$4,$8,$5,$17,$11
  }
EOFx


  # python -m pdb $(which munge_sumstats.py) \
  munge_sumstats.py \
    --sumstats ${project_dir}/outputs/analysis/public_gwas/${gwas_trait}.${gwas_id}.tsv \
    --snp rs_id \
    --a1 effect_allele \
    --a2 other_allele \
    --p p_value \
    --N-cas-col N_case \
    --N-con-col N_ctrl \
    --signed-sumstats beta,0 \
    --N 1392366 \
    --out ${project_dir}/outputs/analysis/public_gwas/${gwas_trait}.${gwas_id}
done


# estiamte SNP heritability and perform LD score regression
# COPD
ldsc.py \
  --h2 ${project_dir}/outputs/analysis/public_gwas/COPD.GCST90399695.sumstats.gz \
  --ref-ld-chr ${project_dir}/inputs/LDSCORE/baselineLD_v2.3/ \
  --w-ld-chr ${project_dir}/inputs/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/ \
  --out ${project_dir}/outputs/analysis/ldscore_regression/COPD.GCST90399695

# IPF
ldsc.py \
  --h2 ${project_dir}/outputs/analysis/public_gwas/IPF.GCST90399721.sumstats.gz \
  --ref-ld-chr ${project_dir}/inputs/LDSCORE/baselineLD_v2.3/ \
  --w-ld-chr ${project_dir}/inputs/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/ \
  --out ${project_dir}/outputs/analysis/ldscore_regression/IPF.GCST90399721

# COPD vs. IPF
ldsc.py \
  --rg ${project_dir}/outputs/analysis/public_gwas/IPF.GCST90399721.sumstats.gz,${project_dir}/outputs/analysis/public_gwas/COPD.GCST90399695.sumstats.gz \
  --ref-ld-chr ${project_dir}/inputs/LDSCORE/baselineLD_v2.3/ \
  --w-ld-chr ${project_dir}/inputs/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/ \
  --out ${project_dir}/outputs/analysis/ldscore_regression/IPF.COPD
