#!/bin/bash

project_dir=~/Documents/projects/wp_ipf
app_container=${project_dir}/scripts/singularity/tre_snp_toolkit.sif
genomic_feature_file=${project_dir}/outputs/references/genomic_feature/gencode.v46.basic.annotation.gtf.gz

# Obtain transcripts
transcript_regions=${project_dir}/outputs/references/genomic_feature/gencode.v46.basic.annotation.transcript_regions.gtf.gz
awk -F $'\t' -f- <<'EOF' <(zcat $genomic_feature_file) | gzip >| $transcript_regions
$3 ~ /transcript/ && $9 ~ /Ensembl_canonical/ {print $0}
EOF


# Obtain exons and UTRs
exon_utr_regions=${project_dir}/outputs/references/genomic_feature/gencode.v46.basic.annotation.exon_and_utr.gtf.gz
awk -F $'\t' -f- <<'EOF' <(zcat $genomic_feature_file) | gzip >| ${exon_utr_regions}
$3 ~ /exon|UTR/ && $9 ~ /Ensembl_canonical/ {print $0}
EOF

# singularity exec ${app_container} bedtools subtract -a ${transcript_regions} -b ${exon_utr_regions} -s | les

# Obtain introns
intron_regions=${project_dir}/outputs/references/genomic_feature/gencode.v46.basic.annotation.intron_regions.gtf
rm -f ${intron_regions}
while read -r line; do
  transcript_id=$(grep -oP 'ENST[0-9]{11}\.[0-9]{1,2}' <<<${line})
  zgrep -w ${transcript_id} ${exon_utr_regions} > cur_exon_utr_regions.gtf
  zgrep -w ${transcript_id} ${transcript_regions} > cur_transcript_regions.gtf
  singularity exec ${app_container} bedtools subtract -a cur_transcript_regions.gtf -b cur_exon_utr_regions.gtf -s \
    | awk -F$'\t' '{OFS="\t"; sub(/transcript/, "intron", $3); intron_id="intron_id \""$1"_"$4"_"$5"\";"; print $0" "intron_id}'
done < <(zcat ${transcript_regions}) >> ${intron_regions}
gzip ${intron_regions}
