#!/usr/bin/env Rscript
# File: collect_tres.r
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Aug 06, 2024
# Updated: Aug 06, 2024

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)

  library(patchwork)
  library(ggsci)
  library(ggalluvial)

  library(biomaRt)
  library(rtracklayer)
  library(GenomicRanges)
})


collect_garbage <- function(...) {
  rm(list = ...); gc()
}


project_dir <- "~/Documents/projects/wp_ipf"

#
## Utilities
#
# Load genomic features
#         biomart                version
# 1         genes      Ensembl Genes 112
# 2 mouse_strains      Mouse strains 112
# 3          snps  Ensembl Variation 112
# 4    regulation Ensembl Regulation 112
# listEnsembl(version = "112")
# listAttributes(mart = ensembl)

# Features fetched from Ensembl database.
save_genomic_features <- file.path(project_dir, "outputs/analysis/tre_identification/selected_genomic_features.txt")
if (file.exists(save_genomic_features)) {
  feature_mat <- fread(save_genomic_features)
} else {
  ensembl <- useEnsembl(biomart = "ensembl", version = "112", dataset = "hsapiens_gene_ensembl")
  selected_attributes <- c(
    "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position",
    "strand", "gene_biotype"
  )
  feature_mat <- getBM(attributes = selected_attributes, mart = ensembl) %>%
    dplyr::filter(chromosome_name %in% as.character(c(1:22, "X", "Y"))) %>%
    dplyr::mutate(strand = dplyr::case_when(strand %in% c(1, "1", "+") ~ "+", strand %in% c(-1, "-1", "-") ~ "-", TRUE ~ "*"))
  fwrite(feature_mat, save_genomic_features)
}

# Biotypes of genes annotated in the database
coding_gene <- c("protein_coding")
long_noncoding_rna <- c("lncRNA")
ig_gene <- c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene", "IG_V_gene", "IG_V_pseudogene", "IG_pseudogene")
tr_gene <- c("TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene")
non_coding_rna <- c("miRNA", "misc_RNA", "rRNA", "rRNA_pseudogene", "scRNA", "snRNA", "snoRNA", "ribozyme", "sRNA", "scaRNA", "vault_RNA")
pseudo_genes <- c("processed_pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene")
others <- c("TEC", "artifact")

# GenomicRanges of features
feature_gr <- makeGRangesFromDataFrame(
  feature_mat, keep.extra.columns = TRUE, seqnames.field = "chromosome_name", start.field = "start_position",
  end.field = "end_position", strand.field = "strand"
) 
seqlevelsStyle(feature_gr) <- "UCSC"
prom_gr <- feature_gr %>% promoters() %>% unique()

sample_info_tab <- file.path(project_dir, 'outputs/overview/sample_information/sample_info.human.txt') %>%
  fread() %>% dplyr::distinct()


#
## TREs from dREG. Unharmonized TREs
#
# Number of TREs
for (which_peak in c("raw", "full")) {
  if (which_peak == "full") {
    col_names <- c("seqnames", "starts", "ends", "dreg_scores", "dreg_p_values", "x")
    file_pattern <- "*.peak.full.bed.gz$"
    fig_save_to <- file.path(project_dir, "outputs/analysis/tre_identification/plots/quantile.significant_peaks.pdf")
  } else if (which_peak == "raw") {
    col_names <- c("seqnames", "starts", "ends", "dreg_scores", "dreg_p_values", "x", "y", "z")
    file_pattern <- "*.peak.raw.bed.gz$"
    fig_save_to <- file.path(project_dir, "outputs/analysis/tre_identification/plots/quantile.all_peaks.pdf")
  }

  tre_tab_list <- file.path(project_dir, "outputs/analysis/tre_identification/tre_results") %>%
    list.files(pattern = file_pattern, recursive = TRUE, full.name = TRUE)
  quantile_tab <- tre_tab_list %>% lapply(
    function(f) fread(f, header = FALSE, col.names = col_names) %>%
      dplyr::mutate(sample_id = stringr::str_extract(f, "/(SRR[0-9]+)/", group = 1)) %>%
      dplyr::reframe(
        sample_id = unique(sample_id), n = dplyr::n(), mean_score = mean(dreg_scores),
        quantile = quantile(dreg_scores, prob = c(0, 0.25, 0.5, 0.75, 1)),
        quantile_name = c("min_score", "q25", "q50", "q75", "max_score")
      ) %>%
      tidyr::pivot_wider(names_from = quantile_name, values_from = quantile)
    ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(sample_id = forcats::fct_reorder(sample_id, q50)) %>%
    tidyr::pivot_longer(-c(sample_id, n), names_to = "quantile", values_to = "value") %>%
    dplyr::mutate(quantile = factor(quantile, levels = c("min_score", "q25", "mean_score", "q50", "q75", "max_score")))

  quantile_plot <- ggplot(data = quantile_tab) +
    geom_line(aes(x = sample_id, y = value, color = quantile, group = quantile)) +
    scale_color_npg() +
    labs(x = "sample", y = "TRE score", color = "quantile") +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(fig_save_to, plot = quantile_plot, width = 8, height = 5)
  # rm(list = c("quantile_tab", "quantile_plot")); gc()
}

n_tre_persample <- quantile_tab %>% dplyr::select(SRR = "sample_id", n_tres = "n") %>%
  dplyr::distinct() %>%
  dplyr::left_join(sample_info_tab, by = "SRR") %>%
  dplyr::mutate(SRR = forcats::fct_reorder(SRR, n_tres)) %>%
  dplyr::mutate(cell_line = dplyr::if_else(cell_line == "", "Others", cell_line)) %>%
  dplyr::mutate(cell_line = forcats::fct_reorder(cell_line, n_tres)) %>%
  dplyr::mutate(tissue = dplyr::if_else(tissue == "", "Others", tissue)) %>%
  dplyr::mutate(tissue = forcats::fct_reorder(tissue, n_tres)) %>%
  dplyr::mutate(cell_type = dplyr::if_else(cell_type == "", "Others", cell_type)) %>%
  dplyr::mutate(cell_type = forcats::fct_reorder(cell_type, n_tres))

p <- ggplot(n_tre_persample, aes(axis1 = tissue, axis2 = SRR, axis3 = cell_type, y = n_tres)) +
  geom_alluvium(aes(fill = tissue)) +
  geom_stratum(width = 0.5) +
  geom_text(aes(label = after_stat(stratum)), stat = "stratum") +
  scale_x_discrete(limits = c("tissue", "SRR", "cell_type")) +
  labs(x = NULL, y = "Number of TREs") +
  theme_minimal()
psaveto <- file.path(project_dir, 'outputs/analysis/tre_identification/plots/sample_by_tissue_and_cell_type.pdf')
ggsave(psaveto, p, width = 10, height = 16)


# Number of TREs per sample grouped by cell type, cell line, or tissue
#p_treatment <- ggplot(n_tre_persample) + geom_boxplot(aes(x = treatment, y = n_tres)) + geom_violin(aes(x = treatment, y = n_tres)) + theme_classic()
p_persample <- ggplot(n_tre_persample) + geom_bar(aes(x = SRR, y = n_tres, color = tissue), stat = "identity") + theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
p_cellline <- ggplot(n_tre_persample) + geom_boxplot(aes(x = cell_line, y = n_tres)) + labs(x = NULL, y = "Nr. TREs by cell line") + theme_classic()
p_celltype <- ggplot(n_tre_persample) + geom_boxplot(aes(x = cell_type, y = n_tres)) + labs(x = NULL, y = "Nr. TREs by cell type") + theme_classic()
p_tissue <- ggplot(n_tre_persample) + geom_boxplot(aes(x = tissue, y = n_tres)) + labs(x = NULL, y = "Nr. TREs by tissue") + theme_classic()

p <- p_persample + ((p_tissue / p_celltype / p_cellline) & theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
psaveto <- file.path(project_dir, "outputs/analysis/tre_identification/plots/number_of_tres_per_sample.by_xxx.pdf")
ggsave(psaveto, p, width = 14, height = 8)



# Enrichment analysis.
# Percentage of promoter TRE, and the biotype percentage of genes whose promoter is overlaped by TRE
col_names <- c("seqnames", "starts", "ends", "dreg_scores", "dreg_p_values", "centroid")
file_pattern <- "*.peak.full.bed.gz$"
tre_tab_list <- file.path(project_dir, "outputs/analysis/tre_identification/tre_results") %>% list.files(pattern = file_pattern, recursive = TRUE, full.name = TRUE)
tre_prom_list <- tre_tab_list %>% lapply(
  function(f) {
     #f <- tre_tab_list[[1]]
    sample_id = stringr::str_extract(f, "/(SRR[0-9]+)/", group = 1)
    sample_gr <- fread(f, header = FALSE, col.names = col_names) %>%
      makeGRangesFromDataFrame(TRUE, FALSE, seqnames.field = "seqnames", start.field = "starts", end.field = "ends")

    tre_biotype <- subsetByOverlaps(prom_gr, sample_gr) %>%
      unique() %>%
      (function(x) x$gene_biotype) %>%
      table() %>%
      dplyr::tibble(gene_biotype = names(.), counts = .) %>%
      dplyr::mutate(gene_biotype = dplyr::case_when(
        gene_biotype %in% ig_gene ~ "IG_gene",
        gene_biotype %in% tr_gene ~ "TR_gene",
        gene_biotype %in% coding_gene ~ "Protein_coding",
        gene_biotype %in% long_noncoding_rna ~ "LncRNA",
        gene_biotype %in% non_coding_rna ~ "Non_coding_RNA",
        gene_biotype %in% pseudo_genes ~ "Pseudo_gene",
        TRUE ~ "Other"
      )) %>%
      dplyr::group_by(gene_biotype) %>%
      dplyr::summarize(counts = sum(counts)) %>%
      tidyr::pivot_wider(names_from = gene_biotype, values_from = counts) %>%
      dplyr::rename_with(~paste0("n_", .x))

    n_tre_prom_gr <- subsetByOverlaps(sample_gr, prom_gr) %>% unique() %>% length()
    dplyr::tibble(sample_id = sample_id, n_tre = length(sample_gr), n_prom = n_tre_prom_gr) %>% dplyr::bind_cols(tre_biotype)
  })

tre_prom_tab <- dplyr::bind_rows(tre_prom_list) %>%
  dplyr::mutate(dplyr::across(-c(sample_id, n_tre), ~ .x / n_tre * 100, .names = "perc_{.col}")) %>%
  dplyr::select(sample_id, n_tre, dplyr::contains("perc_")) %>%
  dplyr::mutate(sample_id = forcats::fct_reorder(sample_id, perc_n_prom)) %>%
  tidyr::pivot_longer(-c(sample_id, n_tre), names_to = "source", values_to = "perc") %>%
  dplyr::mutate(source = stringr::str_remove(source, "perc_n_")) %>%
  dplyr::mutate(perc = dplyr::if_else(is.na(perc), 0, perc)) %>%
  dplyr::mutate(source = factor(source, levels = rev(c("Protein_coding", "LncRNA", "Non_coding_RNA", "Pseudo_genes", "IG_gene", "TR_gene", "Others", "prom"))))

tre_biotype_plot <- ggplot() +
  geom_bar(aes(x = sample_id, y = perc), fill = "grey", stat = "identity", data = dplyr::filter(tre_prom_tab, source == "prom")) +
  scale_y_reverse() +
  scale_fill_npg() +
  labs(x = paste0("Sample (n=", length(tre_tab_list), ")"), y = "Percentage (promoter)")

tre_prom_plot <- ggplot() +
  geom_bar(aes(x = sample_id, y = perc, fill = source), stat = "identity", position = "stack", data = dplyr::filter(tre_prom_tab, source != "prom")) +
  labs(x = NULL, y = "Percentage (biotype)")
tre_fun_plot <- (tre_prom_plot / tre_biotype_plot) & theme_classic() & theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig_save_to <- file.path(project_dir, "outputs/analysis/tre_identification/plots/enrichment_analysis.significant_peaks.pdf")
ggsave(fig_save_to, plot = tre_fun_plot, width = 15, height = 5)


#
## Estimate TRE functions. Harmonized TRE regions.
#
col_names <- c("seqnames", "starts", "ends", "chromatin_states")
tre_tab_list <- file.path(project_dir, "outputs/analysis/tre_identification/harmonization/chromatin_states") %>%
  list.files(pattern = "genome_2_segments.bed", recursive = TRUE, full.name = TRUE) %>%
  purrr::discard(~stringr::str_detect(.x, "all_celltype"))

tre_tab <- tre_tab_list %>%
  lapply(function(fpath) {
    fp_vec <- stringr::str_split(fpath, "/") %>% unlist()
    cell_type <- fp_vec[length(fp_vec) - 2]
    fread(fpath, header = FALSE, col.names = col_names) %>%
      dplyr::filter(chromatin_states == "E2") %>%
      dplyr::select(-chromatin_states) %>%
      dplyr::mutate(cell_type = cell_type)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = table(cell_type) %>% sort(decreasing = T) %>% names())) %>%
  dplyr::mutate(seqnames = factor(seqnames, levels = paste0("chr", c(1:22, "X", "Y"))))


# Basic statistics. number TREs per cell type
n_tre_per_celltype <- tre_tab %>% dplyr::group_by(cell_type) %>% dplyr::summarize(n = n()) %>% dplyr::mutate(cell_type = forcats::fct_reorder(cell_type, n))
p <- ggplot(data = n_tre_per_celltype) +
  geom_bar(aes(x = cell_type, y = n, fill = cell_type), stat = "identity") +
  geom_text(aes(x = cell_type, y = -5, label = n), hjust = 1) +
  ylim(-5, max(n_tre_per_celltype$n) * 1.1) +
  theme_void() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(-.1, -.1, -.1, -.1), "inches")) +
  coord_polar(theta = "y", direction = 1)
plot_saveto <- file.path(project_dir, "outputs/analysis/tre_identification/plots/number_of_tres_per_cell_type.pdf")
ggsave(plot_saveto, plot = p, width = 7.5, height = 5.5)


# Overlaps with genomic features
tre_grl <- tre_tab %>% as.data.frame() %>% makeGRangesFromDataFrame(TRUE, start.field = "starts", end.field = "ends") %>% GenomicRanges::split(f = .$cell_type)

harm_tre_tab <- lapply(names(tre_grl), function(x) {
  x_gr <- tre_grl[[x]]
  sub_gr <- subsetByOverlaps(prom_gr, x_gr)
  sub_gr %>% dplyr::as_tibble() %>%
      dplyr::mutate(gene_biotype = dplyr::case_when(
        gene_biotype %in% ig_gene ~ "IG_gene",
        gene_biotype %in% tr_gene ~ "TR_gene",
        gene_biotype %in% coding_gene ~ "Protein_coding",
        gene_biotype %in% long_noncoding_rna ~ "LncRNA",
        gene_biotype %in% non_coding_rna ~ "Non_coding_RNA",
        gene_biotype %in% pseudo_genes ~ "Pseudo_gene",
        TRUE ~ "Other"
      )) %>%
    dplyr::group_by(gene_biotype) %>%
    dplyr::summarize(n = dplyr::n(), cell_type = x, n_tre_total = length(x_gr), n_prom_total = length(prom_gr), n_overlap_proms = length(sub_gr)) %>%
    tidyr::pivot_wider(names_from = gene_biotype, values_from = n)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.x), 0, .x))) %>%
  dplyr::mutate(dplyr::across(c(n_overlap_proms, LncRNA, Non_coding_RNA, Other, Protein_coding, Pseudo_gene, TR_gene), ~.x / n_prom_total * 100, .names = "perc_{.col}")) %>%
  dplyr::mutate(cell_type = forcats::fct_reorder(cell_type, n_overlap_proms))

harm_tre_ov_prom <- ggplot(harm_tre_tab) +
  geom_bar(aes(y = cell_type, x = n_overlap_proms), fill = "grey95", color = "black", stat = "identity") +
  geom_text(aes(y = cell_type, x = 0, label = n_overlap_proms), hjust = 1, color = "black") +
  labs(y = NULL, x = "Number of TREs overlapping promoters") +
  scale_x_reverse() +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.text.y = element_blank(), legend.position = "none")

harm_tre_ov_biotype_tab <- dplyr::select(harm_tre_tab, cell_type, dplyr::starts_with("perc_"), -perc_n_overlap_proms) %>%
  tidyr::pivot_longer(-cell_type) %>%
  dplyr::mutate(name = stringr::str_remove(name, "perc_") %>% forcats::fct_reorder(value))
harm_tre_ov_biotype <- ggplot(harm_tre_ov_biotype_tab) +
  geom_bar(aes(y = cell_type, x = value, fill = name), stat = "identity", position = "stack") +
  scale_fill_npg() +
  theme_classic() +
  labs(y = NULL, x = "Percentage", fill = "Gene biotype") +
  theme(axis.text.y = element_text(hjust = 0.5), legend.position = "top")

harm_tre_p <- (harm_tre_ov_prom + harm_tre_ov_biotype) + plot_layout(guides = "collect") & theme(legend.position = "top")
plot_saveto <- file.path(project_dir, "outputs/analysis/tre_identification/plots/harmonized_tre.overlap_with_genomic_features.pdf")
ggsave(plot_saveto, plot = harm_tre_p, width = 7.5, height = 5.5)


# TREs overlapped with GWAS loci
loci_in_tre <- file.path(project_dir, "outputs/analysis/functional_annotation/TREs") %>%
  list.files(pattern = "*.overlaps.txt", full.names = TRUE) %>%
  lapply(function(x) {
    cell_type <- basename(x) %>% stringr::str_split("\\.") %>% unlist() %>% .[1]
    fread(x, header = FALSE) %>%
      dplyr::mutate(cell_type = cell_type) %>%
      dplyr::mutate(GWAS = stringr::str_extract(cell_type, "COPD|IPF")) %>%
      dplyr::mutate(cell_type = stringr::str_remove(cell_type, "_(COPD|IPF)")) %>%
      dplyr::select(cell_type, GWAS, CHROM = V1, POS = V2, ID = V4) %>%
      dplyr::distinct()
  }) %>%
  dplyr::bind_rows()

loci_in_tre_per_celltype <- loci_in_tre %>% dplyr::group_by(GWAS, cell_type) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::mutate(cell_type = forcats::fct_reorder(cell_type, n))

p_copd <- ggplot(loci_in_tre_per_celltype %>% dplyr::filter(GWAS == "COPD")) +
  geom_bar(aes(y = cell_type, x = n), fill = "darkred", stat = "identity", position = "stack") +
  labs(x = "Nr. loci in TREs (COPD)", y = NULL) +
  scale_x_reverse() +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.text.y = element_blank())

p_ipf <- ggplot(loci_in_tre_per_celltype %>% dplyr::filter(GWAS == "IPF")) +
  geom_bar(aes(y = cell_type, x = n), fill = "darkblue", stat = "identity", position = "stack") +
  labs(x = "Nr. loci in TREs (IPF)", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_text(hjust = 0.5))

p <- (p_copd + p_ipf) + plot_layout(guides = "collect")
plot_saveto <- file.path(project_dir, "outputs/analysis/tre_identification/plots/gwas_loci_in_tre.by_cell_type.pdf")
ggsave(plot_saveto, plot = p, width = 7.5, height = 5.5)


# Genes example
genes_by_loci_in_tre <- loci_in_tre %>%
  makeGRangesFromDataFrame(TRUE, start.field = "POS", end.field = "POS") %>%
  GenomicRanges::split(f = loci_in_tre$cell_type) %>%
  lapply(function(x) subsetByOverlaps(feature_gr, x) %>% as.data.frame()) %>%
  dplyr::bind_rows() %>%
  dplyr::distinct() %>%
  dplyr::mutate(gene_biotype = dplyr::case_when(
    gene_biotype %in% ig_gene ~ "IG_gene",
    gene_biotype %in% tr_gene ~ "TR_gene",
    gene_biotype %in% coding_gene ~ "Protein_coding",
    gene_biotype %in% long_noncoding_rna ~ "LncRNA",
    gene_biotype %in% non_coding_rna ~ "Non_coding_RNA",
    gene_biotype %in% pseudo_genes ~ "Pseudo_gene",
    TRUE ~ "Other"
  )) %>%
  dplyr::mutate(gene_biotype = factor(gene_biotype, levels = rev(c("IG_gene", "TR_gene", "Protein_coding", "LncRNA", "Pseudo_gene", "Non_coding_RNA", "Other"))))

p <- ggplot(genes_by_loci_in_tre) +
  geom_bar(aes(x = gene_biotype)) +
  labs(x = NULL, y = "Nr. of genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
psaveto <- file.path(project_dir, "outputs/analysis/tre_identification/plots/gene_biotype_in_gwas_loci_by_tre.pdf")
ggsave(psaveto, p, width = 4.5, height = 4.5)
