#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: May 09, 2024
# Updated: May 09, 2024

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(topr)
  library(ggh4x)
  library(ggforce)
  library(rtracklayer)

  # library(coloc)
  # library(TwoSampleMR)
})

rm(list = ls()); gc()

project_dir <- "~/Documents/projects/wp_ipf"

gwas_id <- c("GCST90399721", "GCST90399695")
gwas_traits <- c("IPF", "COPD")


# Manhattan plot
gwas_list <- file.path(project_dir, "inputs/public_gwas", c("GCST90399721.tsv.gz", "GCST90399695.tsv.gz")) %>% purrr::set_names(gwas_traits)
gwas_sumstat_list <- lapply(names(gwas_list), function(x) {
  fread(gwas_list[x]) %>%
    dplyr::select(CHROM = "chromosome", POS = "base_pair_location", ID = "rs_id", effect_allele, other_allele, p_value, beta) %>%
    dplyr::filter(p_value < 0.05, nchar(effect_allele) == 1, nchar(other_allele) == 1, stringr::str_detect(ID, "^rs")) %>%
    dplyr::mutate(ID = dplyr::if_else(ID == "", paste0(CHROM, ":", POS, ":", effect_allele, ":", other_allele), ID))
}) %>% purrr::set_names(names(gwas_list))

p <- topr::manhattan(gwas_sumstat_list, annotate = 5e-8, angle = 90, legend_labels = gwas_traits, color = c("#E64B35B2", "#3C5488B2"), ntop = 1)
plot_saveto <- file.path(project_dir, "outputs/analysis/public_gwas/plots/ipf_func_enrich_manhatten.pdf")
ggsave(plot_saveto, p, width = 16, height = 8)


# Locuszoom plot
plot_saveto <- file.path(project_dir, "outputs/analysis/public_gwas/plots/ipf_func_enrich_locuszoom.MUC5B.pdf")
pdf(plot_saveto, width = 8, height = 5)
topr::regionplot(gwas_sumstat_list, annotate = 5e-8, color = c("#E64B35B2", "#3C5488B2"), legend_labels = gwas_traits, gene = "MUC5B")
dev.off()


# Tissue specific eQTL results
chain_file <- "/home/zzhang/Documents/projects/resources/UCSCGenombrowser/hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain"
chain <- import.chain(chain_file)
eqtl_src <- c("Lung", "Cultured_fibrolbasts")
eqtl_list <- file.path(project_dir, "inputs/public_eqtl", c("ENSG00000138640.Lung.allpairs.txt.gz", "ENSG00000138640.Cells_Transformed_fibroblasts.allpairs.txt.gz")) %>% purrr::set_names(eqtl_src)
eqtl_sumstat_list <- lapply(names(eqtl_list), function(x) {
  ss <- fread(eqtl_list[x]) %>%
    dplyr::select(ID = "variant_id", p_value = "pval_nominal", beta = "slope") %>%
    tidyr::separate(ID, into = c("CHROM", "POS", "other_allele", "effect_allele"), sep = "_", remove = FALSE, extra = "drop") %>%
    dplyr::mutate(POS = as.integer(POS)) %>%
    dplyr::filter(nchar(effect_allele) == 1, nchar(other_allele) == 1)

  ss_cur <- GenomicRanges::makeGRangesFromDataFrame(ss, keep.extra.columns = TRUE, seqnames.field = "CHROM", start.field = "POS", end.field = "POS")
  seqlevelsStyle(ss_cur) <- "UCSC"
  liftOver(ss_cur, chain) %>%
    as.data.frame() %>%
    dplyr::select(CHROM = "seqnames", POS = "start", ID, effect_allele, other_allele, p_value, beta) %>%
    dplyr::mutate(CHROM = stringr::str_remove_all(CHROM, "^chr") %>% as.integer())
}) %>% purrr::set_names(names(eqtl_list))


plot_saveto <- file.path(project_dir, "outputs/analysis/public_gwas/plots/ipf_func_enrich_locuszoom.lung_MUC5B.IPF.pdf")
pdf(plot_saveto, width = 8, height = 5)
topr::regionplot(list(gwas_sumstat_list$IPF, eqtl_sumstat_list$Lung), annotate = 5e-2, color = c("#E64B35B2", "darkblue"), legend_labels = c("IPF GWAS","Lung eQTL"), gene = "MUC5B")
dev.off()

plot_saveto <- file.path(project_dir, "outputs/analysis/public_gwas/plots/ipf_func_enrich_locuszoom.Cultured_fibrolbasts_MUC5B.IPF.pdf")
pdf(plot_saveto, width = 8, height = 5)
topr::regionplot(list(gwas_sumstat_list$IPF, eqtl_sumstat_list$Cultured_fibrolbasts), annotate = 5e-2, color = c("#E64B35B2", "darkgreen"), legend_labels = c("IPF GWAS","Fibroblast"), gene = "MUC5B")
dev.off()


# Overview of annotations around the GWAS SNPs
ann_classes <- c("prom", "enhP", "enhD", "K4m3", "CTCF", "all")
for (per_gwas_id in gwas_id) {
  density_tab <- lapply(ann_classes, function(x) {
    file.path(project_dir, "outputs/analysis/functional_enrichment", paste0("ipf_func_enrich.", per_gwas_id, ".", x, ".txt")) %>%
      fread() %>%
      dplyr::select(Start = V1, End = V2, Count = V3) %>%
      dplyr::mutate(ann_source = x)
  }) %>%
    Reduce(rbind, .)

  p <- ggplot(density_tab) +
    geom_point(aes(x = Start, y = Count, color = ann_source), size = 1) +
    geom_line(aes(x = Start, y = Count, color = ann_source)) +
    facet_zoom(xlim = c(-5e4, 5e4), zoom.size = .5) +
    labs(x = "Position (GRCh38)", y = "Annotation count per bin") +
    theme_classic()
  plot_saveto <- file.path(project_dir, "outputs/analysis/functional_enrichment/plots", paste0("ipf_func_enrich_density.", per_gwas_id, ".zoom_neg5e4_pos5e4.pdf"))
  ggsave(plot_saveto, p, width = 10, height = 5)

  p <- ggplot(density_tab) +
    geom_point(aes(x = Start, y = Count, color = ann_source), size = 1) +
    geom_line(aes(x = Start, y = Count, color = ann_source)) +
    facet_zoom(xlim = c(-6e5, -3e5), y = ann_source %in% c("enhP", "prom"), zoom.size = .5) +
    labs(x = "Position (GRCh38)", y = "Annotation count per bin", color = "Annotation classes") +
    theme_classic()
  plot_saveto <- file.path(project_dir, "outputs/analysis/functional_enrichment/plots", paste0("ipf_func_enrich_density.", per_gwas_id, ".zoom_neg5e5_neg4e5.pdf"))
  ggsave(plot_saveto, p, width = 15, height = 5)
}


# Functional annotations of each GWAS SNPs, pie chart to show percentage of each type of annotation
for (per_gwas_id in gwas_id) {
  anno_tab <- file.path(project_dir, "outputs/analysis/functional_enrichment", paste0("snp_annotation.", per_gwas_id, ".txt")) %>%
    fread() %>%
    dplyr::select(Var_chrom = V1, Var_pos = V3, Var_rsid = V4, Annotation_chrom = V7, Annotation_start = V8, Annotation_end = V9, Annotation_name = V10, Annotation_class = "V19")

  anno_summary <- dplyr::group_by(anno_tab, Annotation_class) %>%
    dplyr::summarise(Count = n()) %>%
    dplyr::mutate(Annotation_class = forcats::fct_reorder(Annotation_class, Count, .desc = TRUE)) %>%
    dplyr::arrange(desc(Annotation_class)) %>%
    dplyr::mutate(prop = Count / sum(Count) *100) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.25 * prop ) %>%
    dplyr::mutate(label = paste0(Count, " (", round(prop, 1), "%)"))

  p <- ggplot(anno_summary) +
    geom_bar(aes(x = "", y = prop, fill = Annotation_class), stat = "identity", position = "stack", width = 1) +
    geom_text(aes(x = 1.25, y = ypos, label = label)) +
    labs(fill = "Annotation") +
    theme_void() +
    coord_polar("y", start = 45)

  plot_saveto <- file.path(project_dir, "outputs/analysis/functional_enrichment/plots", paste0("snp_annotation.", per_gwas_id, ".annotation_classes.pie_chart.pdf"))
  ggsave(plot_saveto, p, width = 5, height = 5)
}
