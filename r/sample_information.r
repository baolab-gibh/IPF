#!/usr/bin/env Rscript
# File: sample_information.r
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: May 30, 2024
# Updated: May 30, 2024

# BASH, collect sample information
# for x in /mnt/UGreen/human/*; do
#   is_pe=$([[ -e $x/${x}_rRNAremoved.fastq.2.gz ]] && echo 1 || echo 0)
#   if [[ $is_pe -eq 1 ]]; then
#     [[ ! -e $x/${x}_rRNAremoved.fastq.gz ]] && echo "$x/${x}_rRNAremoved.fastq.gz not found"
#     echo $x,$x,PE,${x}_rRNAremoved.fastq.gz,${x}_rRNAremoved.fastq.2.gz,$(readlink -f $x/$x/${x}_sort.bam)
#   else
#     echo $x,$x,SE,NA,NA,$(readlink -f $x/$x/${x}_sort.bam)
#   fi
# done > ~/Documents/projects/wp_ipf/misc/sample_information/sample_bam_info.txt

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggalluvial)
  library(ggsci)
})



rm(list=ls()); gc()
project_dir <- '~/Documents/projects/wp_ipf'


selected_sample_info_tab_saveto <- file.path(project_dir, 'outputs/overview/sample_information/sample_info.human.txt')
if (file.exists(selected_sample_info_tab_saveto)) {
  selected_sample_info_tab <- fread(selected_sample_info_tab_saveto)
} else {
  # GSE, GSM, SRR, information per SRR
  library_info_tab <- file.path(project_dir, 'misc/sample_information/all.xls') %>%
    readxl::read_excel(col_names = F, .name_repair = "unique_quiet") %>%
    as.data.table() %>%
    dplyr::select(GSE = `...1`, GSM = `...2`, SRR = `...3`, library = `...4`) %>%
    dplyr::mutate(library = tolower(library) %>% stringr::str_remove('5\'|library strategy: '))
  gsm_lib_info <- library_info_tab %>% dplyr::pull(library, GSM)
  srr_to_gsm <- library_info_tab %>% dplyr::pull(GSM, SRR)
  srr_per_gsm <- library_info_tab %>% dplyr::select(GSM, SRR) %>% dplyr::group_by(GSM) %>%
    dplyr::summarize(SRRs = paste(SRR, collapse = ',')) %>%
    dplyr::pull(SRRs, GSM)

  # GSE, GSM, information per GSM
  sample_info_tab <- file.path(project_dir, 'misc/sample_information/test_uniq.txt') %>% fread()

  # Available alignments, information per alignment
  bam_info_tab <- file.path(project_dir, 'misc/sample_information/sample_bam_info.txt') %>% fread() %>%
    dplyr::rename(GSM = V1, SRR = V2, sequencing = V3, read_r1 = V4, read_r2 = V5, bam_path = V6) %>%
    dplyr::mutate(GSM = dplyr::if_else(GSM %in% names(srr_to_gsm), srr_to_gsm[GSM], GSM)) %>%
    dplyr::mutate(SRR = dplyr::if_else(SRR %in% names(srr_per_gsm), srr_per_gsm[SRR], SRR)) %>%
    dplyr::select(-bam_path) %>%
    dplyr::mutate(library = gsm_lib_info[GSM]) %>%
    tidyr::separate_rows(SRR, sep = ",")  %>%
    dplyr::distinct() %>%
    dplyr::group_by(GSM, SRR) %>%
    dplyr::summarize(info = {
      has_pe = (dplyr::pick(read_r1, read_r2) %>% dplyr::filter(!is.na(read_r1), !is.na(read_r2)) %>% nrow()) > 0
      if (has_pe) {
        dplyr::pick(sequencing, read_r1, read_r2, library) %>% dplyr::filter(sequencing == "PE")
      } else {
        dplyr::pick(sequencing, read_r1, read_r2, library) %>% dplyr::filter(sequencing == "SE")
      }
    }) %>%
    tidyr::unnest(info) %>%
    dplyr::ungroup()

  # GSE, GSM, strand info per GSM
  strand_info_tab <- file.path(project_dir, 'misc/sample_information/info_strand.txt') %>% fread() %>%
    dplyr::select(GSE = V1, GSM = V2, strand = V6) %>%
    dplyr::mutate(strand = dplyr::case_when(strand == 1 ~ '+', strand == 2 ~ '-', TRUE ~ '*'))
  sample_strand <- strand_info_tab %>% dplyr::pull(strand, GSM)

  # Final table
  selected_sample_info_tab <- bam_info_tab %>%
    dplyr::left_join(sample_info_tab, by = c('GSM', "SRR")) %>%
    dplyr::select(GSE, GSM, SRR, species = specie, cell_line = cellLine, cell_type = cellType, tissue = tissue, treatment = treatment, sequencing, library, read_r1, read_r2) %>%
    # dplyr::mutate(tissue = dplyr::if_else(is.na(tissue), "Unknown", tissue)) %>%
    dplyr::mutate(strand_specificity = sample_strand[GSM]) %>%
    dplyr::group_by(GSM) %>%
    tidyr::fill(c(GSE, species, cell_line, cell_type, tissue, treatment))

  selected_sample_info_tab %>% fwrite(selected_sample_info_tab_saveto)
}


# Alluvial plot to show the sequencing technique
technique_tab <- selected_sample_info_tab %>% dplyr::group_by(sequencing, library, strand_specificity) %>% dplyr::summarize(freq = n())
technique_plot <- ggplot(data = technique_tab, aes(axis1 = library, axis2 = sequencing, axis3 = strand_specificity, y = freq)) +
  geom_alluvium(aes(fill = sequencing)) +
  geom_stratum() +
  geom_text(aes(label = after_stat(stratum)), stat = "stratum") +
  scale_x_discrete(limits = c("library", "sequencing", "strand_specificity"), expand = c(.01, .01)) +
  scale_fill_npg() +
  labs(x = NULL, y = "Number of samples") +
  theme_minimal()
technique_plot_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_technique.pdf')
ggsave(technique_plot_saveto, technique_plot, width = 6.5, height = 4)


# Pie chart to show the sequencing technique
sequencing_tab <- selected_sample_info_tab %>% dplyr::group_by(sequencing) %>% dplyr::summarize(freq = n()) %>% dplyr::mutate(sequencing = factor(sequencing, levels = c("SE", "PE")))
sequencing_plot <- ggplot(data = sequencing_tab) +
  geom_bar(aes(x = 1, y = freq, fill = sequencing), stat = "identity", position = "stack") +
  geom_label(aes(x = 1, y = freq * 0.75, label = freq), size = 6, stat = "identity") +
  scale_fill_npg() +
  theme_void() +
  coord_polar(theta = "y", direction = 1)
sequencing_plot_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_sequencing.pdf')
ggsave(sequencing_plot_saveto, sequencing_plot, width = 5, height = 5)


# sample_overview
selected_sample_info_tab$tissue %>% table
selected_sample_info_tab$cell_line %>% table
selected_sample_info_tab$cell_type %>% table

tissue_plot_tab <- selected_sample_info_tab %>%
  dplyr::mutate(tissue_l1 = dplyr::case_when(
    tissue %in% c("Skin", "Skin/Foreskin") ~ "Skin/Foreskin",
    tissue %in% c("Lung", "Lung/Bronchus", "Fetal lung") ~ "Lung",
    tissue %in% c("Heart", "Heart/Aorta") ~ "Heart/Aorta",
    tissue %in% c("Prostate/Lymph node", "Lymph node") ~ "Lymph node",
    tissue %in% c("Kidney", "Embryonic kidney") ~ "Kidney",
    tissue %in% c("Bone", "Bone marrow") ~ "Bone",
    TRUE ~ tissue
  )) %>%
  dplyr::mutate(
    place_holder = "XXX",
    tissue_l1 = factor(tissue_l1, table(tissue_l1) %>% sort(decreasing = T) %>% names()),
    cell_type = factor(cell_type, table(cell_type) %>% sort(decreasing = T) %>% names())
  )

g_pie <- ggplot(tissue_plot_tab) +
  geom_bar(aes(y = place_holder, fill = tissue_l1), color = 'grey', stat = 'count', linewidth = 0.1) +
  coord_polar("x", start = 0) +
  theme_void()
g_pie_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_pie.by_tissue.pdf')
ggsave(g_pie_saveto, g_pie, width = 10, height = 7)

g_pie <- ggplot(tissue_plot_tab) +
  geom_bar(aes(y = place_holder, fill = cell_type), color = 'grey', stat = 'count', linewidth = 0.1) +
  coord_polar("x", start = 0) +
  theme_void()
g_pie_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_pie.by_cell_type.pdf')
ggsave(g_pie_saveto, g_pie, width = 10, height = 7)


tissue_skplot_tab <- tissue_plot_tab %>% dplyr::group_by(tissue_l1, cell_type) %>% dplyr::summarize(freq = n()) %>%
  dplyr::mutate(tissue_l1 = dplyr::if_else(freq < 50, "Others", tissue_l1))
tissue_plot <- ggplot(data = tissue_skplot_tab, aes(axis1 = tissue_l1, axis2 = cell_type, y = freq)) +
  geom_alluvium(aes(fill = tissue_l1)) +
  geom_stratum() +
  geom_text(aes(label = after_stat(stratum)), stat = "stratum") +
  scale_x_discrete(limits = c("tissue_l1", "cell_type"), expand = c(.01, .01)) +
  labs(x = NULL, y = "Number of samples") +
  theme_minimal()
tissue_plot_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_tissue.pdf')
ggsave(tissue_plot_saveto, tissue_plot, width = 10, height = 16)


# Obtain PE samples
selected_pe_sample_info_tab_saveto <- file.path(project_dir, 'outputs/overview/sample_information/sample_info.human.paired_end.txt')
if (file.exists(selected_pe_sample_info_tab_saveto)) {
  selected_pe_sample_info_tab <- fread(selected_pe_sample_info_tab_saveto)
} else {
  selected_pe_sample_info_tab <- selected_sample_info_tab %>% dplyr::filter(sequencing == 'PE') %>%
    dplyr::select(group_id = SRR, sample_id = SRR, library, strand_specificity, read_r1, read_r2) %>%
    dplyr::mutate(read_r1 = file.path("/mnt/UGreen/human", stringr::str_extract(read_r1, "SRR[0-9]+"), read_r1)) %>%
    dplyr::mutate(read_r2 = file.path("/mnt/UGreen/human", stringr::str_extract(read_r2, "SRR[0-9]+"), read_r2))
  selected_pe_sample_info_tab %>% fwrite(selected_pe_sample_info_tab_saveto)
}


# Finished samples
finished_samples <- file.path(project_dir, "outputs/analysis/read_alignment/alignment") %>% list.files()
finished_tab <- selected_sample_info_tab %>% dplyr::filter(SRR %in% finished_samples, sequencing == "PE")

finished_by_celltype <- finished_tab %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(n = dplyr::n()) %>%
  dplyr::mutate(cell_type = dplyr::if_else(cell_type == "", "Others", cell_type)) %>%
  dplyr::mutate(cell_type = forcats::fct_reorder(cell_type, -n)) %>%
  dplyr::arrange(n)

p <- ggplot(data = finished_by_celltype) +
  geom_bar(aes(x = 1, y = n, fill = cell_type), stat = 'identity', linewidth = 0.55) +
  geom_label(aes(x = 1.5, y = n, label = n), alpha = 0.5, position = position_stack(vjust = 0.5)) +
  labs(fill = "Cell types") +
  coord_polar("y", start = 0) +
  theme_void()
p_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_finished.by_cell_type.pdf')
ggsave(p_saveto, p, width = 8, height = 6)


finished_by_cellline <- finished_tab %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarize(n = dplyr::n()) %>%
  dplyr::mutate(cell_line = dplyr::if_else(cell_line == "", "Others", cell_line)) %>%
  dplyr::mutate(cell_line = forcats::fct_reorder(cell_line, -n)) %>%
  dplyr::arrange(n)

p <- ggplot(data = finished_by_cellline) +
  geom_bar(aes(x = 1, y = n, fill = cell_line), stat = 'identity', linewidth = 0.5) +
  geom_label(aes(x = 1.5, y = n, label = n), alpha = 0.5, position = position_stack(vjust = 0.5)) +
  labs(fill = "Cell types") +
  coord_polar("y", start = 0) +
  theme_void()
p_saveto <- file.path(project_dir, 'outputs/overview/plots/sample_finished.by_cell_line.pdf')
ggsave(p_saveto, p, width = 8, height = 6)
