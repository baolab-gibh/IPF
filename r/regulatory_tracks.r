#!/usr/bin/env Rscript
#

project_dir <- "~/Documents/projects/wp_ipf"

score_path <- file.path(project_dir, "outputs/analysis/tre_identification/tre_results") %>%
  list.files("*.dREG.peak.score.bed.gz$", recursive = TRUE, full.name = TRUE) %>%
  purrr::keep(~!stringr::str_detect(.x, pattern = "per_chunk"))

proba_path <- file.path(project_dir, "outputs/analysis/tre_identification/tre_results") %>%
  list.files("*.dREG.peak.prob.bed.gz$", recursive = TRUE, full.name = TRUE) %>%
  purrr::keep(~!stringr::str_detect(.x, pattern = "per_chunk"))

file_path_tab <- data.frame(score_path = score_path, prob_path = proba_path) %>%
  dplyr::mutate(sample_id = stringr::str_extract(score_path, "/(SRR[0-9]+)/", group = 1))

peak_tab_list <- apply(file_path_tab, 1, function(vec) {
  score_tab <- vec[1] %>% fread() %>% dplyr::select(Chrom = "V1", Start = "V2", End = "V3", Score = "V4")
  proba_tab <- vec[2] %>% fread() %>% dplyr::select(Chrom = "V1", Start = "V2", End = "V3", Proba = "V4")
  sample_id <- vec[3]

  dplyr::inner_join(score_tab, proba_tab, by = c("Chrom", "Start", "End")) %>%
    dplyr::mutate(sample_id = sample_id) %>%
    as.data.frame() %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = "Chrom", start.field = "Start", end.field = "End")
}) %>%
  purrr::set_names(file_path_tab$sample_id) %>%
  as("SimpleGRangesList")

selected_samples <- names(peak_tab_list)
meta_info_tab <- file.path(project_dir, "outputs/overview/sample_information/sample_info.human.txt") %>% fread() %>%
  dplyr::select(SRR, cell_line, cell_type, tissue, treatment) %>%
  dplyr::filter(cell_type %in% c("Fibroblast", "fibroblast") | stringr::str_detect(tissue, "[Ll]ung")) %>%
  tidyr::separate_rows(SRR, sep = ",") %>%
  dplyr::filter(SRR %in% selected_samples)

fam13_region <- GRanges("chr4:88725955-89111398") %>% flank(250000, both=TRUE)
fam13_peak_tab <- peak_tab_list %>%
  lapply(function(x) {
    subsetByOverlaps(x, fam13_region) %>%
      as.data.frame() %>%
      dplyr::select(Chrom = "seqnames", Start = "start", End = "end", Score, Proba, Sample_id = "sample_id")
  }) %>%
  Reduce(rbind, .)
  # dplyr::arrange(Chrom, Start) %>%
  # dplyr::mutate(bins = cut(Start, 5000)) %>%
  # dplyr::group_by(Sample_id, bins) %>%
  # dplyr::summarize(Score_mean = mean(Score)) %>%
  # dplyr::mutate(bins = as.character(bins)) %>%
  # dplyr::mutate(Start = stringr::str_extract(bins, "(8.*e\\+07),", group = 1) %>% as.integer()) %>%
  # dplyr::mutate(End = stringr::str_extract(bins, ",(8.*e\\+07)", group = 1) %>% as.integer()) %>%
  # dplyr::select(Sample_id, Start, End, Score = "Score_mean")

# 88625955-89211398
p <- ggplot(fam13_peak_tab) +
  geom_col(aes(x = Start, y = Score)) +
  labs(x = NULL, y = "TRE score") +
  lims(x = c(88625955, 89211398)) +
  theme_classic() +
  theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank())
plot_saveto <- file.path(project_dir, "outputs/analysis/public_gwas/plots/peaks_around_fam13a.pdf")
ggsave(plot_saveto, plot = p, width = 6, height = 2)
