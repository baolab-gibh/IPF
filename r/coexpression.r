#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)

  library(DESeq2)
  library(WGCNA)
})

project_dir <- "~/Documents/projects/wp_ipf/"
sample_info <- file.path(project_dir, "outputs/overview/sample_information/sample_info.human.txt") %>% fread()

selected_chroms <- paste0("chr", seq(1:22))
selected_gene_type <- c(
  "protein_coding",
  "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene", "IG_D_gene",
  "lncRNA", "miRNA", "snRNA", "TEC", "snoRNA", "scaRNA", "scRNA", "sRNA"
)

dds_obj_path <- file.path(project_dir, "outputs/analysis/coexpression/deseq2.dds_object.rds")
norm_expr_tab_path <- file.path(project_dir, "outputs/analysis/coexpression/deseq2.normalized_expression_table.csv")
if (all(file.exists(dds_obj_path, norm_expr_tab_path))) {
  norm_mat <- fread(norm_expr_tab_path) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid") %>% as.matrix()
  dds <- readRDS(dds_obj)
} else {
  readcount_table <- file.path(project_dir, "outputs/analysis/quantify_expression/featureCounts") %>%
    list.files(pattern = "*.transcript_regions.tsv.gz", recursive = TRUE, full = TRUE) %>%
    lapply(function(x) {
      fread(x) %>% dplyr::rename_with(~stringr::str_extract(.x, "SRR[0-9]+"), starts_with("/")) %>%
        dplyr::filter(Chr %in% selected_chroms, gene_type %in% selected_gene_type)
    }) %>%
    Reduce(dplyr::left_join, .)

  readcount_mat <- readcount_table %>%
    dplyr::select(Geneid, starts_with("SRR")) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Geneid") %>%
    as.matrix()

  selected_sample_info <- sample_info %>%
    dplyr::select(GSE, GSM, SRR, species, cell_line, cell_type, tissue, treatment, sequencing, strand_specificity) %>%
    dplyr::filter(SRR %in% colnames(readcount_mat)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("SRR") %>%
    as.matrix()

  # Create DESeq object
  dds <- DESeqDataSetFromMatrix(countData = readcount_mat, colData = selected_sample_info, design = ~1)

  # Filtering
  smallestGroupSize <- 100
  keep <- rowSums(counts(dds) >= 100) >= smallestGroupSize
  dds <- dds[keep,]

  # Basic processing
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  saveRDS(dds, dds_obj_path)

  # Normalization by varianceStablizingTransformation and removing batch effects by limma::removeBatchEffect
  # Note the vsd was not updated, i.e., with batch effects
  model_matrix <- model.matrix(~1, colData(dds))
  vsd <- vst(dds, blind=FALSE)
  norm_mat <- assay(vsd) %>% limma::removeBatchEffect(batch=vsd$GSE, design=model_matrix)
  norm_mat %>% as.data.frame() %>% tibble::rownames_to_column("Geneid") %>% fwrite(norm_expr_tab_path)

  # Check the variance stablized 
  p_tab <- norm_mat %>% as.data.frame() %>% tibble::rownames_to_column("Geneid") %>%
    tidyr::pivot_longer(-Geneid, names_to = "Sample", values_to = "Expression_VST") %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Expr_median = median(Expression_VST)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Sample = forcats::fct_reorder(Sample, Expr_median))
  p <- ggplot(p_tab) +
    geom_violin(aes(y = Expression_VST, x = Sample)) +
    geom_boxplot(aes(y = Expression_VST, x = Sample), outlier.shape=NA, width = 0.5) +
    labs(x = "Exampel samples", y = "Normalized expression per gene") +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  file.path(project_dir, "outputs/analysis/coexpression/plots/deseq2.normalized_expression.pdf") %>%
    ggsave(plot = p, width = 25, height = 5)
}


# Coexpression analysis by WGCNA
allowWGCNAThreads(4)
norm_mat_t <- t(norm_mat)

coexpr_network_path <- file.path(project_dir, "outputs/analysis/coexpression/wgcna.coexpression_network.rds")
coexpr_tom_path <- file.path(project_dir, "outputs/analysis/coexpression/wgcna.coexpression_tom.rds")
tom_path <- file.path(project_dir, "outputs/analysis/coexpression/wgcna.tom")
if (file.exists(coexpr_network_path)) {
  network <- readRDS(coexpr_network_path)
} else {
  # Determine soft thresholds.
  # To identify which genes are in the same modules, WGCNA first creates a weighted network to define which genes are near each other.
  # The measure of “adjacency” it uses is based on the correlation matrix, but requires the definition of a threshold value, which in turn depends on a “power” parameter that defines the exponent used when transforming the correlation values.
  # The choice of power parameter will affect the number of modules identified, and the WGCNA modules provides the pickSoftThreshold() function to help identify good choices for this parameter.
  soft_thresholds <- pickSoftThreshold(norm_mat_t, networkType = "signed", verbose = 5)

  p_tab <- soft_thresholds$fitIndices %>% as.data.frame() %>% dplyr::mutate(signed_rsq = -sign(slope) * SFT.R.sq) %>%
    dplyr::select(Power, `Mean connectivities` = mean.k., `Signed R^2` = signed_rsq) %>%
    tidyr::pivot_longer(-Power)
  p <- ggplot(p_tab) +
    geom_point(aes(x = Power, y = value), size = 4) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    theme_classic() +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "Power", y = NULL)
  file.path(project_dir, "outputs/analysis/coexpression/plots/wgcna.determine_power.pdf") %>%
    ggsave(plot = p, width = 8, height = 4)

  # Here we used power 8
  # power <- p_tab %>% dplyr::filter(name == "Signed R^2", value >= 0.8) %>% dplyr::slice_min(Power, n=1) %>% dplyr::pull("Power")
  power <- soft_thresholds$powerEstimate

  # adjacency <- adjacency(norm_mat_t, power = power, type = "signed")
  # tom_mat <- TOMsimilarity(adjacency)
  # tom_dissam_mat <- 1 - tom_mat
  # tom_dissam_mat <- TOMdist(tom_dissam_mat)
  # gene_tree <- hclust(as.dist(tom_dissam_mat), method = "average")
  # modules <- cutreeDynamic(dendro = gene_tree, distM = tom_dissam_mat, method = "hybrid", minClusterSize = 50, deepSplit = 2, pamRespectsDendro = FALSE)
  # module_colors <- labels2colors(modules)
  # pdf(file.path(project_dir, "outputs/analysis/coexpression/plots/wgcna.dendrogram_and_modules.from_scratch.pdf"), width = 12, height = 7)
  # plotDendroAndColors(gene_tree, module_colors, "Modules", dendroLabels = FALSE, hang = 0.05, addGuide = TRUE, guideHang = 0.05)
  # dev.off()

  network <- blockwiseModules(
    norm_mat_t, power = power, networkType = "signed",
    pamRespectsDendro = FALSE, minModuleSize = 30, maxBlockSize = Inf,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    saveTOMs = TRUE, saveTOMFileBase = tom_path,
    numericLabels = TRUE,
    randomSeed = 31415926, nThreads = 4, verbose = TRUE)

  # Save the network into disk
  saveRDS(network, coexpr_network_path)

  # Check the modules and determine the heigh of tree
  merged_colors <- labels2colors(network$colors)
  for (idx in seq_along(network$dendrograms)) {
    pdf(file.path(project_dir, "outputs/analysis/coexpression/plots", paste0("wgcna.dendrogram_and_modules.", idx, ".pdf")), width = 12, height = 7)
    plotDendroAndColors(
      network$dendrograms[[idx]], merged_colors[network$blockGenes[[idx]]],
      "Module colors", dendroLabels = FALSE, hang = 0.05, addGuide = TRUE, guideHang = 0.05
    )
    dev.off()
  }

  # Check the module
  module_df <- data.frame(gene_id = names(network$colors), colors = labels2colors(network$colors))
  module_eigngenes <- moduleEigengenes(norm_mat_t, merged_colors)$eigengenes %>% orderMEs()

  # module_trait_cor <- WGCNA::cor(module_eigngenes, selected_sample_info, use = "p")
  # module_trait_p_val <- corPvalueStudent(module_trait_p_val, n_samples)

  # Correlation matrix
  tom <- TOMsimilarityFromExpr(norm_mat_t, power = power)
}
