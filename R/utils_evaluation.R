# =============================================================================
# R/utils_evaluation.R -- evaluation helpers for the co-expression benchmarks
# -----------------------------------------------------------------------------
# Sourced by all scripts under S4_spatial_coexpression/.
#   create_ground_truth_clustering()  ground-truth gene modules of the scDesign3
#                                     benchmark data (Supplementary Section S4)
#   evaluate_clustering()             ARI (and NMI) of a clustering vs. truth
#   pairwise_stability()              mean pairwise ARI over repeated runs
#   save_figure()                     ggsave wrapper (300 dpi JPEG by default)
# NOTE: NMI is computed and written to the CSV outputs for completeness, but
# only ARI is reported in the manuscript.
# =============================================================================

suppressPackageStartupMessages({
  library(aricode)
  library(ggplot2)
})

#' Ground truth for the simulated benchmark.
#' Simulated genes are named "<gene>_<replicate>" (see 00b_generate_simulated_data.R);
#' all replicates of the same gene form one true module. The committed file
#' data/simulated/ground_truth.rds (Gene, Gene_Base_Name, True_Cluster) is used
#' when present and is checked against the partition implied by the gene names.
create_ground_truth_clustering <- function(gene_names, verbose = TRUE,
                                           file = if (exists("SIMULATED_DATA")) SIMULATED_DATA$ground_truth else NULL) {
  gene_base_names <- sub("_[0-9]+$", "", gene_names)   # strip the replicate suffix only
  unique_genes    <- unique(gene_base_names)
  result_df <- data.frame(
    Gene           = gene_names,
    Gene_Base_Name = gene_base_names,
    True_Cluster   = match(gene_base_names, unique_genes),
    stringsAsFactors = FALSE
  )
  if (!is.null(file) && file.exists(file)) {
    gt <- readRDS(file)
    if (!all(gene_names %in% gt$Gene)) stop("ground truth file does not cover all genes: ", file)
    gt <- gt[match(gene_names, gt$Gene), ]
    if (aricode::ARI(gt$True_Cluster, result_df$True_Cluster) < 1) {
      stop("ground truth file disagrees with the '<gene>_<replicate>' naming: ", file)
    }
    result_df <- data.frame(Gene = gt$Gene, Gene_Base_Name = gt$Gene_Base_Name,
                            True_Cluster = as.integer(gt$True_Cluster), stringsAsFactors = FALSE)
  }
  result_df <- result_df[order(result_df$True_Cluster, result_df$Gene), ]
  rownames(result_df) <- NULL
  if (verbose) {
    cat("Ground truth: ", length(gene_names), " genes in ",
        length(unique_genes), " true modules\n", sep = "")
    print(table(result_df$True_Cluster))
  }
  result_df
}

#' ARI and NMI of a named label vector against the ground truth.
evaluate_clustering <- function(pred_labels, ground_truth) {
  matched <- intersect(names(pred_labels), ground_truth$Gene)
  if (length(matched) < 2) return(c(ARI = NA_real_, NMI = NA_real_))
  pred  <- as.vector(pred_labels[matched])
  truth <- ground_truth$True_Cluster[match(matched, ground_truth$Gene)]
  c(ARI = aricode::ARI(pred, truth), NMI = aricode::NMI(pred, truth))
}

#' Mean pairwise ARI between all pairs of clusterings in `label_list`
#' (clustering stability; main-text Figure 2b).
pairwise_stability <- function(label_list) {
  n <- length(label_list)
  if (n < 2) return(list(mean = NA_real_, sd = NA_real_, values = numeric(0)))
  vals <- numeric(0)
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      common <- intersect(names(label_list[[i]]), names(label_list[[j]]))
      vals <- c(vals, aricode::ARI(as.vector(label_list[[i]][common]),
                                   as.vector(label_list[[j]][common])))
    }
  }
  list(mean = mean(vals), sd = stats::sd(vals), min = min(vals), max = max(vals), values = vals)
}

#' Save a ggplot as JPEG (creates the directory if needed).
save_figure <- function(plot, path, width = 10, height = 7, dpi = 300, quality = 100) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, plot, width = width, height = height, dpi = dpi,
         device = "jpeg", quality = quality)
  cat("Saved:", path, "\n")
  invisible(path)
}
