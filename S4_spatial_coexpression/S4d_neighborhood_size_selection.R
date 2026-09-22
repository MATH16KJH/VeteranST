# =============================================================================
# S4d_neighborhood_size_selection.R
# Supplementary Section S4.1, Table S4: choice of the neighbourhood size used
# for expression smoothing (7, 19, 37 or 61 spots on the hexagonal grid).
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   With the four selected measures (PM + JH + JVT + JVB) the whole pipeline
#   (smoothing -> similarity matrices -> strict consensus) is run N_REPLICATES
#   times per neighbourhood size on the scDesign3 benchmark with different
#   Leiden seeds, and the mean +/- SD of ARI and NMI and the mean number of
#   recovered clusters (10 true modules) are reported. The two best sizes are
#   compared with Welch's t-test on the per-replicate ARIs (the numbers quoted
#   in the text). Negative Pearson correlations are clipped to 0, as in every
#   other S4 script.
#
# OUTPUTS  results/S4d_neighborhood_size_selection/
#   neighbor_comparison_4measures.csv       Table S4 (mean / SD of ARI and NMI, mean cluster count)
#   neighbor_comparison_replicates.csv      per-replicate values
#
# INPUTS   data/simulated/{count_mat_simul,coordinates_simul}.rds (00b)
# USAGE    Rscript S4_spatial_coexpression/S4d_neighborhood_size_selection.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)
source("R/utils_evaluation.R")
source("R/utils_spatial_similarity.R")
suppressPackageStartupMessages(library(dplyr))

# ---- parameters --------------------------------------------------------------
NEIGHBOR_VALUES  <- c(7, 19, 37, 61)
MEASURE_IDX      <- DEFAULT_MEASURE_IDX     # PM + JH + JVT + JVB
N_REPLICATES     <- 100                     # Leiden runs per neighbourhood size (mean +/- SD)
SEED             <- 2025
OUT_DIR          <- "results/S4d_neighborhood_size_selection"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- data --------------------------------------------------------------------
dat <- load_simulated_data()
gt  <- create_ground_truth_clustering(rownames(dat$expr))

# ---- run ---------------------------------------------------------------------
rows <- list()
for (k in NEIGHBOR_VALUES) {
  cat("\n=== neighbourhood size", k, "===\n")
  mats <- compute_similarity_matrices(dat$expr, dat$coords, n_neighbors = k)
  for (r in seq_len(N_REPLICATES)) {
    labels <- consensus_clustering(mats, MEASURE_IDX, seed = SEED + 1000 * r + k)
    m <- evaluate_clustering(labels, gt)
    rows[[length(rows) + 1]] <- data.frame(N_Neighbors = k, Replicate = r,
                                           N_Clusters = length(unique(labels)),
                                           ARI = unname(m["ARI"]), NMI = unname(m["NMI"]))
    cat(sprintf("  replicate %d: %d clusters, ARI = %.4f, NMI = %.4f\n",
                r, length(unique(labels)), m["ARI"], m["NMI"]))
  }
}
per_rep <- do.call(rbind, rows)
# NB: SDs are computed before the means overwrite ARI / NMI inside summarise().
summary_df <- per_rep %>%
  group_by(N_Neighbors) %>%
  summarise(N_Replicates = n(), Mean_N_Clusters = mean(N_Clusters),
            SD_ARI = sd(ARI), SD_NMI = sd(NMI), ARI = mean(ARI), NMI = mean(NMI), .groups = "drop") %>%
  select(N_Neighbors, N_Replicates, Mean_N_Clusters, ARI, SD_ARI, NMI, SD_NMI)
write.csv(summary_df, file.path(OUT_DIR, "neighbor_comparison_4measures.csv"), row.names = FALSE)
write.csv(per_rep, file.path(OUT_DIR, "neighbor_comparison_replicates.csv"), row.names = FALSE)

cat("\n=== Table S4 (", combo_label(MEASURE_IDX), ") ===\n", sep = "")
print(as.data.frame(summary_df), row.names = FALSE)
# Welch's t-tests on the per-replicate ARIs: the two best sizes against each
# other, and the default (DEFAULT_N_NEIGHBORS) against every other size.
welch <- function(k1, k2) {
  tt <- stats::t.test(per_rep$ARI[per_rep$N_Neighbors == k1], per_rep$ARI[per_rep$N_Neighbors == k2])
  sprintf("k = %d vs k = %d: mean difference %.4f, Welch t = %.2f, p = %.3f", k1, k2,
          diff(rev(tt$estimate)), tt$statistic, tt$p.value)
}
ord <- order(-summary_df$ARI)
cat(welch(summary_df$N_Neighbors[ord[1]], summary_df$N_Neighbors[ord[2]]), " (two best sizes)\n", sep = "")
for (k in setdiff(NEIGHBOR_VALUES, DEFAULT_N_NEIGHBORS)) cat(welch(DEFAULT_N_NEIGHBORS, k), "\n")
cat("Saved:", file.path(OUT_DIR, "neighbor_comparison_4measures.csv"), "\n")
