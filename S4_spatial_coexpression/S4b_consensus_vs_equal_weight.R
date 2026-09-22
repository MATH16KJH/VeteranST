# =============================================================================
# S4b_consensus_vs_equal_weight.R
# Supplementary Section S4.1, Figure S8(b): strict consensus clustering vs.
# equal-weight aggregation for every combination of 2, 3, 4 and 5 measures.
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   Same benchmark and similarity matrices as S4a (neighbourhood N_NEIGHBORS).
#   For each of the 26 measure combinations (10 pairs, 10 triples, 5 quadruples,
#   1 quintuple) accuracy (ARI vs. ground truth) is computed for
#     * Consensus       strict multi-measure consensus (proposed)
#     * Weighted_Equal  equal-weight aggregation + one Leiden run
#   Equal weights are used as the representative aggregation baseline because
#   S4a showed little sensitivity to the weights.
#   As in S4a, every configuration is run N_REPLICATES times with different
#   Leiden seeds and the mean ARI (+/- 1 SD) is reported; a single Leiden run
#   is not a reliable basis for comparing combinations on a 100-gene graph.
#
# OUTPUTS  results/S4b_consensus_vs_equal_weight/
#   all_ari.jpg                    grouped bars, facet = number of measures (Fig S8b)
#   consensus_vs_weighted.csv      mean / SD per combination and method
#   consensus_vs_weighted_replicates.csv   per-replicate values
#
# INPUTS   data/simulated/{count_mat_simul,coordinates_simul}.rds (00b)
# USAGE    Rscript S4_spatial_coexpression/S4b_consensus_vs_equal_weight.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)
source("R/utils_evaluation.R")
source("R/utils_spatial_similarity.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

# ---- parameters --------------------------------------------------------------
N_NEIGHBORS <- 7
N_REPLICATES <- 100         # Leiden runs per configuration (mean is reported)
SEED        <- 2025
OUT_DIR     <- "results/S4b_consensus_vs_equal_weight"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- data --------------------------------------------------------------------
dat <- load_simulated_data()
gt  <- create_ground_truth_clustering(rownames(dat$expr))
cat("Computing similarity matrices (k =", N_NEIGHBORS, ")...\n")
mats <- compute_similarity_matrices(dat$expr, dat$coords, n_neighbors = N_NEIGHBORS)

# ---- all combinations --------------------------------------------------------
combos <- unlist(lapply(2:5, function(n) utils::combn(1:5, n, simplify = FALSE)), recursive = FALSE)
rows <- list()
for (i in seq_along(combos)) {
  combo <- combos[[i]]
  for (r in seq_len(N_REPLICATES)) {
    base <- SEED + 100000 * r + 100 * i
    cl <- consensus_clustering(mats, combo, seed = base)
    wl <- weighted_clustering(mats, combo, seed = base + 50)
    cons <- evaluate_clustering(cl, gt); wgt <- evaluate_clustering(wl, gt)
    rows[[length(rows) + 1]] <- data.frame(
      Combination = combo_label(combo), N_Measures = length(combo), Replicate = r,
      Method = c("Consensus", "Weighted_Equal"),
      N_Clusters = c(length(unique(cl)), length(unique(wl))),
      ARI = unname(c(cons["ARI"], wgt["ARI"])), NMI = unname(c(cons["NMI"], wgt["NMI"])))
  }
  d <- do.call(rbind, rows[(length(rows) - N_REPLICATES + 1):length(rows)])
  cat(sprintf("%-16s consensus ARI = %.4f | weighted ARI = %.4f  (means of %d runs)\n",
              combo_label(combo), mean(d$ARI[d$Method == "Consensus"]),
              mean(d$ARI[d$Method == "Weighted_Equal"]), N_REPLICATES))
}
per_rep <- do.call(rbind, rows)
write.csv(per_rep, file.path(OUT_DIR, "consensus_vs_weighted_replicates.csv"), row.names = FALSE)
# NB: SDs are computed before the means overwrite ARI / NMI inside summarise().
results <- per_rep %>%
  group_by(Combination, N_Measures, Method) %>%
  summarise(N_Replicates = n(), Mean_N_Clusters = mean(N_Clusters),
            SD_ARI = sd(ARI), SD_NMI = sd(NMI), ARI = mean(ARI), NMI = mean(NMI), .groups = "drop") %>%
  select(Combination, N_Measures, Method, N_Replicates, Mean_N_Clusters, ARI, SD_ARI, NMI, SD_NMI) %>%
  as.data.frame()
results$SD_ARI[is.na(results$SD_ARI)] <- 0
results$SD_NMI[is.na(results$SD_NMI)] <- 0
write.csv(results, file.path(OUT_DIR, "consensus_vs_weighted.csv"), row.names = FALSE)

# ---- figures -----------------------------------------------------------------
method_colors <- c(Consensus = "#E74C3C", Weighted_Equal = "#2E86C1")
big_theme <- theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title  = element_text(size = 14, face = "bold"),
        plot.title  = element_text(size = 16, face = "bold"),
        strip.text  = element_text(size = 12, face = "bold"),
        legend.position = "top", legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"))

plot_all <- function(metric) {
  sd_col <- paste0("SD_", metric)
  p <- ggplot(results, aes(x = Combination, y = .data[[metric]], fill = Method)) +
    geom_col(position = position_dodge(0.8), alpha = 0.8)
  if (N_REPLICATES > 1)
    p <- p + geom_errorbar(aes(ymin = .data[[metric]] - .data[[sd_col]],
                               ymax = .data[[metric]] + .data[[sd_col]]),
                           position = position_dodge(0.8), width = 0.25, alpha = 0.7)
  p +
    geom_text(aes(label = sprintf("%.3f", .data[[metric]]), y = .data[[metric]] + .data[[sd_col]]),
              position = position_dodge(0.8), vjust = -0.5, size = 3.5) +
    facet_wrap(~ N_Measures, scales = "free_x", ncol = 1,
               labeller = as_labeller(function(x) paste(x, "Measures"))) +
    scale_fill_manual(values = method_colors) +
    labs(title = paste0(metric, ": Consensus vs Weighted (Equal)"),
         subtitle = if (N_REPLICATES > 1) paste0("Mean of ", N_REPLICATES, " Leiden runs per bar; error bars = +/- 1 SD") else NULL,
         x = "Measure Combination", y = metric) +
    big_theme +
    ylim(0, max(results[[metric]] + results[[sd_col]], na.rm = TRUE) * 1.15)
}
save_figure(plot_all("ARI"), file.path(OUT_DIR, "all_ari.jpg"), 12, 10)
cat("Saved:", file.path(OUT_DIR, "all_ari.jpg"), "\n")
