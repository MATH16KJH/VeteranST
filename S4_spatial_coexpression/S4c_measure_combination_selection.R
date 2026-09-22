# =============================================================================
# S4c_measure_combination_selection.R
# Supplementary Section S4.1, Figure S9: which combination of similarity
# measures gives the most accurate consensus clustering?
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   Same benchmark and similarity matrices as S4a/S4b (neighbourhood
#   N_NEIGHBORS). Strict consensus clustering is run for each of the 26 measure
#   combinations and its ARI against the ground truth is recorded. Because only
#   Leiden is stochastic (the similarity matrices are deterministic), a
#   "replicate" is a re-run of Leiden with a fresh seed; a single run is not a
#   reliable basis for choosing a combination on a 100-gene graph, so
#   N_REPLICATES Leiden runs are averaged and the mean +/- SD is reported,
#   together with the mean number of consensus clusters (10 true modules).
#
# OUTPUTS  results/S4c_measure_combination_selection/
#   consensus_analysis_by_n_measures.jpg   box plot of ARI by number of measures (Fig S9)
#   consensus_combinations.csv             mean / SD of ARI and NMI, cluster count per combination
#   consensus_combinations_replicates.csv  per-replicate values
#
# INPUTS   data/simulated/{count_mat_simul,coordinates_simul}.rds (00b)
# USAGE    Rscript S4_spatial_coexpression/S4c_measure_combination_selection.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)
source("R/utils_evaluation.R")
source("R/utils_spatial_similarity.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

# ---- parameters --------------------------------------------------------------
N_NEIGHBORS  <- 7
N_REPLICATES <- 100        # Leiden runs per combination (mean +/- SD is reported)
TRUE_MODULES <- 10         # ground truth: 10 genes x 10 replicates
SEED         <- 2025
OUT_DIR      <- "results/S4c_measure_combination_selection"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- data --------------------------------------------------------------------
dat <- load_simulated_data()
gt  <- create_ground_truth_clustering(rownames(dat$expr))
cat("Computing similarity matrices (k =", N_NEIGHBORS, ")...\n")
mats <- compute_similarity_matrices(dat$expr, dat$coords, n_neighbors = N_NEIGHBORS)

# ---- consensus for every combination ----------------------------------------
combos <- unlist(lapply(2:5, function(n) utils::combn(1:5, n, simplify = FALSE)), recursive = FALSE)
rows <- list()
for (r in seq_len(N_REPLICATES)) {
  cat("Replicate", r, "/", N_REPLICATES, "\n")
  for (i in seq_along(combos)) {
    combo  <- combos[[i]]
    labels <- consensus_clustering(mats, combo, seed = SEED + 1000 * r + 10 * i)
    m <- evaluate_clustering(labels, gt)
    rows[[length(rows) + 1]] <- data.frame(
      Replicate = r, Combination = combo_label(combo), N_Measures = length(combo),
      N_Clusters = length(unique(labels)), ARI = unname(m["ARI"]), NMI = unname(m["NMI"]))
  }
}
per_rep <- do.call(rbind, rows)

# NB: inside summarise() a new column shadows the input column of the same name,
# so the SDs must be computed before the means are assigned to ARI / NMI.
agg <- per_rep %>%
  group_by(Combination, N_Measures) %>%
  summarise(N_Replicates = n(), Mean_N_Clusters = mean(N_Clusters),
            SD_ARI = sd(ARI), SD_NMI = sd(NMI), ARI = mean(ARI), NMI = mean(NMI), .groups = "drop") %>%
  select(Combination, N_Measures, N_Replicates, Mean_N_Clusters, ARI, SD_ARI, NMI, SD_NMI)
agg <- agg[order(-agg$ARI), ]
write.csv(agg, file.path(OUT_DIR, "consensus_combinations.csv"), row.names = FALSE)
write.csv(per_rep, file.path(OUT_DIR, "consensus_combinations_replicates.csv"), row.names = FALSE)

cat("\n=== Top combinations by mean ARI ===\n")
print(head(as.data.frame(agg), 10), row.names = FALSE)
cat("\nDefault used in the manuscript:", combo_label(DEFAULT_MEASURE_IDX), "\n")
def_row <- agg[agg$Combination == combo_label(DEFAULT_MEASURE_IDX), ]
cat(sprintf("  mean ARI = %.3f (SD %.3f), %.1f clusters on average (rank %d of %d)\n",
            def_row$ARI, def_row$SD_ARI, def_row$Mean_N_Clusters,
            which(agg$Combination == combo_label(DEFAULT_MEASURE_IDX)), nrow(agg)))

# ---- Figure S9: distribution of ARI by number of measures -------------------
plot_by_n <- function(metrics) {
  df <- agg %>% select(N_Measures, all_of(metrics)) %>%
    pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value")
  ggplot(df, aes(x = factor(N_Measures), y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    geom_jitter(aes(color = Metric),
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                alpha = 0.5, size = 2) +
    scale_fill_manual(values = c(ARI = "#E74C3C", NMI = "#3498DB")) +
    scale_color_manual(values = c(ARI = "#C0392B", NMI = "#2874A6")) +
    labs(title = "Consensus Clustering Performance by Number of Measures",
         subtitle = paste0("Mean over ", N_REPLICATES, " Leiden replicates; one point per measure combination"),
         x = "Number of Measures in Consensus", y = "Score") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 14),
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title  = element_text(size = 14, face = "bold"),
          legend.position = "top", legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = "bold")) +
    ylim(0, 1)
}
save_figure(plot_by_n("ARI"), file.path(OUT_DIR, "consensus_analysis_by_n_measures.jpg"), 10, 8)
cat("Saved:", file.path(OUT_DIR, "consensus_analysis_by_n_measures.jpg"), "\n")
