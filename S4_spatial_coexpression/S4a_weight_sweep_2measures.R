# =============================================================================
# S4a_weight_sweep_2measures.R
# Supplementary Section S4.1, Figure S8(a): for every pair of similarity
# measures, strict consensus clustering vs. weighted aggregation over a sweep
# of the pair's weights.
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   On the scDesign3 benchmark (10 true modules x 10 genes) the five similarity
#   matrices are computed once (neighbourhood of N_NEIGHBORS spots). For each of
#   the 10 measure pairs:
#     * consensus:  Leiden on each of the two graphs, strict consensus
#     * weighted:   w * S1 + (1-w) * S2 (min-max normalised) followed by Leiden,
#                   for w = 0, 0.1, ..., 1  (the end points use 0.001 / 0.999 so
#                   that both matrices are always present, as in the original run)
#   Accuracy = ARI against the ground truth. NMI is written to the CSV only.
#
#   Leiden is stochastic and, on a 100-gene graph, a single run can move the
#   ARI of one configuration by about +/- 0.1. Every configuration is
#   therefore evaluated N_REPLICATES times with different seeds and the MEAN
#   ARI is plotted (with +/- 1 SD error bars); the per-replicate values are
#   written to a CSV.
#
# OUTPUTS  results/S4a_weight_sweep_2measures/
#   2m_ari_<A>_<B>.jpg              one panel per pair (Figure S8a)
#   weight_sweep_results.csv        mean and SD of ARI / NMI per configuration
#   weight_sweep_replicates.csv     per-replicate values
#   weight_sweep_summary.csv        consensus vs. best weighted per pair (means)
#
# INPUTS   data/simulated/{count_mat_simul,coordinates_simul}.rds (00b)
# USAGE    Rscript S4_spatial_coexpression/S4a_weight_sweep_2measures.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)
source("R/utils_evaluation.R")
source("R/utils_spatial_similarity.R")
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters --------------------------------------------------------------
N_NEIGHBORS  <- 7                 # neighbourhood used for the measure-selection study (S4.1)
WEIGHT_STEPS <- seq(0, 1, by = 0.1)
N_REPLICATES <- 100               # Leiden runs per configuration (mean is reported)
SEED         <- 2025              # Leiden is stochastic; fix for reproducibility
OUT_DIR      <- "results/S4a_weight_sweep_2measures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- data --------------------------------------------------------------------
dat <- load_simulated_data()
gt  <- create_ground_truth_clustering(rownames(dat$expr))
cat("Computing similarity matrices (k =", N_NEIGHBORS, ")...\n")
mats <- compute_similarity_matrices(dat$expr, dat$coords, n_neighbors = N_NEIGHBORS)

# ---- sweep -------------------------------------------------------------------
pairs <- utils::combn(1:5, 2, simplify = FALSE)
rows  <- list()
for (pi in seq_along(pairs)) {
  pair <- pairs[[pi]]
  lab <- combo_label(pair)
  cat("\n---", lab, "---\n")
  for (r in seq_len(N_REPLICATES)) {
    base <- SEED + 100000 * r + 1000 * pi
    cons <- evaluate_clustering(consensus_clustering(mats, pair, seed = base), gt)
    rows[[length(rows) + 1]] <- data.frame(Measure1 = MEASURE_NAMES[pair[1]],
                                           Measure2 = MEASURE_NAMES[pair[2]], Method = "Consensus",
                                           Weight1 = NA, Weight2 = NA, Replicate = r,
                                           ARI = unname(cons["ARI"]), NMI = unname(cons["NMI"]))
    for (w1 in WEIGHT_STEPS) {
      w <- c(w1, 1 - w1)
      if (w1 == 0) w <- c(0.001, 0.999)
      if (w1 == 1) w <- c(0.999, 0.001)
      wm <- evaluate_clustering(weighted_clustering(mats, pair, weights = w,
                                                    seed = base + 10 + round(w1 * 10)), gt)
      rows[[length(rows) + 1]] <- data.frame(Measure1 = MEASURE_NAMES[pair[1]],
                                             Measure2 = MEASURE_NAMES[pair[2]], Method = "Weighted",
                                             Weight1 = w1, Weight2 = 1 - w1, Replicate = r,
                                             ARI = unname(wm["ARI"]), NMI = unname(wm["NMI"]))
    }
  }
  sub <- do.call(rbind, rows)
  sub <- sub[sub$Measure1 == MEASURE_NAMES[pair[1]] & sub$Measure2 == MEASURE_NAMES[pair[2]], ]
  cat(sprintf("Consensus mean ARI = %.4f | best weighted mean ARI = %.4f  (%d replicates)\n",
              mean(sub$ARI[sub$Method == "Consensus"]),
              max(tapply(sub$ARI[sub$Method == "Weighted"], sub$Weight1[sub$Method == "Weighted"], mean)),
              N_REPLICATES))
}
per_rep <- do.call(rbind, rows)
rownames(per_rep) <- NULL
write.csv(per_rep, file.path(OUT_DIR, "weight_sweep_replicates.csv"), row.names = FALSE)

# mean / SD over replicates (Weight columns are NA for the consensus rows)
key <- with(per_rep, paste(Measure1, Measure2, Method, Weight1))
results <- do.call(rbind, lapply(split(per_rep, key), function(d) data.frame(
  Measure1 = d$Measure1[1], Measure2 = d$Measure2[1], Method = d$Method[1],
  Weight1 = d$Weight1[1], Weight2 = d$Weight2[1], N_Replicates = nrow(d),
  ARI = mean(d$ARI), SD_ARI = sd(d$ARI), NMI = mean(d$NMI), SD_NMI = sd(d$NMI))))
results <- results[order(results$Measure1, results$Measure2, results$Method, results$Weight1), ]
rownames(results) <- NULL
write.csv(results, file.path(OUT_DIR, "weight_sweep_results.csv"), row.names = FALSE)

# ---- figures -----------------------------------------------------------------
plot_sweep <- function(sub, metric) {
  cons_val <- sub[sub$Method == "Consensus", metric]
  wd <- sub[sub$Method == "Weighted", ]
  wd$idx   <- seq_len(nrow(wd))
  wd$label <- sprintf("(%.1f,%.1f)", wd$Weight1, wd$Weight2)
  sd_col <- paste0("SD_", metric)
  cons_sd <- sub[sub$Method == "Consensus", sd_col]
  p <- ggplot()
  if (N_REPLICATES > 1) {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = cons_val - cons_sd, ymax = cons_val + cons_sd,
                      fill = "#E74C3C", alpha = 0.12) +
      geom_errorbar(data = wd, aes(x = idx, ymin = .data[[metric]] - .data[[sd_col]],
                                   ymax = .data[[metric]] + .data[[sd_col]]),
                    width = 0.25, color = "#2E86C1", alpha = 0.7)
  }
  p +
    geom_line(data = wd, aes(x = idx, y = .data[[metric]], color = "Weighted"), linewidth = 1.2) +
    geom_point(data = wd, aes(x = idx, y = .data[[metric]], color = "Weighted"), size = 3.5) +
    geom_hline(aes(yintercept = cons_val, color = "Consensus"), linetype = "dashed", linewidth = 1.2) +
    scale_color_manual(name = NULL, values = c(Consensus = "#E74C3C", Weighted = "#2E86C1"),
                       breaks = c("Consensus", "Weighted")) +
    scale_x_continuous(breaks = wd$idx, labels = wd$label) +
    ylim(0, 1) +
    labs(title = paste0(metric, ": ", sub$Measure1[1], " + ", sub$Measure2[1]),
         subtitle = if (N_REPLICATES > 1) paste0("Mean of ", N_REPLICATES, " Leiden runs; band / bars = +/- 1 SD") else NULL,
         x = "Weight Combination", y = metric) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
          axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12, face = "bold"),
          legend.position = "top", legend.text = element_text(size = 12, face = "bold"),
          legend.key.width = unit(1.5, "cm"),
          panel.grid.major = element_line(color = "gray85", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 1),
          plot.margin = margin(t = 10, r = 15, b = 10, l = 10))
}

summary_rows <- list()
for (pair in pairs) {
  sub <- results[results$Measure1 == MEASURE_NAMES[pair[1]] & results$Measure2 == MEASURE_NAMES[pair[2]], ]
  a <- MEASURE_ABBR[MEASURE_NAMES[pair[1]]]; b <- MEASURE_ABBR[MEASURE_NAMES[pair[2]]]
  save_figure(plot_sweep(sub, "ARI"), file.path(OUT_DIR, sprintf("2m_ari_%s_%s.jpg", a, b)), 10, 7)
  wd <- sub[sub$Method == "Weighted", ]
  best <- which.max(wd$ARI)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Pair = combo_label(pair),
    Consensus_ARI = sub$ARI[sub$Method == "Consensus"],
    Consensus_SD_ARI = sub$SD_ARI[sub$Method == "Consensus"],
    Best_Weighted_ARI = wd$ARI[best],
    Best_Weighted_SD_ARI = wd$SD_ARI[best],
    Best_Weights = sprintf("(%.1f,%.1f)", wd$Weight1[best], wd$Weight2[best]),
    Consensus_NMI = sub$NMI[sub$Method == "Consensus"],
    Best_Weighted_NMI = max(wd$NMI))
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "weight_sweep_summary.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)
cat(sprintf("\nConsensus ARI > best weighted ARI in %d of %d pairs\n",
            sum(summary_df$Consensus_ARI > summary_df$Best_Weighted_ARI), nrow(summary_df)))
