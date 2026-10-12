# =============================================================================
# S4e_benchmark_clustering_methods.R
# Main-text Section 4.1 (Figure 2) and Supplementary Section S4.2 (Figure S10):
# proposed graph-based consensus clustering vs. alternative gene-clustering
# approaches
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   Every method is run N_RUNS times (different seeds) on the scDesign3
#   benchmark (10 true modules x 10 genes) and evaluated by
#     (a) accuracy      ARI vs. ground truth, mean +/- SD over runs  (main-text Fig 2a)
#     (b) stability     mean pairwise ARI between the N_RUNS partitions (main-text Fig 2b)
#     (c) computation   total wall-clock time of the N_RUNS runs      (Fig S10)
#
# METHODS (category)
#   Correlation Graph   (A) single metric: Leiden on (1+Pearson)/2, raw expression
#   Cosine Graph        (A) single metric: Leiden on cosine similarity, raw expression
#   Distance Graph      (A) single metric: Leiden on exp(-d/median d), Euclidean
#   Weighted Aggregated (B) the four smoothed measures of Section S4.1
#                           (PM + JH + JVT + JVB, k = N_NEIGHBORS) min-max
#                           normalised and combined with equal weights, followed
#                           by one Leiden run (weighted_clustering())
#   Classical Consensus (C) Monti et al. (2003): K in 2..15, 50 subsamples of 80%
#                           of genes, k-means, consensus CDF/delta-K rule,
#                           average-linkage cut of the consensus matrix
#   Hierarchical        (D) average linkage on 1 - Pearson, K by silhouette (2..15)
#   Our Method          proposed: smoothing (k = N_NEIGHBORS) -> PM+JH+JVT+JVB
#                           -> strict consensus (R/utils_spatial_similarity.R)
#
# OUTPUTS  results/S4e_benchmark_clustering_methods/
#   ari_comparison.jpg          mean ARI +/- SD per method (main-text Fig 2a)
#   stability_comparison.jpg    mean pairwise ARI +/- SD (main-text Fig 2b)
#   time_comparison.jpg         total run time per method (Fig S10)
#   comparison_results.csv, run_details.csv, stability_results.csv, timing_results.csv
#
# INPUTS   data/simulated/{count_mat_simul,coordinates_simul}.rds (00b)
# RUNTIME  ~10-20 min for N_RUNS = 100 (classical consensus dominates).
# USAGE    Rscript S4_spatial_coexpression/S4e_benchmark_clustering_methods.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)
source("R/utils_evaluation.R")
source("R/utils_spatial_similarity.R")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(cluster); library(igraph)
})

# ---- parameters --------------------------------------------------------------
N_RUNS       <- 100
GLOBAL_SEED  <- 2025
N_NEIGHBORS  <- DEFAULT_N_NEIGHBORS          # 37
MEASURE_IDX  <- DEFAULT_MEASURE_IDX          # PM + JH + JVT + JVB
K_RANGE      <- 2:15                         # candidate K for methods (C) and (D)
OUT_DIR      <- "results/S4e_benchmark_clustering_methods"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Per-method seed offsets (run r of a method uses GLOBAL_SEED + offset + 100 r).
METHODS <- data.frame(
  id       = c("A1_Correlation", "A2_Cosine", "A3_Distance", "B_Weighted",
               "C_Classical", "D_Hierarchical", "Proposed"),
  label    = c("Correlation Graph", "Cosine Graph", "Distance Graph", "Weighted Aggregated",
               "Classical Consensus", "Hierarchical", "Our Method"),
  category = c("(A) Single-metric", "(A) Single-metric", "(A) Single-metric", "(B) Weighted",
               "(C) Classical", "(D) Non-graph", "Proposed"),
  offset   = c(0, 10000, 20000, 30000, 40000, 60000, 70000),
  stringsAsFactors = FALSE)
CATEGORY_COLORS <- c("(A) Single-metric" = "#E74C3C", "(B) Weighted" = "#3498DB",
                     "(C) Classical" = "#9B59B6", "(D) Non-graph" = "#95A5A6", "Proposed" = "#27AE60")

# ---- data --------------------------------------------------------------------
dat  <- load_simulated_data()
expr <- dat$expr; coords <- dat$coords
gt   <- create_ground_truth_clustering(rownames(expr))
gene_names <- rownames(expr)

# ---- alternative methods -----------------------------------------------------
correlation_graph_clustering <- function(expr, seed) {
  cm <- suppressWarnings(stats::cor(t(expr), method = "pearson", use = "pairwise.complete.obs"))
  cm[is.na(cm)] <- 0
  leiden_on_similarity((cm + 1) / 2, seed = seed)
}

cosine_graph_clustering <- function(expr, seed) {
  norms <- sqrt(rowSums(expr^2)); norms[norms == 0] <- 1
  X  <- expr / norms
  cs <- X %*% t(X)
  cs[is.na(cs)] <- 0
  dimnames(cs) <- list(rownames(expr), rownames(expr))
  leiden_on_similarity(cs, seed = seed)                 # negatives are clipped
}

distance_graph_clustering <- function(expr, seed) {
  d <- as.matrix(stats::dist(expr, method = "euclidean"))
  med <- stats::median(d[upper.tri(d)]); if (med == 0) med <- 1
  leiden_on_similarity(exp(-d / med), seed = seed)
}

weighted_aggregated_clustering <- function(expr, coords, seed) {
  mats <- compute_similarity_matrices(expr, coords, n_neighbors = N_NEIGHBORS)
  weighted_clustering(mats, MEASURE_IDX, seed = seed)
}

classical_consensus_clustering <- function(expr, seed, K_range = K_RANGE,
                                           n_iterations = 50, subsample_ratio = 0.8) {
  n <- nrow(expr); size <- floor(n * subsample_ratio)
  area <- numeric(length(K_range)); consensus <- vector("list", length(K_range))
  for (ki in seq_along(K_range)) {
    K <- K_range[ki]; kseed <- seed + K * 1000
    conn <- matrix(0, n, n); ind <- matrix(0, n, n)
    for (h in seq_len(n_iterations)) {
      set.seed(kseed + h * 123)
      s <- sort(sample.int(n, size))
      km <- tryCatch({
        set.seed(kseed + h * 456)
        stats::kmeans(expr[s, , drop = FALSE], centers = K, nstart = 5, iter.max = 100)
      }, error = function(e) NULL)
      if (is.null(km)) next
      conn[s, s] <- conn[s, s] + outer(km$cluster, km$cluster, "==")
      ind[s, s]  <- ind[s, s] + 1
    }
    M <- ifelse(ind > 0, conn / ind, 0); diag(M) <- 1
    dimnames(M) <- list(rownames(expr), rownames(expr))
    consensus[[ki]] <- M
    v <- M[upper.tri(M)]
    xs <- seq(0, 1, by = 0.01); ys <- stats::ecdf(v)(xs)
    area[ki] <- sum(diff(xs) * (head(ys, -1) + tail(ys, -1)) / 2)   # CDF area
  }
  delta <- numeric(length(K_range)); delta[1] <- area[1]
  for (i in seq_along(K_range)[-1]) delta[i] <- if (area[i - 1] > 0) (area[i] - area[i - 1]) / area[i - 1] else 0
  if (length(K_range) >= 3) {
    elbow_idx <- which.min(diff(delta)) + 1
    best_k <- K_range[elbow_idx]
    if (best_k == min(K_range)) best_k <- K_range[which.max(delta[-1]) + 1]
  } else best_k <- K_range[which.max(area)]
  hc <- stats::hclust(stats::as.dist(1 - consensus[[match(best_k, K_range)]]), method = "average")
  labels <- stats::cutree(hc, k = best_k); names(labels) <- rownames(expr)
  labels
}

hierarchical_clustering <- function(expr, K_range = K_RANGE) {
  cm <- suppressWarnings(stats::cor(t(expr), method = "pearson", use = "pairwise.complete.obs"))
  cm[is.na(cm)] <- 0
  d  <- stats::as.dist(1 - cm)
  hc <- stats::hclust(d, method = "average")
  sil <- sapply(K_range, function(k) mean(cluster::silhouette(stats::cutree(hc, k = k), d)[, 3]))
  labels <- stats::cutree(hc, k = K_range[which.max(sil)]); names(labels) <- rownames(expr)
  labels
}

proposed_clustering <- function(expr, coords, seed) {
  mats <- compute_similarity_matrices(expr, coords, n_neighbors = N_NEIGHBORS)
  consensus_clustering(mats, MEASURE_IDX, seed = seed)
}

run_method <- function(id, seed) {
  set.seed(seed)
  switch(id,
         A1_Correlation = correlation_graph_clustering(expr, seed),
         A2_Cosine      = cosine_graph_clustering(expr, seed),
         A3_Distance    = distance_graph_clustering(expr, seed),
         B_Weighted     = weighted_aggregated_clustering(expr, coords, seed),
         C_Classical    = classical_consensus_clustering(expr, seed),
         D_Hierarchical = hierarchical_clustering(expr),
         Proposed       = proposed_clustering(expr, coords, seed))
}

# ---- benchmark ---------------------------------------------------------------
run_details <- list(); summary_rows <- list(); stability_rows <- list(); timing_rows <- list()
for (m in seq_len(nrow(METHODS))) {
  id <- METHODS$id[m]; label <- METHODS$label[m]
  cat("\n=== ", label, " ===\n", sep = "")
  labels_list <- vector("list", N_RUNS)
  t0 <- Sys.time()
  for (r in seq_len(N_RUNS)) {
    labels_list[[r]] <- run_method(id, GLOBAL_SEED + METHODS$offset[m] + 100 * r)
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  metrics <- t(sapply(labels_list, evaluate_clustering, ground_truth = gt))
  n_clust <- sapply(labels_list, function(l) length(unique(l)))
  run_details[[id]] <- data.frame(Method = label, Run = seq_len(N_RUNS),
                                  ARI = metrics[, "ARI"], NMI = metrics[, "NMI"], N_Clusters = n_clust)
  stab <- pairwise_stability(labels_list)
  summary_rows[[id]] <- data.frame(
    Method = label, Category = METHODS$category[m],
    Mean_ARI = mean(metrics[, "ARI"]), SD_ARI = sd(metrics[, "ARI"]),
    Mean_NMI = mean(metrics[, "NMI"]), SD_NMI = sd(metrics[, "NMI"]),
    Mean_N_Clusters = mean(n_clust), N_Runs = N_RUNS, Time_sec = elapsed)
  stability_rows[[id]] <- data.frame(Method = label, Category = METHODS$category[m],
                                     Mean_Pairwise_ARI = stab$mean, SD_Pairwise_ARI = stab$sd,
                                     Min = stab$min, Max = stab$max)
  timing_rows[[id]] <- data.frame(Method = label, Category = METHODS$category[m], Time_sec = elapsed)
  cat(sprintf("mean ARI = %.4f (sd %.4f) | stability = %.4f | %.1f s for %d runs\n",
              mean(metrics[, "ARI"]), sd(metrics[, "ARI"]), stab$mean, elapsed, N_RUNS))
}
results   <- do.call(rbind, summary_rows);   rownames(results) <- NULL
details   <- do.call(rbind, run_details);    rownames(details) <- NULL
stability <- do.call(rbind, stability_rows); rownames(stability) <- NULL
timing    <- do.call(rbind, timing_rows);    rownames(timing) <- NULL
results$ARI_Rank <- rank(-results$Mean_ARI, ties.method = "min")

write.csv(results,   file.path(OUT_DIR, "comparison_results.csv"), row.names = FALSE)
write.csv(details,   file.path(OUT_DIR, "run_details.csv"), row.names = FALSE)
write.csv(stability, file.path(OUT_DIR, "stability_results.csv"), row.names = FALSE)
write.csv(timing,    file.path(OUT_DIR, "timing_results.csv"), row.names = FALSE)
cat("\n=== Summary ===\n"); print(results, row.names = FALSE)
cat("\n=== Stability ===\n"); print(stability, row.names = FALSE)

# ---- figures -----------------------------------------------------------------
bar_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, face = "bold", margin = margin(b = 5)),
        plot.subtitle = element_text(size = 10, color = "gray40", margin = margin(b = 15)),
        legend.position = "bottom", legend.box.margin = margin(t = 15),
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        plot.margin = margin(t = 15, r = 10, b = 10, l = 10))

# (a) accuracy
acc <- results %>% mutate(Method = factor(Method, levels = Method[order(Mean_ARI)]))
p_acc <- ggplot(acc, aes(x = Method, y = Mean_ARI, fill = Category)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(Mean_ARI - SD_ARI, 0), ymax = pmin(Mean_ARI + SD_ARI, 1)),
                width = 0.3, linewidth = 0.7) +
  geom_text(aes(y = pmin(Mean_ARI + SD_ARI, 1), label = sprintf("%.3f", Mean_ARI)),
            hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = CATEGORY_COLORS) +
  coord_flip(ylim = c(0, 1.15)) +
  labs(title = "Clustering Performance Comparison: ARI",
       subtitle = paste0("Mean \u00b1 SD over ", N_RUNS, " runs"),
       x = NULL, y = "Adjusted Rand Index (ARI)", fill = "Method Category") +
  bar_theme
save_figure(p_acc, file.path(OUT_DIR, "ari_comparison.jpg"), 10, 7)

# (b) stability
stb <- stability %>% mutate(Method = factor(Method, levels = Method[order(Mean_Pairwise_ARI)]))
p_stab <- ggplot(stb, aes(x = Method, y = Mean_Pairwise_ARI, fill = Category)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(Mean_Pairwise_ARI - SD_Pairwise_ARI, 0),
                    ymax = pmin(Mean_Pairwise_ARI + SD_Pairwise_ARI, 1)),
                width = 0.3, linewidth = 0.8) +
  geom_text(aes(y = pmin(Mean_Pairwise_ARI + SD_Pairwise_ARI, 1),
                label = sprintf("%.3f", Mean_Pairwise_ARI)),
            hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = CATEGORY_COLORS) +
  coord_flip(ylim = c(0, 1.15)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Clustering Stability: Mean Pairwise ARI",
       subtitle = paste0(N_RUNS, " runs per method; ", N_RUNS * (N_RUNS - 1) / 2,
                         " pairwise comparisons (error bars = \u00b11 SD)"),
       x = NULL, y = "Mean Pairwise ARI", fill = "Category") +
  bar_theme
save_figure(p_stab, file.path(OUT_DIR, "stability_comparison.jpg"), 10, 7)

# (c) computation time
tm <- timing %>% arrange(Time_sec) %>% mutate(Method = factor(Method, levels = Method))
p_time <- ggplot(tm, aes(x = Method, y = Time_sec, fill = Category)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1fs", Time_sec)), hjust = -0.1, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = CATEGORY_COLORS) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Computation Time by Method", subtitle = paste0("Total time for ", N_RUNS, " runs"),
       x = NULL, y = "Time (seconds)", fill = "Category") +
  bar_theme
save_figure(p_time, file.path(OUT_DIR, "time_comparison.jpg"), 10, 7)
cat("\nDone.\n")
