# =============================================================================
# S1_type1_error_gene_level.R
# Supplementary Section S1.2-S1.4: gene-wise type I error of variance-
# homogeneity tests (Figure S1, S3) and of mean-comparison tests (Figure S2, S4)
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   For every gene of each real dataset (P16, BR, D14) the expression values are
#   randomly shuffled across spots N_SIM times, which makes all BANYAN tissue
#   domains samples from one population (null hypothesis of equal variance /
#   equal mean). Each candidate test is applied to every shuffle and the
#   empirical type I error of a gene is the fraction of shuffles rejected at
#   alpha = 0.05 and 0.01.
#
# TESTS
#   family "variance" (Section S1.2, Figure S1; Section S1.3, alpha = 0.05):
#     Brown_Forsythe  one-way F on |x - group median|   (adopted for DVG testing)
#     Fligner         Fligner-Killeen test
#     Bartlett        Bartlett test
#     Levene_mean     one-way F on |x - group mean|
#     Levene_rank     Levene (mean) applied to rank-transformed observations
#     Levene_trim5    one-way F on |x - 5% trimmed group mean|
#     Levene_trim10   one-way F on |x - 10% trimmed group mean|
#   family "mean" (Section S1.2, Figure S2; Section S1.4, alpha = 0.05):
#     ANOVA           one-way ANOVA
#     Kruskal_Wallis  Kruskal-Wallis test
#
# OUTPUTS  results/S1_gene_level_type1_error/
#   <Test>/scatter_comparison_<Test>_<alpha>.jpg   empirical error vs non-zero ratio
#   <Test>/comparison_hist_<Test>_<alpha>.jpg      histogram of empirical errors
#   empirical_error_<family>.rds                   gene-level results (all tests)
#   summary_<family>.csv                           mean / median / variance per test
#
# INPUTS   data/real/count_mat_{p16,br,d14}.rda   (00_data/00a_prepare_real_data.R)
# NOTE     Test definitions follow the development scripts exactly (Levene-type
#          tests = one-way F on |x - centre| with centre = mean, 5%/10% trimmed
#          mean, median; Levene_rank = Levene(mean) on mid-ranks). Those scripts
#          did not set a seed; BASE_SEED + gene index is used here so that
#          results are reproducible.
# RUNTIME  all genes x 100 shuffles x 9 tests: ~ hours on a workstation; set
#          N_GENES to a small number for a smoke test.
# USAGE    Rscript S1_gene_level_type1_error/S1_type1_error_gene_level.R
#          (from the repository root)
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages({
  library(Matrix)
  library(parallel)
  library(ggplot2)
})
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)

# ---- parameters --------------------------------------------------------------
TEST_FAMILIES <- c("variance", "mean")   # which families to run
DATASETS      <- DATASET_ORDER           # c("P16", "BR", "D14")
N_SIM         <- 100                     # shuffles per gene
ALPHAS        <- c(0.05, 0.01)           # nominal levels
N_GENES       <- NULL                    # NULL = all genes; e.g. 200 for a smoke test
BASE_SEED     <- 2024                    # gene i is shuffled with seed BASE_SEED + i
CHUNK_SIZE    <- 500                     # genes per parallel batch
N_CORES       <- default_cores()
OUT_DIR       <- "results/S1_gene_level_type1_error"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- test statistics ---------------------------------------------------------
# One-way ANOVA F test in closed form (fast; avoids lm()).
oneway_F_pvalue <- function(y, group, group_sizes) {
  n <- length(y); m <- length(group_sizes)
  grand  <- mean(y)
  gmeans <- as.vector(rowsum(y, group)) / group_sizes
  ssb    <- sum(group_sizes * (gmeans - grand)^2)
  ssw    <- sum((y - grand)^2) - ssb
  if (!is.finite(ssw) || ssw <= 0) return(NA_real_)
  Fstat <- (ssb / (m - 1)) / (ssw / (n - m))
  stats::pf(Fstat, m - 1, n - m, lower.tail = FALSE)
}

# Levene-type tests: one-way F on absolute deviations from a group centre.
# center_fun = median gives the Brown-Forsythe statistic T_k of Section S1.1.
levene_family_pvalue <- function(y, group, group_sizes, center_fun) {
  centers <- tapply(y, group, center_fun)
  d <- abs(y - centers[group])
  oneway_F_pvalue(d, group, group_sizes)
}

# Every test takes (y, group_int, group_factor, group_sizes) and returns a p-value.
VARIANCE_TESTS <- list(
  Brown_Forsythe = function(y, g, gf, gs) levene_family_pvalue(y, g, gs, stats::median),
  Fligner        = function(y, g, gf, gs) stats::fligner.test(y, gf)$p.value,
  Bartlett       = function(y, g, gf, gs) stats::bartlett.test(y, gf)$p.value,
  Levene_mean    = function(y, g, gf, gs) levene_family_pvalue(y, g, gs, mean),
  Levene_rank    = function(y, g, gf, gs) levene_family_pvalue(rank(y), g, gs, mean),
  Levene_trim5   = function(y, g, gf, gs) levene_family_pvalue(y, g, gs, function(x) mean(x, trim = 0.05)),
  Levene_trim10  = function(y, g, gf, gs) levene_family_pvalue(y, g, gs, function(x) mean(x, trim = 0.10))
)
MEAN_TESTS <- list(
  ANOVA          = function(y, g, gf, gs) oneway_F_pvalue(y, g, gs),
  Kruskal_Wallis = function(y, g, gf, gs) stats::kruskal.test(y, gf)$p.value
)
TEST_SETS <- list(variance = VARIANCE_TESTS, mean = MEAN_TESTS)
FAMILY_TITLE <- c(variance = "Variance Homogeneity Test Comparison",
                  mean     = "Mean Homogeneity Test Comparison")

# ---- null simulation for one gene -------------------------------------------
# Returns a matrix (alpha x test) of empirical type I error rates.
simulate_gene <- function(gene_expr, group, group_factor, group_sizes, tests,
                          n_sim, alphas, seed) {
  set.seed(seed)
  rejections <- matrix(0, nrow = length(alphas), ncol = length(tests),
                       dimnames = list(as.character(alphas), names(tests)))
  for (s in seq_len(n_sim)) {
    y <- sample(gene_expr)                        # shuffle across spots
    p <- vapply(tests, function(f) {
      tryCatch(f(y, group, group_factor, group_sizes), error = function(e) NA_real_)
    }, numeric(1))
    p[is.na(p)] <- 1                              # failed test = no rejection
    for (a in seq_along(alphas)) rejections[a, ] <- rejections[a, ] + (p < alphas[a])
  }
  rejections / n_sim
}
simulate_gene_worker <- function(x, ...) simulate_gene(x$expr, seed = x$seed, ...)

# ---- run one dataset ---------------------------------------------------------
run_dataset <- function(name, tests) {
  cat("\n===", name, "===\n")
  mat   <- load_real_dataset(name)
  gs    <- REAL_DATASETS[[name]]$banyan_sizes
  group <- domain_index(gs)
  gf    <- factor(group)
  genes <- seq_len(if (is.null(N_GENES)) nrow(mat) else min(N_GENES, nrow(mat)))
  cat("Testing", length(genes), "genes with", N_CORES, "cores\n")

  cl <- makeCluster(N_CORES)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("oneway_F_pvalue", "levene_family_pvalue", "simulate_gene"))

  out <- vector("list", length(genes))
  chunks <- split(genes, ceiling(seq_along(genes) / CHUNK_SIZE))
  done <- 0
  for (chunk in chunks) {
    dense <- as.matrix(mat[chunk, , drop = FALSE])
    X <- lapply(seq_along(chunk), function(j) list(expr = dense[j, ], seed = BASE_SEED + chunk[j]))
    out[chunk] <- parLapply(cl, X, simulate_gene_worker,
                            group = group, group_factor = gf, group_sizes = gs,
                            tests = tests, n_sim = N_SIM, alphas = ALPHAS)
    done <- done + length(chunk)
    cat(sprintf("  %d / %d genes (%.1f%%)\n", done, length(genes), 100 * done / length(genes)))
  }

  gene_names <- rownames(mat)[genes]
  res <- lapply(seq_along(ALPHAS), function(a) {
    df <- as.data.frame(do.call(rbind, lapply(out, function(r) r[a, ])))
    rownames(df) <- gene_names
    df
  })
  names(res) <- sprintf("alpha_%.2f", ALPHAS)
  res$nonzero_ratio <- setNames(as.numeric(rowSums(mat[genes, , drop = FALSE] != 0)) / ncol(mat),
                                gene_names)
  res
}

# ---- figures (same layout as in the original scripts) ------------------------
combine_datasets <- function(results, test_name, alpha_key) {
  do.call(rbind, lapply(DATASETS, function(ds) {
    data.frame(dataset = ds,
               nonzero_ratio   = results[[ds]]$nonzero_ratio,
               empirical_error = results[[ds]][[alpha_key]][[test_name]])
  }))
}

plot_scatter <- function(results, test_name, alpha, family) {
  df <- combine_datasets(results, test_name, sprintf("alpha_%.2f", alpha))
  df$dataset <- factor(df$dataset, levels = DATASETS)
  df$error_percentage <- NA_real_
  for (ds in DATASETS) {                      # how frequent each empirical-error value is
    idx <- df$dataset == ds
    tab <- table(df$empirical_error[idx])
    pct <- as.numeric(tab) / sum(idx) * 100
    names(pct) <- names(tab)
    df$error_percentage[idx] <- pct[as.character(df$empirical_error[idx])]
  }
  ggplot(df, aes(x = nonzero_ratio, y = empirical_error)) +
    geom_point(aes(color = error_percentage), alpha = 0.6, size = 2) +
    scale_color_gradient2(low = "blue", mid = "purple", high = "red",
                          midpoint = stats::median(df$error_percentage),
                          name = "Frequency\n(% of genes)") +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "black", linewidth = 0.8) +
    facet_wrap(~ dataset, ncol = 1, scales = "free") +
    labs(title = paste0(FAMILY_TITLE[[family]], ": ", test_name),
         subtitle = paste0("Nominal Error = ", alpha),
         x = "Nonzero Expression Ratio", y = "Empirical Error Rate") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          legend.position = "right", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 12, face = "bold"))
}

plot_histogram <- function(results, test_name, alpha) {
  df <- combine_datasets(results, test_name, sprintf("alpha_%.2f", alpha))
  df$dataset <- factor(df$dataset, levels = DATASETS)
  stats_df <- do.call(rbind, lapply(DATASETS, function(ds) {
    e <- df$empirical_error[df$dataset == ds]
    data.frame(dataset = factor(ds, levels = DATASETS),
               mean_error = mean(e), median_error = stats::median(e), var_error = stats::var(e))
  }))
  ggplot(df, aes(x = empirical_error, fill = dataset)) +
    geom_histogram(bins = 30, color = "white", alpha = 0.7) +
    geom_vline(xintercept = alpha, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(data = stats_df, aes(xintercept = mean_error), color = "darkgreen", linewidth = 1) +
    geom_vline(data = stats_df, aes(xintercept = median_error), linetype = "dotted",
               color = "purple", linewidth = 1) +
    geom_text(data = stats_df, aes(x = Inf, y = Inf, label = sprintf("Mean: %.4f", mean_error)),
              hjust = 1.05, vjust = 1.5, size = 3, color = "darkgreen", fontface = "bold") +
    geom_text(data = stats_df, aes(x = Inf, y = Inf, label = sprintf("Median: %.4f", median_error)),
              hjust = 1.05, vjust = 3, size = 3, color = "purple", fontface = "bold") +
    geom_text(data = stats_df, aes(x = Inf, y = Inf, label = sprintf("Var: %.6f", var_error)),
              hjust = 1.05, vjust = 4.5, size = 3, color = "black", fontface = "bold") +
    facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(P16 = "steelblue", BR = "coral", D14 = "seagreen")) +
    labs(title = paste0("Empirical Error Distribution Comparison: ", test_name),
         subtitle = paste0("Nominal Error = ", alpha, " (red dashed)"),
         x = "Empirical Error Rate", y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.title = element_text(size = 12), legend.position = "none",
          strip.text = element_text(size = 12, face = "bold"))
}

# ---- main --------------------------------------------------------------------
for (family in TEST_FAMILIES) {
  tests <- TEST_SETS[[family]]
  cat("\n################ family:", family, "(", paste(names(tests), collapse = ", "), ")\n")
  results <- lapply(DATASETS, run_dataset, tests = tests)
  names(results) <- DATASETS
  saveRDS(results, file.path(OUT_DIR, paste0("empirical_error_", family, ".rds")))

  # summary table
  summary_df <- do.call(rbind, lapply(DATASETS, function(ds) {
    do.call(rbind, lapply(ALPHAS, function(a) {
      df <- results[[ds]][[sprintf("alpha_%.2f", a)]]
      data.frame(dataset = ds, alpha = a, test = names(df),
                 mean_error = sapply(df, mean), median_error = sapply(df, stats::median),
                 var_error = sapply(df, stats::var), row.names = NULL)
    }))
  }))
  write.csv(summary_df, file.path(OUT_DIR, paste0("summary_", family, ".csv")), row.names = FALSE)
  print(summary_df)

  # figures
  for (test_name in names(tests)) {
    test_dir <- file.path(OUT_DIR, test_name)
    dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
    for (a in ALPHAS) {
      ggsave(file.path(test_dir, sprintf("scatter_comparison_%s_%.2f.jpg", test_name, a)),
             plot_scatter(results, test_name, a, family), width = 10, height = 12, dpi = 300)
      ggsave(file.path(test_dir, sprintf("comparison_hist_%s_%.2f.jpg", test_name, a)),
             plot_histogram(results, test_name, a), width = 10, height = 10, dpi = 300)
    }
    cat("Figures written to", test_dir, "\n")
  }
}
cat("\nDone.\n")
