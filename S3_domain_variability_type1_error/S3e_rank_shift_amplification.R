# =============================================================================
# S3e_rank_shift_amplification.R
# Supplementary Section S3.1, Figure S7 and the numbers of the paragraphs
# "The number of genes amplifies the shift" and "Why the IQR is robust".
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED (for P16, BR, D14; for the sample VARIANCE and the IQR)
#   For N_SHUFFLE per-gene null permutations (permute_genes(), R/utils_data.R)
#   with the BANYAN domain sizes:
#     s_ik   = statistic (variance or IQR) of gene k in domain i
#     R_ik   = descending mid-rank of s_ik among the m domains (Section S2.1)
#     Rbar_i = mean rank of domain i,   delta_i = Rbar_i - (m+1)/2
#     D      = sum_i delta_i^2,   U = 12 N / (m (m+1)) * D   (Friedman statistic,
#              no tie correction) and the tie-corrected statistic / p-value of
#              stats::friedman.test(), which is what the application uses
#   Reported: the distribution of D and U over the permutations against the
#   value needed to reject, chi2_{0.95; m-1} (D_crit = chi2_crit m(m+1)/(12N));
#   the rejection rates; the mean delta_i of every domain; the number of genes
#   whose statistic differs between domains (all other genes are complete ties
#   and contribute nothing); and the decomposition of delta_i of the largest and
#   smallest domain over bins of the nonzero ratio,
#     delta_i = sum_bins (N_bin / N) delta_i^(bin).
#
# OUTPUTS  results/S3e_rank_shift_amplification/
#   fig_S7_sum_delta2.jpg         D of the variance ranks vs D_crit (Figure S7)
#   rank_shift_summary.csv        medians / 95% ranges, rejection rates, mean delta_i,
#                                 tie factor and number of informative genes
#   rank_shift_per_shuffle.csv    D, U, friedman statistic, p, delta_i per permutation
#   delta_by_bin.csv              delta_i^(bin) of the largest and smallest domain
#
# INPUTS   data/real/count_mat_{p16,br,d14}.rda   (00_data/00a_prepare_real_data.R)
# RUNTIME  ~10-12 h with N_SHUFFLE = 1000 on 8 cores (BR dominates; each worker
#          holds a dense copy of the matrix plus its permuted copy, ~1.5 GB for
#          BR). Raise N_CORES on larger machines.
# USAGE    Rscript S3_domain_variability_type1_error/S3e_rank_shift_amplification.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages({ library(Matrix); library(parallel); library(matrixStats)
                                 library(dplyr); library(ggplot2) })
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)

# ---- parameters --------------------------------------------------------------
DATASETS   <- DATASET_ORDER
N_SHUFFLE  <- 1000
ALPHA      <- 0.05
BIN_BREAKS <- c(0, 0.01, 0.03, 0.10, 0.25, 0.30, 1.0001)
SEED       <- 123
N_CORES    <- default_cores(max_cores = 8)
OUT_DIR    <- "results/S3e_rank_shift_amplification"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

bin_label <- function(x) cut(x, BIN_BREAKS, right = FALSE,
                             labels = paste0("[", head(BIN_BREAKS, -1), ", ", pmin(tail(BIN_BREAKS, -1), 1), ")"))

# descending mid-ranks over the m columns, vectorised over genes (Section S2.1)
desc_midranks <- function(S) {
  m <- ncol(S); R <- matrix(1, nrow(S), m)
  for (i in seq_len(m)) for (l in seq_len(m)) if (l != i)
    R[, i] <- R[, i] + (S[, l] > S[, i]) + 0.5 * (S[, l] == S[, i])
  R
}
block_stats <- function(X, n_i) {
  mu <- rowMeans(X)
  q  <- rowQuantiles(X, probs = c(0.25, 0.75), type = 7L)
  list(var = pmax((rowSums(X^2) - n_i * mu^2) / (n_i - 1), 0), iqr = q[, 2] - q[, 1])
}
# tie-correction factor of the Friedman statistic (friedman.test divides by it):
# C = 1 - sum_k sum_ties (t^3 - t) / (N m (m^2 - 1)); a gene tied in all m domains
# contributes exactly 1/N, so with many all-tied genes C ~ N_informative / N
tie_stats <- function(S) {
  m <- ncol(S); N <- nrow(S)
  res <- apply(S, 1, function(x) { t <- tabulate(match(x, unique(x))); c(sum(t^3 - t), length(t) > 1) })
  list(C = 1 - sum(res[1, ]) / (N * m * (m^2 - 1)), n_inf = sum(res[2, ]))
}

one_shuffle <- function(seed, mat, sizes, blocks, bins) {
  set.seed(seed)
  m <- length(sizes); N <- nrow(mat)
  sh <- permute_genes(mat)
  V <- matrix(0, N, m); Q <- matrix(0, N, m)
  for (i in seq_len(m)) {
    st <- block_stats(sh[, blocks[[i]], drop = FALSE], sizes[i])
    V[, i] <- st$var; Q[, i] <- st$iqr
  }
  out <- list()
  for (stat in c("variance", "IQR")) {
    S <- if (stat == "variance") V else Q
    R <- desc_midranks(S)
    delta <- colMeans(R) - (m + 1) / 2
    D <- sum(delta^2); U <- 12 * N / (m * (m + 1)) * D
    ts <- tie_stats(S)
    ft <- tryCatch(friedman.test(R), error = function(e) list(statistic = NA, p.value = NA))
    bin_delta <- sapply(levels(bins), function(b) { sel <- bins == b
      if (any(sel)) colMeans(R[sel, , drop = FALSE]) - (m + 1) / 2 else rep(NA, m) })
    out[[stat]] <- list(seed = seed, stat = stat, D = D, U = U, tie_factor = ts$C, N_informative = ts$n_inf,
                        friedman_stat = unname(ft$statistic), friedman_p = ft$p.value,
                        delta = delta, bin_delta = bin_delta)
  }
  out
}

shuffle_rows <- list(); bin_rows <- list()
for (name in DATASETS) {
  cat("\n=====", name, "=====\n")
  mat   <- load_real_dataset(name, as_dense = TRUE)
  sizes <- REAL_DATASETS[[name]]$banyan_sizes
  n <- ncol(mat); m <- length(sizes); N <- nrow(mat); blocks <- block_columns(sizes)
  bins <- bin_label(rowSums(mat > 0) / n); N_bin <- table(bins)
  chi2_crit <- qchisq(1 - ALPHA, m - 1); D_crit <- chi2_crit * m * (m + 1) / (12 * N)
  cat(sprintf("N = %d genes, m = %d domains (%s); chi2 critical = %.2f, D_crit = %.5f\n",
              N, m, paste(sizes, collapse = "/"), chi2_crit, D_crit))

  seeds <- SEED + seq_len(N_SHUFFLE)
  res <- if (.Platform$OS.type == "unix") {
    mclapply(seeds, one_shuffle, mat = mat, sizes = sizes, blocks = blocks, bins = bins, mc.cores = N_CORES)
  } else lapply(seeds, one_shuffle, mat = mat, sizes = sizes, blocks = blocks, bins = bins)

  for (r in res) for (stat in c("variance", "IQR")) {
    x <- r[[stat]]
    shuffle_rows[[length(shuffle_rows) + 1]] <- data.frame(
      Dataset = name, statistic = stat, seed = x$seed, sum_delta2 = x$D, U = x$U,
      tie_factor = x$tie_factor, N_informative = x$N_informative,
      friedman_stat = x$friedman_stat, friedman_p = x$friedman_p,
      chi2_crit = chi2_crit, D_crit = D_crit, D_crit_tie_corrected = D_crit * x$tie_factor,
      reject = x$U > chi2_crit, reject_tie_corrected = x$friedman_p < ALPHA,
      t(setNames(x$delta, paste0("delta_domain", seq_len(m)))), row.names = NULL)
    for (i in c(which.max(sizes), which.min(sizes))) for (b in colnames(x$bin_delta))
      bin_rows[[length(bin_rows) + 1]] <- data.frame(
        Dataset = name, statistic = stat, seed = x$seed, domain = i, n_i = sizes[i],
        which = ifelse(i == which.max(sizes), "largest", "smallest"),
        bin = b, N_bin = as.integer(N_bin[b]), weight = as.numeric(N_bin[b]) / N,
        delta_bin = x$bin_delta[i, b], contribution = as.numeric(N_bin[b]) / N * x$bin_delta[i, b],
        row.names = NULL)
  }
}
shuffle_all <- bind_rows(shuffle_rows)      # datasets differ in m: delta_domain5 is NA for m = 4
bin_all     <- bind_rows(bin_rows)
write.csv(shuffle_all, file.path(OUT_DIR, "rank_shift_per_shuffle.csv"), row.names = FALSE)

summary_all <- shuffle_all %>% group_by(Dataset, statistic) %>%
  summarise(N_shuffles = n(), chi2_crit = first(chi2_crit), D_crit = first(D_crit),
            D_median = median(sum_delta2), D_q025 = quantile(sum_delta2, 0.025), D_q975 = quantile(sum_delta2, 0.975),
            U_median = median(U), U_q025 = quantile(U, 0.025), U_q975 = quantile(U, 0.975),
            friedman_median = median(friedman_stat, na.rm = TRUE),
            tie_factor_median = median(tie_factor), N_informative_median = median(N_informative),
            rejection_rate_U = mean(reject), rejection_rate_friedman_p = mean(friedman_p < ALPHA, na.rm = TRUE),
            across(starts_with("delta_domain"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
write.csv(summary_all, file.path(OUT_DIR, "rank_shift_summary.csv"), row.names = FALSE)
cat("\n=== sum_i delta_i^2 and Friedman statistic over permutations ===\n")
print(as.data.frame(summary_all %>% mutate(across(where(is.numeric), ~ signif(.x, 4)))), row.names = FALSE)

bin_summary <- bin_all %>% group_by(Dataset, statistic, which, domain, n_i, bin) %>%
  summarise(N_bin = first(N_bin), weight = first(weight),
            delta_bin = mean(delta_bin), contribution = mean(contribution), .groups = "drop")
write.csv(bin_summary, file.path(OUT_DIR, "delta_by_bin.csv"), row.names = FALSE)
cat("\n=== mean-rank shift by nonzero-ratio bin (largest / smallest domain) ===\n")
print(as.data.frame(bin_summary %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# ---- Figure S7: variance ranks only -------------------------------------------
sv <- shuffle_all %>% filter(statistic == "variance") %>%
  mutate(Dataset = factor(Dataset, levels = DATASETS))
crit <- sv %>% distinct(Dataset, D_crit)
p <- ggplot(sv, aes(x = Dataset, y = sum_delta2, colour = Dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  geom_hline(data = crit, aes(yintercept = D_crit), linetype = "dashed") +
  facet_wrap(~ Dataset, scales = "free") +
  scale_y_log10() +
  labs(x = NULL, y = expression(sum((bar(R)[i.] - (m + 1) / 2)^2, i) ~ "(log scale)")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(OUT_DIR, "fig_S7_sum_delta2.jpg"), p, width = 11, height = 5, dpi = 300)
cat("Saved to", OUT_DIR, "\n")
