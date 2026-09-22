# =============================================================================
# S3b_sparse_gene_variance_quantities.R
# Supplementary Section S3.1: the quantities quoted in the text that need no
# simulation (data structure and closed-form hypergeometric expressions), and
# Table S3.
# -----------------------------------------------------------------------------
# WHAT IS COMPUTED (for P16, BR, D14 with the BANYAN domain sizes n_i)
#   (a) Single-count structure of sparse genes ("for genes expressed in fewer
#       than 3% of the spots ... 91-95% of the non-zero values are single
#       counts"): for every gene, the fraction of its non-zero values that sit
#       at its lowest non-zero level v_k (count 1 after SCTransform). Reported
#       for the genes below SPARSE_THRESHOLD and by 1% bins of the nonzero ratio.
#   (b) Exact zero probability P(X_ik = 0) = prod_{l=0}^{c_k-1} (n-n_i-l)/(n-l),
#       X_ik ~ Hypergeom(n, c_k, n_i), averaged over the genes of each sparsity
#       bin for the smallest and the largest domain ("0.27 versus 0.06 for the
#       311- and 609-spot domains of P16 and genes expressed in fewer than
#       0.5% of the spots").
#   (c) Table S3: dataset-dependent coefficients of Var(p_hat) and Skew(p_hat),
#         var coef  = (n/n_i - 1) / (n - 1)
#         skew coef = sqrt(n-1)/(n-2) * (1 - 2 n_i/n) / sqrt((n_i/n)(1 - n_i/n))
#       for the smallest and the largest domain.
#   (d) The value of sum_i (Rbar_i - (m+1)/2)^2 needed to reject the Friedman
#       test at alpha = 0.05: chi2_{0.95; m-1} m(m+1) / (12 N).
#
# OUTPUTS  results/S3b_sparse_gene_variance_quantities/
#   single_count_fraction.csv        (a) below SPARSE_THRESHOLD, per dataset
#   single_count_fraction_by_bin.csv (a) by 1% bins of the nonzero ratio
#   zero_probability_by_bin.csv      (b)
#   table_S3_coefficients.csv        (c)
#   rejection_threshold.csv          (d)
#
# INPUTS   data/real/count_mat_{p16,br,d14}.rda   (00_data/00a_prepare_real_data.R)
# RUNTIME  under a minute.
# USAGE    Rscript S3_domain_variability_type1_error/S3b_sparse_gene_variance_quantities.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages({ library(Matrix); library(dplyr) })
source("R/utils_data.R")

# ---- parameters --------------------------------------------------------------
DATASETS         <- DATASET_ORDER
SPARSE_THRESHOLD <- 0.03                              # "genes expressed in fewer than 3% of the spots"
BIN_BREAKS       <- c(seq(0, 0.30, by = 0.01), 0.40, 0.50, 1.0001)
ZERO_PROB_BINS   <- c(0, 0.005, 0.01, 0.03)           # bins for (b)
ALPHA            <- 0.05
OUT_DIR          <- "results/S3b_sparse_gene_variance_quantities"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

bin_label <- function(x, breaks) cut(x, breaks, right = TRUE, include.lowest = FALSE,
                                     labels = paste0("(", head(breaks, -1), ", ", pmin(tail(breaks, -1), 1), "]"))
p_zero_exact <- function(c_k, n_i, n) dhyper(0, m = c_k, n = n - c_k, k = n_i)   # = prod (n-n_i-l)/(n-l)

single_rows <- list(); single_bin_rows <- list(); zero_rows <- list(); coef_rows <- list(); thr_rows <- list()
for (name in DATASETS) {
  cat("\n=====", name, "=====\n")
  mat   <- load_real_dataset(name, as_dense = TRUE)
  sizes <- REAL_DATASETS[[name]]$banyan_sizes
  n <- ncol(mat); m <- length(sizes); N <- nrow(mat)
  c_k <- rowSums(mat > 0); rho <- c_k / n
  iS <- which.min(sizes); iL <- which.max(sizes); nS <- sizes[iS]; nL <- sizes[iL]

  # ---- (a) single-count structure ------------------------------------------
  expressed <- which(c_k > 0)
  lvl <- t(vapply(expressed, function(k) {
    v <- mat[k, mat[k, ] > 0]; v1 <- min(v)
    c(n_nonzero = length(v), n_at_level1 = sum(abs(v - v1) < 1e-8))
  }, numeric(2)))
  gene_tab <- data.frame(gene = rownames(mat)[expressed], c_k = c_k[expressed], nz_ratio = rho[expressed],
                         n_nonzero = lvl[, "n_nonzero"], n_at_level1 = lvl[, "n_at_level1"],
                         frac_level1 = lvl[, "n_at_level1"] / lvl[, "n_nonzero"])
  sp <- gene_tab[gene_tab$nz_ratio < SPARSE_THRESHOLD, ]
  single_rows[[name]] <- data.frame(
    Dataset = name, threshold = SPARSE_THRESHOLD, n_genes = nrow(sp), share_of_genes = nrow(sp) / N,
    frac_single_count_pooled = sum(sp$n_at_level1) / sum(sp$n_nonzero),   # "of the non-zero values"
    frac_single_count_gene_mean = mean(sp$frac_level1))
  single_bin_rows[[name]] <- gene_tab %>%
    mutate(bin = bin_label(nz_ratio, BIN_BREAKS)) %>% group_by(bin) %>%
    summarise(n_genes = n(), frac_single_count_pooled = sum(n_at_level1) / sum(n_nonzero),
              frac_single_count_gene_mean = mean(frac_level1), .groups = "drop") %>%
    mutate(Dataset = name, .before = 1)

  # ---- (b) exact zero probability -------------------------------------------
  zb <- data.frame(c_k = c_k[expressed], nz_ratio = rho[expressed]) %>%
    filter(nz_ratio <= max(ZERO_PROB_BINS)) %>%
    mutate(bin = bin_label(nz_ratio, ZERO_PROB_BINS),
           P0_smallest = p_zero_exact(c_k, nS, n), P0_largest = p_zero_exact(c_k, nL, n))
  zero_rows[[name]] <- bind_rows(
    zb %>% group_by(bin) %>% summarise(n_genes = n(), P0_smallest = mean(P0_smallest),
                                       P0_largest = mean(P0_largest), .groups = "drop"),
    zb %>% filter(nz_ratio <= 0.01) %>%
      summarise(bin = "(0, 0.01]", n_genes = n(), P0_smallest = mean(P0_smallest), P0_largest = mean(P0_largest))) %>%
    mutate(Dataset = name, n_smallest = nS, n_largest = nL, .before = 1)

  # ---- (c) Table S3 coefficients ----------------------------------------------
  var_coef  <- function(n_i) (n / n_i - 1) / (n - 1)
  skew_coef <- function(n_i) sqrt(n - 1) / (n - 2) * (1 - 2 * n_i / n) / sqrt((n_i / n) * (1 - n_i / n))
  coef_rows[[name]] <- data.frame(
    Dataset = name, n = n, n_smallest = nS, n_largest = nL,
    var_coef_smallest = var_coef(nS), var_coef_largest = var_coef(nL),
    skew_coef_smallest = skew_coef(nS), skew_coef_largest = skew_coef(nL))

  # ---- (d) rejection threshold ------------------------------------------------
  chi2_crit <- qchisq(1 - ALPHA, m - 1)
  thr_rows[[name]] <- data.frame(Dataset = name, N_genes = N, m = m, alpha = ALPHA,
                                 chi2_crit = chi2_crit, sum_delta2_crit = chi2_crit * m * (m + 1) / (12 * N))
}

single_all <- bind_rows(single_rows); single_bin_all <- bind_rows(single_bin_rows)
zero_all <- bind_rows(zero_rows); coef_all <- bind_rows(coef_rows); thr_all <- bind_rows(thr_rows)
write.csv(single_all,     file.path(OUT_DIR, "single_count_fraction.csv"),        row.names = FALSE)
write.csv(single_bin_all, file.path(OUT_DIR, "single_count_fraction_by_bin.csv"), row.names = FALSE)
write.csv(zero_all,       file.path(OUT_DIR, "zero_probability_by_bin.csv"),      row.names = FALSE)
write.csv(coef_all,       file.path(OUT_DIR, "table_S3_coefficients.csv"),        row.names = FALSE)
write.csv(thr_all,        file.path(OUT_DIR, "rejection_threshold.csv"),          row.names = FALSE)

cat(sprintf("\n=== (a) genes expressed in fewer than %.0f%% of the spots ===\n", 100 * SPARSE_THRESHOLD))
print(as.data.frame(single_all %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\n=== (b) exact P(X_ik = 0), mean over genes of each bin ===\n")
print(as.data.frame(zero_all %>% mutate(across(c(P0_smallest, P0_largest), ~ round(.x, 3)))), row.names = FALSE)
cat("\n=== (c) Table S3 ===\n")
print(as.data.frame(coef_all %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
cat("\n=== (d) rejection thresholds ===\n")
print(as.data.frame(thr_all %>% mutate(across(where(is.numeric), ~ signif(.x, 4)))), row.names = FALSE)
cat("Saved to", OUT_DIR, "\n")
