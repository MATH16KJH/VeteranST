# =============================================================================
# S3c_rank_agreement_largest_smallest.R
# Supplementary Section S3.1, Figure S5: do the hypergeometric counts of
# expressing spots alone decide which of the largest (L) and smallest (S)
# BANYAN domain has the larger sample variance?
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED (for P16, BR, D14)
#   N_SHUFFLE per-gene null permutations (permute_genes(), R/utils_data.R);
#   the domains are consecutive column blocks of the permuted matrix with the
#   BANYAN sizes. For every gene k and permutation:
#     observed  : sign(s2_Lk - s2_Sk)          sample variances of the two domains
#     predicted : sign(X_Lk / n_L - X_Sk / n_S) numbers of expressing spots only
#   agreement = fraction of permutations in which the two signs coincide (a tie
#   in one and not in the other counts as disagreement). Summarised by 1% bins
#   of the gene's nonzero ratio (up to 30%, then 30-40% and 40-50%).
#
# OUTPUTS  results/S3c_rank_agreement_largest_smallest/
#   fig_S5_rank_agreement.jpg
#   rank_agreement_by_bin.csv, rank_agreement_by_gene.csv
#
# INPUTS   data/real/count_mat_{p16,br,d14}.rda   (00_data/00a_prepare_real_data.R)
# RUNTIME  ~5 min (P16, D14) to ~15 min (BR) with N_SHUFFLE = 100, single core.
# USAGE    Rscript S3_domain_variability_type1_error/S3c_rank_agreement_largest_smallest.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages({ library(Matrix); library(dplyr); library(ggplot2) })
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)

# ---- parameters --------------------------------------------------------------
DATASETS     <- DATASET_ORDER
N_SHUFFLE    <- 100
MAX_NZ_RATIO <- 0.50
BIN_BREAKS   <- c(seq(0, 0.30, by = 0.01), 0.40, 0.50)
SEED         <- 123
REF_LINE     <- 0.9                          # dotted reference line of the figure
OUT_DIR      <- "results/S3c_rank_agreement_largest_smallest"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

bin_label <- function(x) cut(x, BIN_BREAKS, right = TRUE, include.lowest = FALSE,
                             labels = paste0("(", head(BIN_BREAKS, -1), ", ", tail(BIN_BREAKS, -1), "]"))
bin_upper <- function(b) as.numeric(sub(".*, (.*)\\]", "\\1", as.character(b)))
block_var <- function(X, n_i) { mu <- rowMeans(X); pmax((rowSums(X^2) - n_i * mu^2) / (n_i - 1), 0) }

gene_rows <- list()
for (name in DATASETS) {
  cat("\n=====", name, "=====\n")
  mat   <- load_real_dataset(name, as_dense = TRUE)
  sizes <- REAL_DATASETS[[name]]$banyan_sizes
  n <- ncol(mat); blocks <- block_columns(sizes)
  iL <- which.max(sizes); iS <- which.min(sizes); nL <- sizes[iL]; nS <- sizes[iS]
  nz_ratio <- rowSums(mat > 0) / n
  keep <- which(nz_ratio > 0 & nz_ratio <= MAX_NZ_RATIO)
  X <- mat[keep, , drop = FALSE]; G <- length(keep)
  cat(sprintf("%d genes; largest domain %d (%d spots), smallest domain %d (%d spots)\n", G, iL, nL, iS, nS))

  agree <- numeric(G)
  set.seed(SEED)
  for (s in seq_len(N_SHUFFLE)) {
    sh <- permute_genes(X)
    XL <- sh[, blocks[[iL]], drop = FALSE]; XS <- sh[, blocks[[iS]], drop = FALSE]
    sg_obs  <- sign(block_var(XL, nL) - block_var(XS, nS))
    sg_pred <- sign(rowSums(XL > 0) / nL - rowSums(XS > 0) / nS)
    agree <- agree + (sg_obs == sg_pred)
    if (s %% 25 == 0) cat("  permutation", s, "of", N_SHUFFLE, "\n")
  }
  gene_rows[[name]] <- data.frame(Dataset = name, gene = rownames(mat)[keep], nz_ratio = nz_ratio[keep],
                                  bin = bin_label(nz_ratio[keep]), n_L = nL, n_S = nS,
                                  agreement = agree / N_SHUFFLE, row.names = NULL)
}
gene_all <- bind_rows(gene_rows)
write.csv(gene_all, file.path(OUT_DIR, "rank_agreement_by_gene.csv"), row.names = FALSE)

bin_all <- gene_all %>% group_by(Dataset, bin) %>%
  summarise(n_genes = n(), nz_upper = bin_upper(first(bin)), agreement_mean = mean(agreement),
            agreement_q10 = quantile(agreement, 0.10), .groups = "drop") %>%
  mutate(Dataset = factor(Dataset, levels = DATASETS)) %>% arrange(Dataset, nz_upper)
write.csv(bin_all, file.path(OUT_DIR, "rank_agreement_by_bin.csv"), row.names = FALSE)

cat("\n=== agreement between sign(s2_L - s2_S) and sign(X_L/n_L - X_S/n_S), by nonzero-ratio bin ===\n")
print(as.data.frame(bin_all %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# ---- Figure S5 ---------------------------------------------------------------
p <- ggplot(bin_all, aes(x = nz_upper, y = agreement_mean, colour = Dataset, group = Dataset)) +
  geom_hline(yintercept = REF_LINE, linetype = "dotted", colour = "black") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50), labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
  labs(x = "Nonzero ratio of the gene (bin upper edge)", y = "Agreement", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(OUT_DIR, "fig_S5_rank_agreement.jpg"), p, width = 9, height = 5, dpi = 300)
cat("Saved to", OUT_DIR, "\n")
