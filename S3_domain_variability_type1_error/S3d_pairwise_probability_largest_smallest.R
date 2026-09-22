# =============================================================================
# S3d_pairwise_probability_largest_smallest.R
# Supplementary Section S3.1, Figure S6: P(s2_L > s2_S) for the largest (L)
# and smallest (S) BANYAN domain under the null permutation, by sparsity of
# the gene.
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED (for P16, BR, D14)
#   N_SHUFFLE per-gene null permutations (permute_genes(), R/utils_data.R);
#   the domains are consecutive column blocks of the permuted matrix with the
#   BANYAN sizes. For every gene the sample variances s2_L and s2_S of the two
#   domains are compared in each permutation, giving per gene
#     P(s2_L > s2_S), P(s2_L = s2_S), P(s2_L < s2_S).
#   Under exchangeability (equal sizes) P(>) = P(<); with unequal sizes the
#   larger domain wins more often than it loses (Section S3.1). Summarised by
#   bins of the gene's nonzero ratio.
#
# OUTPUTS  results/S3d_pairwise_probability_largest_smallest/
#   fig_S6_pairwise_probability.jpg
#   pairwise_probability_by_bin.csv, pairwise_probability_by_gene.csv
#
# INPUTS   data/real/count_mat_{p16,br,d14}.rda   (00_data/00a_prepare_real_data.R)
# RUNTIME  ~40 min (P16, D14) to ~2-3 h (BR) with N_SHUFFLE = 1000, single core.
# USAGE    Rscript S3_domain_variability_type1_error/S3d_pairwise_probability_largest_smallest.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages({ library(Matrix); library(dplyr); library(ggplot2) })
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)

# ---- parameters --------------------------------------------------------------
DATASETS   <- DATASET_ORDER
N_SHUFFLE  <- 1000
BIN_BREAKS <- c(0, 0.01, 0.03, 0.05, 0.10, 0.25, 0.50, 1.0001)
SEED       <- 123
OUT_DIR    <- "results/S3d_pairwise_probability_largest_smallest"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

bin_label <- function(x) cut(x, BIN_BREAKS, right = FALSE,
                             labels = paste0("[", head(BIN_BREAKS, -1), ", ", pmin(tail(BIN_BREAKS, -1), 1), ")"))
block_var <- function(X, n_i) { mu <- rowMeans(X); pmax((rowSums(X^2) - n_i * mu^2) / (n_i - 1), 0) }

gene_rows <- list()
for (name in DATASETS) {
  cat("\n=====", name, "=====\n")
  mat   <- load_real_dataset(name, as_dense = TRUE)
  sizes <- REAL_DATASETS[[name]]$banyan_sizes
  n <- ncol(mat); G <- nrow(mat); blocks <- block_columns(sizes)
  iL <- which.max(sizes); iS <- which.min(sizes); nL <- sizes[iL]; nS <- sizes[iS]
  nz_ratio <- rowSums(mat > 0) / n
  cat(sprintf("%d genes, %d spots; largest domain %d (%d spots), smallest domain %d (%d spots)\n",
              G, n, iL, nL, iS, nS))

  gt <- eq <- lt <- numeric(G)
  set.seed(SEED)
  for (s in seq_len(N_SHUFFLE)) {
    sh <- permute_genes(mat)
    vL <- block_var(sh[, blocks[[iL]], drop = FALSE], nL)
    vS <- block_var(sh[, blocks[[iS]], drop = FALSE], nS)
    gt <- gt + (vL > vS); eq <- eq + (vL == vS); lt <- lt + (vL < vS)
    if (s %% 100 == 0) cat("  permutation", s, "of", N_SHUFFLE, "\n")
  }
  gene_rows[[name]] <- data.frame(Dataset = name, gene = rownames(mat), nz_ratio = nz_ratio,
                                  bin = bin_label(nz_ratio), n_L = nL, n_S = nS,
                                  P_gt = gt / N_SHUFFLE, P_eq = eq / N_SHUFFLE, P_lt = lt / N_SHUFFLE,
                                  row.names = NULL)
}
gene_all <- bind_rows(gene_rows)
write.csv(gene_all, file.path(OUT_DIR, "pairwise_probability_by_gene.csv"), row.names = FALSE)

bin_all <- bind_rows(
  gene_all %>% group_by(Dataset, bin, n_L, n_S) %>%
    summarise(n_genes = n(), P_gt = mean(P_gt), P_eq = mean(P_eq), P_lt = mean(P_lt), .groups = "drop") %>%
    mutate(bin = as.character(bin)),
  gene_all %>% group_by(Dataset, n_L, n_S) %>%
    summarise(bin = "all genes", n_genes = n(), P_gt = mean(P_gt), P_eq = mean(P_eq), P_lt = mean(P_lt),
              .groups = "drop")) %>%
  mutate(Dataset = factor(Dataset, levels = DATASETS)) %>% arrange(Dataset)
write.csv(bin_all, file.path(OUT_DIR, "pairwise_probability_by_bin.csv"), row.names = FALSE)

cat("\n=== P(s2_largest > s2_smallest) by nonzero-ratio bin ===\n")
print(as.data.frame(bin_all %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# ---- Figure S6 ---------------------------------------------------------------
plot_df <- bin_all %>% filter(bin != "all genes") %>%
  mutate(bin = factor(bin, levels = levels(gene_all$bin)),
         panel = sprintf("%s * ':' ~ n[L] == %d * ',' ~ n[S] == %d", Dataset, n_L, n_S))
plot_df$panel <- factor(plot_df$panel, levels = unique(plot_df$panel[order(plot_df$Dataset)]))
p <- ggplot(plot_df, aes(x = bin, y = P_gt, colour = Dataset, group = Dataset)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~ panel, labeller = label_parsed) +
  labs(x = "Nonzero ratio of the gene", y = expression(P(s[L]^2 > s[S]^2)), colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        panel.grid.minor = element_blank())
ggsave(file.path(OUT_DIR, "fig_S6_pairwise_probability.jpg"), p, width = 11, height = 4.5, dpi = 300)
cat("Saved to", OUT_DIR, "\n")
