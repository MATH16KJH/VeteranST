# =============================================================================
# S3a_type1_error_domain_variability.R
# Supplementary Section S3: type I error of the spatial-domain variability
# tests -- Table S1 (equal domain sizes) and Table S2 (BANYAN domain sizes)
# -----------------------------------------------------------------------------
# WHAT IS SIMULATED
#   In each run the expression of every gene is shuffled across spots and the
#   first sum(group sizes) shuffled spots are split into consecutive blocks
#   playing the role of spatial domains, so that all domains are drawn from the
#   same population (null hypothesis). For every gene the per-domain
#   variability (IQR or variance) is computed, giving a genes x domains matrix,
#   and each test below is applied once per run. The empirical type I error is
#   the fraction of N_SIM runs rejected at alpha.
#
#   Scenario "equal"  (Table S1): spots are divided evenly among the
#     BANYAN-determined number of domains (P16: 4 x 412, BR: 5 x 759,
#     D14: 4 x 491); spots beyond the last full group are dropped.
#   Scenario "banyan" (Table S2): the actual BANYAN domain sizes
#     (P16 311/609/321/409, BR 861/665/942/741/589, D14 754/392/240/581).
#   Both scenarios are run by default; the number of runs and the seed scheme
#   of each are set in SCENARIO_SETTINGS below and reproduce the tables of the
#   supplement exactly (see SEED SCHEMES).
#
# TESTS (rows of Tables S1 and S2)
#   Friedman Test on IQR                  PROPOSED (genes = blocks, domains = groups)
#   Friedman Test on Variance
#   ANOVA                                 one-way ANOVA on the per-domain variances
#   ANOVA (Variance-Rank Based)           one-way ANOVA on within-gene variance ranks
#   Kruskal-Wallis Test (Variance-Rank Based)  Kruskal-Wallis on the per-domain
#                                         variances (a rank test treating genes as
#                                         independent samples)
#
# SEED SCHEMES
#   "submitted" (Table S1, scenario "equal", 1,000 runs) replays the two
#     development runs that produced Table S1, so the table is reproduced bit
#     for bit given the committed input matrices:
#       * IQR stream      -- run s is shuffled with seed 123 + s, for every
#                            dataset and both alpha levels (one stream);
#       * variance stream -- run s of block k is shuffled with seed
#                            12345 + 10000 k + s, where k numbers the
#                            (dataset, alpha) combinations in the order
#                            P16/0.01, P16/0.05, D14/0.01, D14/0.05, BR/0.01,
#                            BR/0.05 (k = 1..6 for "equal", 7..12 for "banyan").
#     The four variance-based tests are applied to the same shuffled matrix of
#     a run; the IQR test uses its own shuffles. 3 x n_sim shuffles per dataset.
#   "unified" (Table S2, scenario "banyan", 1,000 runs) applies all five tests
#     to one shuffled matrix per run (seed 123 + s), i.e. n_sim shuffles per
#     dataset.
#   Either scheme can be used for either scenario; the two differ only within
#   Monte Carlo error (SE = sqrt(a(1-a)/n_sim), about 0.007 at alpha = 0.05
#   for 1,000 runs).
#
# OUTPUTS  results/S3_domain_variability_type1_error/
#   table_S1_type1_error_equal.csv    Table S1
#   table_S2_type1_error_banyan.csv   Table S2
#   details_<scenario>.rds            per-run rejection indicators and seeds
#
# INPUTS   data/real/count_mat_{p16,br,d14}.rda   (00_data/00a_prepare_real_data.R)
# RUNTIME  Table S1 (3,000 shuffles x 3 datasets) about 1 h with 20 cores;
#          Table S2 (1,000 shuffles x 3 datasets) about 20 min. Each worker
#          holds a dense copy of the matrix (BR ~ 0.6 GB); reduce N_CORES on
#          small machines.
# USAGE    Rscript S3_domain_variability_type1_error/S3a_type1_error_domain_variability.R
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages({
  library(Matrix)
  library(parallel)
})
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)

# ---- parameters --------------------------------------------------------------
# One entry per scenario: table of the supplement, seed scheme and number of
# runs (the settings of the runs reported in the supplement).
SCENARIO_SETTINGS <- list(
  equal  = list(table = "S1", seed_scheme = "submitted", n_sim = 1000),
  banyan = list(table = "S2", seed_scheme = "unified",   n_sim = 1000)
)
DATASETS       <- DATASET_ORDER            # c("P16", "BR", "D14") -- output order only
ALPHAS         <- c(0.01, 0.05)
IQR_BASE_SEED  <- 123                      # IQR stream:      seed = 123 + s
VAR_BASE_SEED  <- 12345                    # variance stream: seed = 12345 + 10000 k + s
N_CORES        <- default_cores(max_cores = 20)
OUT_DIR        <- "results/S3_domain_variability_type1_error"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Block index k of the variance stream in the original run (dataset order
# p16, d14, br; alphas 0.01, 0.05; "equal" blocks first, then "banyan").
VARIANCE_SEED_BLOCK <- function(dataset, alpha, scenario) {
  ds_pos <- match(dataset, c("P16", "D14", "BR"))
  a_pos  <- match(alpha, c(0.01, 0.05))
  if (is.na(ds_pos) || is.na(a_pos)) stop("No original seed block for ", dataset, " / ", alpha)
  (if (scenario == "equal") 0 else 6) + 2 * (ds_pos - 1) + a_pos
}

# ---- group definition --------------------------------------------------------
make_group_sizes <- function(banyan_sizes, n_spots, scenario) {
  if (scenario == "equal") {
    n_groups <- length(banyan_sizes)
    rep(floor(n_spots / n_groups), n_groups)
  } else if (scenario == "banyan") {
    banyan_sizes
  } else stop("Unknown scenario: ", scenario)
}

# ---- one shuffle -------------------------------------------------------------
# Identical to the development code: set.seed(seed), then row i of the dense
# matrix is permuted with sample(n_spots) in row order.
shuffle_genes <- function(mat_dense, seed) {
  set.seed(seed)
  n_spots <- ncol(mat_dense)
  out <- matrix(0, nrow(mat_dense), n_spots)
  for (i in seq_len(nrow(mat_dense))) out[i, ] <- mat_dense[i, sample(n_spots)]
  out
}

# genes x groups matrix of a per-domain statistic (stats::var or stats::IQR).
group_statistic <- function(shuffled, group_sizes, fun) {
  idx_list <- split(seq_len(sum(group_sizes)), domain_index(group_sizes))
  vapply(idx_list, function(idx) apply(shuffled[, idx, drop = FALSE], 1, fun),
         numeric(nrow(shuffled)))
}

# ---- tests (p-values) --------------------------------------------------------
long_format <- function(m) {
  data.frame(value = as.vector(m), group = factor(rep(seq_len(ncol(m)), each = nrow(m))))
}
safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
friedman_p     <- function(m) safe(stats::friedman.test(m)$p.value)     # rows = blocks (genes)
anova_p        <- function(m) safe({ df <- long_format(m)
                                     summary(stats::aov(value ~ group, data = df))[[1]][["Pr(>F)"]][1] })
anova_rank_p   <- function(m) anova_p(t(apply(m, 1, function(x) rank(-x, ties.method = "average"))))
kruskal_p      <- function(m) safe({ df <- long_format(m)
                                     stats::kruskal.test(value ~ group, data = df)$p.value })

TEST_NAMES <- c("Friedman Test on IQR", "Friedman Test on Variance", "ANOVA",
                "ANOVA (Variance-Rank Based)", "Kruskal-Wallis Test (Variance-Rank Based)")

# p-values of the four variance-based tests for one shuffle
variance_tests <- function(seed, mat_dense, group_sizes) {
  v <- group_statistic(shuffle_genes(mat_dense, seed), group_sizes, stats::var)
  c(friedman_p(v), anova_p(v), anova_rank_p(v), kruskal_p(v))
}
# p-value of the Friedman test on IQR for one shuffle
iqr_test <- function(seed, mat_dense, group_sizes) {
  friedman_p(group_statistic(shuffle_genes(mat_dense, seed), group_sizes, stats::IQR))
}
# all five tests on one shuffle (unified scheme)
all_tests <- function(seed, mat_dense, group_sizes) {
  sh <- shuffle_genes(mat_dense, seed)
  v  <- group_statistic(sh, group_sizes, stats::var)
  q  <- group_statistic(sh, group_sizes, stats::IQR)
  c(friedman_p(q), friedman_p(v), anova_p(v), anova_rank_p(v), kruskal_p(v))
}

run_parallel <- function(seeds, FUN, mat_dense, group_sizes) {
  if (.Platform$OS.type == "unix") {
    mclapply(seeds, FUN, mat_dense = mat_dense, group_sizes = group_sizes, mc.cores = N_CORES)
  } else {
    cl <- makeCluster(N_CORES); on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("shuffle_genes", "group_statistic", "domain_index", "long_format", "safe",
                        "friedman_p", "anova_p", "anova_rank_p", "kruskal_p",
                        "variance_tests", "iqr_test", "all_tests"))
    parLapply(cl, seeds, FUN, mat_dense = mat_dense, group_sizes = group_sizes)
  }
}

# ---- main --------------------------------------------------------------------
for (GROUP_SIZE_SCENARIO in names(SCENARIO_SETTINGS)) {
SEED_SCHEME <- SCENARIO_SETTINGS[[GROUP_SIZE_SCENARIO]]$seed_scheme
N_SIM       <- SCENARIO_SETTINGS[[GROUP_SIZE_SCENARIO]]$n_sim
TABLE       <- SCENARIO_SETTINGS[[GROUP_SIZE_SCENARIO]]$table
all_rows <- list()
details  <- list()
for (name in DATASETS) {
  cat("\n===", name, "| scenario:", GROUP_SIZE_SCENARIO, "| runs:", N_SIM, "| seed scheme:", SEED_SCHEME, "===\n")
  mat_dense   <- load_real_dataset(name, as_dense = TRUE)
  group_sizes <- make_group_sizes(REAL_DATASETS[[name]]$banyan_sizes, ncol(mat_dense),
                                  GROUP_SIZE_SCENARIO)
  cat("Genes:", nrow(mat_dense), "| spots:", ncol(mat_dense),
      "| group sizes:", paste(group_sizes, collapse = ", "), "\n")
  t0 <- Sys.time()

  # rejection indicators: test x alpha x run
  rej <- array(NA_integer_, dim = c(length(TEST_NAMES), length(ALPHAS), N_SIM),
               dimnames = list(TEST_NAMES, as.character(ALPHAS), NULL))
  seeds_used <- list()

  if (SEED_SCHEME == "submitted") {
    seeds_iqr <- IQR_BASE_SEED + seq_len(N_SIM)
    p_iqr <- unlist(run_parallel(seeds_iqr, iqr_test, mat_dense, group_sizes))
    for (a in seq_along(ALPHAS)) rej[1, a, ] <- as.integer(p_iqr < ALPHAS[a])
    seeds_used$iqr <- seeds_iqr
    cat("  IQR stream done\n")
    for (a in seq_along(ALPHAS)) {
      k <- VARIANCE_SEED_BLOCK(name, ALPHAS[a], GROUP_SIZE_SCENARIO)
      seeds_var <- VAR_BASE_SEED + 10000 * k + seq_len(N_SIM)
      p_var <- simplify2array(run_parallel(seeds_var, variance_tests, mat_dense, group_sizes))
      rej[2:5, a, ] <- (p_var < ALPHAS[a]) * 1L
      seeds_used[[sprintf("variance_alpha_%.2f", ALPHAS[a])]] <- seeds_var
      cat(sprintf("  variance stream, alpha = %.2f (seed block %d) done\n", ALPHAS[a], k))
    }
  } else if (SEED_SCHEME == "unified") {
    seeds <- IQR_BASE_SEED + seq_len(N_SIM)
    p_all <- simplify2array(run_parallel(seeds, all_tests, mat_dense, group_sizes))   # 5 x N_SIM
    for (a in seq_along(ALPHAS)) rej[, a, ] <- (p_all < ALPHAS[a]) * 1L
    seeds_used$all <- seeds
  } else stop("Unknown SEED_SCHEME: ", SEED_SCHEME)

  cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  details[[name]] <- list(rejections = rej, seeds = seeds_used, group_sizes = group_sizes)
  for (a in seq_along(ALPHAS)) {
    all_rows[[length(all_rows) + 1]] <- data.frame(
      Method  = TEST_NAMES,
      Dataset = name,
      alpha   = ALPHAS[a],
      empirical_type1_error = apply(rej[, a, , drop = FALSE], 1, mean, na.rm = TRUE),
      n_valid_runs = apply(rej[, a, , drop = FALSE], 1, function(x) sum(!is.na(x))),
      row.names = NULL)
  }
}

table_long <- do.call(rbind, all_rows)
table_wide <- stats::reshape(table_long[, c("Method", "Dataset", "alpha", "empirical_type1_error")],
                             idvar = c("Method", "Dataset"), timevar = "alpha", direction = "wide")
names(table_wide) <- sub("empirical_type1_error\\.", "alpha_", names(table_wide))
table_wide <- table_wide[order(match(table_wide$Method, TEST_NAMES),
                               match(table_wide$Dataset, DATASETS)), ]

out_csv <- file.path(OUT_DIR, sprintf("table_%s_type1_error_%s.csv", TABLE, GROUP_SIZE_SCENARIO))
write.csv(table_wide, out_csv, row.names = FALSE)
saveRDS(list(details = details, seed_scheme = SEED_SCHEME, n_sim = N_SIM, alphas = ALPHAS),
        file.path(OUT_DIR, sprintf("details_%s.rds", GROUP_SIZE_SCENARIO)))

cat("\n=== Empirical type I error (Table ", TABLE, ": ", GROUP_SIZE_SCENARIO, " group sizes, ",
    N_SIM, " runs, ", SEED_SCHEME, " seeds) ===\n", sep = "")
print(table_wide, row.names = FALSE)
cat("Saved:", out_csv, "\n")
}   # end of scenario loop
