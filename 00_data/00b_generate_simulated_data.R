# =============================================================================
# 00b_generate_simulated_data.R -- ground-truth benchmark data for Section S4
# -----------------------------------------------------------------------------
# PURPOSE
#   Generate the synthetic spatial transcriptomics data used to evaluate the
#   spatial gene co-expression (clustering) component (Supplementary Section S4:
#   Figures S8-S10 and Table S4; main-text Figure 2).
#
# DESIGN (Supplementary Section S4)
#   * Reference: scDesign3 (Song et al., 2024) fitted to the 10x Visium mouse
#     brain example dataset distributed with scDesign3 (figshare 40582019).
#   * N_GENES genes are taken from the reference and N_REPLICATES independent
#     scDesign3 realisations are generated for each gene, each with its own seed.
#     Every replicate keeps the spatial expression pattern of its parent gene,
#     so the ground-truth module of a simulated gene is its parent gene:
#     N_GENES modules x N_REPLICATES genes (10 x 10 = 100 genes by default).
#   * Simulated genes are named "<gene>_<replicate>"; R/utils_evaluation.R
#     recovers the truth from that suffix.
#   * Expression is log1p-transformed; coordinates are the reference spot
#     positions (identical across replicates).
#
# OUTPUT (inputs of every S4 script; the files used for the manuscript are
#         COMMITTED in data/simulated/, so running this script is optional)
#   data/simulated/count_mat_simul.rds     100 genes x 2,696 spots, log1p counts
#   data/simulated/coordinates_simul.rds   spots x 2  (spatial1, spatial2; image pixels)
#   data/simulated/ground_truth.rds        Gene, Gene_Base_Name, True_Cluster
#   The reference object is the scDesign3 example SingleCellExperiment
#   (figshare 40582019: 10x Visium mouse brain, anterior section, 2,696 spots),
#   which the scDesign3 authors built from the Seurat `stxBrain` "anterior1"
#   object; that Seurat object is kept as data/simulated/source/brain_data.rda
#   (object `brain`) for provenance. The first 10 genes of the reference are
#   Calb2, Gng4, Ttr, S100a5, Nrgn, Doc2g, Ppp1r1b, Pcbp3, Fabp7, Gpsm1.
#
# RUNTIME  ~ minutes per replicate with a GP smoother (k = 400); uses all but
#           two cores by default.
# REPRODUCIBILITY  Replicate i is drawn after set.seed(BASE_SEED + i), which
#           is how the manuscript data were generated. Whether a re-run gives
#           bit-identical counts also depends on the scDesign3 version and on
#           its parallel back-end, so the files actually used for the manuscript
#           are committed in data/simulated/ and every S4 script reads those
#           files; regenerating them is optional. If you regenerate, compare
#           with the committed files before replacing them.
# =============================================================================

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(parallel)
})

suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))   # pin the RNG

# ---- parameters --------------------------------------------------------------
N_GENES      <- 10        # number of reference genes (= number of true modules)
N_REPLICATES <- 10        # scDesign3 realisations per gene
BASE_SEED    <- 123       # replicate i uses seed BASE_SEED + i
OUT_DIR      <- "data/simulated"
N_CORES      <- max(1, detectCores() - 2)
REFERENCE_URL <- "https://figshare.com/ndownloader/files/40582019"  # scDesign3 Visium example

cat("Using", N_CORES, "cores\n")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- reference data ----------------------------------------------------------
example_sce <- readRDS(url(REFERENCE_URL))
print(example_sce)
# Gene selection: the first N_GENES rows of the scDesign3 Visium example object,
# exactly as in the scDesign3 spatial vignette ("we subset the top 10 genes").
# The example object holds 1,000 genes that the scDesign3 authors pre-selected
# as the top spatially variable genes of the section, so the ten reference genes
# of the benchmark are its first ten spatially variable genes (Calb2, Gng4, Ttr,
# S100a5, Nrgn, Doc2g, Ppp1r1b, Pcbp3, Fabp7, Gpsm1).
example_sce <- example_sce[seq_len(N_GENES), ]

# ---- simulate ----------------------------------------------------------------
simulation_results <- vector("list", N_REPLICATES)
last_simu_sce <- NULL

for (i in seq_len(N_REPLICATES)) {
  cat("Running scDesign3 replicate", i, "of", N_REPLICATES, "\n")
  set.seed(BASE_SEED + i)

  example_simu <- scdesign3(
    sce              = example_sce,
    assay_use        = "counts",
    celltype         = "cell_type",
    pseudotime       = NULL,
    spatial          = c("spatial1", "spatial2"),
    other_covariates = NULL,
    mu_formula       = "s(spatial1, spatial2, bs = 'gp', k = 400)",
    sigma_formula    = "1",
    family_use       = "nb",
    n_cores          = N_CORES,
    usebam           = FALSE,
    corr_formula     = "1",
    copula           = "gaussian",
    DT               = TRUE,
    pseudo_obs       = FALSE,
    return_model     = FALSE,
    nonzerovar       = FALSE,
    parallelization  = "pbmcapply"
  )

  simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count),
                                   colData = example_simu$new_covariate)
  log_counts <- log1p(counts(simu_sce))
  rownames(log_counts) <- paste0(rownames(log_counts), "_", i)   # "<gene>_<replicate>"
  simulation_results[[i]] <- log_counts
  last_simu_sce <- simu_sce
}

# ---- merge and save ----------------------------------------------------------
count_mat_simul <- as.matrix(do.call(rbind, simulation_results))
coordinates_simul <- as.data.frame(colData(last_simu_sce)[, c("spatial1", "spatial2")])
rownames(coordinates_simul) <- colnames(count_mat_simul)

cat("Reference dimensions :", dim(example_sce), "\n")
cat("Simulated dimensions :", dim(count_mat_simul), "\n")
cat("Coordinates          :", dim(coordinates_simul), "\n")
cat("First gene names     :", head(rownames(count_mat_simul), 12), "\n")

gene_base <- sub("_[0-9]+$", "", rownames(count_mat_simul))
ground_truth <- data.frame(Gene = rownames(count_mat_simul),
                           Gene_Base_Name = gene_base,
                           True_Cluster = match(gene_base, unique(gene_base)),
                           stringsAsFactors = FALSE)
ground_truth <- ground_truth[order(ground_truth$Gene), ]
rownames(ground_truth) <- NULL

# Do not overwrite the committed files silently: compare first.
for (f in c("count_mat_simul.rds", "coordinates_simul.rds", "ground_truth.rds")) {
  path <- file.path(OUT_DIR, f)
  new  <- switch(f, count_mat_simul.rds = count_mat_simul, coordinates_simul.rds = coordinates_simul,
                 ground_truth.rds = ground_truth)
  if (file.exists(path)) {
    same <- isTRUE(all.equal(readRDS(path), new, check.attributes = TRUE))
    cat(sprintf("%-22s exists: regenerated object %s the committed one\n", f,
                if (same) "EQUALS" else "DIFFERS FROM"))
    saveRDS(new, sub("\\.rds$", "_regenerated.rds", path))
  } else saveRDS(new, path)
}
cat("Saved to", OUT_DIR, "\n")
