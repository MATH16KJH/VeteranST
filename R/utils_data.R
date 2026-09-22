# =============================================================================
# R/utils_data.R -- dataset registry and loaders
# -----------------------------------------------------------------------------
# Sourced by every simulation script. Defines where the input data live and
# how the tissue domains (groups) are encoded.
#
# REAL DATA (Supplementary Sections S1 and S3):
#   data/real/count_mat_<key>.rda   the exact expression objects used for the
#       submission: `Matrix::dgCMatrix`, genes x spots, SCTransform log-
#       normalised values (Seurat `@assays$SCT@data`), same row / column order
#       as in the development scripts. Loaded with load_real_dataset().
#   data/real/clusters_<key>.csv    BANYAN tissue-domain label (1..K) of every
#       spot, same spot order as the matrix columns (`banyan_fit$z`).
#   data/real/coordinates_<key>.csv spot coordinates from Seurat
#       GetTissueCoordinates(), same spot order.
#   The null simulations permute expression across spots and use the domain
#   labels only through their sizes (`banyan_sizes`, = table(labels)); the
#   spot-level labels and coordinates are used for validation (00a) and for
#   exporting web-application inputs (00c).
#   See data/real/README.md for provenance and 00_data/00a_prepare_real_data.R
#   for the derivation from the Seurat / BANYAN objects.
#
# SIMULATED DATA (Supplementary Section S4):
#   Produced by 00_data/00b_generate_simulated_data.R (scDesign3). Gene names
#   are "<gene>_<replicate>", which is what defines the ground-truth modules.
# =============================================================================

REAL_DATASETS <- list(
  P16 = list(
    file         = "data/real/count_mat_p16.rda",
    object       = "count_mat_p16",
    labels       = "data/real/clusters_p16.csv",
    coords       = "data/real/coordinates_p16.csv",
    banyan_sizes = c(311, 609, 321, 409),                    # 4 domains, 1,650 spots
    n_genes      = 16074,
    description  = "Human melanoma brain metastasis, 10x Visium (Sudmeier et al., 2022; GEO)"
  ),
  BR = list(
    file         = "data/real/count_mat_br.rda",
    object       = "count_mat_br",
    labels       = "data/real/clusters_br.csv",
    coords       = "data/real/coordinates_br.csv",
    banyan_sizes = c(861, 665, 942, 741, 589),               # 5 domains, 3,798 spots
    n_genes      = 21129,
    description  = "Human breast cancer (invasive ductal carcinoma), 10x Visium (10x Genomics, 2019)"
  ),
  D14 = list(
    file         = "data/real/count_mat_d14.rda",
    object       = "count_mat_d14",
    labels       = "data/real/clusters_d14.csv",
    coords       = "data/real/coordinates_d14.csv",
    banyan_sizes = c(754, 392, 240, 581),                    # 4 domains, 1,967 spots
    n_genes      = 11356,
    description  = "Chicken heart, day 14, 10x Visium (Mantri et al., 2021; GEO)"
  )
)

# Order in which the three datasets are shown in the supplementary figures.
DATASET_ORDER <- c("P16", "BR", "D14")

# Simulated benchmark data for Section S4 (100 genes = 10 gene programs x 10
# scDesign3 replicates; see 00b_generate_simulated_data.R).
SIMULATED_DATA <- list(
  expr         = "data/simulated/count_mat_simul.rds",    # 100 genes x 2,696 spots, log1p counts
  coords       = "data/simulated/coordinates_simul.rds",  # spots x 2 (spatial1, spatial2)
  ground_truth = "data/simulated/ground_truth.rds"        # Gene, Gene_Base_Name, True_Cluster
)

#' Expand per-domain spot counts into a spot-level group label vector.
domain_index <- function(group_sizes) {
  rep(seq_along(group_sizes), times = group_sizes)
}

#' Column index ranges of consecutive blocks of the given sizes, i.e. the
#' spots of every domain after the spots of a permuted matrix are cut into
#' blocks (list of integer vectors, one per domain).
block_columns <- function(group_sizes) {
  starts <- cumsum(c(0, group_sizes))
  lapply(seq_along(group_sizes), function(i) (starts[i] + 1):starts[i + 1])
}

#' Independent random permutation of the columns of every row (gene) of a
#' dense matrix -- the null model of Supplementary Section S3: every gene is
#' shuffled across all spots on its own, so a single permutation shared by all
#' genes (which would preserve spot-level effects such as sequencing depth) is
#' deliberately NOT used. Draws from R's RNG: call set.seed() before it.
permute_genes <- function(mat) {
  n <- ncol(mat)
  for (k in seq_len(nrow(mat))) mat[k, ] <- mat[k, sample.int(n)]
  mat
}

#' Read an R object from .rds or .rda/.RData (the latter must contain `object`).
read_r_object <- function(path, object = NULL) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (is.null(object)) {
    if (length(loaded) != 1) stop(path, " contains several objects: ", paste(loaded, collapse = ", "))
    object <- loaded
  }
  if (!object %in% loaded) stop(path, " does not contain an object named ", object)
  get(object, envir = env)
}

#' Load one real dataset and verify that it matches the registry.
#' @param name      one of names(REAL_DATASETS)
#' @param as_dense  convert the sparse matrix to a dense base matrix
load_real_dataset <- function(name, as_dense = FALSE) {
  cfg <- REAL_DATASETS[[name]]
  if (is.null(cfg)) stop("Unknown dataset: ", name)
  if (!file.exists(cfg$file)) {
    stop("File not found: ", cfg$file, "\nSee data/real/README.md.")
  }
  mat <- read_r_object(cfg$file, cfg$object)
  if (ncol(mat) != sum(cfg$banyan_sizes)) {
    stop(sprintf("%s: matrix has %d columns but banyan_sizes sum to %d.",
                 name, ncol(mat), sum(cfg$banyan_sizes)))
  }
  if (!is.null(cfg$n_genes) && nrow(mat) != cfg$n_genes) {
    warning(sprintf("%s: %d genes, but the manuscript analyses used %d.", name, nrow(mat), cfg$n_genes))
  }
  if (as_dense) mat <- as.matrix(mat)
  mat
}

#' Domain label of every spot (integer 1..K) of a real dataset, in matrix
#' column order. `spot_names` (e.g. colnames of the matrix) is checked when given.
load_domain_labels <- function(name, spot_names = NULL) {
  cfg <- REAL_DATASETS[[name]]
  if (!file.exists(cfg$labels)) stop("File not found: ", cfg$labels)
  tab <- utils::read.csv(cfg$labels, stringsAsFactors = FALSE)
  lab <- as.integer(tab[[2]])
  if (!identical(as.integer(table(lab)), as.integer(cfg$banyan_sizes))) {
    stop(sprintf("%s: domain label counts (%s) differ from banyan_sizes (%s).", name,
                 paste(table(lab), collapse = ", "), paste(cfg$banyan_sizes, collapse = ", ")))
  }
  if (!is.null(spot_names) && !identical(as.character(tab[[1]]), as.character(spot_names))) {
    stop(name, ": spot order of ", cfg$labels, " differs from the expression matrix.")
  }
  lab
}

#' Spot coordinates of a real dataset (spots x 2 data.frame, spot names as
#' rownames, same order as the matrix columns).
load_real_coordinates <- function(name, spot_names = NULL) {
  cfg <- REAL_DATASETS[[name]]
  if (!file.exists(cfg$coords)) stop("File not found: ", cfg$coords)
  coords <- utils::read.csv(cfg$coords, row.names = 1, check.names = FALSE)
  if (nrow(coords) != sum(cfg$banyan_sizes)) {
    stop(sprintf("%s: coordinates have %d rows, expected %d.", name, nrow(coords),
                 sum(cfg$banyan_sizes)))
  }
  if (!is.null(spot_names) && !identical(rownames(coords), as.character(spot_names))) {
    stop(name, ": spot order of ", cfg$coords, " differs from the expression matrix.")
  }
  coords[, 1:2]
}

#' Load the scDesign3 benchmark data (expression + coordinates, aligned).
load_simulated_data <- function(paths = SIMULATED_DATA) {
  for (p in c(paths$expr, paths$coords)) {
    if (!file.exists(p)) {
      stop("File not found: ", p,
           "\nThe committed benchmark data are missing; see README, section 'Data'.")
    }
  }
  expr   <- as.matrix(readRDS(paths$expr))
  coords <- as.data.frame(readRDS(paths$coords))
  if (nrow(coords) != ncol(expr)) {
    stop("Simulated data: number of coordinate rows differs from number of expression columns.")
  }
  if (!is.null(rownames(coords)) && !is.null(colnames(expr)) &&
      all(colnames(expr) %in% rownames(coords))) {
    coords <- coords[colnames(expr), , drop = FALSE]
  }
  list(expr = expr, coords = as.matrix(coords[, 1:2]))
}

#' Pin the random number generator so that seeds give the same streams on every
#' R version >= 3.6 (Mersenne-Twister, inversion, rejection sampling). Called at
#' the top of every script; the Leiden implementation of igraph (>= 1.3) draws
#' from R's RNG, so set.seed() also fixes the community detection.
set_rng <- function() {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  invisible(NULL)
}

#' Number of physical cores to use (leaves `reserve` cores free).
default_cores <- function(reserve = 1, max_cores = Inf) {
  n <- parallel::detectCores()
  if (is.na(n)) n <- 1
  max(1, min(n - reserve, max_cores))
}
