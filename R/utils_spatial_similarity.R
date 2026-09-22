# =============================================================================
# R/utils_spatial_similarity.R -- spatial gene co-expression machinery
# -----------------------------------------------------------------------------
# Shared implementation used by every script in S4_spatial_coexpression/:
#   compute_smoothed_features()    local mean / local variance over the k nearest
#                                  spots (k includes the spot itself; Section 2.3)
#   compute_similarity_matrices()  the five similarity measures of Section 2.3
#   leiden_on_similarity()         Leiden (modularity) on one similarity graph
#   consensus_clustering()         PROPOSED: strict multi-measure consensus
#                                  (genes are co-clustered only if co-clustered
#                                  in every measure-specific Leiden partition)
#   weighted_clustering()          BASELINE: weighted aggregation of the
#                                  min-max-normalised similarity matrices,
#                                  followed by a single Leiden run (Section S4.1)
#
# This is a compact, dense-matrix re-implementation of the code originally used
# for the simulation studies (all gene pairs are compared explicitly; suitable
# for the 100-gene benchmark). The VeteranST web application applies the same
# measures on an HNSW/disparity-filter sparsified graph for genome-scale
# gene sets (main text, Section 2.3).
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
  library(RANN)
})

MEASURE_NAMES <- c("Jaccard-Zero", "Pearson-Mean", "Jaccard-Hotspot",
                   "Jaccard-Var-Top", "Jaccard-Var-Bot")
MEASURE_ABBR  <- c("JZ", "PM", "JH", "JVT", "JVB")
names(MEASURE_ABBR) <- MEASURE_NAMES

# Defaults of the method (manuscript Section 2.3 and the web application),
# selected in Supplementary Section S4.1 (S4c: measure set; S4d: neighbourhood
# size). S4e and the application use them.
DEFAULT_MEASURE_IDX <- 2:5   # PM + JH + JVT + JVB (Zero Jaccard excluded)
DEFAULT_N_NEIGHBORS <- 37    # neighbourhood size for smoothing (7 / 19 / 37 / 61 compared in S4d)

#' Short label for a measure combination, e.g. c(2,3,4,5) -> "PM+JH+JVT+JVB".
combo_label <- function(measure_idx, sep = "+") {
  paste(MEASURE_ABBR[MEASURE_NAMES[measure_idx]], collapse = sep)
}

#' Make sure `expr` (genes x spots) and `coords` (spots x 2) are aligned.
align_expr_coords <- function(expr, coords) {
  expr   <- as.matrix(expr)
  coords <- as.matrix(as.data.frame(coords)[, 1:2])
  if (nrow(coords) != ncol(expr)) {
    stop("coords has ", nrow(coords), " rows but expr has ", ncol(expr), " columns.")
  }
  if (!is.null(rownames(coords)) && !is.null(colnames(expr)) &&
      all(colnames(expr) %in% rownames(coords))) {
    coords <- coords[colnames(expr), , drop = FALSE]
  }
  if (is.null(rownames(expr))) rownames(expr) <- paste0("Gene_", seq_len(nrow(expr)))
  list(expr = expr, coords = coords)
}

#' Local (neighbourhood-smoothed) mean and variance of every gene at every spot.
#' The neighbourhood of a spot is its `n_neighbors` nearest spots including
#' itself (hexagonal Visium grid: 7 / 19 / 37 / 61 spots). The local variance is
#' the sample variance (denominator n-1), as in the original simulation code.
compute_smoothed_features <- function(expr, coords, n_neighbors = DEFAULT_N_NEIGHBORS) {
  n_spots <- ncol(expr)
  if (n_neighbors < 2) stop("n_neighbors must be >= 2")
  n_neighbors <- min(n_neighbors, n_spots)
  nn_idx <- RANN::nn2(coords, coords, k = n_neighbors)$nn.idx   # column 1 = the spot itself
  # Sparse averaging operator: W[neighbour, spot] = 1/k
  W <- Matrix::sparseMatrix(i = as.vector(nn_idx),
                            j = rep(seq_len(n_spots), times = n_neighbors),
                            x = 1 / n_neighbors,
                            dims = c(n_spots, n_spots))
  local_mean <- as.matrix(expr %*% W)
  local_sq   <- as.matrix((expr^2) %*% W)
  local_var  <- (local_sq - local_mean^2) * n_neighbors / (n_neighbors - 1)
  local_var[local_var < 0] <- 0                                  # numerical guard
  dimnames(local_mean) <- dimnames(local_var) <- dimnames(expr)
  list(local_mean = local_mean, local_var = local_var, nn_idx = nn_idx)
}

#' Jaccard similarity between the row-sets of a genes x spots 0/1 matrix.
#' Empty sets give similarity 0 (the diagonal is set to 1).
jaccard_matrix <- function(B) {
  B <- as.matrix(B) * 1
  inter <- tcrossprod(B)
  sizes <- rowSums(B)
  uni   <- outer(sizes, sizes, "+") - inter
  J <- ifelse(uni > 0, inter / uni, 0)
  diag(J) <- 1
  dimnames(J) <- list(rownames(B), rownames(B))
  J
}

#' The five similarity measures of Section 2.3, computed on smoothed expression.
#' Returns a named list of genes x genes matrices (names = MEASURE_NAMES).
#' @param hotspot_quantile   spots with local mean >= this quantile are "hotspots"
#' @param high_var_quantile  spots with local variance >= this quantile
#' @param low_var_quantile   spots with local variance <= this quantile
compute_similarity_matrices <- function(expr, coords,
                                        n_neighbors       = DEFAULT_N_NEIGHBORS,
                                        hotspot_quantile  = 0.90,
                                        high_var_quantile = 0.90,
                                        low_var_quantile  = 0.10) {
  al     <- align_expr_coords(expr, coords)
  feats  <- compute_smoothed_features(al$expr, al$coords, n_neighbors)
  M <- feats$local_mean
  V <- feats$local_var
  gene_names <- rownames(al$expr)

  row_q <- function(X, p) apply(X, 1, stats::quantile, probs = p, na.rm = TRUE, names = FALSE)
  Z  <- (M == 0)                        # 1. Zero Jaccard: spots with zero smoothed expression
  H  <- (M >= row_q(M, hotspot_quantile))   # 2. Hotspot Jaccard: top-10% smoothed expression
  VT <- (V >= row_q(V, high_var_quantile))  # 3. High-variability Jaccard: top-10% local variance
  VB <- (V <= row_q(V, low_var_quantile))   # 4. Low-variability Jaccard: bottom-10% local variance
  rownames(Z) <- rownames(H) <- rownames(VT) <- rownames(VB) <- gene_names

  PM <- suppressWarnings(stats::cor(t(M), method = "pearson"))   # 5. Pearson on smoothed mean
  PM[is.na(PM)] <- 0
  diag(PM) <- 1
  dimnames(PM) <- list(gene_names, gene_names)

  mats <- list(jaccard_matrix(Z), PM, jaccard_matrix(H), jaccard_matrix(VT), jaccard_matrix(VB))
  names(mats) <- MEASURE_NAMES
  mats
}

#' Leiden community detection (modularity objective) on one similarity matrix.
#' Negative similarities (only Pearson can be negative) are set to 0 (no edge).
#' @param seed  optional RNG seed for reproducibility
leiden_on_similarity <- function(sim_mat, resolution = 1, seed = NULL) {
  adj <- as.matrix(sim_mat)
  adj[adj < 0] <- 0
  diag(adj) <- 0
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  if (!is.null(seed)) set.seed(seed)
  cl <- igraph::cluster_leiden(g, objective_function = "modularity", resolution = resolution)
  m <- as.integer(igraph::membership(cl))
  names(m) <- rownames(sim_mat)
  m
}

#' PROPOSED METHOD: strict consensus over measure-specific Leiden partitions.
#' Two genes end up in the same module only if they are co-clustered in ALL of
#' the selected measures; modules are the connected components of that
#' agreement graph (equivalently, the meet of the partitions).
#' @param matrices          output of compute_similarity_matrices()
#' @param measure_idx       indices into MEASURE_NAMES (default PM+JH+JVT+JVB)
#' @param seed              if given, measure k uses seed + k
#' @return named integer vector of module labels (attribute "memberships" holds
#'         the per-measure partitions when return_memberships = TRUE)
consensus_clustering <- function(matrices, measure_idx = DEFAULT_MEASURE_IDX,
                                 resolution = 1, seed = NULL, return_memberships = FALSE) {
  gene_names <- rownames(matrices[[1]])
  n_genes    <- length(gene_names)
  memberships <- vector("list", length(measure_idx))
  for (k in seq_along(measure_idx)) {
    memberships[[k]] <- leiden_on_similarity(matrices[[measure_idx[k]]], resolution = resolution,
                                             seed = if (is.null(seed)) NULL else seed + k)
  }
  names(memberships) <- names(matrices)[measure_idx]
  agree <- matrix(TRUE, n_genes, n_genes)
  for (m in memberships) agree <- agree & outer(m, m, "==")
  g <- igraph::graph_from_adjacency_matrix(agree * 1, mode = "undirected", diag = FALSE)
  labels <- igraph::components(g)$membership
  names(labels) <- gene_names
  if (return_memberships) attr(labels, "memberships") <- memberships
  labels
}

#' BASELINE: weighted aggregation of similarity matrices + one Leiden run.
#' Each matrix is min-max normalised to [0,1] (over all entries, as in the
#' original code) and combined as sum_k w_k * S_k (weights normalised to 1).
#' @param weights  numeric vector aligned with measure_idx; NULL = equal weights
weighted_clustering <- function(matrices, measure_idx = DEFAULT_MEASURE_IDX,
                                weights = NULL, resolution = 1, seed = NULL) {
  if (is.null(weights)) weights <- rep(1, length(measure_idx))
  stopifnot(length(weights) == length(measure_idx))
  weights <- weights / sum(weights)
  combined <- 0
  for (k in seq_along(measure_idx)) {
    mat <- as.matrix(matrices[[measure_idx[k]]])
    rng <- range(mat, na.rm = TRUE)
    if (rng[2] > rng[1]) mat <- (mat - rng[1]) / (rng[2] - rng[1])
    combined <- combined + weights[k] * mat
  }
  dimnames(combined) <- dimnames(matrices[[measure_idx[1]]])
  leiden_on_similarity(combined, resolution = resolution, seed = seed)
}
