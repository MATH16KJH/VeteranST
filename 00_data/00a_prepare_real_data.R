# =============================================================================
# 00a_prepare_real_data.R -- provenance of the real-data inputs (data/real/)
# -----------------------------------------------------------------------------
# The files in data/real/ are COMMITTED to this repository; nothing has to be
# run before the simulations. This script documents how they were derived from
# the processed Seurat objects and the BANYAN fits, and can regenerate them
# (and check that the regenerated files equal the committed ones) if those
# objects are downloaded from the Zenodo archive.
#
# SOURCE OBJECTS (Zenodo, not committed because of their size)
#   pt16_norm.rda                  Seurat object `p16_norm`   (P16, 213 MB)
#   heart_D14_seurat_processed.RData  Seurat object `heart_D14` (D14, 78 MB)
#   breast_processed_seurat.RData  Seurat object `breast`     (BR, 284 MB)
#   pt16_fit4.rda, D14_banyan_fit4.rda, breast_fit_banyan.RData
#       BANYAN fits `pt16_fit4`, `D14_fit4`, `breast_fit` (committed in
#       data/real/banyan/, 0.5-1.8 MB each)
#
# DERIVATION (one line per file; identical for the three datasets)
#   count_mat_<key>.rda    <- <seurat>@assays[["SCT"]]@data      (dgCMatrix, genes x spots,
#                              SCTransform log-normalised values; 16,074 / 11,356 /
#                              21,129 genes x 1,650 / 1,967 / 3,798 spots)
#   coordinates_<key>.csv  <- Seurat::GetTissueCoordinates(<seurat>)
#                              (P16: imagerow/imagecol pixels; D14, BR: array x/y)
#   clusters_<key>.csv     <- <banyan_fit>[["z"]]   (posterior domain label, K = 4 / 4 / 5)
#   All three share the spot order of the Seurat object (checked below).
#
# USAGE  Rscript 00_data/00a_prepare_real_data.R [path/to/source_objects]
#        Without the source objects the script only validates data/real/.
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
suppressPackageStartupMessages(library(Matrix))
source("R/utils_data.R")

SOURCE_DIR <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else NULL
BANYAN_DIR <- "data/real/banyan"          # committed fit objects

SOURCES <- list(
  P16 = list(seurat_file = "pt16_norm.rda",                    seurat_object = "p16_norm",
             banyan_file = "pt16_fit4.rda",                    banyan_object = "pt16_fit4"),
  D14 = list(seurat_file = "heart_D14_seurat_processed.RData", seurat_object = "heart_D14",
             banyan_file = "D14_banyan_fit4.rda",              banyan_object = "D14_fit4"),
  BR  = list(seurat_file = "breast_processed_seurat.RData",    seurat_object = "breast",
             banyan_file = "breast_fit_banyan.RData",          banyan_object = "breast_fit")
)

# ---- derivation --------------------------------------------------------------
derive_dataset <- function(name, source_dir) {
  src <- SOURCES[[name]]
  seurat <- read_r_object(file.path(source_dir, src$seurat_file), src$seurat_object)
  fit    <- read_r_object(file.path(BANYAN_DIR, src$banyan_file), src$banyan_object)
  expr   <- seurat@assays[["SCT"]]@data                         # dgCMatrix, genes x spots
  coords <- Seurat::GetTissueCoordinates(seurat)
  labels <- as.integer(fit[["z"]])
  stopifnot(ncol(expr) == nrow(coords), ncol(expr) == length(labels))
  # BANYAN was run on the same spot order as the Seurat object
  # (fit$coords equals the Seurat image coordinates row by row).
  list(expr = expr, coords = coords[, 1:2], labels = labels, spots = colnames(expr))
}

write_dataset <- function(name, d) {
  cfg <- REAL_DATASETS[[name]]
  assign(cfg$object, d$expr)
  save(list = cfg$object, file = cfg$file)
  write.csv(data.frame(spot = d$spots, clusters = d$labels), cfg$labels, row.names = FALSE)
  write.csv(d$coords, cfg$coords, row.names = TRUE)
  cat("Wrote", cfg$file, cfg$labels, cfg$coords, "\n")
}

# ---- validation of the committed files ---------------------------------------
validate_dataset <- function(name) {
  cfg    <- REAL_DATASETS[[name]]
  expr   <- load_real_dataset(name)
  labels <- load_domain_labels(name, spot_names = colnames(expr))
  coords <- load_real_coordinates(name, spot_names = colnames(expr))
  cat(sprintf("%-3s %s x %s  | domains: %s | class %s\n", name, nrow(expr), ncol(expr),
              paste(table(labels), collapse = ", "), class(expr)[1]))
  invisible(list(expr = expr, labels = labels, coords = coords))
}

for (name in DATASET_ORDER) {
  if (!is.null(SOURCE_DIR)) {
    src <- SOURCES[[name]]
    if (file.exists(file.path(SOURCE_DIR, src$seurat_file))) {
      d <- derive_dataset(name, SOURCE_DIR)
      committed <- tryCatch(validate_dataset(name), error = function(e) NULL)
      if (!is.null(committed)) {
        same <- isTRUE(all.equal(as.matrix(committed$expr), as.matrix(d$expr), check.attributes = TRUE)) &&
                identical(committed$labels, d$labels)
        cat(name, ": regenerated files", if (same) "EQUAL" else "DIFFER FROM", "the committed ones\n")
      } else write_dataset(name, d)
    } else cat(name, ": Seurat object not found in", SOURCE_DIR, "- skipping derivation\n")
  }
  validate_dataset(name)
  # the committed labels must equal the committed BANYAN fit
  fit <- read_r_object(file.path(BANYAN_DIR, SOURCES[[name]]$banyan_file), SOURCES[[name]]$banyan_object)
  stopifnot(identical(as.integer(fit[["z"]]), load_domain_labels(name)))
}
cat("data/real/ is consistent with R/utils_data.R\n")
