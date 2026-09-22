# =============================================================================
# 00c_export_app_test_input.R -- export a real dataset as VeteranST app input
# -----------------------------------------------------------------------------
# Writes the expression matrix and spot coordinates of one real dataset as CSV
# (and optionally XLSX) files in the layout accepted by the VeteranST web
# application (genes x spots matrix with gene names in the first column; spots x
# 2 coordinates with spot barcodes in the first column). Used to create the
# example / test inputs distributed with the application.
#
# USAGE   Rscript 00_data/00c_export_app_test_input.R [dataset] [outdir]
#         e.g.  Rscript 00_data/00c_export_app_test_input.R D14 data/app_test_input
# =============================================================================

if (!file.exists("R/utils_data.R")) stop("Run this script from the repository root.")
source("R/utils_data.R")
set_rng()                                   # pin the RNG (see R/utils_data.R)

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else "D14"
OUT_DIR <- if (length(args) >= 2) args[2] else "data/app_test_input"
WRITE_XLSX <- FALSE                                   # requires the openxlsx package
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

expr   <- load_real_dataset(DATASET, as_dense = TRUE)
coords <- load_real_coordinates(DATASET)
if (is.null(coords)) stop("No coordinates file for ", DATASET, " (see R/utils_data.R).")
key <- tolower(DATASET)

write.csv(as.data.frame(expr), file.path(OUT_DIR, paste0("count_mat_", key, ".csv")), row.names = TRUE)
write.csv(coords,              file.path(OUT_DIR, paste0("coordinates_", key, ".csv")), row.names = TRUE)
if (WRITE_XLSX) {
  openxlsx::write.xlsx(coords, file.path(OUT_DIR, paste0("coordinates_", key, ".xlsx")), rowNames = TRUE)
}
cat("Exported", DATASET, "to", OUT_DIR, "\n")
