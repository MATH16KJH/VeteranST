# =============================================================================
# install_dependencies.R -- install every R package used by the simulation
# scripts (the web application has its own installer, app/install_app_dependencies.R)
# Run once:  Rscript install_dependencies.R
# =============================================================================

cran_pkgs <- c(
  "Matrix", "igraph", "RANN", "aricode", "cluster", "matrixStats",   # methods / evaluation
  "ggplot2", "dplyr", "tidyr", "scales"                              # figures and tables
)
missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

# Only needed to re-derive the committed data from the source objects
# (00_data/00a: Seurat; 00_data/00b: SingleCellExperiment + scDesign3).
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc_pkgs <- c("SingleCellExperiment")
bioc_missing <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(bioc_missing)) BiocManager::install(bioc_missing, update = FALSE, ask = FALSE)
if (!requireNamespace("scDesign3", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools", repos = "https://cloud.r-project.org")
  devtools::install_github("SONGDONGYUAN1994/scDesign3")
}
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat", repos = "https://cloud.r-project.org")

# igraph >= 1.3.0 is required (cluster_leiden() with the `resolution` argument).
stopifnot(utils::packageVersion("igraph") >= "1.3.0")
cat("All dependencies are installed.\n")
