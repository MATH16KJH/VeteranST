# Installs the packages needed to run the VeteranST web application locally.
#   Rscript app/install_app_dependencies.R
pkgs <- c(
  "shiny", "shinythemes", "shinyjs", "dplyr", "DT", "car", "ggplot2", "Matrix",
  "foreach", "doParallel", "webshot", "cowplot", "Cairo", "openxlsx", "zip",
  "igraph", "visNetwork", "RANN", "RcppHNSW", "matrixStats", "viridis"
)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

# igraph >= 1.3.0 is required for cluster_leiden().
stopifnot(utils::packageVersion("igraph") >= "1.3.0")

# webshot needs PhantomJS once for the network-image export (Download Module Information).
if (!webshot::is_phantomjs_installed()) webshot::install_phantomjs()

cat("All application dependencies are installed.\n")
