# VeteranST web application

Source of the R/Shiny application described in the main text (Section 3,
*Implementation*) and served at <https://chunglab.bmi.osumc.edu/VeteranST>.
The public server requires no installation; this folder lets you run the same
application locally or deploy it on your own Shiny Server.

## Files

| File | Content |
|---|---|
| `app.R` | The application: UI and server (data upload, tissue-domain assignment, DVG test, spatial domain variability test, spatial gene co-expression analysis). |
| `help_content.R` | Help & Tutorial tab and the per-step "?" pop-ups; sourced by `app.R`. The values still marked `★` (Zenodo DOI, runtime, version) must be filled in before release. |
| `test_data/count_mat_p16.rda` | Example expression matrix (`count_mat_p16`, 16,074 genes × 1,650 spots, SCTransform-normalised log values = `data` slot of the Seurat `SCT` assay); identical to `data/real/count_mat_p16.rda`. |
| `test_data/coordinates_p16.rda` | Example spot coordinates (`coordinates_p16`, `imagerow`/`imagecol` from `Seurat::GetTissueCoordinates()`). |
| `test_data/clusters_p16.rds` | Example BANYAN tissue-domain labels (named integer vector, 4 domains of 311/609/321/409 spots). |
| `install_app_dependencies.R` | Installs the R packages used by the application. |

The three `test_data/` files are what the **Load Example Data** button loads;
they are the melanoma brain metastasis section (P16) analysed in the paper.

## Running locally

```r
# once
source("app/install_app_dependencies.R")

# every time (from the repository root, or from within app/)
shiny::runApp("app")
```

The application opens in the system browser (`options(shiny.launch.browser = TRUE)`
at the end of `app.R`). The DVG test and the similarity computation run in
parallel on all but one core (`doParallel`); a full-transcriptome DVG test on the
example data takes a few minutes.

## Method settings implemented in the application

The statistical components follow the main text and Supplementary Sections
S1–S2 exactly; the simulation scripts in the parent folders re-implement the
same definitions (see `R/utils_spatial_similarity.R`).

| Component | Implementation in `app.R` |
|---|---|
| DVG test | `car::leveneTest(center = median)` = Brown–Forsythe; Bonferroni; one-way ANOVA at the same FWER defines DVG-only genes |
| Spatial domain variability test | IQR per domain (`stats::IQR`); descending mid-ranks (`rank(-x, ties.method = "average")`); `friedman.test`; Nemenyi CD = `qtukey(1 − α, k, Inf) · sqrt(k(k+1)/(12N))` |
| Spatial smoothing | `RANN::nn2`, neighbourhoods of 7 / 19 / 37 / 61 spots (spot itself included); local variance = mean of squares − squared mean (1/\|N\| denominator) |
| Sparse gene graph | 19 feature spaces → HNSW (`RcppHNSW`, cosine, M = 32, ef = 200, k = ⌈√N⌉, ef_search = max(2k, 200)) → support-weighted union → disparity filter (α = 0.05) |
| Similarity measures | on the 37-spot smoothed data: Pearson of local means (negative values clipped to 0), Hotspot Jaccard (top 10 % local mean), High-/Low-variability Jaccard (top / bottom 10 % local variance) |
| Consensus | Leiden (`cluster_leiden(objective_function = "modularity", resolution = 1)`, igraph ≥ 1.3) on each of the four graphs; genes linked only if co-clustered in all four; modules = connected components; minimum module size filter |

## Deployment notes

* Upload limit is set to 1 GB (`options(shiny.maxRequestSize)`); adjust to the
  server's memory.
* `webshot` + PhantomJS are used only for the PNG export of the module network.
* The `rsconnect/` folder created by `rsconnect::deployApp()` is git-ignored.
