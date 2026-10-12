# VeteranST

Source code, test data and simulation studies of

> Kim J, Jung K, Shin J, Chung D, Shin S. **VeteranST: Transcriptional Variability and Gene
> Co-Expression Analysis in Spatial Transcriptomics.** *Bioinformatics* (submitted).

The repository has two parts:

* **`app/`** — the VeteranST web application (R/Shiny) that is served at
  <https://chunglab.bmi.osumc.edu/VeteranST>, with the example data it ships
  with. It runs locally with `shiny::runApp("app")`; see `app/README.md`.
* **Simulation studies** (everything else) — stand-alone R scripts that
  reproduce every figure and table of the simulation studies (main-text
  Section 4.1 and Figure 2; Supplementary Sections S1, S3 and S4), together
  with the exact input data they used.

Archived version of this repository used for the submission: **Zenodo DOI TODO**.

## Web application (`app/`)

| File | Content |
|---|---|
| `app/app.R` | UI and server of the application (five analysis tabs + Help & Tutorial). |
| `app/help_content.R` | Help texts and per-step "?" pop-ups (values marked `★` are release placeholders). |
| `app/test_data/` | Example data loaded by the *Load Example Data* button: melanoma brain metastasis section P16 (`count_mat_p16.rda`, `coordinates_p16.rda`, `clusters_p16.rds` = BANYAN domains). |
| `app/install_app_dependencies.R` | Installs the packages the application needs. |

```r
source("app/install_app_dependencies.R")   # once
shiny::runApp("app")
```

The statistical settings implemented in the application (Brown–Forsythe DVG
test, IQR-based Friedman/Nemenyi test with descending ranks, 19-feature HNSW
graph with disparity filter, four similarity measures on 37-spot smoothed
data, strict four-way Leiden consensus) are summarised in `app/README.md`;
the simulation code below re-implements the same definitions.

# Simulation studies

## What reproduces what

Figure and table numbers refer to the supplementary file unless marked
"main text". All outputs are written under `results/`.

| Manuscript item | Script | Main output |
|---|---|---|
| Fig. S1, S3 – type I error of variance-homogeneity tests, α = 0.01 / 0.05 (Sec. S1.2–S1.3) | `S1_gene_level_type1_error/S1_type1_error_gene_level.R` (family `variance`) | `S1_gene_level_type1_error/<Test>/scatter_comparison_<Test>_<α>.jpg`, `comparison_hist_<Test>_<α>.jpg` |
| Fig. S2, S4 – type I error of mean-comparison tests, α = 0.01 / 0.05 (Sec. S1.2, S1.4) | same script (family `mean`) | `S1_gene_level_type1_error/ANOVA/…`, `Kruskal_Wallis/…` |
| Table S1 – domain variability tests, equal domain sizes (Sec. S3) | `S3_domain_variability_type1_error/S3a_type1_error_domain_variability.R` | `S3_domain_variability_type1_error/table_S1_type1_error_equal.csv` |
| Table S2 – same with the BANYAN domain sizes (Sec. S3) | same run | `…/table_S2_type1_error_banyan.csv` |
| Sec. S3.1 text: single-count fractions, exact P(X = 0), rejection thresholds; **Table S3** | `S3_domain_variability_type1_error/S3b_sparse_gene_variance_quantities.R` | `S3b_sparse_gene_variance_quantities/*.csv` (`table_S3_coefficients.csv`) |
| Fig. S5 – ordering by variance vs. ordering by expressing-spot proportions (Sec. S3.1) | `S3_domain_variability_type1_error/S3c_rank_agreement_largest_smallest.R` | `S3c_rank_agreement_largest_smallest/fig_S5_rank_agreement.jpg` |
| Fig. S6 – P(s²_L > s²_S) by sparsity (Sec. S3.1) | `S3_domain_variability_type1_error/S3d_pairwise_probability_largest_smallest.R` | `S3d_pairwise_probability_largest_smallest/fig_S6_pairwise_probability.jpg` |
| Fig. S7 – Σᵢ(R̄ᵢ − (m+1)/2)² under the null; mean-rank shifts, rejection rates, share of informative genes for the IQR (Sec. S3.1) | `S3_domain_variability_type1_error/S3e_rank_shift_amplification.R` | `S3e_rank_shift_amplification/fig_S7_sum_delta2.jpg`, `rank_shift_summary.csv`, `delta_by_bin.csv` |
| Fig. S8(a) – consensus vs. weighted aggregation, 2-measure weight sweep (Sec. S4.1) | `S4_spatial_coexpression/S4a_weight_sweep_2measures.R` | `S4a_weight_sweep_2measures/2m_ari_<A>_<B>.jpg` |
| Fig. S8(b) – consensus vs. equal-weight aggregation, 2–5 measures (Sec. S4.1) | `S4_spatial_coexpression/S4b_consensus_vs_equal_weight.R` | `S4b_consensus_vs_equal_weight/all_ari.jpg` |
| Fig. S9 – ARI of consensus clustering by number of measures (Sec. S4.1) | `S4_spatial_coexpression/S4c_measure_combination_selection.R` | `S4c_measure_combination_selection/consensus_analysis_by_n_measures.jpg`, `consensus_combinations.csv` |
| Table S4 – neighbourhood size (Sec. S4.1) | `S4_spatial_coexpression/S4d_neighborhood_size_selection.R` | `S4d_neighborhood_size_selection/neighbor_comparison_4measures.csv` (Welch tests printed to the console) |
| Main-text Fig. 2(a, b) – accuracy and stability of clustering methods; Fig. S10 – run time (main-text Sec. 4.1, Sec. S4.2) | `S4_spatial_coexpression/S4e_benchmark_clustering_methods.R` | `S4e_benchmark_clustering_methods/ari_comparison.jpg`, `stability_comparison.jpg`, `time_comparison.jpg` |
| Fig. S11 – BANYAN tissue domains (Sec. S5) | not a simulation; the fits are in `data/real/banyan/` (see `00_data/00a_prepare_real_data.R`) | — |

## Repository layout

```
.
├── app/                                 VeteranST web application (see app/README.md)
│   ├── app.R, help_content.R
│   ├── install_app_dependencies.R
│   └── test_data/                       example data of the application (P16)
├── R/                                   shared functions (sourced by the scripts)
│   ├── utils_data.R                     dataset registry, loaders, BANYAN domain sizes,
│   │                                    per-gene permutation (null model of S3)
│   ├── utils_spatial_similarity.R       smoothing, 5 similarity measures, Leiden,
│   │                                    strict consensus (proposed), weighted aggregation
│   └── utils_evaluation.R               ground truth, ARI/NMI, stability, figure saving
├── 00_data/
│   ├── 00a_prepare_real_data.R          provenance + validation of data/real/
│   ├── 00b_generate_simulated_data.R    scDesign3 benchmark -> data/simulated/*.rds
│   └── 00c_export_app_test_input.R      real dataset -> CSV input for the web app
├── S1_gene_level_type1_error/           Section S1  (Figs S1–S4)
├── S3_domain_variability_type1_error/   Section S3  (Tables S1–S3, Figs S5–S7)
├── S4_spatial_coexpression/             Section S4  (Figs S8–S10, Table S4; main-text Fig. 2)
├── data/
│   ├── real/                            count_mat_*.rda, clusters_*.csv, coordinates_*.csv,
│   │                                    banyan/, annotations/ (see data/real/README.md)
│   └── simulated/                       count_mat_simul.rds, coordinates_simul.rds,
│                                        ground_truth.rds, source/brain_data.rda (see data/simulated/README.md)
├── results/                             all outputs are written here (git-ignored)
├── install_dependencies.R
└── LICENSE
```

## Requirements (simulation scripts)

* R ≥ 4.1
* CRAN: `Matrix`, `igraph` (≥ 1.3.0 for `cluster_leiden(resolution=)`), `RANN`,
  `aricode`, `cluster`, `matrixStats`, `ggplot2`, `dplyr`, `tidyr`, `scales`
* Only for re-deriving the committed data: `Seurat` (00a), Bioconductor
  `SingleCellExperiment` and `scDesign3` (GitHub `SONGDONGYUAN1994/scDesign3`) (00b)

`Rscript install_dependencies.R` installs everything.
Package versions used for the manuscript: **TODO – paste `sessionInfo()`**.

## Data

### Real datasets (Sections S1 and S3)

Three public 10x Visium datasets (see *Data availability* in the main text).
The exact objects used for the analyses are **committed in `data/real/`**
(≈ 65 MB), so Sections S1 and S3 run out of the box:

| Key | Dataset | File (object) | Genes × spots | BANYAN domains (spots per domain) |
|---|---|---|---|---|
| `P16` | melanoma brain metastasis, Sudmeier et al. 2022 (GEO) | `count_mat_p16.rda` (`count_mat_p16`) | 16,074 × 1,650 | 311, 609, 321, 409 |
| `BR`  | human breast cancer (IDC), 10x Genomics 2019 | `count_mat_br.rda` (`count_mat_br`) | 21,129 × 3,798 | 861, 665, 942, 741, 589 |
| `D14` | chicken heart day 14, Mantri et al. 2021 (GEO) | `count_mat_d14.rda` (`count_mat_d14`) | 11,356 × 1,967 | 754, 392, 240, 581 |

Each matrix is a `dgCMatrix` of SCTransform-normalised log values
(`<seurat>@assays$SCT@data`); `clusters_<key>.csv` holds the BANYAN domain label
of every spot (`<fit>$z`) and `coordinates_<key>.csv` the tissue coordinates
(`Seurat::GetTissueCoordinates()`), all in the spot order of the matrix
columns. The null simulations permute expression across these columns with
fixed seeds and use the labels only through their sizes, so the row/column
order of the committed objects must be preserved. Provenance, the BANYAN fit
objects and the breast-cancer region annotation used in the main text are
described in `data/real/README.md`; `00_data/00a_prepare_real_data.R` re-derives
and verifies the files from the processed Seurat objects, which are archived on
Zenodo (**TODO: DOI**) because of their size. `00_data/00c_export_app_test_input.R`
exports a dataset as CSV in the layout accepted by the web application.

### Simulated benchmark (Section S4)

The exact benchmark data used for the manuscript are **committed in
`data/simulated/`** (see `data/simulated/README.md`): `count_mat_simul.rds`
(100 genes × 2,696 spots, log1p counts; genes named `<gene>_<replicate>`),
`coordinates_simul.rds` (image coordinates of the spots) and `ground_truth.rds`
(the true module of every simulated gene = its parent gene; ten modules of ten
genes). `00_data/00b_generate_simulated_data.R` documents how they were
generated with scDesign3 (Song et al., 2024) from the first ten genes of the
scDesign3 Visium mouse-brain example object (figshare 40582019; the underlying
Seurat object is kept in `data/simulated/source/brain_data.rda` for
provenance) and can regenerate them for comparison, but the S4 scripts run out
of the box.

## Running

All scripts are run **from the repository root**, e.g.

```bash
Rscript S3_domain_variability_type1_error/S3a_type1_error_domain_variability.R
```

All input data are committed, so the scripts are independent of each other and
can be run in any order (`00a`/`00b` only validate or regenerate the inputs).
Every script has a parameter block at the top (number of runs, cores, seeds,
output folder) documented in its header comment; nothing needs to be edited
for a default run. Approximate run times on a 20-core workstation:

| Script | Runs | Time |
|---|---|---|
| S1 | 100 shuffles × all genes × 9 tests | hours (set `N_GENES` for a smoke test) |
| S3a | Table S1: 1,000 runs (`submitted` seeds); Table S2: 1,000 runs (`unified` seeds) | ≈ 1 h + 20 min |
| S3b | closed-form | < 1 min |
| S3c | 100 per-gene permutations | ≈ 25 min (single core) |
| S3d | 1,000 per-gene permutations | ≈ 4 h (single core; datasets can be run in parallel sessions) |
| S3e | 1,000 per-gene permutations | ≈ 10 h on 8 cores (BR dominates) |
| S4a–S4d | 100 Leiden runs per configuration | 10–30 min each |
| S4e | 100 runs × 7 methods | 10–20 min |

## Method summary (what the shared code does)

* `permute_genes()` – the null model of Section S3: every gene is permuted
  across all spots independently, and the domains are consecutive blocks of the
  permuted spots (a single permutation shared by all genes would preserve
  spot-level effects such as sequencing depth and is not used).
* `compute_smoothed_features()` – local mean / local variance of each gene over
  the *k* nearest spots (*k* = 7, 19, 37 or 61 including the spot itself).
* `compute_similarity_matrices()` – the five similarity measures of Section 2.3:
  Zero Jaccard, Hotspot Jaccard, High-/Low-variability Jaccard, Pearson on the
  smoothed mean (negative correlations are clipped to 0).
* `consensus_clustering()` – **proposed method**: Leiden (modularity) on each
  selected similarity graph, then strict consensus (genes are grouped only if
  co-clustered in every measure). Default: PM + JH + JVT + JVB, *k* = 37.
* `weighted_clustering()` – baseline: min-max-normalised matrices combined with
  weights, single Leiden run.

The web application applies the same measures on an HNSW/disparity-filter
sparsified graph so that genome-scale gene sets are tractable; the scripts here
use dense matrices, which is sufficient for the 100-gene benchmark.

## Reproducibility notes

* Random-number generation is pinned in every script (`set_rng()` in
  `R/utils_data.R`: Mersenne-Twister / inversion / rejection sampling), and
  every stochastic step (permutations, Leiden, k-means, subsampling) is seeded
  from the parameter block of the script. Re-running a script therefore gives
  identical numbers on the same package versions.
* Leiden community detection is stochastic, and on the 100-gene benchmark a
  single run can move the ARI of one configuration by about ±0.1. `S4a`–`S4d`
  therefore average `N_REPLICATES = 100` seeded Leiden runs per configuration
  and report the mean ± SD (per-replicate values are written to
  `*_replicates.csv`); the selections of Section S4.1 (measure set in `S4c`,
  neighbourhood size in `S4d`) are made on the means, and `S4d` prints the
  Welch tests quoted in the text. The method defaults (`DEFAULT_MEASURE_IDX`,
  `DEFAULT_N_NEIGHBORS` in `R/utils_spatial_similarity.R`) are the outcome of
  these two experiments; `S4e` and the web application use them.
* Table S1 combines two development runs (Friedman on IQR with seeds
  `123 + s`; the four variance-based tests with seeds `12345 + 10000 k + s`);
  Table S2 applies all five tests to the same 1,000 shuffles (seeds
  `123 + s`). `S3a` replays exactly these settings
  (`SCENARIO_SETTINGS`), so both tables are reproduced bit for bit provided
  the input matrices have the row/column order of the committed objects.
  Either seed scheme can be used for either scenario; they differ only
  within Monte Carlo error.
* `S3c`–`S3e` seed the per-gene permutations with `123` (`S3c`, `S3d`: one
  stream per dataset) and `123 + s` (`S3e`: one seed per permutation), so the
  figures and the numbers quoted in Section S3.1 are reproduced exactly.
* NMI is computed alongside ARI and stored in the CSV outputs; only ARI is
  shown in the figures.

## Provenance

The scripts consolidate the development code used for the submission (one
script per figure/table instead of several overlapping files). The statistical
definitions, parameter values, seeds where they existed, and output file names
were kept unchanged; only the code organisation, parallelisation and plotting
boilerplate differ. Exploratory analyses that preceded the final design and
are not reported in the manuscript are not included.

## License and citation

The application and the simulation code are released under the MIT License (see
`LICENSE`). If you use this code, please cite the VeteranST manuscript above and
the Zenodo archive (**DOI TODO**).

Contact: chung.911@osu.edu, sunyoungshin@postech.ac.kr
