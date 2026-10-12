# Real datasets (Supplementary Sections S1 and S3)

Three public 10x Visium datasets. The files here are the exact objects used
for the submitted analyses; the simulation scripts read them directly and no
preprocessing has to be re-run.

| Key | Dataset | File | Object | Genes × spots | Values |
|---|---|---|---|---|---|
| `P16` | human melanoma brain metastasis, patient 16 (Sudmeier et al., 2022; GEO) | `count_mat_p16.rda` | `count_mat_p16` | 16,074 × 1,650 | SCT log-normalised |
| `D14` | chicken heart, day 14 (Mantri et al., 2021; GEO) | `count_mat_d14.rda` | `count_mat_d14` | 11,356 × 1,967 | SCT log-normalised |
| `BR`  | human breast cancer, IDC (10x Genomics, 2019) | `count_mat_br.rda` | `count_mat_br` | 21,129 × 3,798 | SCT log-normalised |

* `count_mat_<key>.rda` — `Matrix::dgCMatrix` (genes × spots) taken from the
  processed Seurat object as `<seurat>@assays[["SCT"]]@data`, i.e. the
  log-normalised values of `SCTransform()`. Row and column order are those of
  the Seurat object and must not be changed: the null simulations permute
  these columns with fixed seeds, so the order is part of the reproducibility.
* `clusters_<key>.csv` — BANYAN tissue-domain label of every spot
  (`<banyan_fit>[["z"]]`, K = 4 / 4 / 5), in the same spot order. Domain sizes:
  P16 311 / 609 / 321 / 409; D14 754 / 392 / 240 / 581;
  BR 861 / 665 / 942 / 741 / 589 (these are the `banyan_sizes` in
  `R/utils_data.R`). The null simulations use only these sizes.
* `coordinates_<key>.csv` — `Seurat::GetTissueCoordinates()` of the same
  object (P16: image pixel coordinates `imagerow`/`imagecol`; D14 and BR:
  array coordinates `x`/`y`), same spot order.
* `banyan/pt16_fit4.rda`, `banyan/D14_banyan_fit4.rda`,
  `banyan/breast_fit_banyan.RData` — the BANYAN fit objects (`pt16_fit4`,
  `D14_fit4`, `breast_fit`; class `banyan`, K = 4 / 4 / 5, 900 posterior
  samples; `$z` is the label vector above, `$coords` the spot coordinates used
  by BANYAN, `$bic` the BIC of the fit).
* `annotations/clusters_br_annotation.csv` — the 20 manually annotated tissue
  regions of the breast-cancer section (Xu et al., 2024), transferred to the
  Visium spots with the region-selection tool of the VeteranST application.
  Used for the application in the main text, not for the simulations. Spot
  identifiers are R-mangled (`-1` written as `.1`).

Not committed (size): the processed Seurat objects (`pt16_norm.rda`,
`heart_D14_seurat_processed.RData`, `breast_processed_seurat.RData`) and the
CSV exports of the matrices. They are archived on Zenodo (**TODO: DOI**).
`00_data/00a_prepare_real_data.R` re-derives and verifies the committed files
from those objects.
