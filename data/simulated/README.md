# Simulated benchmark data (Supplementary Section S4)

The exact files used for the manuscript. Every script under
`S4_spatial_coexpression/` reads them directly; they never have to be
regenerated.

| File | Content |
|---|---|
| `count_mat_simul.rds` | numeric matrix, 100 genes × 2,696 spots, `log1p` of scDesign3 counts. Rows are named `<gene>_<replicate>` (`Calb2_1`, `Gng4_1`, …, `Gpsm1_10`); columns are the Visium barcodes of the reference section. |
| `coordinates_simul.rds` | data.frame, 2,696 × 2 (`spatial1`, `spatial2`; image pixel coordinates), row names = barcodes, same order as the matrix columns. |
| `ground_truth.rds` | data.frame with `Gene`, `Gene_Base_Name`, `True_Cluster` (1–10): the true module of a simulated gene is its parent gene, i.e. the ten replicates of `Calb2` form module 1, and so on. |
| `source/brain_data.rda` | Seurat object `brain` (10x Visium mouse brain, anterior section; `stxBrain` "anterior1"), the tissue section from which the scDesign3 example object (figshare 40582019) used as the simulation reference was built. Provenance only. |

Generation: `00_data/00b_generate_simulated_data.R` — scDesign3 (negative-
binomial marginals, Gaussian-process mean `s(spatial1, spatial2, bs = "gp",
k = 400)`, Gaussian copula) fitted to the first 10 genes of the reference
object (Calb2, Gng4, Ttr, S100a5, Nrgn, Doc2g, Ppp1r1b, Pcbp3, Fabp7, Gpsm1);
replicate *i* generated after `set.seed(123 + i)`, *i* = 1…10. Re-running the
script writes `*_regenerated.rds` next to the committed files and reports
whether they are equal instead of overwriting them.
