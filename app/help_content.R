## =============================================================================
## help_content.R -- VeteranST Help & Tutorial 탭 + STEP 별 "?" 팝업 콘텐츠
##
##   app.R 과 같은 폴더에 두고, app.R 에서  source("help_content.R", local = TRUE, encoding = "UTF-8")
##   로 불러오기.
##
##   정의되는 것
##     help_css          스타일
##     step_header()     사이드바의 "STEP n." 제목 + ? 링크  (기존 tags$h4 대체)
##     help_content      섹션별 title / body  (Help 탭과 ? 모달이 공유)
##     help_tab_ui()     Help & Tutorial tabPanel
##     help_modal()      ? 클릭 시 열리는 모달
##
##   용어: tissue domain = "Domain", gene module = "Module" (app.R 의 UI 문자열과 동일)
##   ★ 표시된 곳은 제출 전에 채워야 하는 값 (URL, DOI, runtime, 버전).
## =============================================================================

help_version   <- "1.0.0"                 # ★ 앱 버전 (Zenodo release 와 맞추기)
help_github    <- "https://github.com/MATH16KJH/VeteranST"
help_zenodo    <- "https://doi.org/10.5281/zenodo.XXXXX" # ★
help_runtime   <- "approximately XX minutes on the public server (16,074 genes x 1,650 spots)"  # ★

## -----------------------------------------------------------------------------
## CSS
## -----------------------------------------------------------------------------
help_css <- "
.help-link { margin-left: 8px; font-size: 15px; color: #95a5a6; vertical-align: middle; }
.help-link:hover { color: #2980b9; text-decoration: none; }
.help-container { padding: 10px 30px 40px 30px; background-color: #fff; border-radius: 8px;
                  box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
.help-toc { position: sticky; top: 15px; font-size: 13px; padding-right: 10px; }
.help-toc h4 { color: #2980b9; font-size: 15px; margin-bottom: 8px; }
.help-toc a { display: block; padding: 3px 0; color: #2c3e50; }
.help-toc a:hover { color: #2980b9; text-decoration: none; }
.help-section { margin-bottom: 30px; }
.help-section h3 { color: #2980b9; font-size: 20px; border-bottom: 2px solid #e8f4f8;
                   padding-bottom: 6px; margin-top: 34px; }
.help-section h4 { color: #2c3e50; font-size: 15px; margin-top: 20px; font-weight: bold; }
.help-section p, .help-section li, .help-section td { font-size: 14px; line-height: 1.55; color: #2c3e50; }
.help-table { font-size: 13px; margin-top: 8px; }
.help-table th { background-color: #2980b9; color: #fff; font-weight: bold; }
.help-note { background-color: #fff8e1; border-left: 4px solid #f39c12; padding: 8px 12px;
             margin: 12px 0; font-size: 13.5px; }
.help-tip  { background-color: #e8f4f8; border-left: 4px solid #2980b9; padding: 8px 12px;
             margin: 12px 0; font-size: 13.5px; }
.help-mono { font-family: Menlo, Consolas, monospace; font-size: 12.5px; background: #f4f6f7;
             padding: 8px 12px; border-radius: 4px; white-space: pre; overflow-x: auto; margin: 6px 0 12px 0; }
.modal-lg { width: 860px; }
.modal-body p, .modal-body li, .modal-body td { font-size: 14px; line-height: 1.5; }
.modal-body h4 { color: #2c3e50; font-size: 15px; font-weight: bold; margin-top: 18px; }
"

## -----------------------------------------------------------------------------
## 작은 헬퍼
## -----------------------------------------------------------------------------
step_header <- function(text, step_id) {
  tags$h4(
    style = "color: #2980b9; font-size: 18px; margin-bottom: 10px;",
    text,
    actionLink(paste0("help_", step_id), label = NULL,
               icon = icon("question-circle"), class = "help-link",
               title = "Show help for this step")
  )
}

help_table <- function(df) {
  tags$table(
    class = "table table-condensed table-bordered help-table",
    tags$thead(tags$tr(lapply(names(df), function(n) tags$th(n)))),
    tags$tbody(lapply(seq_len(nrow(df)), function(i) {
      tags$tr(lapply(seq_along(df), function(j) tags$td(HTML(as.character(df[i, j])))))
    }))
  )
}

help_note <- function(...) div(class = "help-note", ...)
help_tip  <- function(...) div(class = "help-tip", ...)
help_mono <- function(text) div(class = "help-mono", text)

## -----------------------------------------------------------------------------
## 콘텐츠
## -----------------------------------------------------------------------------
help_content <- list()

## ---- Overview ---------------------------------------------------------------
help_content$overview <- list(
  title = "Overview",
  body = tagList(
    p("VeteranST (", strong("V"), "ariability H", strong("eter"), "ogeneity ", strong("An"), "alysis in ",
      strong("S"), "patial ", strong("T"), "ranscriptomics) tests whether transcriptional ", em("variability"),
      " - rather than mean expression - differs between tissue domains, and detects spatial gene modules ",
      "organised by that variability. The whole pipeline runs in this browser session: no installation, ",
      "registration or scripting is needed."),
    h4("Workflow (one tab per step)"),
    tags$ol(
      tags$li(strong("Data Upload"), " - expression matrix, spot coordinates and, optionally, tissue-domain labels; ",
              "or the built-in example data set."),
      tags$li(strong("Tissue Domain Assignment"), " - paint domains directly on the spot map, or keep the uploaded labels."),
      tags$li(strong("Test for DVGs"), " - Brown-Forsythe test for unequal expression variability across domains ",
              "(family-wise error rate controlled by Bonferroni correction)."),
      tags$li(strong("Spatial Domain Variability Testing"), " - IQR-based Friedman test and Nemenyi post-hoc ",
              "comparison that rank domains by transcriptional variability."),
      tags$li(strong("Spatial Gene Co-expression Analysis"), " - graph-based consensus clustering of DVGs into ",
              "spatially co-expressed gene modules.")
    ),
    p("Two words are used consistently throughout the app: a ", strong("domain"), " is a region of tissue (STEP 2-4), ",
      "and a ", strong("module"), " is a group of genes with a shared spatial expression pattern (STEP 5)."),
    p("Results propagate downstream: changing the domain assignment, the domain criterion or the FWER threshold ",
      "clears the results of STEP 3-5 so that they can be rerun consistently."),
    help_tip(strong("Quick start. "), "Data Upload → ", em("Load Example Data"), " → Test for DVGs → ",
             em("Test Start!"), " → Spatial Domain Variability Testing → ", em("Run Analysis"),
             " → Spatial Gene Co-expression Analysis → ", em("Run Clustering"), ". ",
             "The example is the melanoma brain metastasis section (P16) analysed in the paper, with its BANYAN ",
             "tissue domains preselected as the domain criterion.")
  )
)

## ---- Input files ------------------------------------------------------------
help_content$inputs <- list(
  title = "Input files",
  body = tagList(
    p("Three comma-separated (.csv) files are accepted. Uploads of up to 1 GB are allowed; a 100 MB expression ",
      "matrix takes roughly half a minute to load."),
    help_table(data.frame(
      File = c("Expression matrix", "Spot coordinates", "Tissue domain labels"),
      Required = c("Yes", "Yes", "No (optional)"),
      Layout = c("Rows = genes, columns = spots. First column = gene symbol, header row = spot IDs (barcodes).",
                 "Three columns: spot ID, X, Y. Only the first three columns are used; column names are ignored.",
                 "Two columns: spot ID, domain label (integer). Column names are ignored. Label <code>0</code> = spot excluded from analysis."),
      Notes = c("Any expression scale can be uploaded (raw counts, log-normalised counts, SCTransform values, ...); VeteranST does not normalise. The paper used SCTransform-normalised log values. Gene symbols must be unique.",
                "Numeric pixel or array coordinates. Rows with missing X/Y are dropped, which then breaks the spot count - remove them beforehand.",
                "Same number of rows as spots. Labels are used as-is; any integer coding is fine."),
      stringsAsFactors = FALSE)),
    help_note(strong("Order matters. "), "Rows of the coordinate file and of the label file must be in the same ",
              "order as the columns of the expression matrix (the app matches spots by position, and checks only ",
              "that the numbers agree)."),
    h4("Example: expression matrix"),
    help_mono('"",AAACACCAATAACTGC-1,AAACAGCTTTCAGAAG-1,AAACAGGGTCTATATT-1, ...\nLINC01409,0.693,0,0, ...\nLINC01128,0,0,0.511, ...\nTAP1,1.204,0,0.842, ...'),
    h4("Example: spot coordinates"),
    help_mono("Spot,X,Y\nAAACACCAATAACTGC-1,4812.3,2075.9\nAAACAGCTTTCAGAAG-1,3327.8,6412.0"),
    h4("Example: tissue domain labels"),
    help_mono("Spot,Domain\nAAACACCAATAACTGC-1,2\nAAACAGCTTTCAGAAG-1,2\nAAACAGGGTCTATATT-1,0"),
    p("Labels can come from any domain-detection method (BANYAN, BayesSpace, Seurat clusters, pathologist ",
      "annotation, ...). The file exported by ", em("Download Domain Assignments"), " in STEP 2 has exactly this format ",
      "and can be re-uploaded later.")
  )
)

## ---- STEP 1 ---------------------------------------------------------------
help_content$step1 <- list(
  title = "STEP 1. Data Upload",
  body = tagList(
    p("Upload the expression matrix and the coordinate file (formats in ", em("Input files"), "). ",
      "Once both are loaded and their sizes agree, the ", strong("Tissue Spots"), " map is drawn and the optional ",
      "label upload (", em("Predefined Tissue Domain"), ") is enabled. Spots are coloured by label if a label file is provided."),
    p(strong("Load Example Data"), " loads the melanoma brain metastasis section (P16; Sudmeier ", em("et al."),
      ", 2022; 16,074 genes x 1,650 spots) together with its four BANYAN tissue domains, and switches the domain ",
      "criterion in STEP 3 to ", em("Predefined Tissue Domain"), ". It replaces anything uploaded before."),
    help_note("If the number of coordinate rows differs from the number of expression columns, an error is shown ",
              "and all downstream buttons stay disabled (hover a disabled button to see why).")
  )
)

## ---- STEP 2 ---------------------------------------------------------------
help_content$step2 <- list(
  title = "STEP 2. Tissue Domain Assignment",
  body = tagList(
    p("Use this tab to define tissue domains by hand. If you uploaded a label file and want to use it as-is, ",
      "skip this tab and choose ", em("Predefined Tissue Domain"), " in STEP 3."),
    tags$ol(
      tags$li(strong("Start!"), " creates Domain 1. Drag on the map to paint spots; only unassigned spots (black outline) are painted."),
      tags$li(strong("New Domain"), " creates the next domain. A domain must contain at least one spot before a new one can be added."),
      tags$li(strong("Erase / Brush"), " toggles eraser mode, which removes spots from the ", em("current"), " domain only."),
      tags$li(strong("Brush/Eraser Size"), " sets the radius of the brush (1-50 px)."),
      tags$li(strong("Set the rest as a new domain"), " assigns every still-unassigned spot to a new domain in one click."),
      tags$li(strong("Initialize"), " clears all assignments."),
      tags$li(strong("Download Domain Assignments"), " exports the assignment as a two-column CSV (Spot, Domain) that can be uploaded as a label file later.")
    ),
    p("The ", strong("Domain Status"), " cards show the spot count per domain; clicking a card briefly highlights that domain on the map, ",
      "and hovering a spot shows its domain. ",
      "Unassigned spots (domain 0) are excluded from all tests. At least two domains and at least 10 assigned spots are required."),
    help_note("Editing the assignment clears the results of STEP 3-5.")
  )
)

## ---- STEP 3 ---------------------------------------------------------------
help_content$step3 <- list(
  title = "STEP 3. Test for Differentially Variable Genes (DVGs)",
  body = tagList(
    h4("Settings"),
    tags$ul(
      tags$li(strong("FWER Threshold"), " (default 0.01; any real number strictly between 0 and 1): a gene is called a DVG when its Bonferroni-adjusted p-value is at or below this value."),
      tags$li(strong("Tissue Domain Criteria"), ": ", em("Brush-assigned Tissue Domain"), " uses the domains painted in STEP 2; ",
              em("Predefined Tissue Domain"), " uses the uploaded label file. This choice also governs STEP 4 and 5.")
    ),
    h4("What is computed"),
    p("For every gene, the ", strong("Brown-Forsythe test"), " compares expression variability across domains using ",
      "absolute deviations from the domain medians (one-way F statistic on those deviations). It is robust to the ",
      "zero-inflated, non-normal distributions typical of spatial transcriptomics. In parallel, a one-way ANOVA tests ",
      "differences in mean expression at the same FWER."),
    tags$ul(
      tags$li(strong("DVGs"), " tab: all genes with significantly different variability, ranked by p-value (DVG Rank 1 = smallest p)."),
      tags$li(strong("DVG-only genes"), " tab: DVGs whose ", em("mean"), " expression is not significantly different between domains. ",
              "These carry heterogeneity information that mean-based differential expression would miss.")
    ),
    h4("Gene detail panel"),
    p("Click a gene card to see the domain map and its spatial expression map (viridis scale: purple = low, yellow = high), ",
      "the raw and adjusted Brown-Forsythe p-values, the mean and variance per domain, and a link to its GeneCards entry. ",
      "Point size and transparency only affect the display. ", em("Download Gene Analysis File"), " packs the maps (PNG) ",
      "and the table (XLSX) into a ZIP."),
    h4("Downloads"),
    p("The download button under each list exports the corresponding gene table as CSV with columns ",
      code("Gene"), ", ", code("P_value"), ", ", code("Adjusted P-value"), " and ", code("DVG_only"), " (TRUE/FALSE)."),
    help_note("Genes with zero variance across all spots are skipped (p-value NA). The test can take a few minutes for ",
              "a full transcriptome; the progress bar reports the stage.")
  )
)

## ---- STEP 4 ---------------------------------------------------------------
help_content$step4 <- list(
  title = "STEP 4. Spatial Domain Variability Testing",
  body = tagList(
    p("This step asks which tissue domains are the most transcriptionally variable, using the DVGs from STEP 3. ",
      "Press ", strong("Run Analysis"), " (available once STEP 3 has produced at least one DVG)."),
    h4("Settings"),
    tags$ul(
      tags$li(strong("Significance level (Friedman / Nemenyi)"), " (default 0.01; any real number strictly between 0 and 1): the level α used both to decide ",
              "whether the Friedman test is significant and, if it is, to compute the Nemenyi critical difference. ",
              "After changing it, press ", em("Run Analysis"), " again.")
    ),
    h4("What is computed"),
    tags$ol(
      tags$li("For each DVG and each domain, the interquartile range (IQR) of expression over the spots of that domain."),
      tags$li("Within each gene, domains are ranked by IQR (rank 1 = highest variability; ties receive the average rank)."),
      tags$li(strong("Friedman test"), " on the rank matrix (genes = blocks, domains = treatments) tests whether variability ranks differ systematically across domains. The chi-squared statistic, degrees of freedom and p-value are shown in the sidebar; the p-value is printed in red when it is below α."),
      tags$li("If p < α, the ", strong("Nemenyi post-hoc test"), " compares all pairs of domains: two domains differ significantly when the difference of their mean ranks exceeds the critical difference CD = q(1−α; k, ∞) · √(k(k+1)/(12N)), with k domains and N genes.")
    ),
    h4("Reading the output"),
    tags$ul(
      tags$li(strong("Tissue Domain Variability map"), ": spots are coloured by the mean rank of their domain (yellow = most variable, purple = least). Hover a spot to see its domain and rank; click a card in ", em("Variability Ranking"), " to highlight that domain."),
      tags$li(strong("Variability Ranking"), " cards: overall rank and mean rank of each domain (lower mean rank = higher variability)."),
      tags$li(strong("Critical Difference Diagram"), ": domains are placed on the ", em("Mean Rank"), " axis, most variable on the left. The short bar labelled CD shows the critical difference. ",
              strong("Red horizontal bars connect domains that are not significantly different"), " from each other; domains not joined by a red bar differ significantly at the chosen α. The diagram is drawn only when the Friedman test is significant.")
    ),
    help_note(strong("Interpretation caveat. "), "Because the test is run on genes that were themselves selected for ",
              "differential variability, the Friedman p-value is a conditional summary of how systematically variability ",
              "ranks differ within the selected gene set, not an unconditional transcriptome-wide test.")
  )
)

## ---- STEP 5 ---------------------------------------------------------------
help_content$step5 <- list(
  title = "STEP 5. Spatial Gene Co-expression Analysis",
  body = tagList(
    h4("Settings"),
    tags$ul(
      tags$li(strong("Number of Top DVGs"), " (default 100, or all DVGs if fewer): how many DVGs, in order of DVG Rank, enter the clustering. The paper used 1,000 for the P16 data and 5,000 for the breast cancer data."),
      tags$li(strong("Min Module Size"), " (default 5): modules with fewer genes are discarded.")
    ),
    h4("What is computed"),
    p("Genes with similar spatial expression patterns are grouped into modules by graph-based consensus clustering:"),
    tags$ol(
      tags$li("Expression of each gene is smoothed over spatial neighbourhoods of 7, 19, 37 and 61 nearest spots, giving local mean, local variance, local zero proportion and a local top-10% indicator at each size; with raw expression, a zero mask and a top-10% mask this yields 19 feature spaces."),
      tags$li("In each feature space an approximate k-nearest-neighbour graph over genes is built with HNSW (cosine distance). The 19 graphs are summed into a support-weighted gene graph and sparsified with a disparity filter (α = 0.05)."),
      tags$li("On the remaining edges four similarity measures are computed from the 37-spot smoothed data: Pearson correlation of local means, and Jaccard similarities of hotspot spots (top 10% local mean), high-variability spots (top 10% local variance) and low-variability spots (bottom 10% local variance)."),
      tags$li("The Leiden algorithm is run on each of the four similarity graphs; two genes end up in the same module only if they are co-clustered in ", strong("all four"), " runs. Modules smaller than ", em("Min Module Size"), " are removed.")
    ),
    h4("Reading the output"),
    tags$ul(
      tags$li(strong("Module network"), ": one node per module, size proportional to the number of genes. Click a node (or pick it from the drop-down) to list its genes on the left under ", em("Genes in Module n"), ", with their DVG rank and DVG-only rank."),
      tags$li(strong("Find Gene's Module"), ": type a gene symbol to find the module it belongs to; click the module name to open it."),
      tags$li("Click a gene card to show its spatial expression map and per-domain statistics."),
      tags$li(strong("Download Module Information"), " exports one Excel sheet per module (gene, p-values, DVG-only flag); ",
              em("Download Gene Analysis"), " exports the maps and table of the selected gene.")
    ),
    help_note("HNSW and Leiden contain random steps, so module boundaries can differ slightly between runs; the four-way ",
              "consensus keeps the core of each module stable. If no module survives the size filter, lower ",
              em("Min Module Size"), " or raise ", em("Number of Top DVGs"), ".")
  )
)

## ---- Parameters -------------------------------------------------------------
help_content$params <- list(
  title = "Parameter reference",
  body = tagList(
    h4("User-adjustable parameters"),
    help_table(data.frame(
      Parameter = c("FWER Threshold", "Tissue Domain Criteria", "Brush/Eraser Size",
                    "Significance level (Friedman / Nemenyi)",
                    "Number of Top DVGs", "Min Module Size", "Point Size / Point Transparency"),
      Tab = c("STEP 3", "STEP 3", "STEP 2", "STEP 4", "STEP 5", "STEP 5", "STEP 3, 5"),
      Default = c("0.01", "Brush-assigned", "20", "0.01", "100 (all if fewer DVGs)", "5", "0.2 / 0.8"),
      Range = c("0 &lt; value &lt; 1", "Brush-assigned / Predefined", "1 - 50", "0 &lt; value &lt; 1", "0 - number of DVGs", "1 - 100", "0.1 - 1"),
      Meaning = c("Family-wise error rate for calling DVGs (Bonferroni); also used for the ANOVA that defines DVG-only genes.",
                  "Which domain labels are used in STEP 3-5.",
                  "Radius of the painting brush, in pixels.",
                  "&alpha; for the Friedman test and for the Nemenyi critical difference; the post-hoc test runs only when the Friedman p-value is below &alpha;.",
                  "DVGs entering the co-expression analysis, taken in order of DVG Rank.",
                  "Modules with fewer genes are discarded.",
                  "Display only; do not affect any result."),
      stringsAsFactors = FALSE)),
    h4("Fixed settings"),
    help_table(data.frame(
      Component = c("DVG test", "DVG test", "DVG test", "DVG test",
                    "Domain variability", "Domain variability",
                    "Co-expression", "Co-expression", "Co-expression", "Co-expression", "Co-expression", "Co-expression", "Co-expression"),
      Setting = c("Variability statistic", "Multiple testing", "Mean-expression test (DVG-only definition)", "Minimum data",
                  "Variability statistic", "Nemenyi critical difference",
                  "Spatial neighbourhoods", "Feature spaces for the kNN graph", "HNSW index", "Disparity filter",
                  "Similarity measures", "Community detection", "Consensus rule"),
      Value = c("Brown-Forsythe test (Levene-type test on absolute deviations from the domain median)",
                "Bonferroni",
                "One-way ANOVA, Bonferroni, same FWER threshold",
                "&ge; 10 assigned spots and &ge; 2 domains; zero-variance genes skipped",
                "IQR of expression per domain; rank 1 = largest IQR; ties = average rank",
                "CD = q(1&minus;&alpha;; k, &infin;) &middot; &radic;(k(k+1) / (12N)), k = domains, N = DVGs, &alpha; = the significance level chosen in STEP 4",
                "7, 19, 37, 61 nearest spots (including the spot itself)",
                "19 = (local mean, local variance, local zero proportion, local top-10% indicator) &times; 4 sizes + raw expression, zero mask, top-10% mask",
                "cosine distance, M = 32, ef<sub>construction</sub> = 200, k = &lceil;&radic;N&rceil;, ef<sub>search</sub> = max(2k, 200)",
                "&alpha; = 0.05 on kNN-support weights (0-19)",
                "Pearson (37-spot local mean); Hotspot Jaccard (top 10% local mean); High-variability Jaccard (top 10% local variance); Low-variability Jaccard (bottom 10% local variance)",
                "Leiden, modularity objective, resolution 1, on each similarity graph",
                "Two genes are linked only if co-clustered in all four graphs; modules = connected components"),
      stringsAsFactors = FALSE))
  )
)

## ---- Reproducing the paper --------------------------------------------------
help_content$reproduce <- list(
  title = "Reproducing the paper results with the example data",
  body = tagList(
    tags$ol(
      tags$li("Data Upload → ", strong("Load Example Data"), ". The four BANYAN domains of the P16 section are loaded: Domain 2 = tumour, Domains 3 and 4 = peritumoral inflamed tissue, Domain 1 = blood."),
      tags$li("Test for DVGs → keep FWER Threshold = 0.01 and ", em("Predefined Tissue Domain"), " → ", strong("Test Start!"),
              ". The DVG-only list contains 25 genes, including ", em("HLA-C"), ", ", em("TAP1"), " and ", em("PSME2"),
              " (Figure 3c of the paper)."),
      tags$li("Spatial Domain Variability Testing → keep the significance level at 0.01 → ", strong("Run Analysis"),
              ". The tumour domain (Domain 2) ranks first, followed by the inflamed domains (3, 4) and blood (1); the critical difference diagram corresponds to Figure 3b."),
      tags$li("Spatial Gene Co-expression Analysis → set ", strong("Number of Top DVGs = 1,000"), " and ", strong("Min Module Size = 5"), " → ",
              strong("Run Clustering"), ". Search ", em("TAP1"), " in ", em("Find Gene's Module"), ": its module contains ",
              em("PSME2"), ", ", em("HLA-DRA"), ", ", em("HLA-DRB1"), ", ", em("HLA-DPA1"), ", ", em("HLA-DPB1"), ", ",
              em("CD74"), ", ", em("IFI30"), ", ", em("TAPBP"), " and ", em("PSMB9"), " (Figure 4).")
    ),
    p("Expected running time of the full pipeline: ", help_runtime, "."),
    help_note("Module membership can differ marginally between runs because of the random steps in HNSW and Leiden; ",
              "the antigen-presentation core of the module is stable.")
  )
)

## ---- FAQ / troubleshooting --------------------------------------------------
help_content$faq <- list(
  title = "FAQ and troubleshooting",
  body = tagList(
    help_table(data.frame(
      Message = c("Number of spots (...) does not match number of samples in expression matrix (...)",
                  "Coordinates file must have at least 3 columns (Spot, X, Y)",
                  "Domain label file must contain exactly 2 columns",
                  "Domain label file rows (...) do not match the number of spots",
                  "Spot IDs in the domain label file do not match those in the coordinate file",
                  "Too few spots selected",
                  "At least 2 domains (excluding unassigned spots) are required / At least 2 domains are required for analysis",
                  "No DVGs available",
                  "Insufficient genes for clustering",
                  "No modules meet the minimum size requirement",
                  "No edges after disparity filter"),
      Cause_and_fix = c("Rows of the coordinate file &ne; columns of the expression matrix. Check for missing X/Y values, an extra index column, or spots filtered from one file only.",
                        "The coordinate file needs spot ID, X and Y as its first three columns.",
                        "The label file needs exactly spot ID and label.",
                        "The label file must list every spot, in the same order as the expression columns.",
                        "Warning only: the labels are still applied by position. Make sure the label file was made for this coordinate file.",
                        "Fewer than 10 spots are assigned to a domain. Paint more spots or use a label file.",
                        "Only one domain has spots. Create at least two domains (or check that the label file has more than one non-zero label).",
                        "STEP 3 has not been run, or no gene passed the FWER threshold. Rerun STEP 3, possibly with a larger threshold.",
                        "Fewer than two DVGs were selected. Increase <i>Number of Top DVGs</i>.",
                        "Every module is smaller than <i>Min Module Size</i>. Lower it, or increase <i>Number of Top DVGs</i>.",
                        "The gene set is too small for a stable graph. Increase <i>Number of Top DVGs</i>."),
      stringsAsFactors = FALSE)),
    h4("Are results deterministic?"),
    p("STEP 3 and 4 are fully deterministic. STEP 5 uses HNSW (approximate nearest neighbours) and Leiden, both of which ",
      "contain random steps; the four-way consensus makes modules stable, but small differences at module boundaries between runs are possible."),
    h4("Which normalisation should I use?"),
    p("VeteranST does not require a particular normalisation: the workflow runs on whatever expression values you upload, ",
      "whether raw counts, log-normalised counts (Seurat's ", code("NormalizeData"), ", scanpy's ", code("normalize_total"),
      " + ", code("log1p"), ") or model-based normalised values. All analyses in the paper, including the example data, used ",
      strong("SCTransform"), "-normalised values (Seurat's ", code("SCTransform()"), "; the ", code("data"), " slot of the ",
      code("SCT"), " assay, i.e. log1p of the corrected counts). Note that the choice affects what variability means - ",
      "on raw counts, differences in variability partly reflect differences in sequencing depth between spots - so choose ",
      "the scale that matches your question and report it with the results."),
    h4("Can I use my own domain detection?"),
    p("Yes - export the labels of any method to a two-column CSV (spot ID, label) and upload it as the tissue domain file, ",
      "then select ", em("Predefined Tissue Domain"), " in STEP 3."),
    h4("The app seems frozen"),
    p("Long steps (the DVG test on a full transcriptome, or co-expression on thousands of genes) run for several minutes; ",
      "the progress bar at the bottom right shows the current stage. If the browser loses the connection, reload the page and repeat the steps.")
  )
)

## ---- About ------------------------------------------------------------------
help_content$about <- list(
  title = "Citation, source code and contact",
  body = tagList(
    p(strong("Citation. "), "Kim J, Jung K, Shin J, Chung D, Shin S. VeteranST: Transcriptional Variability and Gene Co-Expression Analysis ",
      "in Spatial Transcriptomics. ", em("Bioinformatics"), " (submitted)."),     # ★ 게재 후 갱신
    p(strong("Source code and test data. "), tags$a(href = help_github, target = "_blank", help_github),
      " (open-source licence: MIT). The version described in the paper is archived at ",
      tags$a(href = help_zenodo, target = "_blank", help_zenodo), "."),
    p(strong("Example data. "), "Human melanoma brain metastasis, 10x Visium, sample P16 (Sudmeier ", em("et al."),
      ", 2022; GEO accession GSE179572). Tissue domains were detected with BANYAN (Xie ", em("et al."), ", 2024)."),
    p(strong("Contact and support. "), "Dongjun Chung (", tags$a(href = "mailto:chung.911@osu.edu", "chung.911@osu.edu"), ") and Sunyoung Shin (",
      tags$a(href = "mailto:sunyoungshin@postech.ac.kr", "sunyoungshin@postech.ac.kr"), "). ",
      "Bug reports and feature requests are welcome through the GitHub issue tracker. ",
      "The web server will be maintained for at least three years after publication."),
    p(style = "color: #7f8c8d; font-size: 12px;", paste0("VeteranST web application, version ", help_version, "."))
  )
)

## -----------------------------------------------------------------------------
## Help & Tutorial 탭
## -----------------------------------------------------------------------------
help_tab_ui <- function() {
  sections <- c(overview = "Overview", inputs = "Input files",
                step1 = "STEP 1. Data Upload", step2 = "STEP 2. Tissue Domain Assignment",
                step3 = "STEP 3. Test for DVGs", step4 = "STEP 4. Spatial Domain Variability",
                step5 = "STEP 5. Spatial Co-expression", params = "Parameter reference",
                reproduce = "Reproducing the paper", faq = "FAQ / troubleshooting",
                about = "Citation & contact")
  
  tabPanel(
    tags$span(icon("book"), " Help & Tutorial", class = "tab-panel-title"),
    value = "help_tab",
    div(
      class = "help-container",
      fluidRow(
        column(
          3,
          div(class = "help-toc",
              h4("Contents"),
              lapply(names(sections), function(s) tags$a(href = paste0("#help-", s), sections[[s]])))
        ),
        column(
          9,
          lapply(names(sections), function(s) {
            div(id = paste0("help-", s), class = "help-section",
                h3(help_content[[s]]$title),
                help_content[[s]]$body)
          })
        )
      )
    )
  )
}

## -----------------------------------------------------------------------------
## "?" 모달
## -----------------------------------------------------------------------------
help_modal <- function(section) {
  sec <- help_content[[section]]
  showModal(modalDialog(
    title = sec$title,
    sec$body,
    footer = tagList(
      tags$a(href = "#", class = "btn btn-default btn-sm",
             onclick = "Shiny.setInputValue('goto_help_tab', Date.now(), {priority: 'event'}); return false;",
             icon("book"), " Open full tutorial"),
      modalButton("Close")
    ),
    size = "l", easyClose = TRUE
  ))
}
