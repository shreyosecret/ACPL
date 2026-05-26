# ==============================================================================
# polar_vs_mst_buettner.R
#
# Third-dataset benchmark using Buettner et al. (2015) mouse embryonic stem
# cell data: N = 288 cells FACS-sorted into G1, S, and G2M, profiled by
# Smart-Seq. Available via Bioconductor scRNAseq::BuettnerESCData().
#
# Reference: Buettner F, Natarajan KN, et al. Computational analysis of
#   cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals
#   hidden subpopulations of cells. Nat Biotechnol. 2015;33(2):155-160.
#
# Run: Rscript polar_vs_mst_buettner.R
# Expected runtime: 5-10 minutes (after first-time package install).
# ==============================================================================

cat("\n=== Polar / ACPL / MST on Buettner mESC ===\n\n")

required_cran <- c("Seurat", "dplyr", "igraph", "boot", "Matrix", "BiocManager")
required_bioc <- c("scRNAseq", "SingleCellExperiment", "SummarizedExperiment")

missing_cran <- required_cran[!required_cran %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}
missing_bioc <- required_bioc[!required_bioc %in% installed.packages()[, "Package"]]
if (length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

suppressMessages({
  library(Seurat); library(dplyr); library(igraph); library(boot); library(Matrix)
  library(scRNAseq); library(SingleCellExperiment); library(SummarizedExperiment)
})

set.seed(42)

# ── Load Buettner data ───────────────────────────────────────────────────────
cat("Loading Buettner mESC data from Bioconductor...\n")
sce <- BuettnerESCData()
cat(sprintf("  Dimensions: %d genes x %d cells\n", nrow(sce), ncol(sce)))
cat("  Assays available:", paste(assayNames(sce), collapse = ", "), "\n")
cat("  colData columns:", paste(colnames(colData(sce)), collapse = ", "), "\n")

# Find phase column. Buettner uses 'phase' or 'Phase' in scRNAseq metadata.
cd <- as.data.frame(colData(sce))
phase_col <- NULL
for (col in c("phase", "Phase", "cell_cycle_phase", "stage")) {
  if (col %in% colnames(cd)) {
    phase_col <- col
    break
  }
}
if (is.null(phase_col)) {
  cat("  No phase column found. Available columns:\n")
  for (col in colnames(cd)) {
    vals <- head(unique(cd[[col]]), 10)
    cat(sprintf("    %s: %s\n", col, paste(vals, collapse = ", ")))
  }
  stop("Could not find phase column.")
}
cat("  Using phase column:", phase_col, "\n")
phases_raw <- as.character(cd[[phase_col]])
cat("  Phase distribution:\n")
print(table(phases_raw))

# Standardise phase labels to G1 / S / G2M
phase_map_in <- c("G1" = "G1", "S" = "S", "G2M" = "G2M", "G2/M" = "G2M",
                  "G2" = "G2M", "M" = "G2M")
phases <- phase_map_in[phases_raw]
if (any(is.na(phases))) {
  cat("  Warning: some phases unmapped. Setting to NA.\n")
}

# ── Build Seurat object ──────────────────────────────────────────────────────
cat("\nBuilding Seurat object...\n")
counts_assay <- if ("counts" %in% assayNames(sce)) "counts" else assayNames(sce)[1]
counts <- as.matrix(assay(sce, counts_assay))
# Buettner data from scRNAseq sometimes uses Ensembl rownames; that's fine since
# we don't need gene symbols for this analysis (no Seurat cell cycle scoring).
keep <- !is.na(phases)
counts <- counts[, keep]
phases <- phases[keep]
cat(sprintf("  Final: %d genes x %d cells with valid phases\n",
            nrow(counts), ncol(counts)))

obj <- CreateSeuratObject(counts = counts, project = "buettner_mesc")
# Detect whether the data is already log-normalised
max_val <- max(counts)
if (max_val < 50) {
  cat("  Data already normalised; copying to data slot.\n")
  obj[["RNA"]]$data <- counts
} else {
  obj <- NormalizeData(obj, verbose = FALSE)
}
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# Assign Phase metadata positionally (Seurat may rename cells, so we can't
# rely on name matching; the order matches since we filtered counts and
# phases together above).
phase_vec <- as.character(phases)
names(phase_vec) <- colnames(obj)
obj <- AddMetaData(obj, metadata = phase_vec, col.name = "Phase")
cat("  Phase assigned. Distribution in Seurat object:\n")
print(table(obj$Phase))

# ── Engine functions ─────────────────────────────────────────────────────────
phase_map_num <- c("G1" = 1L, "S" = 2L, "G2M" = 3L)

calc_swsa <- function(time_proxy, phases, window = 10) {
  df <- data.frame(t = time_proxy, Phase = as.character(phases)) %>%
    dplyr::filter(!is.na(Phase), Phase %in% c("G1", "S", "G2M")) %>%
    dplyr::mutate(PhaseNum = phase_map_num[Phase]) %>%
    dplyr::arrange(t)
  n <- nrow(df)
  if (n <= window) return(NA_real_)
  100 * mean(df$PhaseNum[1:(n - window)] <=
               df$PhaseNum[(window + 1):n], na.rm = TRUE)
}

find_adaptive_origin <- function(x, y) {
  optim(c(median(x), median(y)),
        function(p) var(sqrt((x - p[1])^2 + (y - p[2])^2)))$par
}
unwrap_theta <- function(theta) {
  d <- diff(theta); d <- d - 2 * pi * round(d / (2 * pi))
  c(theta[1], theta[1] + cumsum(d))
}

polar_alone <- function(coords, phases) {
  uv <- find_adaptive_origin(coords[, 1], coords[, 2])
  th <- atan2(coords[, 2] - uv[2], coords[, 1] - uv[1])
  max(calc_swsa(th, phases), calc_swsa(-th, phases), na.rm = TRUE)
}
acpl_full <- function(coords, phases) {
  uv <- find_adaptive_origin(coords[, 1], coords[, 2])
  th <- atan2(coords[, 2] - uv[2], coords[, 1] - uv[1])
  ord <- order(th)
  th_uw <- unwrap_theta(th[ord])
  arc <- c(0, cumsum(sqrt(diff(cos(th_uw))^2 + diff(sin(th_uw))^2)))
  arc_sm <- predict(loess(arc ~ th_uw, span = 0.25))
  final <- arc_sm[order(seq_along(th)[ord])]
  max(calc_swsa(final, phases), calc_swsa(-final, phases), na.rm = TRUE)
}
mst_method <- function(coords, phases) {
  dmat <- as.matrix(dist(coords))
  g <- graph_from_adjacency_matrix(dmat, weighted = TRUE, mode = "undirected")
  mst <- igraph::mst(g)
  root <- which(phases == "G1")[1]
  if (is.na(root)) root <- 1
  pt <- as.numeric(distances(mst, v = root, to = V(mst),
                             weights = E(mst)$weight))
  max(calc_swsa(pt, phases), calc_swsa(-pt, phases), na.rm = TRUE)
}

# ── 10 UMAP seeds ────────────────────────────────────────────────────────────
cat("\nRunning 10 UMAP seeds (Buettner is small; should be fast)...\n")
n_seeds <- 10
polar_res <- numeric(n_seeds)
acpl_res  <- numeric(n_seeds)
mst_res   <- numeric(n_seeds)
for (i in seq_len(n_seeds)) {
  cat(sprintf("  Seed %d/%d...\n", i, n_seeds))
  set.seed(40 + i)
  obj <- RunUMAP(obj,
                 dims = 1:min(20, ncol(Embeddings(obj, "pca"))),
                 n.neighbors = min(30, ncol(obj) - 1),
                 seed.use = 40 + i, verbose = FALSE)
  um <- Embeddings(obj, "umap")
  ph <- obj$Phase
  polar_res[i] <- polar_alone(um, ph)
  acpl_res[i]  <- acpl_full(um, ph)
  mst_res[i]   <- mst_method(um, ph)
}

# ── Statistics ───────────────────────────────────────────────────────────────
d_pm <- polar_res - mst_res
d_ap <- acpl_res - polar_res
wt_pm <- wilcox.test(polar_res, mst_res, paired = TRUE, exact = FALSE)
wt_ap <- wilcox.test(acpl_res, polar_res, paired = TRUE, exact = FALSE)
boot_pm <- boot(d_pm, function(d, i) mean(d[i]), R = 2000)
boot_ap <- boot(d_ap, function(d, i) mean(d[i]), R = 2000)
ci_pm <- tryCatch(boot.ci(boot_pm, type = "bca")$bca, error = function(e) NULL)
ci_ap <- tryCatch(boot.ci(boot_ap, type = "bca")$bca, error = function(e) NULL)

# ── Report ───────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("=== BUETTNER mESC RESULTS ===\n")
cat("========================================\n\n")
cat(sprintf("Cells: %d (FACS-sorted)\n", ncol(obj)))
cat(sprintf("\nPolar-alone: %.2f%% +/- %.2f\n", mean(polar_res), sd(polar_res)))
cat(sprintf("Full ACPL:   %.2f%% +/- %.2f\n", mean(acpl_res),  sd(acpl_res)))
cat(sprintf("MST:         %.2f%% +/- %.2f\n", mean(mst_res),   sd(mst_res)))
cat(sprintf("\nPolar vs MST: delta %+.2f +/- %.2f pp, Wilcoxon p=%.4g\n",
            mean(d_pm), sd(d_pm), wt_pm$p.value))
if (!is.null(ci_pm)) cat(sprintf("  95%% BCa CI: [%.2f, %.2f] pp\n",
                                 ci_pm[4], ci_pm[5]))
cat(sprintf("\nACPL vs Polar: delta %+.2f +/- %.2f pp, Wilcoxon p=%.4g\n",
            mean(d_ap), sd(d_ap), wt_ap$p.value))
if (!is.null(ci_ap)) cat(sprintf("  95%% BCa CI: [%.2f, %.2f] pp\n",
                                 ci_ap[4], ci_ap[5]))

write.csv(data.frame(seed = 41:50, polar = polar_res, acpl = acpl_res,
                     mst = mst_res),
          "buettner_per_seed.csv", row.names = FALSE)
cat("\nSaved buettner_per_seed.csv. Paste the printed summary to chat.\n")
