# ==============================================================================
# polar_vs_mst_comparison.R
#
# Three analyses that strengthen the new paper:
#   1. Polar-alone vs MST paired comparison on Nestorowa (10 UMAP seeds,
#      paired Wilcoxon, BCa bootstrap CI on the delta)
#   2. Same comparison on Spellman
#   3. Runtime scaling on synthetic cyclic data at N = 100 to 20,000
#
# Run: Rscript polar_vs_mst_comparison.R
# Expected runtime: 15-25 minutes total.
# ==============================================================================

cat("\n=== Polar-alone vs MST: paired comparison and runtime scaling ===\n\n")

required_cran <- c("Seurat", "dplyr", "igraph", "boot", "BiocManager")
required_bioc <- c("scRNAseq", "yeastCC", "SingleCellExperiment",
                   "SummarizedExperiment", "AnnotationDbi", "org.Mm.eg.db")

missing_cran <- required_cran[!required_cran %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) install.packages(missing_cran, repos = "https://cloud.r-project.org")
missing_bioc <- required_bioc[!required_bioc %in% installed.packages()[, "Package"]]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

suppressMessages({
  library(Seurat); library(dplyr); library(igraph); library(boot)
  library(scRNAseq); library(yeastCC); library(SingleCellExperiment)
  library(SummarizedExperiment); library(AnnotationDbi); library(org.Mm.eg.db)
})

set.seed(42)

# ── Helpers ──────────────────────────────────────────────────────────────────

phase_map <- c("G1" = 1L, "S" = 2L, "G2M" = 3L)

calc_swsa <- function(time_proxy, phases, window = 10) {
  df <- data.frame(t = time_proxy, Phase = as.character(phases)) %>%
    dplyr::filter(!is.na(Phase), Phase %in% c("G1", "S", "G2M")) %>%
    dplyr::mutate(PhaseNum = phase_map[Phase]) %>%
    dplyr::arrange(t)
  n <- nrow(df)
  if (n <= window) return(NA_real_)
  100 * mean(df$PhaseNum[1:(n - window)] <= df$PhaseNum[(window + 1):n], na.rm = TRUE)
}

find_adaptive_origin <- function(x, y) {
  optim(c(median(x), median(y)),
        function(p) var(sqrt((x - p[1])^2 + (y - p[2])^2)))$par
}

polar_alone_swsa <- function(umap_coords, phases) {
  uv <- find_adaptive_origin(umap_coords[, 1], umap_coords[, 2])
  u  <- umap_coords[, 1] - uv[1]
  v  <- umap_coords[, 2] - uv[2]
  theta <- atan2(v, u)
  max(calc_swsa(theta, phases), calc_swsa(-theta, phases), na.rm = TRUE)
}

mst_swsa <- function(umap_coords, phases) {
  dmat <- as.matrix(dist(umap_coords))
  g    <- igraph::graph_from_adjacency_matrix(dmat, weighted = TRUE, mode = "undirected")
  mst  <- igraph::mst(g)
  root <- which(as.character(phases) == "G1")[1]
  if (is.na(root)) root <- 1
  pt <- as.numeric(igraph::distances(mst, v = root,
                                     to = V(mst),
                                     weights = E(mst)$weight))
  max(calc_swsa(pt, phases), calc_swsa(-pt, phases), na.rm = TRUE)
}

# ── Build Nestorowa Seurat object (with proper gene mapping) ─────────────────

build_nestorowa <- function() {
  cat("Building Nestorowa Seurat object...\n")
  sce <- NestorowaHSCData()
  ensembl_clean <- sub("\\..*$", "", rownames(sce))
  symbols <- mapIds(org.Mm.eg.db, keys = ensembl_clean,
                    column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  keep <- !is.na(symbols) & !duplicated(symbols)
  counts <- as.matrix(assay(sce, assayNames(sce)[1])[keep, ])
  rownames(counts) <- symbols[keep]
  colnames(counts) <- colnames(sce)
  
  obj <- CreateSeuratObject(counts = counts, project = "nestorowa")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  
  to_title <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
  s_g   <- to_title(cc.genes$s.genes)
  g2m_g <- to_title(cc.genes$g2m.genes)
  obj <- CellCycleScoring(obj, s.features = s_g, g2m.features = g2m_g,
                          set.ident = FALSE)
  obj
}

build_spellman <- function() {
  cat("Building Spellman matrix...\n")
  e <- new.env(); data("yeastCC", package = "yeastCC", envir = e)
  yc <- get("yeastCC", envir = e)
  expr <- Biobase::exprs(yc)
  alpha_cols <- grep("^alpha", colnames(expr), value = TRUE, ignore.case = TRUE)
  if (length(alpha_cols) < 5) alpha_cols <- colnames(expr)[seq_len(min(18, ncol(expr)))]
  expr <- expr[, alpha_cols]
  expr <- expr[complete.cases(expr), ]
  list(expr = t(expr), n_timepoints = length(alpha_cols))
}

# ── Analysis 1: Nestorowa paired comparison across UMAP seeds ────────────────

cat("\n=== Analysis 1: Nestorowa paired comparison (polar vs MST) ===\n")
nest_obj <- build_nestorowa()
cat("Phase distribution:\n"); print(table(nest_obj$Phase))

n_seeds <- 10
polar_results <- numeric(n_seeds)
mst_results   <- numeric(n_seeds)

for (i in seq_len(n_seeds)) {
  cat(sprintf("  Seed %d/%d...\n", i, n_seeds))
  set.seed(40 + i)
  nest_obj <- RunUMAP(nest_obj, dims = 1:20, verbose = FALSE)
  um <- Embeddings(nest_obj, "umap")
  ph <- nest_obj$Phase
  polar_results[i] <- polar_alone_swsa(um, ph)
  mst_results[i]   <- mst_swsa(um, ph)
}

deltas <- polar_results - mst_results
wt <- wilcox.test(polar_results, mst_results, paired = TRUE,
                  exact = FALSE)
boot_delta <- boot(deltas, function(d, i) mean(d[i]), R = 2000)
boot_ci <- boot.ci(boot_delta, type = "bca")

cat("\nNestorowa results (10 UMAP seeds):\n")
cat(sprintf("  polar-alone: %.1f%% +/- %.2f\n", mean(polar_results), sd(polar_results)))
cat(sprintf("  MST:         %.1f%% +/- %.2f\n", mean(mst_results),   sd(mst_results)))
cat(sprintf("  delta:       %+.2f +/- %.2f pp\n", mean(deltas), sd(deltas)))
cat(sprintf("  Wilcoxon paired: V=%.0f, p=%.4g\n", wt$statistic, wt$p.value))
if (!is.null(boot_ci$bca)) {
  cat(sprintf("  95%% BCa CI on delta: [%.2f, %.2f] pp\n",
              boot_ci$bca[4], boot_ci$bca[5]))
}

nest_summary <- data.frame(
  dataset = "Nestorowa", N = ncol(nest_obj),
  polar_mean = round(mean(polar_results), 1),
  polar_sd   = round(sd(polar_results), 2),
  mst_mean   = round(mean(mst_results), 1),
  mst_sd     = round(sd(mst_results), 2),
  delta_mean = round(mean(deltas), 2),
  delta_sd   = round(sd(deltas), 2),
  wilcox_p   = signif(wt$p.value, 3),
  ci_low     = if (!is.null(boot_ci$bca)) round(boot_ci$bca[4], 2) else NA,
  ci_high    = if (!is.null(boot_ci$bca)) round(boot_ci$bca[5], 2) else NA
)

# ── Analysis 2: Spellman paired comparison ───────────────────────────────────

cat("\n=== Analysis 2: Spellman paired comparison ===\n")
spell <- build_spellman()
# Spellman truth = time, treated as phase G1->S->G2M tertiles
# Build phase labels from timepoint number
n_t <- spell$n_timepoints
phase_spell <- rep(NA_character_, n_t)
phase_spell[1:floor(n_t/3)]                          <- "G1"
phase_spell[(floor(n_t/3)+1):floor(2*n_t/3)]         <- "S"
phase_spell[(floor(2*n_t/3)+1):n_t]                  <- "G2M"

# Need uwot for Spellman (no Seurat pipeline)
suppressMessages(library(uwot))

polar_s <- numeric(n_seeds); mst_s <- numeric(n_seeds)
for (i in seq_len(n_seeds)) {
  cat(sprintf("  Seed %d/%d...\n", i, n_seeds))
  set.seed(40 + i)
  um <- umap(spell$expr,
             n_neighbors = min(n_t - 1, 14),
             min_dist = 0.3, n_components = 2)
  polar_s[i] <- polar_alone_swsa(um, phase_spell)
  mst_s[i]   <- mst_swsa(um, phase_spell)
}
deltas_s <- polar_s - mst_s

wt_s <- wilcox.test(polar_s, mst_s, paired = TRUE, exact = FALSE)
boot_s <- boot(deltas_s, function(d, i) mean(d[i]), R = 2000)
boot_ci_s <- tryCatch(boot.ci(boot_s, type = "bca"), error = function(e) NULL)

cat("\nSpellman results (10 UMAP seeds):\n")
cat(sprintf("  polar-alone: %.1f%% +/- %.2f\n", mean(polar_s), sd(polar_s)))
cat(sprintf("  MST:         %.1f%% +/- %.2f\n", mean(mst_s),   sd(mst_s)))
cat(sprintf("  delta:       %+.2f +/- %.2f pp\n", mean(deltas_s), sd(deltas_s)))
cat(sprintf("  Wilcoxon paired: V=%.0f, p=%.4g\n", wt_s$statistic, wt_s$p.value))
ci_low_s <- if (!is.null(boot_ci_s) && !is.null(boot_ci_s$bca)) boot_ci_s$bca[4] else NA
ci_high_s <- if (!is.null(boot_ci_s) && !is.null(boot_ci_s$bca)) boot_ci_s$bca[5] else NA
if (!is.na(ci_low_s)) cat(sprintf("  95%% BCa CI on delta: [%.2f, %.2f] pp\n",
                                  ci_low_s, ci_high_s))

spell_summary <- data.frame(
  dataset = "Spellman", N = n_t,
  polar_mean = round(mean(polar_s), 1),
  polar_sd   = round(sd(polar_s), 2),
  mst_mean   = round(mean(mst_s), 1),
  mst_sd     = round(sd(mst_s), 2),
  delta_mean = round(mean(deltas_s), 2),
  delta_sd   = round(sd(deltas_s), 2),
  wilcox_p   = signif(wt_s$p.value, 3),
  ci_low     = if (!is.na(ci_low_s)) round(ci_low_s, 2) else NA,
  ci_high    = if (!is.na(ci_high_s)) round(ci_high_s, 2) else NA
)

# ── Analysis 3: Runtime scaling on synthetic cyclic data ─────────────────────

cat("\n=== Analysis 3: Runtime scaling (synthetic cyclic data) ===\n")

generate_cyclic <- function(N, sigma = 0.2, dim = 50) {
  t <- seq(0, 2 * pi, length.out = N)
  base_x <- cos(t); base_y <- sin(t)
  # Embed in `dim`-dimensional space with noise
  signal <- cbind(base_x, base_y, matrix(0, N, dim - 2))
  rotation <- qr.Q(qr(matrix(rnorm(dim^2), dim, dim)))
  signal <- signal %*% rotation
  signal + matrix(rnorm(N * dim, sd = sigma), N, dim)
}

n_levels <- c(100, 500, 1000, 5000, 10000, 20000)
runtime_polar <- numeric(length(n_levels))
runtime_mst   <- numeric(length(n_levels))

for (j in seq_along(n_levels)) {
  N <- n_levels[j]
  cat(sprintf("  N = %d...\n", N))
  set.seed(N)
  X <- generate_cyclic(N)
  um <- umap(X, n_neighbors = min(N - 1, 30), min_dist = 0.3,
             n_components = 2, verbose = FALSE)
  
  # Time polar-alone
  t0 <- Sys.time()
  uv <- find_adaptive_origin(um[, 1], um[, 2])
  u  <- um[, 1] - uv[1]; v <- um[, 2] - uv[2]
  theta <- atan2(v, u)
  ord_polar <- order(theta)
  runtime_polar[j] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # Time MST (skip if N too large to avoid memory issues)
  if (N <= 10000) {
    t0 <- Sys.time()
    dmat <- as.matrix(dist(um))
    g    <- graph_from_adjacency_matrix(dmat, weighted = TRUE, mode = "undirected")
    mst  <- igraph::mst(g)
    pt   <- as.numeric(distances(mst, v = 1, to = V(mst),
                                 weights = E(mst)$weight))
    ord_mst <- order(pt)
    runtime_mst[j] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  } else {
    runtime_mst[j] <- NA  # too memory-intensive to attempt
  }
  
  cat(sprintf("    polar: %.3f s   MST: %s s\n",
              runtime_polar[j],
              if (is.na(runtime_mst[j])) "memory-skipped" else sprintf("%.3f", runtime_mst[j])))
}

runtime_summary <- data.frame(
  N = n_levels,
  polar_seconds = round(runtime_polar, 3),
  mst_seconds   = round(runtime_mst, 3),
  speedup       = round(runtime_mst / runtime_polar, 1)
)

# ── Final summary ────────────────────────────────────────────────────────────

cat("\n\n========================================\n")
cat("=== FINAL RESULTS ===\n")
cat("========================================\n\n")

cat("Accuracy comparison (polar-alone vs MST):\n")
acc_summary <- rbind(nest_summary, spell_summary)
print(acc_summary, row.names = FALSE)

cat("\nRuntime scaling:\n")
print(runtime_summary, row.names = FALSE)

write.csv(acc_summary, "polar_vs_mst_accuracy.csv", row.names = FALSE)
write.csv(runtime_summary, "polar_vs_mst_runtime.csv", row.names = FALSE)
cat(sprintf("\nSaved to: polar_vs_mst_accuracy.csv, polar_vs_mst_runtime.csv\n"))
cat("Paste both tables back to chat.\n")
