# ==============================================================================
# Analysis 05: Pairwise bootstrap CIs across all method/dataset combinations
# Section 4.5 of the paper. The single most important statistical analysis.
# Outputs: results/bootstrap_all_pairs.csv, figures/fig_bootstrap_forest.png
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("R/statistics.R")
source("R/plotting.R")
source("data-raw/load_datasets.R")

library(Seurat); library(ggplot2); library(boot); library(dplyr)

results <- list()
window  <- 10

# ── Nestorowa: ACPL vs MST ──
cat("Nestorowa: ACPL vs MST\n")
seu <- load_nestorowa()
if (!is.null(seu)) {
  coords <- Embeddings(seu, "umap")
  phases <- as.character(seu$Phase)
  wa <- window_correctness(acpl_arc(coords), phases, window)
  wm <- window_correctness(mst_pseudotime(coords, "G1", phases), phases,
                            window, mirror = TRUE)
  results[["nest_mst"]] <- bootstrap_ci(wa, wm, "ACPL", "MST", "Nestorowa HSC")
}

# ── Nestorowa: ACPL vs Slingshot ──
cat("Nestorowa: ACPL vs Slingshot\n")
if (!is.null(seu)) {
  coords <- Embeddings(seu, "umap")
  phases <- as.character(seu$Phase)
  ws <- tryCatch({
    sling_pt <- slingshot_pseudotime(coords, phases, "G1")
    window_correctness(sling_pt, phases, window, mirror = TRUE)
  }, error = function(e) NULL)
  if (!is.null(ws)) {
    wa2 <- window_correctness(acpl_arc(coords), phases, window)
    n_use <- min(length(wa2), length(ws))
    results[["nest_sling"]] <- bootstrap_ci(wa2[1:n_use], ws[1:n_use],
                                             "ACPL", "Slingshot",
                                             "Nestorowa HSC")
  }
}

# ── Buettner: ACPL vs MST ──
cat("Buettner: ACPL vs MST\n")
buet <- load_buettner()
if (!is.null(buet)) {
  cb <- Embeddings(buet, "umap")
  pb <- as.character(buet$Phase)
  wa <- window_correctness(acpl_arc(cb), pb, window)
  wm <- window_correctness(mst_pseudotime(cb, "G1", pb), pb, window,
                            mirror = TRUE)
  results[["buet_mst"]] <- bootstrap_ci(wa, wm, "ACPL", "MST",
                                         "Buettner mESC")
}

# ── Leng: ACPL vs MST ──
cat("Leng: ACPL vs MST\n")
leng <- load_leng()
if (!is.null(leng)) {
  cl <- Embeddings(leng, "umap")
  pl <- as.character(leng$Phase)
  wa <- window_correctness(acpl_arc(cl), pl, window)
  wm <- window_correctness(mst_pseudotime(cl, "G1", pl), pl, window,
                            mirror = TRUE)
  results[["leng_mst"]] <- bootstrap_ci(wa, wm, "ACPL", "MST",
                                         "Leng mESC")
}

all_res <- do.call(rbind, results)
print(all_res[, c("Dataset", "Method_A", "Method_B",
                  "Obs_pp", "CI_lo_bca", "CI_hi_bca",
                  "Wilcoxon_p", "Sig")], row.names = FALSE)

write.csv(all_res, "results/bootstrap_all_pairs.csv", row.names = FALSE)

if (nrow(all_res) > 0) {
  p <- bootstrap_forest(all_res)
  ggsave("figures/fig_bootstrap_forest.png", p, width = 9, height = 5, dpi = 300)
  cat("\nFigure saved: figures/fig_bootstrap_forest.png\n")
}
cat("Results saved: results/bootstrap_all_pairs.csv\n")
