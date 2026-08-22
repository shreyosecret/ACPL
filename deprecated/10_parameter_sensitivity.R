# ==============================================================================
# Analysis 10: Parameter sensitivity (LOESS span and UMAP hyperparameters)
# Section 4.10. Demonstrates practical robustness of ACPL.
# Outputs: results/loess_span_sensitivity.csv, results/umap_sensitivity.csv
# ==============================================================================

source("R/acpl_engine.R")
source("data-raw/load_datasets.R")

library(Seurat); library(dplyr)

seu <- load_nestorowa()
if (is.null(seu)) stop("Nestorowa load failed.")

# ── LOESS span sensitivity ──
cat("LOESS span sensitivity\n")
spans <- c(0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 0.95)
Ks    <- c(2, 5, 10, 20)

span_grid <- expand.grid(span = spans, K = Ks, SW_SA = NA_real_)
coords <- Embeddings(seu, "umap")
phases <- as.character(seu$Phase)

for (i in seq_len(nrow(span_grid))) {
  arc_i <- acpl_arc(coords, span = span_grid$span[i])
  span_grid$SW_SA[i] <- acpl_swsa(arc_i, phases)
}
cat(sprintf("Span variance: %.4f pp^2 (range %.2f%%-%.2f%%)\n",
            var(span_grid$SW_SA, na.rm = TRUE),
            min(span_grid$SW_SA), max(span_grid$SW_SA)))
write.csv(span_grid, "results/loess_span_sensitivity.csv", row.names = FALSE)

# ── UMAP sensitivity ──
cat("\nUMAP sensitivity (25 configurations)\n")
nn_vals <- c(10, 20, 30, 40, 60)
md_vals <- c(0.1, 0.2, 0.3, 0.5, 0.8)
sens_grid <- expand.grid(n_neighbors = nn_vals, min_dist = md_vals,
                          SW_SA = NA_real_)

for (i in seq_len(nrow(sens_grid))) {
  seu_tmp <- RunUMAP(seu, dims = 1:10,
                      n.neighbors = sens_grid$n_neighbors[i],
                      min.dist    = sens_grid$min_dist[i],
                      verbose = FALSE, seed.use = 42)
  coords_i <- Embeddings(seu_tmp, "umap")
  arc_i    <- acpl_arc(coords_i)
  sens_grid$SW_SA[i] <- acpl_swsa(arc_i, phases)
  if (i %% 5 == 0) cat(sprintf("  %d/%d\n", i, nrow(sens_grid)))
}

cat(sprintf("\nUMAP SW-SA range: %.2f%% to %.2f%% (SD = %.3f pp)\n",
            min(sens_grid$SW_SA), max(sens_grid$SW_SA),
            sd(sens_grid$SW_SA)))
write.csv(sens_grid, "results/umap_sensitivity.csv", row.names = FALSE)
cat("\nResults saved.\n")
