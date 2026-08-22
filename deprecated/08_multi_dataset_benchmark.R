# ==============================================================================
# Analysis 08: Multi-dataset benchmark
# Section 4.8. ACPL vs MST across all five biological datasets.
# Outputs: results/multi_dataset_benchmark.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("data-raw/load_datasets.R")

library(Seurat); library(uwot); library(igraph)

results <- list()

# ── Spellman ──
sp <- load_spellman()
if (!is.null(sp)) {
  pca_s <- prcomp(t(sp$counts), scale. = TRUE)$x[, 1:2]
  set.seed(42)
  umap_s <- uwot::umap(pca_s, n_neighbors = 5, min_dist = 0.5, verbose = FALSE)
  true_idx <- seq_len(nrow(umap_s))

  arc <- acpl_arc(umap_s)
  g_a <- acpl_graph(umap_s, arc)
  acpl_sa <- 100 * mean(true_idx[g_a$edges[, "to"]] > true_idx[g_a$edges[, "from"]])

  dm <- as.matrix(dist(umap_s))
  g  <- igraph::graph_from_adjacency_matrix(dm, weighted = TRUE, mode = "undirected")
  el <- igraph::as_edgelist(igraph::mst(g))
  fwd <- mean(true_idx[as.integer(el[,2])] > true_idx[as.integer(el[,1])])
  mst_sa <- 100 * max(fwd, 1 - fwd)

  results[["spellman"]] <- data.frame(Dataset = "Spellman yeast", N = ncol(sp$counts),
                                       ACPL = round(acpl_sa,1), MST = round(mst_sa,1),
                                       Diff = round(acpl_sa - mst_sa, 1))
}

# ── Helper function for Seurat-based datasets ──
benchmark_seurat <- function(seu, name) {
  if (is.null(seu)) return(NULL)
  coords <- Embeddings(seu, "umap")
  phases <- as.character(seu$Phase)
  acpl_sa <- acpl_swsa(acpl_arc(coords), phases)
  mst_sa  <- acpl_swsa(mst_pseudotime(coords, "G1", phases), phases, mirror = TRUE)
  data.frame(Dataset = name, N = ncol(seu),
             ACPL = acpl_sa, MST = mst_sa, Diff = round(acpl_sa - mst_sa, 1))
}

results[["nestorowa"]] <- benchmark_seurat(load_nestorowa(), "Nestorowa HSC")
results[["buettner"]]  <- benchmark_seurat(load_buettner(),  "Buettner mESC")
results[["leng"]]      <- benchmark_seurat(load_leng(),      "Leng mESC")
results[["richard"]]   <- benchmark_seurat(load_richard(),   "Richard T Cells")

out <- do.call(rbind, results)
print(out, row.names = FALSE)
write.csv(out, "results/multi_dataset_benchmark.csv", row.names = FALSE)
cat("\nResults saved: results/multi_dataset_benchmark.csv\n")
