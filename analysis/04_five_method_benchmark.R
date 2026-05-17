# ==============================================================================
# Analysis 04: Five-method benchmark on Nestorowa HSC
# Section 4.4 of the paper. ACPL vs MST, PAGA+DPT, Slingshot, DPT.
# Outputs: results/five_method_benchmark.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("data-raw/load_datasets.R")

library(Seurat)

seu     <- load_nestorowa()
if (is.null(seu)) stop("Nestorowa load failed.")
coords  <- Embeddings(seu, "umap")
phases  <- as.character(seu$Phase)

# ACPL
arc_acpl    <- acpl_arc(coords)
sw_acpl     <- acpl_swsa(arc_acpl, phases)

# MST
mst_pt      <- mst_pseudotime(coords, "G1", phases)
sw_mst      <- acpl_swsa(mst_pt, phases, mirror = TRUE)

# Slingshot (requires slingshot package)
sw_sling <- tryCatch({
  sling_pt <- slingshot_pseudotime(coords, phases, "G1")
  acpl_swsa(sling_pt, phases, mirror = TRUE)
}, error = function(e) {
  message("Slingshot failed: ", conditionMessage(e))
  NA
})

# PAGA+DPT (requires reticulate + scanpy)
sw_paga <- tryCatch({
  paga_pt <- paga_dpt_pseudotime(Embeddings(seu, "pca")[, 1:20], phases)
  acpl_swsa(paga_pt, phases, mirror = TRUE)
}, error = function(e) {
  message("PAGA failed: ", conditionMessage(e))
  NA
})

results <- data.frame(
  Method = c("ACPL", "MST (Monocle-style)", "PAGA+DPT", "Slingshot"),
  SW_SA  = c(sw_acpl, sw_mst, sw_paga, sw_sling),
  Notes  = c("Polar prior; O(N log N); no root",
             "Greedy MST; G1 root cell",
             "Cluster graph; Leiden+DPT",
             "Principal curves; G1 start cluster")
)
print(results, row.names = FALSE)
write.csv(results, "results/five_method_benchmark.csv", row.names = FALSE)
cat("\nResults saved: results/five_method_benchmark.csv\n")
