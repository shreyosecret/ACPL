# ==============================================================================
# Analysis 03: Nestorowa ablation at scale
# Section 4.3 of the paper. Same ablation as Spellman but on N=1,920 cells.
# Key finding: polar transform alone matches full ACPL on this manifold.
# Outputs: results/nestorowa_ablation.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("data-raw/load_datasets.R")

library(Seurat); library(dplyr)

seu <- load_nestorowa()
if (is.null(seu)) stop("Nestorowa load failed.")

coords  <- Embeddings(seu, "umap")
phases  <- as.character(seu$Phase)

# All five pipeline steps
ablation <- list()

# Cartesian MST baseline
mst_pt <- mst_pseudotime(coords, "G1", phases)
ablation[["Cartesian MST"]] <- list(
  sa = acpl_swsa(mst_pt, phases, mirror = TRUE),
  note = "Greedy Euclidean MST")

# Polar transform only
uv <- find_adaptive_origin(coords[,1], coords[,2])
th <- atan2(coords[,2] - uv[2], coords[,1] - uv[1])
ablation[["Polar only"]] <- list(
  sa = acpl_swsa(th, phases, mirror = TRUE),
  note = "theta as time proxy")

# Polar + LOESS on r
r_raw <- sqrt((coords[,1] - uv[1])^2 + (coords[,2] - uv[2])^2)
r_smooth <- predict(loess(r_raw ~ th, span = 0.25))
ablation[["Polar+LOESS on r"]] <- list(
  sa = acpl_swsa(r_smooth, phases, mirror = TRUE),
  note = "r uninformative in loop")

# Raw arc length
ord <- order(th)
theta_uw <- unwrap_theta(th[ord])
un <- cos(theta_uw); vn <- sin(theta_uw)
arc_raw <- c(0, cumsum(sqrt(diff(un)^2 + diff(vn)^2)))
arc_back <- arc_raw[order(seq_along(th)[ord])]
ablation[["Arc, no LOESS (deg.)"]] <- list(
  sa = acpl_swsa(arc_back, phases, mirror = TRUE),
  note = "Degenerate")

# Full ACPL
arc <- acpl_arc(coords)
ablation[["ACPL (Arc+LOESS)"]] <- list(
  sa = acpl_swsa(arc, phases, mirror = TRUE),
  note = "Principled")

out <- data.frame(
  pipeline_step = names(ablation),
  SW_SA_pct     = sapply(ablation, `[[`, "sa"),
  note          = sapply(ablation, `[[`, "note")
)
print(out, row.names = FALSE)
write.csv(out, "results/nestorowa_ablation.csv", row.names = FALSE)
cat("\nResults saved: results/nestorowa_ablation.csv\n")
