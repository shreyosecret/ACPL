# ==============================================================================
# Analysis 02: Spellman yeast ablation
# Section 4.2 of the paper. Each pipeline component evaluated separately.
# Outputs: results/spellman_ablation.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("data-raw/load_datasets.R")

library(igraph); library(uwot); library(dplyr)

dat <- load_spellman()
if (is.null(dat)) stop("Spellman load failed; cannot run this analysis.")

# Reduce dimensions
pca_s <- prcomp(t(dat$counts), scale. = TRUE)$x[, 1:2]
set.seed(42)
umap_s <- uwot::umap(pca_s, n_neighbors = 5, min_dist = 0.5, verbose = FALSE)

# True ordering from timepoint labels (already sorted in Spellman experiment)
true_idx <- seq_len(nrow(umap_s))

calc_sa <- function(graph, true_idx) {
  if (is.null(graph$edges) || nrow(graph$edges) == 0) return(NA)
  fwd <- mean(true_idx[graph$edges[, "to"]] > true_idx[graph$edges[, "from"]])
  100 * max(fwd, 1 - fwd)
}

ablation <- list()

# Cartesian MST baseline
dm <- as.matrix(dist(umap_s))
g  <- igraph::graph_from_adjacency_matrix(dm, weighted = TRUE, mode = "undirected")
mt <- igraph::mst(g)
el <- igraph::as_edgelist(mt)
sa_cart <- 100 * max(mean(true_idx[as.integer(el[, 2])] >
                          true_idx[as.integer(el[, 1])]),
                     mean(true_idx[as.integer(el[, 1])] >
                          true_idx[as.integer(el[, 2])]))
ablation[["Cartesian MST"]] <- list(sa = sa_cart, note = "Greedy Euclidean MST")

# Polar transform only
uv  <- find_adaptive_origin(umap_s[,1], umap_s[,2])
th  <- atan2(umap_s[,2] - uv[2], umap_s[,1] - uv[1])
g_t <- acpl_graph(umap_s, th, K = 2)
ablation[["Polar transform only"]] <- list(sa = calc_sa(g_t, true_idx),
                                           note = "theta as time proxy")

# Polar + LOESS on r
r_raw    <- sqrt((umap_s[,1] - uv[1])^2 + (umap_s[,2] - uv[2])^2)
r_smooth <- predict(loess(r_raw ~ th, span = 0.25))
g_r <- acpl_graph(umap_s, r_smooth, K = 2)
ablation[["Polar+LOESS on r"]] <- list(sa = calc_sa(g_r, true_idx),
                                       note = "r uninformative in loop")

# Raw arc length (degenerate)
ord       <- order(th)
theta_uw  <- unwrap_theta(th[ord])
un <- cos(theta_uw); vn <- sin(theta_uw)
arc_raw   <- c(0, cumsum(sqrt(diff(un)^2 + diff(vn)^2)))
arc_orig  <- arc_raw[order(seq_along(th)[ord])]
g_a <- acpl_graph(umap_s, arc_orig, K = 2)
ablation[["Arc, no LOESS (deg.)"]] <- list(sa = calc_sa(g_a, true_idx),
                                           note = "Degenerate, Theorem 2")

# Full ACPL
arc_acpl <- acpl_arc(umap_s)
g_acpl   <- acpl_graph(umap_s, arc_acpl, K = 2)
ablation[["ACPL (Arc+LOESS)"]] <- list(sa = calc_sa(g_acpl, true_idx),
                                       note = "Principled prior")

out <- data.frame(
  pipeline_step = names(ablation),
  SA_pct        = round(sapply(ablation, `[[`, "sa"), 1),
  note          = sapply(ablation, `[[`, "note")
)
print(out, row.names = FALSE)
write.csv(out, "results/spellman_ablation.csv", row.names = FALSE)
cat("\nResults saved: results/spellman_ablation.csv\n")
