# ==============================================================================
# R/methods.R -- ACPL and the comparison methods
# ==============================================================================

# ---- ACPL: Arc plus LOESS Polar-Linearisation --------------------------------

#' Adaptive hub minimising radial variance
find_adaptive_origin <- function(x, y)
  optim(c(median(x), median(y)),
        function(p) var(sqrt((x - p[1])^2 + (y - p[2])^2)))$par

unwrap_theta <- function(theta) {
  d <- diff(theta); d <- d - 2 * pi * round(d / (2 * pi))
  c(theta[1], theta[1] + cumsum(d))
}

#' ACPL arc length
#'
#' LOESS smoothing is not cosmetic. Raw cumulative arc length is strictly
#' increasing in traversal index by positive definiteness of the norm, so a
#' monotonicity filter on it admits every edge and reports perfect accuracy as
#' an arithmetic identity. Any implementation reporting 100% is measuring a
#' tautology.
acpl_arc <- function(coords, span = 0.25) {
  if (ncol(coords) != 2) stop("acpl_arc requires 2D coordinates")
  uv <- find_adaptive_origin(coords[, 1], coords[, 2])
  th <- atan2(coords[, 2] - uv[2], coords[, 1] - uv[1])
  o  <- order(th); tu <- unwrap_theta(th[o])
  arc <- c(0, cumsum(sqrt(diff(cos(tu))^2 + diff(sin(tu))^2)))
  predict(loess(arc ~ tu, span = span))[order(seq_along(coords[, 1])[o])]
}

#' Raw polar angle around the adaptive hub
acpl_theta <- function(coords) {
  uv <- find_adaptive_origin(coords[, 1], coords[, 2])
  atan2(coords[, 2] - uv[2], coords[, 1] - uv[1])
}

#' Directed acyclic graph from the arc-length ordering
acpl_graph <- function(coords, tau, K = 2, lambda = 1.0) {
  N <- nrow(coords)
  th <- acpl_theta(coords)
  edges <- vector("list", N)
  for (i in seq_len(N)) {
    valid <- which(tau > tau[i])
    if (!length(valid)) next
    dth <- abs(th[valid] - th[i]); dth <- pmin(dth, 2 * pi - dth)
    d <- abs(tau[valid] - tau[i]) + lambda * dth
    edges[[i]] <- cbind(from = i, to = valid[order(d)[1:min(K, length(valid))]])
  }
  list(edges = do.call(rbind, edges), tau = tau)
}

acpl_is_acyclic <- function(graph) {
  if (is.null(graph$edges) || !nrow(graph$edges)) return(TRUE)
  all(graph$tau[graph$edges[, "to"]] > graph$tau[graph$edges[, "from"]])
}

# ---- baselines ---------------------------------------------------------------
# Both accept any coordinate matrix, so they can be run on their NATIVE input
# (principal components) as well as on the 2D embedding. v1 ran them only on the
# embedding, which was the stated reason for a desk rejection.

mst_pseudotime <- function(coords, phases, root_phase = "G1") {
  g  <- igraph::graph_from_adjacency_matrix(as.matrix(dist(coords)),
          weighted = TRUE, mode = "undirected")
  tr <- igraph::mst(g)
  ri <- which(phases == root_phase)[1]
  if (is.na(ri)) stop("No cell with root_phase found")
  as.numeric(igraph::distances(tr, v = ri, weights = igraph::E(tr)$weight))
}

slingshot_pseudotime <- function(coords, phases, start_phase = "G1",
                                 reduced_dim = "PCA") {
  keep <- as.character(phases) %in% names(PHASE_MAP)
  sds <- slingshot::slingshot(data = coords[keep, , drop = FALSE],
           clusterLabels = as.character(phases)[keep],
           start.clus = start_phase, reducedDim = reduced_dim)
  pt <- rep(NA_real_, nrow(coords))
  pt[keep] <- rowMeans(slingshot::slingPseudotime(sds), na.rm = TRUE)
  pt
}

#' tricycle cell-cycle position (reference-based)
#'
#' Requires the expression matrix, not the cached embedding. Its reference is
#' mammalian, so it cannot be applied outside the organisms that reference
#' covers; on yeast it returns no projection genes and fails.
tricycle_position <- function(sce, species = c("mouse", "human")) {
  species <- match.arg(species)
  SummarizedExperiment::assay(sce, "counts") <-
    as.matrix(SummarizedExperiment::assay(sce, pick_assay(sce)))
  sce <- scuttle::logNormCounts(sce)
  out <- tricycle::estimate_cycle_position(sce, gname.type = "SYMBOL",
                                           species = species)
  as.numeric(SummarizedExperiment::colData(out)$tricyclePosition)
}
