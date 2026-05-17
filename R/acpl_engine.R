# ==============================================================================
# ACPL Engine: Arc + LOESS Polar-Linearization
# Core algorithm functions. Source this file before running any analysis.
# ==============================================================================

#' Find adaptive radial origin
#'
#' Computes hub coordinates that minimise radial variance across the cell
#' population. For circular manifolds this recovers the centroid; for branching
#' manifolds it finds a hub balancing radial spread across lineages.
#'
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#' @return Numeric vector c(cx, cy) of hub coordinates
#' @export
find_adaptive_origin <- function(x, y) {
  fn <- function(p) var(sqrt((x - p[1])^2 + (y - p[2])^2))
  optim(c(median(x), median(y)), fn)$par
}

#' Unwrap angular coordinate
#'
#' Converts a wrapped angle (in [-pi, pi]) into a continuous monotonically
#' increasing variable to handle multi-cycle datasets.
#'
#' @param theta Numeric vector of angles in radians
#' @return Numeric vector of unwrapped angles
#' @export
unwrap_theta <- function(theta) {
  d <- diff(theta)
  d <- d - 2 * pi * round(d / (2 * pi))
  c(theta[1], theta[1] + cumsum(d))
}

#' Compute ACPL arc length
#'
#' Main ACPL pipeline. Takes 2D coordinates, transforms to polar, computes
#' cumulative arc length along the angular ordering, and applies LOESS smoothing.
#'
#' @param coords 2D coordinate matrix (rows = cells, columns = dimensions)
#' @param span LOESS smoothing span (default 0.25; results are flat across [0.05, 0.95])
#' @return Numeric vector of arc length values, one per cell, in original order
#' @examples
#' set.seed(1)
#' theta <- seq(0, 2*pi, length.out = 100)
#' coords <- cbind(cos(theta) + rnorm(100, 0, 0.05),
#'                 sin(theta) + rnorm(100, 0, 0.05))
#' arc <- acpl_arc(coords)
#' @export
acpl_arc <- function(coords, span = 0.25) {
  if (ncol(coords) != 2) stop("acpl_arc requires 2D coordinates")
  uv         <- find_adaptive_origin(coords[, 1], coords[, 2])
  u          <- coords[, 1] - uv[1]
  v          <- coords[, 2] - uv[2]
  theta_raw  <- atan2(v, u)
  ord        <- order(theta_raw)
  theta_uw   <- unwrap_theta(theta_raw[ord])
  un         <- cos(theta_uw)
  vn         <- sin(theta_uw)
  arc_raw    <- c(0, cumsum(sqrt(diff(un)^2 + diff(vn)^2)))
  arc_smooth <- predict(loess(arc_raw ~ theta_uw, span = span))
  arc_smooth[order(seq_along(coords[, 1])[ord])]
}

#' Build directed acyclic graph from arc length ordering
#'
#' Constructs the directed graph where edges are admitted only in the direction
#' of increasing arc length. By Theorem 1 of the paper, this graph is acyclic
#' by construction.
#'
#' @param coords 2D coordinate matrix
#' @param tau Numeric vector of arc length values (output of acpl_arc)
#' @param K Number of nearest valid neighbours per node (default 2)
#' @param lambda Angular weight in distance metric (default 1.0)
#' @return List with components:
#'   \item{edges}{Two-column matrix of (from, to) edge indices}
#'   \item{tau}{The arc length vector}
#' @export
acpl_graph <- function(coords, tau, K = 2, lambda = 1.0) {
  N <- nrow(coords)
  if (length(tau) != N) stop("tau must have length nrow(coords)")
  uv     <- find_adaptive_origin(coords[, 1], coords[, 2])
  theta  <- atan2(coords[, 2] - uv[2], coords[, 1] - uv[1])
  edges  <- vector("list", N)
  for (i in seq_len(N)) {
    valid <- which(tau > tau[i])
    if (length(valid) == 0) next
    d_tau   <- abs(tau[valid] - tau[i])
    d_theta <- abs(theta[valid] - theta[i])
    d_theta <- pmin(d_theta, 2 * pi - d_theta)  # angular wrap
    d       <- d_tau + lambda * d_theta
    nbr_idx <- valid[order(d)[1:min(K, length(valid))]]
    edges[[i]] <- cbind(from = i, to = nbr_idx)
  }
  list(edges = do.call(rbind, edges), tau = tau)
}

#' Compute sliding-window structural accuracy
#'
#' Primary evaluation metric. Fraction of consecutive cell pairs at fixed lag
#' for which the known phase label is non-decreasing along the arc-length
#' ordering. Random baseline approximately 50%.
#'
#' @param tau Numeric vector of arc length values
#' @param phases Character or factor vector of phase labels (G1 / S / G2M)
#' @param window Window spacing for the sliding comparison (default 10)
#' @param mirror Logical; if TRUE, automatically flip ordering if reverse SW-SA exceeds forward
#' @return Numeric SW-SA value in [0, 100]
#' @export
acpl_swsa <- function(tau, phases, window = 10, mirror = FALSE) {
  phase_map <- c(G1 = 1L, S = 2L, G2M = 3L)
  df <- data.frame(t = tau, phase = as.character(phases))
  df <- df[!is.na(df$phase) & df$phase %in% c("G1", "S", "G2M"), , drop = FALSE]
  df$phase_num <- phase_map[df$phase]
  df <- df[order(df$t), , drop = FALSE]
  n <- nrow(df)
  if (n <= window) return(NA_real_)
  fwd <- mean(df$phase_num[1:(n - window)] <=
              df$phase_num[(window + 1):n], na.rm = TRUE)
  if (mirror) {
    rev_ <- mean(df$phase_num[1:(n - window)] >=
                 df$phase_num[(window + 1):n], na.rm = TRUE)
    if (!is.na(rev_) && rev_ > fwd) fwd <- rev_
  }
  round(100 * fwd, 2)
}

#' Verify acyclicity of an ACPL graph
#'
#' Sanity check that Theorem 1 holds in practice. Should always return TRUE
#' for a properly constructed ACPL graph. Returns FALSE only if there is a bug.
#'
#' @param graph Output of acpl_graph()
#' @return Logical TRUE if the graph is acyclic
#' @export
acpl_is_acyclic <- function(graph) {
  if (is.null(graph$edges) || nrow(graph$edges) == 0) return(TRUE)
  tau <- graph$tau
  all(tau[graph$edges[, "to"]] > tau[graph$edges[, "from"]])
}

#' Compute loop circularity score
#'
#' Diagnostic measuring how circular a 2D embedding is. Loop score 1.0 = perfect
#' circle; higher values indicate elongated or fragmented manifolds. Use as a
#' pre-check before applying ACPL, though see paper for the limitations of this
#' metric (it measures circularity, not angular monotonicity).
#'
#' @param coords 2D coordinate matrix
#' @return List with loop_score, rmsd_norm, and circle parameters
#' @export
loop_score <- function(coords) {
  x <- coords[, 1]
  y <- coords[, 2]
  n  <- length(x)
  mx <- mean(x); my <- mean(y)
  xi <- x - mx; yi <- y - my
  Sxx <- sum(xi^2); Syy <- sum(yi^2); Sxy <- sum(xi * yi)
  Sxxx <- sum(xi^3); Syyy <- sum(yi^3)
  Sxyy <- sum(xi * yi^2); Sxxy <- sum(xi^2 * yi)
  A <- matrix(c(Sxx, Sxy, sum(xi),
                Sxy, Syy, sum(yi),
                sum(xi), sum(yi), n), 3, 3, byrow = TRUE)
  b <- c(-(Sxxx + Sxyy), -(Syyy + Sxxy), -(Sxx + Syy))
  sol <- tryCatch(solve(A, b), error = function(e) NULL)
  if (is.null(sol)) return(list(loop_score = NA, rmsd_norm = NA))
  cx <- -sol[1] / 2 + mx; cy <- -sol[2] / 2 + my
  r  <- sqrt(sol[1]^2 / 4 + sol[2]^2 / 4 - sol[3])
  dist_to_circle <- sqrt((x - cx)^2 + (y - cy)^2)
  rmsd <- sqrt(mean((dist_to_circle - r)^2))
  theta_p <- atan2(y - cy, x - cx)
  ord     <- order(theta_p)
  arc_len <- sum(sqrt(diff(x[ord])^2 + diff(y[ord])^2))
  list(
    loop_score = arc_len / (2 * pi * r),
    rmsd_norm  = rmsd / r,
    circle     = list(cx = cx, cy = cy, r = r)
  )
}
