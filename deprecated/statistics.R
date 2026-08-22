# ==============================================================================
# Statistical helpers: bootstrap CIs, Moran's I, mixed model wrappers
# ==============================================================================

#' Bootstrap CI for SW-SA difference between two methods
#'
#' Resamples per-window binary correctness vectors (the right level for
#' assessing window-to-window dependence) and computes BCa CI.
#'
#' @param windows_a Integer 0/1 correctness vector for method A
#' @param windows_b Integer 0/1 correctness vector for method B
#' @param label_a Method A name (string)
#' @param label_b Method B name (string)
#' @param dataset Dataset name (string)
#' @param R Number of bootstrap iterations (default 10000)
#' @param seed Random seed for reproducibility
#' @return Data frame with observed difference, BCa CI bounds, p-value,
#'   and a Sig flag for whether the CI excludes zero
#' @export
bootstrap_ci <- function(windows_a, windows_b, label_a, label_b, dataset,
                         R = 10000, seed = 2025) {
  if (!requireNamespace("boot", quietly = TRUE))
    stop("Install boot: install.packages('boot')")
  n <- min(length(windows_a), length(windows_b))
  windows_a <- windows_a[1:n]
  windows_b <- windows_b[1:n]
  df <- data.frame(a = windows_a, b = windows_b)
  set.seed(seed)
  br <- boot::boot(df, function(d, i) mean(d$a[i]) - mean(d$b[i]), R = R)
  ci <- boot::boot.ci(br, type = c("perc", "bca"), conf = 0.95)
  wt <- wilcox.test(windows_a, windows_b, paired = TRUE,
                    alternative = "two.sided", exact = FALSE)
  data.frame(
    Dataset    = dataset,
    Method_A   = label_a,
    Method_B   = label_b,
    Obs_pp     = round(100 * (mean(windows_a) - mean(windows_b)), 3),
    CI_lo_bca  = round(100 * ci$bca[4], 3),
    CI_hi_bca  = round(100 * ci$bca[5], 3),
    Wilcoxon_p = signif(wt$p.value, 3),
    N          = n,
    Sig        = (ci$bca[4] > 0 | ci$bca[5] < 0)
  )
}

#' Convert pseudotime + phases to per-window correctness vector
#'
#' Used as input to bootstrap_ci. Sorts cells by pseudotime, optionally
#' applying a mirror flip if reverse ordering scores higher.
#'
#' @param pseudotime Numeric vector of pseudotime values
#' @param phases Phase labels
#' @param window Window spacing
#' @param mirror Logical: auto-flip ordering if reversed scores higher
#' @return Integer 0/1 vector of per-window correctness
#' @export
window_correctness <- function(pseudotime, phases, window = 10, mirror = FALSE) {
  phase_map <- c(G1 = 1L, S = 2L, G2M = 3L)
  df <- data.frame(t = pseudotime, p = as.character(phases))
  df <- df[!is.na(df$p) & df$p %in% c("G1", "S", "G2M"), , drop = FALSE]
  df$pn <- phase_map[df$p]
  df <- df[order(df$t), , drop = FALSE]
  n <- nrow(df)
  if (mirror) {
    fwd <- mean(df$pn[1:(n - window)] <= df$pn[(window + 1):n], na.rm = TRUE)
    rev_<- mean(df$pn[1:(n - window)] >= df$pn[(window + 1):n], na.rm = TRUE)
    if (!is.na(rev_) && rev_ > fwd) df <- df[order(-df$t), , drop = FALSE]
  }
  as.integer(df$pn[1:(n - window)] <= df$pn[(window + 1):n])
}

#' Moran's I for spatial autocorrelation along a 1D ordering
#'
#' Used for the biological anchor validation. Tests whether marker expression
#' is autocorrelated along the arc length axis, a prerequisite for using mixed
#' models rather than naive ANOVA.
#'
#' @param values Numeric vector of marker expression values
#' @param order_var Numeric vector defining the 1D ordering (e.g. arc length)
#' @param n_neighbours Number of neighbours per cell in the spatial weight matrix
#' @return List with Moran's I statistic and p-value
#' @export
moran_i_1d <- function(values, order_var, n_neighbours = 5) {
  if (!requireNamespace("ape", quietly = TRUE))
    stop("Install ape: install.packages('ape')")
  ord <- order(order_var)
  v   <- values[ord]
  n   <- length(v)
  W   <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- setdiff(max(1, i - n_neighbours):min(n, i + n_neighbours), i)
    W[i, nbrs] <- 1
  }
  W <- W / rowSums(W)
  res <- ape::Moran.I(v, W)
  list(I = res$observed, p = res$p.value)
}

#' Three-point local curvature
#'
#' Computes the reciprocal radius of the circle through three coordinate
#' points, used in the curvature analysis. Returns 0 for nearly collinear
#' points (curvature undefined).
#'
#' @param A First point (numeric vector of length 2)
#' @param B Middle point
#' @param C Third point
#' @return Numeric curvature value (1/R), or Inf if collinear
#' @export
three_point_curvature <- function(A, B, C) {
  a <- sqrt(sum((B - C)^2))
  b <- sqrt(sum((A - C)^2))
  c <- sqrt(sum((A - B)^2))
  s <- (a + b + c) / 2
  area <- sqrt(max(s * (s - a) * (s - b) * (s - c), 0))
  if (area < 1e-10) return(Inf)
  R <- (a * b * c) / (4 * area)
  1 / R
}
