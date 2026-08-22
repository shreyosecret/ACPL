# ==============================================================================
# R/metrics.R -- evaluation statistics
#
# NEW in v2. The v1 release scored everything with sliding-window structural
# accuracy against a stated 50% random baseline. That baseline is wrong; see
# chance_level() below and CHANGELOG.md.
# ==============================================================================

PHASE_MAP <- c(G1 = 1L, S = 2L, G2M = 3L)

clean_phases <- function(p) { p <- as.character(p); p[p %in% names(PHASE_MAP)] }

#' Chance level of tie-counting windowed concordance
#'
#' SW-SA counts ties as successes, so under a random ordering of the same label
#' multiset E[SW-SA] = (1 + sum p^2)/2, bounded below by 2/3 for three phase
#' levels. It is never 50 percent.
chance_level <- function(phases) {
  p <- prop.table(table(factor(clean_phases(phases), levels = names(PHASE_MAP))))
  100 * (1 + sum(p^2)) / 2
}

#' Chance-corrected agreement (Cohen 1960, applied to a cyclic ordering)
kappa_like <- function(obs, chance) 100 * (obs - chance) / (100 - chance)

#' Sliding-window structural accuracy
#'
#' Retained for comparison with v1 only. Do not use as a primary metric: it is
#' dominated by the tie rate and is insensitive to differences below about
#' 0.16 in tau_cyc. `mirror` is deliberately absent; resolving orientation by
#' whichever direction scores higher against the labels uses the answer key to
#' orient the answer.
swsa <- function(tau, phases, window = 10) {
  d <- data.frame(t = tau, p = as.character(phases))
  d <- d[!is.na(d$t) & d$p %in% names(PHASE_MAP), , drop = FALSE]
  d$pn <- PHASE_MAP[d$p]
  d <- d[order(d$t), , drop = FALSE]
  n <- nrow(d); if (n <= window) return(NA_real_)
  100 * mean(d$pn[1:(n - window)] <= d$pn[(window + 1):n])
}

#' Decompose SW-SA into its tie component and its directional component
swsa_decompose <- function(tau, phases, window = 10) {
  d <- data.frame(t = tau, p = as.character(phases))
  d <- d[!is.na(d$t) & d$p %in% names(PHASE_MAP), , drop = FALSE]
  d$pn <- PHASE_MAP[d$p]; d <- d[order(d$t), , drop = FALSE]; n <- nrow(d)
  a <- d$pn[1:(n - window)]; b <- d$pn[(window + 1):n]
  tie <- mean(a == b); s <- mean(a <= b)
  c(swsa = 100 * s, ties = 100 * tie, no_direction = 100 * (1 + tie) / 2,
    directional = 100 * (s - (1 + tie) / 2))
}

# ---- origin-free rank statistic ---------------------------------------------
# Concordant minus discordant pairs, O(n). Positions are distinct so there are
# no ties on the x side.
.cd_stat <- function(pn) {
  b1 <- c(0, head(cumsum(pn == 1L), -1))
  b2 <- c(0, head(cumsum(pn == 2L), -1))
  b3 <- c(0, head(cumsum(pn == 3L), -1))
  sum(ifelse(pn == 1L, 0, ifelse(pn == 2L, b1, b1 + b2))) -
  sum(ifelse(pn == 1L, b2 + b3, ifelse(pn == 2L, b3, 0)))
}

#' Origin-free cyclic Kendall tau
#'
#' Maximum tau_b over a grid of origin rotations and both traversal directions.
#' Required because a method that takes no starting cluster is penalised by a
#' fixed G1 < S < G2M scale: the same perfect cyclic ordering scores +0.82 or
#' -0.27 depending only on where it is taken to begin. Apply the identical
#' maximisation under permutation so the selection is priced into the null.
tau_cyc <- function(pn_ordered, nrot = 24) {
  n <- length(pn_ordered); n0 <- n * (n - 1) / 2
  n2 <- sum(vapply(1:3, function(k) { m <- sum(pn_ordered == k); m * (m - 1) / 2 },
                   numeric(1)))
  den <- sqrt(n0 * (n0 - n2))
  if (!is.finite(den) || den <= 0) return(NA_real_)
  cuts <- unique(round(seq(0, n - 1, length.out = nrot + 1)[1:nrot]))
  best <- -Inf
  for (r in cuts) {
    rolled <- if (r == 0) pn_ordered else c(pn_ordered[(r + 1):n], pn_ordered[1:r])
    for (s in list(rolled, rev(rolled))) {
      v <- .cd_stat(s) / den
      if (v > best) best <- v
    }
  }
  best
}

#' tau_cyc from a pseudotime vector and phase labels
tau_of <- function(pseudotime, phases, nrot = 24) {
  ok <- phases %in% names(PHASE_MAP) & !is.na(pseudotime)
  tau_cyc(PHASE_MAP[phases[ok]][order(pseudotime[ok])], nrot = nrot)
}

# ---- circular helpers --------------------------------------------------------
circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))
circ_R    <- function(a) sqrt(mean(sin(a))^2 + mean(cos(a))^2)
rank_angle <- function(v) 2 * pi * (rank(v, ties.method = "first") - 0.5) / length(v)

#' Jammalamadaka-Sarma circular correlation
#'
#' WARNING: degenerate when one variable takes three equally spaced values with
#' near-equal weights, because the circular mean resultant approaches zero.
#' Do not apply to three-level phase labels. Valid for continuous ground truth,
#' e.g. experimental timepoints.
circ_cor <- function(a, b) {
  da <- a - circ_mean(a); db <- b - circ_mean(b)
  sum(sin(da) * sin(db)) / sqrt(sum(sin(da)^2) * sum(sin(db)^2))
}

#' Per-phase concentration and cyclic order verdict, both origin-free
phase_geometry <- function(pseudotime, phases) {
  ok <- phases %in% names(PHASE_MAP) & !is.na(pseudotime)
  pn <- PHASE_MAP[phases[ok]]; ang <- rank_angle(pseudotime[ok])
  mu <- vapply(1:3, function(k) circ_mean(ang[pn == k]), numeric(1))
  Rk <- vapply(1:3, function(k) circ_R(ang[pn == k]), numeric(1))
  g <- c((mu[2] - mu[1]) %% (2*pi), (mu[3] - mu[2]) %% (2*pi), (mu[1] - mu[3]) %% (2*pi))
  verdict <- if (abs(sum(g) - 2*pi) < 1e-6) "correct" else
             if (abs(sum(g) - 4*pi) < 1e-6) "reversed-handedness" else "indeterminate"
  list(mean_angle = mu, resultant = Rk, gaps = g,
       concentration = mean(Rk), cyclic_order = verdict)
}
