# ==============================================================================
# demo_circ_cor_tradeoff.R
#
# Circular correlation cannot be both origin-free and direction-sensitive.
#
# The recommended protocol asks for two things that turn out to be incompatible
# for this statistic:
#
#   (i)  use a direction-sensitive statistic, and circular correlation is the
#        natural choice for a cyclic process;
#   (ii) when a root-free method is in the comparison, use an origin-free
#        statistic, since one that assumes a fixed start favours the methods
#        that are handed one.
#
# This script shows, with no data and in a few seconds, that no version of the
# Jammalamadaka-Sarma coefficient satisfies both.
#
#   Rscript demo_circ_cor_tradeoff.R
#
# RESULT 1. Against a rank-based pseudotime angle the coefficient is not merely
#   noisy, it is arbitrary. The angle is uniform on the circle by construction,
#   so its resultant is zero, circ_mean() returns floating-point noise, and the
#   reference angle is undetermined. A PERFECT cyclic ordering scores anywhere
#   from about -0.35 to +0.83 depending only on where rank one happens to land.
#
# RESULT 2. The obvious repair, maximising over the reference angle, provably
#   destroys direction sensitivity. Rotating the origin by pi negates the
#   coefficient exactly, so reversal and a half-turn of the origin are the same
#   operation. Maximising over the origin therefore scores a trajectory and its
#   reverse identically, to machine precision.
#
# CONSEQUENCE. The direction problem cannot be solved from inside the statistic
#   when the method supplies no origin. It has to be solved before the statistic
#   is applied, by fixing orientation against an external anchor such as marker
#   gene expression. For methods that are given a start cluster the origin is
#   already fixed and the coefficient is fine as it stands.
# ==============================================================================

set.seed(2026)
NROT <- 720

circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))
circ_cor  <- function(a, b) {
  da <- a - circ_mean(a); db <- b - circ_mean(b)
  sum(sin(da) * sin(db)) / sqrt(sum(sin(da)^2) * sum(sin(db)^2))
}
# the coefficient at an explicit reference angle c
cc_at <- function(a, b, c0) {
  da <- a - c0; db <- b - circ_mean(b)
  sum(sin(da) * db) / sqrt(sum(sin(da)^2) * sum(db^2))
}
max_over_origin <- function(a, b, nrot = NROT) {
  cs <- seq(0, 2 * pi, length.out = nrot + 1)[-(nrot + 1)]
  max(vapply(cs, function(c0) cc_at(a, b, c0), numeric(1)))
}

p_ang <- function(pn) c(0, 2 * pi / 3, 4 * pi / 3)[pn]

cat("== RESULT 1: a perfect ordering scores arbitrarily ==\n")
cat("   n     naive circ_cor on a PERFECT G1->S->G2M ordering\n")
for (n in c(247, 288, 572, 1920)) {
  pn <- sort(sample(rep(1:3, length.out = n)))
  b  <- p_ang(pn)
  a  <- 2 * pi * (seq_len(n) - 0.5) / n
  cat(sprintf("  %5d   %+.4f\n", n, circ_cor(a, b)))
}
cat("  The ordering is identical in every row. Only n differs.\n")

cat("\n== RESULT 2: a pi rotation of the origin negates the coefficient ==\n")
n  <- 572; pn <- sample(rep(1:3, length.out = n)); b <- p_ang(pn)
a  <- 2 * pi * (rank(runif(n), ties.method = "average") - 0.5) / n
for (c0 in c(0.0, 0.7, 2.3)) {
  v1 <- cc_at(a, b, c0); v2 <- cc_at(a, b, c0 + pi)
  cat(sprintf("  c = %.2f -> %+.8f    c + pi -> %+.8f    sum = %.1e\n",
              c0, v1, v2, v1 + v2))
}
cat("  The sums are zero to machine precision.\n")

cat("\n== CONSEQUENCE: origin-maximised circular correlation is direction-blind ==\n")
cat("   n     forward    reversed      gap\n")
worst <- 0
for (n in c(247, 288, 572, 1920)) {
  pn <- sort(sample(rep(1:3, length.out = n)))
  b  <- p_ang(pn)
  a  <- 2 * pi * (seq_len(n) - 0.5) / n
  ar <- 2 * pi * (rev(seq_len(n)) - 0.5) / n
  f  <- max_over_origin(a, b); r <- max_over_origin(ar, b)
  worst <- max(worst, abs(f - r))
  cat(sprintf("  %5d   %+.6f  %+.6f   %.2e\n", n, f, r, abs(f - r)))
}
cat(sprintf("\n  largest gap across all n: %.2e\n", worst))

cat("\n---------------------------------------------------------------\n")
if (worst < 1e-6) {
  cat("VERDICT: origin-free and direction-sensitive are incompatible for this\n")
  cat("statistic. Orientation must be fixed by an external anchor before the\n")
  cat("coefficient is applied; no maximisation inside it can recover direction.\n")
} else {
  cat("VERDICT: the gap is not negligible. Re-examine the argument.\n")
}
