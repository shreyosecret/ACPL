# ==============================================================================
# test_methods.R -- unit tests for R/methods.R and R/metrics.R
#
# Run with:  testthat::test_dir("tests/testthat")  from the repo root
# ==============================================================================

library(testthat)
source("../../R/methods.R")
source("../../R/metrics.R")

# ── find_adaptive_origin ─────────────────────────────────────────────────────
test_that("find_adaptive_origin recovers centroid on a perfect circle", {
  set.seed(1)
  theta <- seq(0, 2 * pi, length.out = 100)
  x <- 5 + cos(theta)
  y <- 3 + sin(theta)
  origin <- find_adaptive_origin(x, y)
  expect_equal(origin[1], 5, tolerance = 0.05)
  expect_equal(origin[2], 3, tolerance = 0.05)
})

# ── unwrap_theta ─────────────────────────────────────────────────────────────
test_that("unwrap_theta produces monotonically increasing output", {
  theta <- seq(-pi, pi, length.out = 50)
  unwrapped <- unwrap_theta(theta)
  expect_true(all(diff(unwrapped) > 0))
})

test_that("unwrap_theta handles wrap-around correctly", {
  # Constructed sequence with a 2*pi jump
  theta <- c(0, 1, 2, -3, -2, -1)  # actually wraps from +pi area to -pi area
  unwrapped <- unwrap_theta(theta)
  # The unwrap should not contain large jumps in the unwrapped diff
  expect_true(all(abs(diff(unwrapped)) < pi))
})

# ── acpl_arc ─────────────────────────────────────────────────────────────────
test_that("acpl_arc returns numeric vector of correct length", {
  set.seed(1)
  theta <- seq(0, 2 * pi, length.out = 100)
  coords <- cbind(cos(theta) + rnorm(100, 0, 0.05),
                  sin(theta) + rnorm(100, 0, 0.05))
  arc <- acpl_arc(coords)
  expect_type(arc, "double")
  expect_length(arc, 100)
  expect_true(all(is.finite(arc)))
})

test_that("acpl_arc rejects non-2D input", {
  expect_error(acpl_arc(matrix(1:30, ncol = 3)), "2D")
  expect_error(acpl_arc(matrix(1:10, ncol = 1)), "2D")
})

# ── acpl_graph ───────────────────────────────────────────────────────────────
test_that("acpl_graph produces edges only forward in tau", {
  set.seed(1)
  theta <- seq(0, 2 * pi, length.out = 50)
  coords <- cbind(cos(theta), sin(theta))
  arc <- acpl_arc(coords)
  graph <- acpl_graph(coords, arc, K = 2)
  expect_true(!is.null(graph$edges))
  expect_true(all(arc[graph$edges[, "to"]] > arc[graph$edges[, "from"]]))
})

test_that("acpl_graph rejects mismatched tau length", {
  coords <- cbind(rnorm(50), rnorm(50))
  expect_error(acpl_graph(coords, rnorm(40)), "length")
})

# ── acpl_is_acyclic (Theorem 1) ──────────────────────────────────────────────
test_that("acpl_is_acyclic returns TRUE on properly constructed graph", {
  set.seed(1)
  theta <- seq(0, 2 * pi, length.out = 50)
  coords <- cbind(cos(theta) + rnorm(50, 0, 0.1),
                  sin(theta) + rnorm(50, 0, 0.1))
  arc <- acpl_arc(coords)
  graph <- acpl_graph(coords, arc, K = 3)
  expect_true(acpl_is_acyclic(graph))
})

# ── acpl_swsa ────────────────────────────────────────────────────────────────
test_that("acpl_swsa returns 100 for perfectly ordered phases", {
  pseudotime <- 1:90
  phases <- rep(c("G1", "S", "G2M"), each = 30)
  expect_equal(acpl_swsa(pseudotime, phases, window = 10), 100)
})

test_that("acpl_swsa returns NA for very small N", {
  expect_true(is.na(acpl_swsa(1:5, rep(c("G1", "S", "G2M"), length.out = 5),
                               window = 10)))
})

# ── loop_score ───────────────────────────────────────────────────────────────
test_that("loop_score returns approximately 1.0 for a perfect circle", {
  theta <- seq(0, 2 * pi, length.out = 200)
  coords <- cbind(cos(theta), sin(theta))
  ls <- loop_score(coords)
  expect_equal(ls$loop_score, 1.0, tolerance = 0.05)
  expect_lt(ls$rmsd_norm, 0.01)
})

test_that("loop_score returns sensible result for noisy circle", {
  set.seed(1)
  theta <- seq(0, 2 * pi, length.out = 100)
  coords <- cbind(cos(theta) + rnorm(100, 0, 0.1),
                  sin(theta) + rnorm(100, 0, 0.1))
  ls <- loop_score(coords)
  # Loop score will be modestly above 1.0 due to noise
  expect_gt(ls$loop_score, 1.0)
  expect_lt(ls$loop_score, 1.5)
})


# -- chance level, new in v2 ---------------------------------------------------
test_that("chance level is (1 + sum p^2)/2 and never 50 percent", {
  balanced <- rep(c("G1", "S", "G2M"), each = 100)
  expect_equal(chance_level(balanced), 100 * (1 + 1/3) / 2, tolerance = 1e-8)
  expect_gt(chance_level(balanced), 66)
  skewed <- c(rep("G1", 700), rep("S", 150), rep("G2M", 150))
  expect_gt(chance_level(skewed), chance_level(balanced))
})

test_that("a random ordering scores near the analytic chance level", {
  set.seed(42)
  ph <- sample(rep(c("G1", "S", "G2M"), each = 400))
  obs <- mean(replicate(200, swsa(runif(length(ph)), ph)))
  expect_equal(obs, chance_level(ph), tolerance = 1.5)
})

# -- tau_cyc is origin-free ----------------------------------------------------
test_that("tau_cyc is invariant to where the cycle is taken to begin", {
  pn <- rep(1:3, each = 60)
  rolled <- c(pn[61:180], pn[1:60])
  expect_equal(tau_cyc(pn), tau_cyc(rolled), tolerance = 1e-6)
})

test_that("raw arc length is strictly increasing, i.e. degenerate", {
  set.seed(7)
  th <- sort(runif(200, -pi, pi))
  cd <- cbind(cos(th), sin(th))
  raw <- c(0, cumsum(sqrt(diff(cd[order(th), 1])^2 + diff(cd[order(th), 2])^2)))
  expect_true(all(diff(raw) > 0))
})

# -- acyclicity ----------------------------------------------------------------
test_that("the constructed graph is acyclic", {
  set.seed(3)
  th <- seq(0, 2 * pi, length.out = 150)
  cd <- cbind(cos(th) + rnorm(150, 0, .03), sin(th) + rnorm(150, 0, .03))
  g <- acpl_graph(cd, acpl_arc(cd))
  expect_true(acpl_is_acyclic(g))
})
