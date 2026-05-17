# ==============================================================================
# minimal_example.R — 50-line reproducible ACPL demonstration
#
# Generates a noisy synthetic cyclic manifold, applies ACPL, and verifies
# the two formal theorems hold. Run from repository root:
#
#   Rscript examples/minimal_example.R
# ==============================================================================

source("R/acpl_engine.R")
library(igraph); library(ggplot2)

# 1. Generate a noisy synthetic cell cycle
set.seed(42)
N <- 50
theta_true <- seq(0, 2 * pi * (1 - 1/N), length.out = N)
r <- 1 + rnorm(N, 0, 0.05)
x <- r * cos(theta_true) + rnorm(N, 0, 0.15)
y <- r * sin(theta_true) + rnorm(N, 0, 0.15)
coords <- cbind(x, y)

# 2. Apply ACPL
arc <- acpl_arc(coords, span = 0.25)
graph <- acpl_graph(coords, arc, K = 2)

# 3. Verify Theorem 1: graph is acyclic
cat("Theorem 1 (acyclicity by construction):",
    if (acpl_is_acyclic(graph)) "PASS" else "FAIL", "\n")

# 4. Verify Theorem 2: raw arc length is degenerate
ord <- order(atan2(y - mean(y), x - mean(x)))
arc_raw <- c(0, cumsum(sqrt(diff(x[ord])^2 + diff(y[ord])^2)))
cat("Theorem 2 (raw arc length is strictly increasing):",
    if (all(diff(arc_raw) > 0)) "PASS" else "FAIL", "\n")
cat("  -> Raw arc length filter would admit",
    sum(diff(arc_raw) > 0), "of", N - 1,
    "edges (100% by mathematical necessity, not inference)\n")

# 5. Score against ground truth
true_idx <- seq_len(N)
forward_frac <- mean(true_idx[graph$edges[, "to"]] >
                     true_idx[graph$edges[, "from"]])
cat(sprintf("ACPL structural accuracy: %.1f%%\n",
            100 * max(forward_frac, 1 - forward_frac)))

# 6. Visualise
df_pts <- data.frame(x = x, y = y, idx = true_idx)
df_edges <- data.frame(
  x1 = x[graph$edges[, "from"]], y1 = y[graph$edges[, "from"]],
  x2 = x[graph$edges[, "to"]],   y2 = y[graph$edges[, "to"]])
p <- ggplot() +
  geom_segment(data = df_edges,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "#7b5ea7", alpha = 0.6, linewidth = 0.5,
               arrow = arrow(length = unit(0.12, "cm"))) +
  geom_point(data = df_pts, aes(x = x, y = y, colour = idx), size = 2.5) +
  scale_colour_viridis_c(option = "plasma", name = "True\norder") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "ACPL on a noisy synthetic cycle",
       subtitle = sprintf("SA = %.1f%%, %d directed edges",
                          100 * max(forward_frac, 1 - forward_frac),
                          nrow(graph$edges)))

if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/minimal_example.png", p, width = 6, height = 5, dpi = 150)
cat("\nFigure: figures/minimal_example.png\n")
