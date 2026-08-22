# ==============================================================================
# Analysis 01: Synthetic stress test
# Section 4.1 of the paper. Generates 50 cyclic manifolds at three noise levels
# and compares ACPL to Cartesian MST on each seed.
# Outputs: results/synthetic_stress.csv, figures/fig_stress_test.png
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("R/plotting.R")

library(igraph); library(ggplot2); library(dplyr); library(tidyr)

set.seed(42)
N_seeds <- 50
N_nodes <- 30
noise_levels <- c(0.3, 0.8, 1.5)

results <- list()

for (noise in noise_levels) {
  acpl_acc <- numeric(N_seeds); cart_acc <- numeric(N_seeds)
  for (seed in seq_len(N_seeds)) {
    set.seed(seed * 11 + as.integer(noise * 100))
    theta <- seq(0, 2 * pi, length.out = N_nodes + 1)[-1]
    r     <- 1 + rnorm(N_nodes, 0, 0.05)
    x <- r * cos(theta) + rnorm(N_nodes, 0, noise)
    y <- r * sin(theta) + rnorm(N_nodes, 0, noise)
    coords <- cbind(x, y)
    true_idx <- seq_len(N_nodes)

    arc <- acpl_arc(coords)
    g_acpl <- acpl_graph(coords, arc, K = 2)
    if (!is.null(g_acpl$edges)) {
      forward <- mean(true_idx[g_acpl$edges[, "to"]] >
                      true_idx[g_acpl$edges[, "from"]])
      acpl_acc[seed] <- 100 * forward
    }

    dm <- as.matrix(dist(coords))
    g  <- igraph::graph_from_adjacency_matrix(dm, weighted = TRUE,
                                              mode = "undirected")
    mst_t <- igraph::mst(g)
    el <- igraph::as_edgelist(mst_t)
    forward_c <- mean(true_idx[as.integer(el[, 2])] >
                      true_idx[as.integer(el[, 1])])
    cart_acc[seed] <- 100 * max(forward_c, 1 - forward_c)
  }
  wt <- wilcox.test(acpl_acc, cart_acc, paired = TRUE,
                    alternative = "greater", exact = FALSE)
  results[[as.character(noise)]] <- data.frame(
    noise        = noise,
    cartesian    = sprintf("%.1f%%%s%.1f%%", mean(cart_acc), "\u00b1", sd(cart_acc)),
    acpl         = sprintf("%.1f%%%s%.1f%%", mean(acpl_acc), "\u00b1", sd(acpl_acc)),
    W            = wt$statistic,
    p            = signif(wt$p.value, 3)
  )
  cat(sprintf("Noise %.1f: ACPL %.1f%% vs Cartesian %.1f%% (W=%.0f, p=%.2e)\n",
              noise, mean(acpl_acc), mean(cart_acc), wt$statistic, wt$p.value))
}

out <- do.call(rbind, results)
write.csv(out, "results/synthetic_stress.csv", row.names = FALSE)
cat("\nResults saved: results/synthetic_stress.csv\n")
