# ==============================================================================
# Analysis 12: 3D helix synthetic validation
# Section 4.12. Tests whether the 3D helix lift resolves cross-cycle confusion.
# Self-contained synthetic data; no real datasets required.
# Outputs: results/helix_validation.csv
# ==============================================================================

source("R/acpl_engine.R")
library(dplyr)

set.seed(42)
N_per_cycle <- 100
n_cycles    <- 2
noise       <- 0.08
k_helix     <- 0.5
window      <- 10
N_total     <- N_per_cycle * n_cycles

# Two-cycle ground truth
theta <- seq(0, 2 * pi * n_cycles * (1 - 1/N_total), length.out = N_total)
true_time <- seq_len(N_total)
phase_frac <- (theta %% (2 * pi)) / (2 * pi)
phases <- cut(phase_frac, breaks = c(0, 0.33, 0.66, 1.0),
              labels = c("G1", "S", "G2M"), include.lowest = TRUE)
r <- 1 + rnorm(N_total, 0, noise)

# 2D collapse (both cycles overlap)
x_2d <- r * cos(theta %% (2 * pi))
y_2d <- r * sin(theta %% (2 * pi))

# 3D helix (z = k * theta separates cycles)
x_3d <- r * cos(theta)
y_3d <- r * sin(theta)
z_3d <- k_helix * theta

# 2D control (single cycle)
N_single <- N_per_cycle
theta_s  <- seq(0, 2*pi*(1-1/N_single), length.out = N_single)
r_s      <- 1 + rnorm(N_single, 0, noise)
x_s      <- r_s * cos(theta_s); y_s <- r_s * sin(theta_s)
phase_s  <- cut((theta_s %% (2*pi)) / (2*pi),
                breaks = c(0, 0.33, 0.66, 1.0),
                labels = c("G1", "S", "G2M"), include.lowest = TRUE)

# 3D ACPL (use x,y for angle, include z in arc length)
acpl_3d_arc <- function(x, y, z, span = 0.25) {
  uv <- find_adaptive_origin(x, y)
  u <- x - uv[1]; v <- y - uv[2]
  th <- atan2(v, u); ord <- order(th); tu <- unwrap_theta(th[ord])
  arc_raw <- c(0, cumsum(sqrt(diff(x[ord])^2 + diff(y[ord])^2 + diff(z[ord])^2)))
  arc_smooth <- predict(loess(arc_raw ~ tu, span = span))
  arc_smooth[order(seq_along(x)[ord])]
}

# Run all three conditions
res <- list()

# 2D single cycle (control)
arc_s <- acpl_arc(cbind(x_s, y_s))
sw_s  <- acpl_swsa(arc_s, phase_s)
sp_s  <- cor(arc_s, seq_len(N_single), method = "spearman")
res[["2D control"]] <- list(N = N_single, sw = sw_s, sp = sp_s)

# 2D two cycles collapsed
arc_2d <- acpl_arc(cbind(x_2d, y_2d))
sw_2d  <- acpl_swsa(arc_2d, phases)
sp_2d  <- cor(arc_2d, true_time, method = "spearman")
res[["2D 2 cycles"]] <- list(N = N_total, sw = sw_2d, sp = sp_2d)

# 3D helix
arc_3d <- acpl_3d_arc(x_3d, y_3d, z_3d)
sw_3d  <- acpl_swsa(arc_3d, phases)
sp_3d  <- cor(arc_3d, true_time, method = "spearman")
res[["3D helix"]] <- list(N = N_total, sw = sw_3d, sp = sp_3d)

out <- do.call(rbind, lapply(names(res), function(n)
  data.frame(condition = n, N = res[[n]]$N,
             SW_SA = round(res[[n]]$sw, 1),
             Spearman_rho = round(res[[n]]$sp, 4))))
print(out, row.names = FALSE)
write.csv(out, "results/helix_validation.csv", row.names = FALSE)
cat("\nResults saved: results/helix_validation.csv\n")
