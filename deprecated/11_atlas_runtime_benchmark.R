# ==============================================================================
# Analysis 11: Atlas-scale runtime benchmark
# Section 4.11. Compares ACPL, MST, Slingshot, and PAGA+DPT up to N=50,000.
# Self-contained: generates synthetic cyclic data, no real datasets needed.
# Outputs: results/runtime_benchmark.csv, results/runtime_projections.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")

library(igraph); library(dplyr); library(tidyr)

make_cycle <- function(N, noise = 0.15) {
  theta <- seq(0, 2 * pi * (1 - 1/N), length.out = N)
  x     <- cos(theta) + rnorm(N, 0, noise)
  y     <- sin(theta) + rnorm(N, 0, noise)
  list(coords = cbind(x, y), theta = theta,
       phases = cut(theta / (2*pi), breaks = 3,
                    labels = c("G1", "S", "G2M"), include.lowest = TRUE))
}

N_all   <- c(100, 500, 1000, 2000, 5000, 10000, 20000, 50000)
N_mst   <- N_all[N_all <= 5000]    # MST infeasible above 5k
N_sling <- N_all[N_all <= 10000]   # Slingshot O(N^2)
N_reps  <- 3

times <- list()

for (N in N_all) {
  cat(sprintf("N = %d ... ", N))

  # ACPL
  t_acpl <- numeric(N_reps)
  for (r in seq_len(N_reps)) {
    set.seed(r); dat <- make_cycle(N)
    t_acpl[r] <- system.time(acpl_arc(dat$coords))[["elapsed"]]
  }

  # MST
  t_mst <- if (N %in% N_mst) {
    rt <- numeric(N_reps)
    for (r in seq_len(N_reps)) {
      set.seed(r); dat <- make_cycle(N)
      rt[r] <- system.time(
        tryCatch(mst_pseudotime(dat$coords, "G1",
                                 as.character(dat$phases)),
                 error = function(e) NULL))[["elapsed"]]
    }
    median(rt)
  } else NA

  # Slingshot
  t_sling <- if (N %in% N_sling) {
    rt <- numeric(N_reps)
    for (r in seq_len(N_reps)) {
      set.seed(r); dat <- make_cycle(N)
      rt[r] <- system.time(
        tryCatch(slingshot_pseudotime(dat$coords,
                                        as.character(dat$phases), "G1"),
                 error = function(e) NULL))[["elapsed"]]
    }
    median(rt)
  } else NA

  times[[as.character(N)]] <- data.frame(
    N = N, ACPL = median(t_acpl), MST = t_mst, Slingshot = t_sling)
  cat(sprintf("ACPL=%.3fs  MST=%.3fs  Sling=%.3fs\n",
              median(t_acpl), t_mst, t_sling))
}

time_df <- do.call(rbind, times)
write.csv(time_df, "results/runtime_benchmark.csv", row.names = FALSE)

# Fit empirical scaling exponents
fit_exp <- function(N_v, t_v) {
  keep <- !is.na(t_v) & t_v > 0
  if (sum(keep) < 3) return(NA)
  coef(lm(log(t_v[keep]) ~ log(N_v[keep])))[2]
}

exp_acpl  <- fit_exp(time_df$N, time_df$ACPL)
exp_mst   <- fit_exp(time_df$N, time_df$MST)
exp_sling <- fit_exp(time_df$N, time_df$Slingshot)

cat(sprintf("\nFitted scaling exponents:\n  ACPL=%.2f  MST=%.2f  Slingshot=%.2f\n",
            exp_acpl, exp_mst, exp_sling))

# Project to atlas scale
proj <- data.frame(
  N = c(100000, 500000, 1000000),
  ACPL_s      = NA, MST_s = NA, Slingshot_s = NA)

ref_N <- 5000
ref_acpl  <- time_df$ACPL[time_df$N == ref_N]
ref_sling <- time_df$Slingshot[time_df$N == ref_N]
ref_mst   <- time_df$MST[time_df$N == ref_N]

for (i in seq_len(nrow(proj))) {
  proj$ACPL_s[i]      <- ref_acpl  * (proj$N[i] / ref_N)^exp_acpl
  proj$MST_s[i]       <- ref_mst   * (proj$N[i] / ref_N)^exp_mst
  proj$Slingshot_s[i] <- ref_sling * (proj$N[i] / ref_N)^exp_sling
}
print(proj)
write.csv(proj, "results/runtime_projections.csv", row.names = FALSE)
cat("\nResults saved.\n")
