# ==============================================================================
# run_all.R — master script to reproduce every analysis in the paper
#
# Usage from repository root:
#   Rscript run_all.R
#
# This script reproduces every figure and table reported in:
#   Ghosh, S. (2026). ACPL: A Coordinate-System Prior for Topology-Aware
#   Directed Trajectory Inference on Cyclic Biological Manifolds.
#   Bioinformatics Advances (submitted).
#
# Total runtime: ~90 min on a modern laptop (4-core, 16GB RAM).
# Bottlenecks: atlas runtime benchmark (~45 min) and UMAP sensitivity (~30 min).
#
# All required packages are installed automatically on first run.
# Each analysis writes its own files to results/ and figures/.
# Failures in any single analysis are reported but do not halt the rest.
# ==============================================================================

# ── Step 1: Ensure all dependencies are installed ────────────────────────────
required_cran <- c("igraph", "uwot", "boot", "lme4", "ape",
                   "ggplot2", "ggridges", "ggrepel",
                   "dplyr", "tidyr", "stringr", "Seurat", "BiocManager")
required_bioc <- c("scRNAseq", "yeastCC",
                   "SingleCellExperiment", "SummarizedExperiment",
                   "slingshot")

missing_cran <- required_cran[!required_cran %in% installed.packages()[, "Package"]]
missing_bioc <- required_bioc[!required_bioc %in% installed.packages()[, "Package"]]

if (length(missing_cran) > 0 || length(missing_bioc) > 0) {
  cat("Missing dependencies detected. Running install_dependencies.R...\n\n")
  source("install_dependencies.R")
  cat("\nResuming run_all.R\n\n")
}

# ── Step 2: Create output directories ────────────────────────────────────────
suppressMessages({
  if (!dir.exists("results")) dir.create("results")
  if (!dir.exists("figures")) dir.create("figures")
})

start_time <- Sys.time()

# ── Step 3: Run all twelve analyses in scheduled order ───────────────────────
# Order matters: cheap and self-contained first, expensive last.
# Each script self-sources its dependencies and writes its own outputs.
analyses <- c(
  "analysis/01_synthetic_stress_test.R",          # ~2 min  (synthetic)
  "analysis/12_helix_synthetic_validation.R",     # ~1 min  (synthetic)
  "analysis/02_spellman_ablation.R",              # ~2 min  (Spellman)
  "analysis/03_nestorowa_ablation.R",             # ~5 min  (Nestorowa)
  "analysis/04_five_method_benchmark.R",          # ~10 min (Nestorowa, all 5 methods)
  "analysis/08_multi_dataset_benchmark.R",        # ~15 min (all 5 datasets)
  "analysis/09_marker_validation.R",              # ~5 min  (Nestorowa, mixed models)
  "analysis/06_curvature_analysis.R",             # ~5 min  (Nestorowa)
  "analysis/07_loop_circularity_diagnostic.R",    # ~10 min (4 datasets)
  "analysis/05_pairwise_bootstrap_cis.R",         # ~15 min (3 datasets, bootstrap)
  "analysis/10_parameter_sensitivity.R",          # ~30 min (25 UMAP configs)
  "analysis/11_atlas_runtime_benchmark.R"         # ~45 min (synthetic, up to N=50k)
)

results <- list()
for (script in analyses) {
  cat(sprintf("\n%s\n%s %s\n%s\n",
              strrep("=", 70), "Running:", script, strrep("=", 70)))
  t0 <- Sys.time()
  ok <- tryCatch({
    source(script, echo = FALSE)
    TRUE
  }, error = function(e) {
    cat(sprintf("\n!!! FAILED: %s\n  %s\n", script, conditionMessage(e)))
    FALSE
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  results[[script]] <- list(ok = ok, mins = round(dt, 1))
  cat(sprintf("\n  >>> %s in %.1f min\n", if (ok) "DONE" else "FAILED", dt))
}

# ── Step 4: Summary report ───────────────────────────────────────────────────
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n%s\nSUMMARY (total %.1f min)\n%s\n", strrep("=", 70),
            elapsed, strrep("=", 70)))
for (script in names(results)) {
  cat(sprintf("  [%s] %5.1f min  %s\n",
              if (results[[script]]$ok) "OK " else "ERR",
              results[[script]]$mins, script))
}
n_ok  <- sum(sapply(results, `[[`, "ok"))
n_err <- length(results) - n_ok
cat(sprintf("\n%d/%d analyses completed successfully.\n",
            n_ok, length(results)))
if (n_err > 0) cat(sprintf("%d failed; check messages above.\n", n_err))
cat("\nResults in results/, figures in figures/.\n")

# ── Step 5: Write session info for reproducibility ───────────────────────────
sink("results/session_info_run.txt")
cat("Session info from run_all.R execution at:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()
cat("\nSession info saved to results/session_info_run.txt\n")
