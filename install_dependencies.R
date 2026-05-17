# ==============================================================================
# install_dependencies.R
#
# Run this once before run_all.R. Installs every CRAN and Bioconductor package
# needed by any analysis script in the repository. Safe to re-run; only
# packages not already installed will be downloaded.
#
# Usage:
#   Rscript install_dependencies.R
#
# Expected install time on a fresh R: 15-25 minutes depending on connection.
# ==============================================================================

cat("\n=== ACPL repository: installing dependencies ===\n\n")

# Core CRAN packages
cran_pkgs <- c(
  "igraph", "uwot", "boot", "lme4", "ape",
  "ggplot2", "ggridges", "ggrepel",
  "dplyr", "tidyr", "stringr",
  "Seurat", "testthat",
  "BiocManager", "reticulate"
)

missing_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) {
  cat(sprintf("Installing %d CRAN packages: %s\n",
              length(missing_cran), paste(missing_cran, collapse = ", ")))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages already installed.\n")
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "scRNAseq", "yeastCC",
  "SingleCellExperiment", "SummarizedExperiment",
  "slingshot"
)

missing_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_bioc) > 0) {
  cat(sprintf("\nInstalling %d Bioconductor packages: %s\n",
              length(missing_bioc), paste(missing_bioc, collapse = ", ")))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

# Optional Python dependencies for PAGA + DPT comparison via reticulate
cat("\n=== Optional: Python dependencies for PAGA+DPT ===\n")
cat("To enable PAGA+DPT in analysis/04_five_method_benchmark.R:\n")
cat("  reticulate::py_install(c('scanpy','anndata','leidenalg'), pip=TRUE)\n")
cat("\nWithout this, the PAGA+DPT comparison will be skipped gracefully.\n")

cat("\n=== Dependency installation complete ===\n")
cat("You can now run: Rscript run_all.R\n\n")
