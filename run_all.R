# ==============================================================================
# run_all.R -- reproduce every result in engrXiv 10.31224/7087 version 2
#
# Run from the repository root:   Rscript run_all.R
#
# First run downloads four datasets through ExperimentHub and caches embeddings
# to cache/. Budget 45-60 minutes end to end. Later runs are much faster.
# ==============================================================================

stopifnot(file.exists("R/metrics.R"))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

need_cran <- c("igraph", "boot", "Seurat", "stringr", "dplyr", "uwot")
need_bioc <- c("scRNAseq", "SummarizedExperiment", "SingleCellExperiment",
               "org.Mm.eg.db", "AnnotationDbi", "slingshot", "DelayedMatrixStats",
               "tricycle", "scuttle", "yeastCC")
miss <- c(need_cran, need_bioc)[
  !vapply(c(need_cran, need_bioc), requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing packages. See install_dependencies.R\n  ",
                       paste(miss, collapse = ", "))

# ExperimentHub prompts on first use and a prompt inside a sourced script eats
# the rest of the file. Create the caches up front.
for (h in c("AnnotationHub", "ExperimentHub", "BiocFileCache")) {
  d <- tools::R_user_dir(h, which = "cache")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

dir.create("cache",   showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

scripts <- sort(list.files("analysis", pattern = "^0[1-8]_.*\\.R$",
                           full.names = TRUE))
for (s in scripts) {
  message("\n========== ", s, " ==========")
  source(s, echo = FALSE)
}

message("\nAll analyses complete. Tables are in results/.")
writeLines(capture.output(sessionInfo()), "session_info.txt")
message("session_info.txt regenerated.")
