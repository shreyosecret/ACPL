# ==============================================================================
# run_all.R -- reproduce every result in engrXiv 10.31224/7087 version 2
#
# Run from the repository root:   Rscript run_all.R
#
# First run downloads four datasets through ExperimentHub, builds the embedding
# cache in cache/, and writes tables to results/. Budget 60-90 minutes end to
# end. Later runs reuse the cache and are much faster.
#
# This script fails loudly. If an input is missing or an analysis produces no
# rows, it stops rather than reporting success over an empty table. Version
# 2.0.0 of this repository did the latter, and every table it wrote on a clean
# clone was empty.
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

# ---- build the dataset cache -------------------------------------------------
# Every analysis script reads cache/cache_<tag>.rds and has no way to create it.
# Without this block a clean clone finds nothing, skips silently, and writes
# empty tables. cache_dataset() returns immediately when the file already exists.
LOADERS <- list(nestorowa = load_nestorowa, buettner = load_buettner,
                leng      = load_leng,      richard  = load_richard)

for (tag in names(LOADERS)) {
  f <- file.path("cache", paste0("cache_", tag, ".rds"))
  if (file.exists(f)) { message("[cache] ", tag, " present"); next }
  message("[build] ", tag, " -- downloading and embedding, this is the slow part")
  ok <- tryCatch({ cache_dataset(LOADERS[[tag]], tag); TRUE },
                 error = function(e) { message("   FAILED: ", conditionMessage(e)); FALSE })
  if (!ok) message("   ", tag, " could not be built")
}

absent <- names(LOADERS)[!file.exists(
  file.path("cache", paste0("cache_", names(LOADERS), ".rds")))]
if (length(absent))
  stop("Dataset cache could not be built for: ", paste(absent, collapse = ", "),
       "\nThe analyses below read that cache and cannot run without it.",
       "\nDo not treat any table in results/ as a result of this run.")

# ---- run the analyses --------------------------------------------------------
scripts <- sort(list.files("analysis", pattern = "^0[1-9]_.*\\.R$",
                           full.names = TRUE))
for (s in scripts) {
  message("\n========== ", s, " ==========")
  source(s, echo = FALSE)
}

# ---- refuse to report success over empty output ------------------------------
csvs  <- list.files("results", pattern = "\\.csv$", full.names = TRUE)
empty <- csvs[file.size(csvs) < 32]
if (length(empty))
  stop("These tables are empty, so the run did not actually produce them:\n  ",
       paste(basename(empty), collapse = ", "),
       "\nTreat this run as failed.")

message("\nAll analyses complete. ", length(csvs), " tables in results/.")
writeLines(capture.output(sessionInfo()), "session_info.txt")
message("session_info.txt regenerated.")
