# ==============================================================================
# check_cache_reproducible.R
#
# Does a clean rebuild of the embedding cache reproduce the cache the published
# numbers were computed from?
#
# This matters because every table in results/ is a function of cache/*.rds, and
# those files were built once, months before the repository was assembled. If
# Seurat's PCA and UMAP are seeded end to end the rebuild is identical and the
# published point estimates are exactly reproducible. If they are not, the
# published point estimates are specific to one embedding, the conclusions rest
# on analysis/06_seed_replication.R instead, and the paper must say so.
#
# Do not guess which. Run this.
#
#   Rscript check_cache_reproducible.R
#
# Runtime is a full rebuild of four datasets, roughly 20-40 minutes. Nothing in
# cache/ is modified; the rebuild goes to cache_check/.
# ==============================================================================

stopifnot(file.exists("R/loaders.R"))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

REF <- "cache"; NEW <- "cache_check"
dir.create(NEW, showWarnings = FALSE)

LOADERS <- list(nestorowa = load_nestorowa, buettner = load_buettner,
                leng      = load_leng,      richard  = load_richard)

report <- list()
for (tag in names(LOADERS)) {
  rf <- file.path(REF, paste0("cache_", tag, ".rds"))
  if (!file.exists(rf)) { message("[skip] no reference cache for ", tag); next }
  message("\n[rebuild] ", tag)
  ok <- tryCatch({ cache_dataset(LOADERS[[tag]], tag, cache_dir = NEW); TRUE },
                 error = function(e) { message("  FAILED: ", conditionMessage(e)); FALSE })
  if (!ok) next

  a <- readRDS(rf); b <- readRDS(file.path(NEW, paste0("cache_", tag, ".rds")))

  same_n     <- length(a$phases) == length(b$phases)
  same_phase <- same_n && identical(as.character(a$phases), as.character(b$phases))
  # UMAP axes can come out sign-flipped or swapped without changing the geometry,
  # so compare the pairwise distance matrix, which is invariant to both.
  geom <- NA_real_
  if (same_n && nrow(a$coords) == nrow(b$coords)) {
    da <- as.numeric(dist(a$coords)); db <- as.numeric(dist(b$coords))
    geom <- suppressWarnings(cor(da, db, method = "spearman"))
  }
  ident <- isTRUE(all.equal(a$coords, b$coords, tolerance = 1e-8))

  report[[length(report) + 1]] <- data.frame(
    dataset = tag, n_ref = length(a$phases), n_new = length(b$phases),
    phases_identical = same_phase, coords_identical = ident,
    geometry_rho = round(geom, 4), row.names = NULL)
  message(sprintf("  n %d/%d | phases identical %s | coords identical %s | geometry rho %.4f",
                  length(a$phases), length(b$phases), same_phase, ident, geom))
}

out <- do.call(rbind, report)
if (is.null(out) || !nrow(out)) stop("Nothing was compared. Check that cache/ holds the reference files.")
print(out, row.names = FALSE)

cat("\n---------------------------------------------------------------\n")
if (all(out$coords_identical)) {
  cat("VERDICT: the rebuild is bit-identical. Published point estimates are\n")
  cat("exactly reproducible from a clean clone. No change to the papers.\n")
} else if (all(out$geometry_rho > 0.95, na.rm = TRUE)) {
  cat("VERDICT: the rebuild is NOT identical but the geometry is preserved.\n")
  cat("Point estimates in the tables will move. The conclusions rest on\n")
  cat("analysis/06_seed_replication.R, and the availability statement in both\n")
  cat("papers must say that exact values are embedding-specific.\n")
} else {
  cat("VERDICT: the rebuild differs materially. Stop and investigate before\n")
  cat("making any reproducibility claim.\n")
}
