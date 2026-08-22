# ==============================================================================
# R/loaders.R -- dataset loaders
#
# REPLACES data-raw/load_datasets.R from the v1 release, which no longer runs.
# The scRNAseq package now returns Ensembl gene identifiers, so four of the five
# v1 loaders produced no phase labels or failed outright. See CHANGELOG.md.
# ==============================================================================

pick_assay <- function(sce) {
  an <- SummarizedExperiment::assayNames(sce)
  for (a in c("counts", "normcounts", "logcounts", "tpm")) if (a %in% an) return(a)
  an[1]
}

#' Map Ensembl gene IDs to symbols, dropping unmapped and duplicated symbols
ens_to_symbol <- function(sce, db_name = "org.Mm.eg.db") {
  if (!requireNamespace(db_name, quietly = TRUE))
    stop("Need annotation DB. BiocManager::install('", db_name, "')")
  db  <- getExportedValue(db_name, db_name)
  ids <- sub("\\..*$", "", rownames(sce))
  sym <- suppressMessages(AnnotationDbi::mapIds(
    db, keys = ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  keep <- !is.na(sym) & !duplicated(sym)
  message(sprintf("  mapped %d/%d Ensembl IDs to unique symbols", sum(keep), length(ids)))
  sce <- sce[keep, ]
  rownames(sce) <- as.character(sym[keep])
  sce
}

seurat_pipeline <- function(sce, dims, nn, md) {
  seu <- Seurat::CreateSeuratObject(
    counts = SummarizedExperiment::assay(sce, pick_assay(sce)))
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, verbose = FALSE)
  Seurat::RunUMAP(seu, dims = dims, n.neighbors = nn, min.dist = md, verbose = FALSE)
}

score_phases <- function(seu, species = c("mouse", "human")) {
  species <- match.arg(species)
  sg <- Seurat::cc.genes.updated.2019$s.genes
  gg <- Seurat::cc.genes.updated.2019$g2m.genes
  if (species == "mouse") {
    sg <- stringr::str_to_title(sg); gg <- stringr::str_to_title(gg)
  }
  si <- intersect(sg, rownames(seu)); gi <- intersect(gg, rownames(seu))
  message(sprintf("  cell-cycle genes matched: S %d/%d, G2M %d/%d",
                  length(si), length(sg), length(gi), length(gg)))
  if (length(si) < 5 || length(gi) < 5)
    stop("Too few cell-cycle genes matched; gene identifier format is wrong.")
  seu <- Seurat::CellCycleScoring(seu, s.features = si, g2m.features = gi,
                                  set.ident = FALSE)
  seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))
  seu
}

# ---- the four single-cell datasets ------------------------------------------
# Species and cell line corrected relative to v1: Leng is HUMAN H1 ESC (v1 said
# mouse); Richard is MOUSE (v1 said human).

load_nestorowa <- function() {                       # mouse HSC, N = 1920
  sce <- ens_to_symbol(scRNAseq::NestorowaHSCData(), "org.Mm.eg.db")
  score_phases(seurat_pipeline(sce, 1:10, 30, 0.3), "mouse")
}

load_buettner <- function() {                        # mouse ESC, N = 288, FACS
  sce  <- scRNAseq::BuettnerESCData()
  keep <- !is.na(SummarizedExperiment::colData(sce)$phase)
  seu  <- seurat_pipeline(sce[, keep], 1:10, 15, 0.4)
  seu$Phase <- factor(
    as.character(SummarizedExperiment::colData(sce)$phase[keep]),
    levels = c("G1", "S", "G2M"))
  seu
}

load_leng <- function() {                            # HUMAN H1 ESC, N = 247, FACS
  sce <- scRNAseq::LengESCData()                     # assay is normcounts, col is Phase
  raw <- as.character(SummarizedExperiment::colData(sce)$Phase)
  phs <- ifelse(raw == "G1", "G1",
         ifelse(raw == "S",  "S",
         ifelse(raw %in% c("G2", "G2M"), "G2M", NA_character_)))
  keep <- !is.na(phs)
  seu  <- seurat_pipeline(sce[, keep], 1:15, 20, 0.4)
  seu$Phase <- factor(phs[keep], levels = c("G1", "S", "G2M"))
  seu
}

load_richard <- function() {                         # MOUSE T cells, N = 572
  sce <- ens_to_symbol(scRNAseq::RichardTCellData(), "org.Mm.eg.db")
  score_phases(seurat_pipeline(sce, 1:15, 30, 0.3), "mouse")
}

# ---- yeast ------------------------------------------------------------------
# v1 called data(spellman); the object in current yeastCC is spYCCES.
load_spellman_alpha <- function() {
  avail <- as.character(utils::data(package = "yeastCC")$results[, "Item"])
  eset <- NULL
  for (nm in c("spYCCES", avail)) {
    ok <- tryCatch({ utils::data(list = nm, package = "yeastCC",
                                 envir = environment()); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) next
    obj <- get(nm, envir = environment())
    if (methods::is(obj, "ExpressionSet")) { eset <- obj; break }
  }
  if (is.null(eset)) stop("No ExpressionSet found in yeastCC.")
  X <- Biobase::exprs(eset)
  cols <- grep("^alpha", colnames(X), value = TRUE)
  X <- X[, cols, drop = FALSE]
  tmin <- as.numeric(gsub("[^0-9]", "", colnames(X)))
  o <- order(tmin); X <- X[, o, drop = FALSE]; tmin <- tmin[o]
  X <- X[rowMeans(is.na(X)) < 0.2, , drop = FALSE]
  for (i in seq_len(nrow(X))) {
    na <- is.na(X[i, ]); if (any(na)) X[i, na] <- mean(X[i, ], na.rm = TRUE)
  }
  list(expr = X, time = tmin, period = 66)
}

#' Cache a dataset's embedding, PCA and phase labels to disk
cache_dataset <- function(loader, tag, cache_dir = "cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir, paste0("cache_", tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  seu <- loader()
  obj <- list(coords = Seurat::Embeddings(seu, "umap"),
              pca    = Seurat::Embeddings(seu, "pca"),
              phases = as.character(seu$Phase))
  saveRDS(obj, f)
  obj
}
