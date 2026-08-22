# ==============================================================================
# Dataset loaders for the five biological benchmarks
# Each function returns a list with: counts (or coords), phases, name.
# Defensive: each loader has a tryCatch with informative error messages.
# ==============================================================================

#' Load Spellman yeast cell cycle dataset
#'
#' Microarray data from synchronised S. cerevisiae cells across 18 timepoints.
#' Returns expression matrix and timepoint labels (which serve as phase proxies).
#'
#' @return List(counts, timepoints, name) or NULL if loading fails
#' @export
load_spellman <- function() {
  tryCatch({
    if (!requireNamespace("yeastCC", quietly = TRUE))
      BiocManager::install("yeastCC", update = FALSE, ask = FALSE)
    library(yeastCC)
    data(spellman)
    sm <- as.matrix(spellman)
    sm <- sm[, colMeans(is.na(sm)) < 0.5]
    for (j in seq_len(ncol(sm))) {
      nas <- is.na(sm[, j])
      if (any(nas)) sm[nas, j] <- mean(sm[, j], na.rm = TRUE)
    }
    list(counts = sm, timepoints = colnames(sm), name = "Spellman yeast")
  }, error = function(e) {
    message("load_spellman failed: ", conditionMessage(e))
    NULL
  })
}

#' Load Nestorowa HSC dataset
#'
#' Mouse haematopoietic stem cell SMART-seq2 data (1,920 cells).
#' Phase labels are assigned by Seurat::CellCycleScoring on murine
#' homologues of canonical human markers.
#'
#' @return Seurat object with Phase metadata, or NULL on failure
#' @export
load_nestorowa <- function() {
  tryCatch({
    if (!requireNamespace("scRNAseq", quietly = TRUE))
      stop("Install scRNAseq: BiocManager::install('scRNAseq')")
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Install Seurat: install.packages('Seurat')")
    library(scRNAseq); library(Seurat)
    nest <- NestorowaHSCData()
    seu  <- CreateSeuratObject(counts = SummarizedExperiment::assay(nest, "counts"))
    seu  <- NormalizeData(seu, verbose = FALSE)
    seu  <- FindVariableFeatures(seu, verbose = FALSE)
    seu  <- ScaleData(seu, verbose = FALSE)
    seu  <- RunPCA(seu, verbose = FALSE)
    seu  <- RunUMAP(seu, dims = 1:10, n.neighbors = 30,
                     min.dist = 0.3, verbose = FALSE)
    s_genes   <- stringr::str_to_title(cc.genes.updated.2019$s.genes)
    g2m_genes <- stringr::str_to_title(cc.genes.updated.2019$g2m.genes)
    seu <- CellCycleScoring(
      seu,
      s.features   = intersect(s_genes,   rownames(seu)),
      g2m.features = intersect(g2m_genes, rownames(seu)),
      set.ident    = FALSE)
    seu
  }, error = function(e) {
    message("load_nestorowa failed: ", conditionMessage(e))
    NULL
  })
}

#' Load Buettner mESC dataset
#'
#' Mouse embryonic stem cells with FACS-validated phase labels (gold standard).
#'
#' @return Seurat object with Phase metadata, or NULL on failure
#' @export
load_buettner <- function() {
  tryCatch({
    library(scRNAseq); library(Seurat); library(SingleCellExperiment)
    sce <- BuettnerESCData()
    keep <- !is.na(SummarizedExperiment::colData(sce)$phase)
    seu  <- CreateSeuratObject(counts = SummarizedExperiment::assay(sce[, keep], "counts"))
    seu  <- NormalizeData(seu, verbose = FALSE) |>
            FindVariableFeatures(verbose = FALSE) |>
            ScaleData(verbose = FALSE) |>
            RunPCA(verbose = FALSE) |>
            RunUMAP(dims = 1:10, n.neighbors = 15, min.dist = 0.4, verbose = FALSE)
    seu$Phase <- factor(
      as.character(SummarizedExperiment::colData(sce)$phase[keep]),
      levels = c("G1", "S", "G2M"))
    seu
  }, error = function(e) {
    message("load_buettner failed: ", conditionMessage(e))
    NULL
  })
}

#' Load Leng mESC dataset
#'
#' FACS-sorted mouse embryonic stem cells. Note: phase labels use 'G2'
#' rather than 'G2M'; we recode for consistency.
#'
#' @return Seurat object with Phase metadata, or NULL on failure
#' @export
load_leng <- function() {
  tryCatch({
    library(scRNAseq); library(Seurat); library(SingleCellExperiment)
    sce  <- LengESCData()
    raw  <- as.character(SummarizedExperiment::colData(sce)$phase)
    phs  <- dplyr::recode(raw, "G1" = "G1", "S" = "S", "G2" = "G2M",
                          .default = NA_character_)
    keep <- !is.na(phs)
    seu  <- CreateSeuratObject(counts = SummarizedExperiment::assay(sce[, keep], "counts"))
    seu  <- NormalizeData(seu, verbose = FALSE) |>
            FindVariableFeatures(verbose = FALSE) |>
            ScaleData(verbose = FALSE) |>
            RunPCA(verbose = FALSE) |>
            RunUMAP(dims = 1:15, n.neighbors = 20, min.dist = 0.4, verbose = FALSE)
    seu$Phase <- factor(phs[keep], levels = c("G1", "S", "G2M"))
    seu
  }, error = function(e) {
    message("load_leng failed: ", conditionMessage(e))
    NULL
  })
}

#' Load Richard cytotoxic T-cell dataset
#'
#' Human T cells; phases assigned by Seurat::CellCycleScoring.
#'
#' @return Seurat object with Phase metadata, or NULL on failure
#' @export
load_richard <- function() {
  tryCatch({
    library(scRNAseq); library(Seurat); library(SingleCellExperiment)
    sce <- RichardTCellData()
    seu <- CreateSeuratObject(counts = SummarizedExperiment::assay(sce, "counts"))
    seu <- NormalizeData(seu, verbose = FALSE) |>
           FindVariableFeatures(verbose = FALSE) |>
           ScaleData(verbose = FALSE) |>
           RunPCA(verbose = FALSE) |>
           RunUMAP(dims = 1:15, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
    seu <- CellCycleScoring(
      seu,
      s.features   = intersect(cc.genes.updated.2019$s.genes,   rownames(seu)),
      g2m.features = intersect(cc.genes.updated.2019$g2m.genes, rownames(seu)),
      set.ident    = FALSE)
    seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))
    seu
  }, error = function(e) {
    message("load_richard failed: ", conditionMessage(e))
    NULL
  })
}
