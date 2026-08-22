# ==============================================================================
# Comparison methods: MST, Slingshot wrapper, PAGA wrapper
# Used in benchmarking against ACPL.
# ==============================================================================

#' Monocle-style MST pseudotime
#'
#' Builds a minimum spanning tree over Euclidean distances and assigns
#' pseudotime as the shortest-path distance from a user-specified root cell.
#'
#' @param coords Numeric coordinate matrix
#' @param root_phase Phase label to use as root (e.g. "G1")
#' @param phases Vector of phase labels matching rows of coords
#' @return Numeric vector of pseudotime values
#' @export
mst_pseudotime <- function(coords, root_phase, phases) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Install igraph: install.packages('igraph')")
  dist_mat <- as.matrix(dist(coords))
  g <- igraph::graph_from_adjacency_matrix(dist_mat, weighted = TRUE,
                                            mode = "undirected")
  mst_tree <- igraph::mst(g)
  root_idx <- which(phases == root_phase)[1]
  if (is.na(root_idx)) stop("No cell with root_phase found in phases")
  as.numeric(igraph::distances(mst_tree, v = root_idx,
                                weights = igraph::E(mst_tree)$weight))
}

#' Slingshot pseudotime
#'
#' Wrapper around the slingshot package, with a reasonable default of taking
#' the mean pseudotime across all fitted lineages.
#'
#' @param coords 2D coordinate matrix (UMAP recommended)
#' @param phases Phase labels for clustering
#' @param start_phase Starting phase cluster (default "G1")
#' @return Numeric vector of mean pseudotime across lineages
#' @export
slingshot_pseudotime <- function(coords, phases, start_phase = "G1") {
  if (!requireNamespace("slingshot", quietly = TRUE))
    stop("Install slingshot: BiocManager::install('slingshot')")
  keep <- as.character(phases) %in% c("G1", "S", "G2M")
  sds <- slingshot::slingshot(
    data = coords[keep, , drop = FALSE],
    clusterLabels = as.character(phases)[keep],
    start.clus = start_phase,
    reducedDim = "UMAP"
  )
  pt_full <- rep(NA_real_, nrow(coords))
  pt_full[keep] <- rowMeans(slingshot::slingPseudotime(sds), na.rm = TRUE)
  pt_full
}

#' PAGA + diffusion pseudotime via reticulate
#'
#' Bridges to scanpy's PAGA implementation. Requires reticulate and scanpy
#' installed in the active Python environment.
#'
#' @param coords_pca PCA coordinate matrix
#' @param phases Phase labels for root cluster identification
#' @param resolution Leiden clustering resolution (default 0.5)
#' @return Numeric vector of DPT pseudotime
#' @export
paga_dpt_pseudotime <- function(coords_pca, phases, resolution = 0.5) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Install reticulate: install.packages('reticulate')")
  sc  <- reticulate::import("scanpy")
  ad  <- reticulate::import("anndata")
  spy <- reticulate::import("scipy.sparse")

  phase_chr <- as.character(phases)
  obs_df    <- data.frame(phase = phase_chr,
                          row.names = paste0("cell", seq_len(nrow(coords_pca))))
  X_sparse  <- spy$csr_matrix(coords_pca)
  adata     <- ad$AnnData(X = X_sparse, obs = obs_df)
  adata$obsm[["X_pca"]] <- coords_pca

  sc$pp$neighbors(adata, n_neighbors = 40L, n_pcs = ncol(coords_pca),
                  use_rep = "X_pca")
  sc$tl$leiden(adata, resolution = resolution)
  sc$tl$paga(adata, groups = "leiden")

  leiden_lab <- as.character(adata$obs[["leiden"]])
  if ("G1" %in% phase_chr) {
    cross_tab <- table(leiden_lab, phase_chr)
    g1_frac   <- if ("G1" %in% colnames(cross_tab))
                   cross_tab[, "G1"] / rowSums(cross_tab)
                 else rep(0, nrow(cross_tab))
    root_clus <- names(which.max(g1_frac))
  } else {
    root_clus <- names(table(leiden_lab))[1]
  }

  sc$tl$diffmap(adata)
  adata$uns[["iroot"]] <- as.integer(which(leiden_lab == root_clus)[1] - 1L)
  sc$tl$dpt(adata)
  as.numeric(adata$obs[["dpt_pseudotime"]])
}
