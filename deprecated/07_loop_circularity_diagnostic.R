# ==============================================================================
# Analysis 07: Loop circularity diagnostic across datasets
# Section 4.7 of the paper. Computes loop scores and tests whether they
# correlate with ACPL advantage. Note: counterintuitive result (Nestorowa
# scores high but ACPL wins) led to the angular-monotonicity reframing.
# Outputs: results/loop_scores.csv, figures/fig_loop_diagnostic.png
# ==============================================================================

source("R/acpl_engine.R")
source("R/plotting.R")
source("data-raw/load_datasets.R")

library(Seurat); library(ggplot2); library(dplyr); library(uwot)

datasets <- list()

# Nestorowa
seu_n <- load_nestorowa()
if (!is.null(seu_n)) {
  ls <- loop_score(Embeddings(seu_n, "umap"))
  datasets[["nestorowa"]] <- list(
    name = "Nestorowa HSC", N = ncol(seu_n),
    loop_score = ls$loop_score, rmsd_norm = ls$rmsd_norm,
    acpl_sa = 75.9, mst_sa = 74.3, diff_pp = 1.6)
}

# Buettner
seu_b <- load_buettner()
if (!is.null(seu_b)) {
  ls <- loop_score(Embeddings(seu_b, "umap"))
  datasets[["buettner"]] <- list(
    name = "Buettner mESC", N = ncol(seu_b),
    loop_score = ls$loop_score, rmsd_norm = ls$rmsd_norm,
    acpl_sa = 75.5, mst_sa = 77.0, diff_pp = -1.5)
}

# Leng
seu_l <- load_leng()
if (!is.null(seu_l)) {
  ls <- loop_score(Embeddings(seu_l, "umap"))
  datasets[["leng"]] <- list(
    name = "Leng mESC", N = ncol(seu_l),
    loop_score = ls$loop_score, rmsd_norm = ls$rmsd_norm,
    acpl_sa = 89.5, mst_sa = 96.2, diff_pp = -6.7)
}

# Spellman (build UMAP on the fly)
sp <- load_spellman()
if (!is.null(sp)) {
  pca_s <- prcomp(t(sp$counts), scale. = TRUE)$x[, 1:2]
  set.seed(42)
  umap_s <- uwot::umap(pca_s, n_neighbors = 5, min_dist = 0.5, verbose = FALSE)
  ls <- loop_score(umap_s)
  datasets[["spellman"]] <- list(
    name = "Spellman yeast", N = ncol(sp$counts),
    loop_score = ls$loop_score, rmsd_norm = ls$rmsd_norm,
    acpl_sa = 97.0, mst_sa = 52.8, diff_pp = 44.2)
}

# Compile
loop_df <- do.call(rbind, lapply(datasets, function(d)
  data.frame(Dataset = d$name, N = d$N,
             Loop_score = round(d$loop_score, 3),
             RMSD_norm  = round(d$rmsd_norm, 3),
             ACPL_SA = d$acpl_sa, MST_SA = d$mst_sa,
             Diff_pp = d$diff_pp)))
print(loop_df, row.names = FALSE)

# Correlation if 3+ datasets loaded
n_ok <- sum(!is.na(loop_df$Loop_score))
if (n_ok >= 3) {
  cor_loop <- cor.test(loop_df$Loop_score, loop_df$Diff_pp,
                        method = "spearman", exact = FALSE)
  cor_rmsd <- cor.test(loop_df$RMSD_norm, loop_df$Diff_pp,
                        method = "spearman", exact = FALSE)
  cat(sprintf("\nSpearman rho (Loop ~ Diff): %.3f, p=%.4f\n",
              cor_loop$estimate, cor_loop$p.value))
  cat(sprintf("Spearman rho (RMSD ~ Diff): %.3f, p=%.4f\n",
              cor_rmsd$estimate, cor_rmsd$p.value))
} else {
  cat(sprintf("\nOnly %d datasets loaded; correlation requires n>=3.\n", n_ok))
}

write.csv(loop_df, "results/loop_scores.csv", row.names = FALSE)
cat("Results saved: results/loop_scores.csv\n")
