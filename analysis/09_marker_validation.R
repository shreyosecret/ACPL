# ==============================================================================
# Analysis 09: Biological anchor validation
# Section 4.9 of the paper. Tests whether ACPL arc length recovers biological
# progression as measured by canonical cell-cycle markers, with Moran's I
# correction for spatial autocorrelation.
# Outputs: results/marker_validation.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/statistics.R")
source("data-raw/load_datasets.R")

library(Seurat); library(lme4); library(dplyr)

seu <- load_nestorowa()
if (is.null(seu)) stop("Nestorowa load failed.")

coords <- Embeddings(seu, "umap")
phases <- as.character(seu$Phase)
arc    <- acpl_arc(coords)

# Canonical markers (use mouse gene symbols — capitalised first letter)
markers <- c("Top2a", "Mki67", "Pcna", "Mcm2")

# Pull expression
expr <- GetAssayData(seu, layer = "data")
present <- markers[markers %in% rownames(expr)]
if (length(present) == 0) stop("No markers found in dataset.")

results <- list()
for (m in present) {
  vals <- as.numeric(expr[m, ])
  # Moran's I along arc-length ordering
  mor <- moran_i_1d(vals, arc, n_neighbours = 5)

  # Naive ANOVA
  df <- data.frame(expr = vals, Phase = factor(phases))
  df_phased <- df[df$Phase %in% c("G1", "S", "G2M"), ]
  aov_p <- summary(aov(expr ~ Phase, data = df_phased))[[1]][1, "Pr(>F)"]

  # Mixed model with arc-length bin random effect
  bin <- cut(arc, breaks = 50, labels = FALSE)
  df_mm <- data.frame(expr = vals, Phase = factor(phases), bin = factor(bin))
  df_mm <- df_mm[df_mm$Phase %in% c("G1", "S", "G2M"), ]
  m1 <- tryCatch(lme4::lmer(expr ~ Phase + (1 | bin), data = df_mm,
                             REML = FALSE), error = function(e) NULL)
  m0 <- tryCatch(lme4::lmer(expr ~ 1 + (1 | bin), data = df_mm,
                             REML = FALSE), error = function(e) NULL)
  mm_p <- if (!is.null(m1) && !is.null(m0)) anova(m0, m1)$`Pr(>Chisq)`[2] else NA

  # Mean expression per phase
  means <- df_phased %>% dplyr::group_by(Phase) %>%
    dplyr::summarise(mean = mean(expr, na.rm = TRUE), .groups = "drop")

  results[[m]] <- data.frame(
    marker     = m,
    morans_I   = round(mor$I, 3),
    morans_p   = signif(mor$p, 3),
    anova_p    = signif(aov_p, 3),
    mixed_p    = signif(mm_p, 3),
    G1  = round(means$mean[means$Phase == "G1"],  3),
    S   = round(means$mean[means$Phase == "S"],   3),
    G2M = round(means$mean[means$Phase == "G2M"], 3)
  )
}
out <- do.call(rbind, results)
print(out, row.names = FALSE)
write.csv(out, "results/marker_validation.csv", row.names = FALSE)
cat("\nResults saved: results/marker_validation.csv\n")
