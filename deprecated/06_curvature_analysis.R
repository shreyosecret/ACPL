# ==============================================================================
# Analysis 06: Local manifold curvature vs ACPL/MST advantage
# Section 4.6 of the paper. Tests whether ACPL corrections concentrate at
# high-curvature regions. NULL RESULT: this analysis returned p=0.32 and
# rho=-0.005, which led to withdrawal of the spatial-localisation claim
# from earlier drafts.
# Outputs: results/curvature_summary.csv
# ==============================================================================

source("R/acpl_engine.R")
source("R/comparison_methods.R")
source("R/statistics.R")
source("data-raw/load_datasets.R")

library(Seurat); library(dplyr)

seu <- load_nestorowa()
if (is.null(seu)) stop("Nestorowa load failed.")
coords <- Embeddings(seu, "umap")
phases <- as.character(seu$Phase)

arc    <- acpl_arc(coords)
ord    <- order(arc)
co     <- coords[ord, ]
window <- 8

# Local curvature at each window position
n <- nrow(co)
curv <- rep(NA_real_, n)
for (i in (window + 1):(n - window)) {
  curv[i] <- three_point_curvature(co[i - window, ], co[i, ], co[i + window, ])
}

# ACPL correctness vs MST correctness at each window
wa <- window_correctness(acpl_arc(coords), phases, 10)
wm <- window_correctness(mst_pseudotime(coords, "G1", phases), phases,
                          10, mirror = TRUE)
n_use <- min(length(wa), length(wm), n - 10)

cat_vec <- ifelse(wa[1:n_use] == 1 & wm[1:n_use] == 0, "ACPL wins",
           ifelse(wa[1:n_use] == 0 & wm[1:n_use] == 1, "MST wins",
           ifelse(wa[1:n_use] == 1 & wm[1:n_use] == 1, "Both correct",
                  "Both wrong")))

mid_idx <- pmin(round((seq_len(n_use) + seq_len(n_use) + 10) / 2), n)
df <- data.frame(category = cat_vec, curvature = curv[mid_idx])
df <- df[!is.na(df$curvature) & is.finite(df$curvature), ]

summary_tbl <- df %>%
  group_by(category) %>%
  summarise(n = n(),
            mean_curv = round(mean(curvature, na.rm = TRUE), 3),
            median_curv = round(median(curvature, na.rm = TRUE), 3),
            sd = round(sd(curvature, na.rm = TRUE), 3),
            .groups = "drop")
print(summary_tbl)

# Wilcoxon test: ACPL-win curvature vs MST-win curvature
acpl_win <- df$curvature[df$category == "ACPL wins"]
mst_win  <- df$curvature[df$category == "MST wins"]
wt <- wilcox.test(acpl_win, mst_win, alternative = "two.sided", exact = FALSE)
sp <- cor.test(df$curvature,
               as.integer(df$category == "ACPL wins"),
               method = "spearman", exact = FALSE)

cat(sprintf("\nWilcoxon W=%.0f, p=%.4f\n", wt$statistic, wt$p.value))
cat(sprintf("Spearman rho=%.4f, p=%.4f\n", sp$estimate, sp$p.value))
cat("\nNULL RESULT: spatial-localisation claim WITHDRAWN.\n")

write.csv(summary_tbl, "results/curvature_summary.csv", row.names = FALSE)
cat("Results saved: results/curvature_summary.csv\n")
