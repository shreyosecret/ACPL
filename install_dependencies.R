# Run once before run_all.R
install.packages(c("igraph", "boot", "Seurat", "stringr", "dplyr", "uwot"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "scRNAseq", "SummarizedExperiment", "SingleCellExperiment",
  "org.Mm.eg.db", "AnnotationDbi",
  "slingshot", "DelayedMatrixStats",   # DelayedMatrixStats is a slingshot dep
  "tricycle", "scuttle",               # reference-based baseline
  "yeastCC"                            # Spellman yeast series (object: spYCCES)
))
