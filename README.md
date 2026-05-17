# ACPL: Polar-Linearization for Trajectory Inference

R implementation of ACPL (Arc + LOESS Polar-Linearization), a coordinate-system prior for directed graph construction on cyclic biological manifolds.

**Repository:** https://github.com/shreyosecret/ACPL

This is the official code repository for:

> Ghosh, S. (2026). *ACPL: A Coordinate-System Prior for Topology-Aware Directed Trajectory Inference on Cyclic Biological Manifolds.* Submitted to *Bioinformatics Advances* (Oxford University Press).

**Submission status:** Under review at *Bioinformatics Advances*.

## Quick reproduction

```bash
git clone https://github.com/shreyosecret/ACPL.git
cd ACPL
Rscript run_all.R
```

This installs all required CRAN and Bioconductor packages on first run (via `install_dependencies.R`), then executes all twelve analysis scripts in scheduled order, writing CSV results to `results/` and PNG figures to `figures/`. Total runtime: approximately 90 minutes on a modern laptop.

## What ACPL Does

Standard trajectory inference methods (Monocle, Slingshot, PAGA) construct an undirected graph from Euclidean proximity, then attempt to assign direction afterward. On cyclic processes like the cell cycle, this produces backward edges where late G2M cells are connected to early G1 cells across the loop closure — a structural failure we call the **Euclidean shortcut problem**.

ACPL takes a different approach: it transforms the 2D UMAP embedding into polar coordinates aligned with the manifold's angular structure, computes LOESS-smoothed cumulative arc length as a topology-aware temporal proxy, and admits edges only in the direction of increasing arc length. The result is a directed acyclic graph by construction at O(N log N) cost, with no root cell required.

Two formal theorems hold regardless of dataset:

1. **Acyclicity by construction.** The ACPL graph contains no directed cycles.
2. **LOESS necessity.** Without LOESS smoothing, raw arc length produces 100% structural accuracy as a trivial mathematical artefact rather than a geometric result.

## Installation

```r
# Core dependencies
install.packages(c("igraph", "uwot", "boot", "lme4", "ape", "ggplot2",
                   "ggridges", "ggrepel", "dplyr", "tidyr", "stringr"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("scRNAseq", "yeastCC", "SingleCellExperiment",
                       "SummarizedExperiment", "slingshot"))

# Seurat
install.packages("Seurat")
```

For PAGA+DPT comparisons, you also need scanpy:

```r
install.packages("reticulate")
reticulate::py_install(c("scanpy", "anndata", "leidenalg"), pip = TRUE)
```

## Quick Start

```r
source("R/acpl_engine.R")

# Synthetic cycle
set.seed(1)
theta  <- seq(0, 2 * pi, length.out = 100)
coords <- cbind(cos(theta) + rnorm(100, 0, 0.05),
                sin(theta) + rnorm(100, 0, 0.05))

# Compute ACPL arc length
arc <- acpl_arc(coords, span = 0.25)

# Build directed acyclic graph
graph <- acpl_graph(coords, arc, K = 2)

# Verify acyclicity (Theorem 1)
stopifnot(acpl_is_acyclic(graph))
```

A more substantial worked example is in `examples/minimal_example.R`.

## Repository Structure

```
acpl-repo/
├── R/                          # Core ACPL package functions
│   ├── acpl_engine.R           #   Main algorithm + theorems
│   ├── comparison_methods.R    #   MST, Slingshot, PAGA wrappers
│   ├── statistics.R            #   Bootstrap CIs, Moran's I
│   └── plotting.R              #   ACPL colour palette, ridge/forest plots
├── data-raw/                   # Dataset loaders
│   └── load_datasets.R         #   All five biological datasets
├── analysis/                   # Paper analysis scripts (numbered by section)
│   ├── 01_synthetic_stress_test.R
│   ├── 02_spellman_ablation.R
│   ├── 03_nestorowa_ablation.R
│   ├── 04_five_method_benchmark.R
│   ├── 05_pairwise_bootstrap_cis.R
│   ├── 06_curvature_analysis.R
│   ├── 07_loop_circularity_diagnostic.R
│   ├── 08_multi_dataset_benchmark.R
│   ├── 09_marker_validation.R
│   ├── 10_parameter_sensitivity.R
│   ├── 11_atlas_runtime_benchmark.R
│   └── 12_helix_synthetic_validation.R
├── examples/
│   └── minimal_example.R       # 50-line reproducible example
├── tests/testthat/
│   └── test_acpl_engine.R      # Unit tests for engine functions
├── figures/                    # Output figures land here
├── results/                    # Output CSVs land here
├── run_all.R                   # Master script: runs every analysis
├── DESCRIPTION                 # R package metadata
├── LICENSE                     # MIT
├── CITATION                    # How to cite the paper
└── .gitignore
```

## Reproducing the Paper

```bash
# Run from the repository root
Rscript run_all.R
```

This runs all twelve analyses end-to-end and writes results to `results/` and figures to `figures/`. Total runtime: approximately 90 minutes on a modern laptop, dominated by the atlas runtime benchmark (45 min) and the UMAP sensitivity analysis (30 min).

To run a single analysis:

```bash
Rscript analysis/05_pairwise_bootstrap_cis.R
```

## Datasets

| Dataset         | Source                          | N     | Phase labels         |
|-----------------|----------------------------------|-------|----------------------|
| Spellman yeast  | `yeastCC` (Bioconductor)         | 18    | Microarray timepoints|
| Nestorowa HSC   | `scRNAseq::NestorowaHSCData()`   | 1,920 | Seurat scoring       |
| Buettner mESC   | `scRNAseq::BuettnerESCData()`    | 288   | FACS-validated       |
| Leng mESC       | `scRNAseq::LengESCData()`        | 247   | FACS-sorted          |
| Richard T Cells | `scRNAseq::RichardTCellData()`   | 572   | Seurat scoring       |
| HeLa            | GEO `GSE64016`                   | 349   | Asynchronous         |

All datasets except HeLa are downloaded automatically on first use through the loaders in `data-raw/load_datasets.R`. HeLa is referenced in the paper as the boundary condition (no method scores above the 50% random baseline) and is not required for the main reproduction.

## Headline Results

| Dataset         | ACPL  | MST   | Slingshot | Diff (ACPL-MST) | Bootstrap 95% BCa CI |
|-----------------|-------|-------|-----------|-----------------|----------------------|
| Spellman yeast  | 97.0% | 52.8% | n/a       | +44.2 pp        | not computed (N=18)  |
| Nestorowa HSC   | 75.9% | 74.3% | 77.3%     | +1.6 pp         | [-1.31, +4.24]       |
| Buettner mESC   | 75.5% | 77.0% | n/a       | -1.5 pp         | [-8.63, +5.04]       |
| Leng mESC       | 89.5% | 96.2% | n/a       | -6.7 pp         | not computed         |
| Richard T Cells | 67.8% | 69.2% | n/a       | -1.4 pp         | not computed         |

No pairwise bootstrap comparison reaches conventional significance. The Nestorowa advantage is directional evidence; ACPL's claims are formal (acyclicity by construction, LOESS necessity, parameter invariance) rather than aggregate accuracy claims.

## Citation

See `CITATION` for the BibTeX entry. The short form is:

```
Ghosh, S. (2026). Polar-Linearization: A Coordinate-System Prior for
Topology-Aware Directed Trajectory Inference on Cyclic Biological Manifolds.
Preprint.
```

## License

MIT License. See `LICENSE` for details.

## Contact

Shreyo Ghosh — `ghoshshreyo2007@gmail.com`

## Acknowledgements

This work benefited from extensive peer review feedback, particularly on the curvature quantification (which produced a null result and led to withdrawal of the spatial-localisation claim), the loop circularity diagnostic (which produced a counterintuitive result and led to the angular monotonicity reframing), and the atlas-scale runtime benchmark (which revised an earlier claim about ACPL's scaling advantage).
