# ACPL — a root-free polar ordering for cyclic single-cell trajectories

Code accompanying **engrXiv [10.31224/7087](https://doi.org/10.31224/7087),
version 2**.

> **This release supersedes v1.** The v1 code does not run against current
> Bioconductor and its benchmark scripts contain a scoring asymmetry and a wrong
> chance level. Every error is documented in [CHANGELOG.md](CHANGELOG.md). The
> old scripts are kept in [`deprecated/`](deprecated/) for the record and should
> not be used.

## What ACPL is

Graph constructions that minimise Euclidean cost route across the interior of a
closed loop rather than around it, placing a backward edge at the G2M-to-G1
closure. ACPL (Arc plus LOESS Polar-Linearisation) prevents this by ordering
cells along the angular coordinate of a planar embedding *before* any graph is
built, then admitting edges only in the direction of increasing smoothed arc
length. The result is acyclic by construction, needs no starting cell, and is
built in O(N log N).

Angular ordering of cell-cycle data is **established practice** and is not a
contribution of this work — see tricycle, Cyclum, reCAT, CCPE and peco. The
contribution is the directed-graph construction and its evaluation.

## What the evaluation shows

| | Result |
|---|---|
| vs MST and Slingshot, native inputs | ACPL higher on all four datasets, 40/40 embedding seeds; 7 of 16 intervals exclude zero, none favouring a baseline |
| vs tricycle | tricycle higher on two of three datasets, significantly |
| Yeast (tricycle cannot run) | ACPL recovers no significant signal, n = 10, p = 0.81 |
| Handedness | ACPL cannot determine it; on two of four datasets the cycle runs backwards |

Root-freedom is a trade, not an advantage: independence from a starting cluster
costs the ability to orient the result, so an external anchor is required.

## Running it

```r
source("install_dependencies.R")   # once
```
```bash
Rscript run_all.R                  # from the repository root
```

First run downloads four datasets through ExperimentHub and builds the embedding
cache in `cache/`. Budget 60–90 minutes. Tables land in `results/`. The run
aborts rather than reporting success over an empty table.

## Layout

```
R/loaders.R      corrected dataset loaders (Ensembl to symbol, right assays)
R/metrics.R      chance levels, tau_cyc, circular helpers
R/methods.R      ACPL, MST, Slingshot, tricycle wrappers
analysis/01..08  the analyses, in run order
analysis/09      Holm and Bonferroni adjustment over the sixteen comparisons
check_cache_reproducible.R  does a clean rebuild match the published cache?
deprecated/      superseded v1 scripts, kept for the record
CHANGELOG.md     every correction, in full
```

Run `run_all.R` and nothing else. It builds the cache first, and every script
downstream reads that cache and cannot create it. Running a single analysis
script against an empty `cache/` used to produce an empty table and no error;
it now stops.

## Why tau_cyc differs between scripts

`tau_cyc` is the maximum Kendall tau against ordinal phase rank taken over a
grid of origin rotations and both traversal directions. Its value therefore
depends on how fine that grid is, and the three scripts that compute it use
different grids for different reasons:

| script | `NROT` | why |
|---|---|---|
| `02_cyclic_ordering.R` | 36 | no bootstrap, so a fine grid is cheap |
| `03_tau_pairwise_ci.R` | 12 | 500 paired resamples on a finer grid is not tractable |
| `05_native_tau_ci.R` | 24 | **the values quoted in the manuscripts** |

So ACPL on Nestorowa reads 0.4401 in script 02, 0.4314 in script 03 and 0.4499
in script 05. These are the same ordering measured on three grids, not three
different results. Quote script 05.

Within a script the grid is identical for every method, so it shifts both arms
of a comparison equally. Differences and confidence intervals are unaffected by
the choice; only the per-method point estimates move.

## Two things worth knowing if you build on this

**The chance level of tie-counting windowed concordance is not 50 percent.** It
is `(1 + sum p^2)/2`, bounded below by 2/3 for three phase levels, and 66.7 to
69.3 percent on these datasets. `chance_level()` in `R/metrics.R` computes it.

**Raw cumulative arc length reports perfect accuracy as an arithmetic identity.**
It is strictly increasing in traversal index by positive definiteness of the
norm, so a monotonicity filter on it admits every edge. Any implementation of
this construction reporting 100 percent is measuring a tautology, not
performance. This is why the LOESS step exists.

## Citation

Ghosh, S. *ACPL: a root-free polar ordering for cyclic single-cell
trajectories.* engrXiv 10.31224/7087, version 2.

## Licence

MIT, see [LICENSE](LICENSE).
