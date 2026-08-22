# Changelog

## v2 — corrections to the v1 release

The v1 code and the results it produced contained errors. They are listed here
in full. The corrected manuscript is engrXiv 10.31224/7087 version 2.

### Errors in the evaluation

**Wrong chance level.** v1 stated that sliding-window structural accuracy has "a
random baseline of approximately 50%". The statistic counts ties as successes,
so under a random ordering of the same label multiset its expectation is
`(1 + sum_k p_k^2) / 2`, bounded below by 2/3 for three phase levels. The true
values are 66.67 (Buettner), 66.77 (Leng), 67.59 (Richard) and 69.31
(Nestorowa) percent. Every absolute accuracy figure in v1 was compared against a
baseline wrong by 16.7 to 19.3 percentage points. Differences between methods
are unaffected; absolute claims are not.

**Asymmetric orientation handling.** In `04_five_method_benchmark.R` and
`08_multi_dataset_benchmark.R`, baselines were scored with `mirror = TRUE` and
ACPL with the default `mirror = FALSE`. `02_spellman_ablation.R` applied
`max(f, 1-f)` to MST only. The asymmetry ran against ACPL in every table it
touched. Separately, `mirror = TRUE` selects orientation using the ground-truth
labels and then scores against those same labels; neither setting is defensible.
The `mirror` argument has been removed in v2.

**Wrong dataset designations.**
- "Leng mESC" is *human* H1 embryonic stem cells, not mouse.
- "Richard T Cells" is *mouse*, not human; the v1 loader scored it with un-cased
  human gene symbols, which match nothing.
- Section 4.9 of v1 analysed data attributed to GEO GSE64016 and described as
  HeLa. That accession is the Leng et al. Oscope series (213 unsorted H1 hESC
  plus 247 H1-Fucci sorted cells, 460 total). It contains no HeLa cells and no
  subset of size 349. That analysis is withdrawn and its underlying data cannot
  presently be identified.

**Overstated formal results.** The three "theorems" are direct consequences of
totally ordering cells by a scalar. They are stated as properties in v2.

**Uncited prior art.** Angular ordering of cells for cell-cycle inference is
established practice (tricycle, Cyclum, reCAT, CCPE, peco). v1 cited none of it
and claimed to introduce a polar coordinate prior. That claim is withdrawn.

### Effect of the corrections

Correcting the orientation asymmetry and the identifier faults reverses the sign
of three of four reported ACPL minus MST comparisons, one from -6.7 to +0.8
percentage points, and removes v1's central claim that MST outperformed ACPL on
three of five datasets.

### What v2 adds

- Corrected loaders for all five datasets (`R/loaders.R`).
- An origin-free rank statistic `tau_cyc` (`R/metrics.R`), needed because ACPL
  takes no starting cluster while its baselines are given one.
- Baselines run on their native inputs (principal components), not only on
  ACPL's embedding.
- Replication across ten embedding seeds.
- tricycle as a baseline. It exceeds ACPL on two of three datasets tested.
- A yeast experiment, where tricycle cannot run and ACPL recovers no significant
  signal. Reported as a negative result.

### Reproducibility note

The v1 release cannot regenerate the v1 figures on current Bioconductor. This is
most likely version drift rather than error at the time of writing, but the
consequence stands.
