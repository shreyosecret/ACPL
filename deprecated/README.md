# Superseded code from the v1 release

These scripts produced the results in engrXiv 10.31224/7087 version 1. They are
retained for the record and **should not be used**.

Three reasons.

1. **They do not run.** `load_datasets_v1.R` fails on four of five datasets
   because `scRNAseq` now returns Ensembl gene identifiers. Nestorowa and
   Richard intersect gene symbols against `ENSMUSG` IDs and produce no phase
   labels; Leng requests an assay `counts` and a column `phase` where the object
   provides `normcounts` and `Phase`; Spellman calls `data(spellman)`, which no
   longer resolves.

2. **They score asymmetrically.** `04_five_method_benchmark.R` and
   `08_multi_dataset_benchmark.R` pass `mirror = TRUE` for the baselines and the
   default `mirror = FALSE` for ACPL, giving the baselines a free maximum over
   two orientations. `02_spellman_ablation.R` does the same via `max(f, 1-f)`
   for MST only.

3. **They score against a wrong baseline.** All of them compare sliding-window
   structural accuracy to a stated 50 percent random baseline. The true chance
   level is `(1 + sum p^2)/2`, which is 66.7 to 69.3 percent on these datasets.

See `../CHANGELOG.md`.
