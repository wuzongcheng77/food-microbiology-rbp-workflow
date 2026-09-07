# Food Microbiology RBP workflow v1.0.0

Paper: **An auditable sparse-evidence computational–experimental workflow for prioritizing phage receptor-binding proteins and Gp17-mediated magnetic capture of Yersinia enterocolitica O:3 in pork- and seafood-derived matrices**.

This standalone, offline package reproduces the frozen downstream computational ranking. Release v1.0.0 is publicly available on GitHub and archived in Zenodo. Author metadata and the Apache-2.0 selection were supplied by the author for this export. The PATENT ADVISORY NOTICE in LICENSE and NOTICE concerns only patent claims outside the Apache-2.0 Section 3 grant and does not restrict commercial use granted under that license. The author has approved public release. Historical review records are retained for provenance.

Run from PowerShell, with the pinned dependencies available:

```powershell
.\reproduce.ps1
# Or choose a fresh output directory:
.\reproduce.ps1 -Python python -OutputDirectory ./outputs/new_rbp_replay
python -B -m pytest -p no:cacheprovider tests -q
```

`environment.yml` and `requirements.lock.txt` pin Python and the observed direct/transitive Python packages. Creating a fresh environment requires a separate dependency installation; this task used existing packages without network access. This is a version lock, not a platform-specific wheel or Conda build lock. Reproduction has been checked on Windows/Python 3.12.3; portability to another platform requires rerunning the tests. All replay outputs go into a new directory; existing files are not overwritten. `scripts/assemble_inputs.py` is a provenance extraction helper for the original repository, and is **not** required for reproduction from this package.

The primary analysis uses **all 827 frozen candidates**, with an additional Top50 analysis for comparison. Frozen version: `v5.4.8z2_G_pre1_target_engineering_composite_rank`. Weights are 0.55 target context, 0.20 docking support, 0.20 solubility/no-chaperone, and 0.05 E7 residual. The original Top50 bytes have the specified SHA-256. No upstream docking, structure prediction, sequence fetching, or experimental assays are rerun.

The full baseline reproduces 827/827 scores and ranks; Top50 reproduces 50/50. Q9T0Z9 has score 0.856763 and rank 1. The baseline reconstructs the unrounded solubility subscore from its frozen subcomponents before final six-decimal formatting. Using the already rounded solubility column instead causes 41 one-unit differences in the sixth decimal across the full table; `precision_differences.csv` records each one. This is a precision correction to reproduction, not a weight change. Frozen tie-breaks are score, context, organism, phage, docking, solubility, E7 descending, then lexical accession.

All ablations, single-component scoring-source baselines and perturbations use score followed only by a seeded anonymous ordinal for exact ties. Omitted evidence is not reintroduced through tie-breaks. Spearman rho and Kendall tau compare the resulting unique ordinal ranks with frozen ranks; overlaps are counts out of 5 or 10. Tie-rich single-component results depend on this declared convention and are not classifier performance or prediction accuracy.

Identifier masking accepts frozen numeric components without metadata text; it cannot erase literature or target context already encoded in those numbers. Strict context-blind additionally removes target context and renormalizes the remaining weights. It does not remove literature-derived solubility evidence. For the full universe the target ranks are 1, 1 and 8 for full, identifier-masked and strict context-blind. The anonymous tie allocation uses NumPy PCG64, seed 20260906, and the frozen row order.

Weight sensitivity uses 10,000 independent four-factor Uniform(0.8,1.2) draws, seed 20260906, followed by normalization. Q9T0Z9 is Top1 and Top5 in 100% of these draws. Its margin over the best other candidate has mean 0.248211772 and central empirical 95% interval [0.217899277, 0.274700172]. These are perturbation intervals, not experimental confidence intervals. Every candidate has inclusion frequencies; original rank 51–827 candidates have separate Top50/10/5 results. None entered Top10 or Top5 in these draws; some entered Top50.

G_pro1 is **deterministic feedback, not machine-learning retraining**. Its scope remains the frozen Top50, with one directly verified candidate and 49 `not_tested` candidates. The latter retain the archived neutral status score of 0.5, not a failure label. The packaged implementation derives verification from the existing evidence registry and computes similarities uniformly. Original accession-specific status assignment and self-similarity shortcuts were removed; all 50 archived G_pro1 scores/ranks reproduce. Evidence registry entries are copied historical assertions about Figures 2–5, not independent verification of raw experiments.

`src/scoring.py` reads an exact numeric-only schema and cannot consume outcome columns. `src/feedback.py` is a separate downstream consumer. The historical composite generator also opened G_pre2 integrity tables for auditing; that runtime dependency is not included in this pre-feedback scorer. Tests cover schema rejection, identifier relabeling, frozen reproduction, not-tested semantics, expected hashes, and byte-identical repeated runs.

Key artifacts:

- `outputs/reference_full827/`: full analysis, complete rankings, trial-level weights/ranks/margins, S4 and Figure 6 source CSVs.
- `outputs/reference_top50/`: separately computed Top50 comparison.
- `figures/Figure6/`: 600 dpi PNG, LZW TIFF, editable text SVG and review preview; editable plotting source is `scripts/figure6.py`.
- `manifests/`: input provenance/hashes, expected output hashes, runtime and manuscript comparison.
- `IP_REVIEW_PACKAGE/`: proposed file inventory and unresolved review questions. This is not authorization to distribute.

The current manuscript's B-stage numerical table was not supplied or located. Its absent values remain marked `AUTHOR_COMPARISON_PENDING`; no manuscript numerical values were invented. The independent v1.0.0 release points to commit `604ee08829ac5ee628f4aa0c98fbb5a21e3c3e17`. Source provenance is recorded in SOURCE_PROVENANCE.md.

Corresponding author: Shiying Lu (lushiying1129@163.com). Source commit and export changes: see SOURCE_PROVENANCE.md.

## Code availability

The paper-specific source code, frozen computational inputs, and reproducibility scripts are publicly available in the GitHub repository https://github.com/wuzongcheng77/food-microbiology-rbp-workflow (release v1.0.0) and permanently archived in Zenodo at https://doi.org/10.5281/zenodo.22625899. The archive includes the complete 827-candidate G_pre1 analysis, robustness analyses, and deterministic C_verified/G_pro1 feedback workflow.

Version DOI: https://doi.org/10.5281/zenodo.22625899

Concept DOI (all versions): https://doi.org/10.5281/zenodo.22625898

## Software reference

Wu, Z., & Lu, S. (2026). Food Microbiology RBP Prioritization and Magnetic Bioseparation Workflow (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22625899

The DOI metadata was added on main after archiving; the immutable v1.0.0 tag remains unchanged.
