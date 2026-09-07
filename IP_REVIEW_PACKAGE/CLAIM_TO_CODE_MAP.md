# Claims and implementation boundaries

| Claim | Implementation/evidence | Boundary |
|---|---|---|
| Frozen pre-feedback ranking | src/scoring.py; baseline.csv; input hashes | Downstream frozen numeric reproduction, not upstream tool rerun |
| Weight and component robustness | scripts/run_analysis.py; S4; trial CSV | Scoring sensitivity, not experimental discrimination |
| Identifier/context masking | numeric schema and seeded tie keys | Cannot erase context encoded numerically |
| Deterministic G_pro1 | src/feedback.py; registry; archived Top50 | No training; one verified anchor; others not_tested |
| Magnetic capture / pork / seafood claims | Historical registry refers to Figures 2–5 | Raw assay and matrix data not supplied; no new experimental verification |
| Absence of candidate forcing in packaged scorer | identifier relabeling and source tests | Does not retroactively prove every historical curation choice was blind |

Authors must reconcile the manuscript wording with these boundaries before submission.
