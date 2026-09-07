> Historical source record. Export path references sanitized; current metadata is in SOURCE_PROVENANCE.md and CITATION.cff.

# Repository and scientific audit

Root: `./`.
Branch: `v5-clean-rebuild-mainline-abf`.
HEAD: `1dd4323d19d3718b12988cf53b8e55633210aa01`.
Initial `git status --short --branch`: `## v5-clean-rebuild-mainline-abf` with no modified/untracked paths.
No AGENTS.md was found by recursive search including ignored paths. Checked applicable ancestor locations; none found.
Read root README, ROOT_LAYOUT.md, PIPELINE_RUNBOOK.md, the v5.4.8z2 task/rulebook, E6/E7 guardrails and v5.4.9g task/policy. Historical task instructions are recorded as provenance, not new authorization to commit or overwrite.

Located original implementations:

- E6: `tools/yrbp_v5_4_8u_e6_scoring_after_policy_repair.py`.
- E7: `tools/yrbp_v5_4_8v_e7_unlocked_rank.py`; score equals frozen E6 total with declared docking/contact tie-breaks.
- Composite: `tools/yrbp_v5_4_8z2_g_pre1_target_engineering_composite_rank.py`.
- C_verified registration and G_pro1 generation: `data/derived/YRBP_branch_runs/v5_4_9g_G_pro1_C_verified_anchor_feedback_sequence_complete_generate.py`.
- Original E6/E7/full827/Top50 and feedback outputs: corresponding versioned directories under `data/derived/YRBP_branch_runs/`.
- Environment: root requirements.txt had broad lower bounds, not a lock; this package records exact installed versions.
- Tests: repository tests exist, but no YRBP-named unit-test files were found. Historical YRBP validation was largely inline script/report assertions. New tests are dedicated to this package.

The original G_pro1 report asserted no hard coding, while its code assigned `status_score = 1.0 if acc == ANCHOR else 0.5` and shortcut anchor similarities to 1.0. The new implementation derives status from the registry and uses the same similarity calculations for every record. This removes those implementation shortcuts without changing the frozen G_pro1 numbers. It does not independently authenticate evidence claims.

The original composite used literature/context text and opened G_pre2 audit inputs. Packaged pre scoring has no C_verified/G_pro1/wet-lab dependency; it consumes only frozen numeric upstream evidence. Therefore the package tests support downstream computational separation, not a claim that all historical literature curation was blinded. Upstream E6/E7 pipelines are not rerun or represented as independently rebuilt.

The original domain-evidence path no longer exists. The real local file is `data/B3_interproscan/v5_2/B3_candidate_domain_evidence.csv`; row extraction and hashes are recorded. The frozen Top50 hash matches the specified published hash. The complete 827 table exists, has unique ranks 1–827, and all scores/ranks reproduce after preserving original intermediate precision.

No dirty-worktree branch/worktree confirmation was required because the initial worktree was clean. Only the new `paper_release/food_microbiology_rbp_workflow_v1.0/` tree was authored. No original generator was executed.
