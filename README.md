# RBP prioritization workflow version 1.1.0

This release accompanies the manuscript **Computational prioritization of phage receptor binding proteins under incomplete evidence with a Gp17 application case**. The contribution is an application-oriented candidate decision method with a preliminary Gp17 functional case. Publication of the manuscript is not implied.

## Scientific question and contribution

Selecting phage receptor-binding proteins for bacterial recognition reagents requires combining target-related support, structural interaction hypotheses and engineering evidence. The workflow makes these choices inspectable through identity checks, declared scoring policies, missing-evidence bounds and checks on receptor/structure declarations. Gp17 is an established recognition protein; the experimental case assesses the tagged preparation and preliminary magnetic capture in the reported systems. No new wet experiments were added for this release.

Historical prioritization and later retrospective extensions are distinguished. The frozen 827-record ranking represents 744 exact unique sequences. In the same-pool internal comparisons, Gp17 ranks first under the integrated policy, 250 with docking alone and eighth when target context is removed. With target context retained as an unknown quantity under the original weights, all 827 records can potentially enter the Top5 and none is guaranteed. These are decision comparisons, not a test of external predictive superiority.

The input extension generates sequence descriptors and raw E6 annotations and issues evidence requests for incomplete inputs. It does not regenerate all E6/E7 predictors. The external I7LEH0 record demonstrates intake and abstention; it is not a validated independent biological success. Receptor declarations and sequence correspondence do not establish a physiological receptor, a correct fold, assembly or docking grid.

## Contents

- The original `src`, `tests`, `configs`, `inputs`, `outputs` and `figures` retain the historical v1.0.0 scientific baseline without changes.
- `extension_v04/` contains the previously verified version 0.4.0 input/evidence code, source data, standard-library tests and portable reproduction entry.
- `submission_evidence_v5/` contains the paper's internal comparisons, five full-pool masking scenarios (4,135 record-scenario rows), ELISA derivatives, claim-to-source map and source hashes. Counts refer to accession records, not independent proteins or positive biological labels.
- `EXTENSION_RELEASE_MANIFEST_SHA256.csv` inventories the additions; `VERSION_1_1_0_MANIFEST_SHA256.csv` inventories the full version. Older manifests and publication reports describe v1.0.0 and are retained for provenance.

## Reproduction and validation

The portable extension instructions are in `extension_v04/reproduce_extension.py`. The standard-library mode verifies archived identity/missingness outputs; its optional scientific-dependency mode additionally runs the existing annotation consumer and ELISA reanalysis. It does not install environments or databases. The previously completed 24-test suite passed both locally and in the portable copy; these are the same 24 tests, not 48 distinct tests. Release preparation checks file identity and provenance and does not add biological predictions.

The historical baseline can be reproduced with `reproduce.ps1`, using the dependencies in `requirements.lock.txt`. New outputs should be kept in a separate output directory. Source-code and data availability does not establish prospective experimental chronology or the ability to predict every newly supplied phage end to end.

## Citation and licensing

Creators: Zongcheng Wu and Shiying Lu. Software: Apache-2.0, as in the prior release. Retained third-party data remain subject to their source terms and attribution. See `LICENSE`, `NOTICE`, `THIRD_PARTY_NOTICES.md` and the dated UniProt source records. The old `LICENSE_PENDING.txt` explicitly states that its historical placeholder has been superseded.

Historical software v1.0.0: https://doi.org/10.5281/zenodo.22625899

Historical research data: https://doi.org/10.17632/49smydnjjn.3

The old software DOI does not identify this update. Version 1.1.0 archival citation is to be added after a real repository receipt exists; no DOI has been invented or reused. See `CITATION.cff` for the current software metadata and `RELEASE_STATUS.json` for the packaging status.
