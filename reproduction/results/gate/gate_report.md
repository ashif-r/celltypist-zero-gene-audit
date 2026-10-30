# Feasibility gate report — 23 October 2026

**Decision: technical pilot passed; eligible input route recovered for two studies. Main inference has not started.** Gate review and the prospective public freeze remain outstanding. This package is private working material; it is not a submission-ready preprint.

## Checkpoint and environment

Released Immune_All_Low v2, official source https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl. SHA-256: `290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502`. The scaler has `with_mean=True`, `with_std=True`; there are 98 classes and 6639 features. Exact feature/class ordering digests, finite-parameter checks and serialization warnings are retained in `checkpoint.json`. Saved metadata says 18 studies and 20 tissues.

Python 3.12.13; CellTypist 1.7.1; NumPy 1.26.4; SciPy 1.13.1; scikit-learn 1.6.1. The checkpoint was serialized by scikit-learn 0.24.1. We checked native prediction behavior against independent float64 calculations; this does not establish universal cross-version compatibility. `requirements.lock.txt` records all installed packages. `logs/environment_repairs.md` records the initial incompatible wheel and implementation corrections.

## Model-only screen

`model_only.npz` contains all class/gene zero-deletion contributions. `model_gene_screen.csv` contains R_g and M_g for every feature. The screen timestamp is 2026-10-23T11:43:33.336704+00:00, before expression predictions. It did not select cells, genes, masks or settings. These are conditional score-shift diagnostics, not a new biological importance measure.

## PBMC3k pilot

Original 10x archive: 2700 cells and 32738 recorded features. The fixed duplicate rule excludes 186 rows representing 91 symbols, leaving 32552 features and 5256 model genes. Counts were normalized once over all original recorded features; excluded counts remain in the denominator. Full identifier and per-cell denominator records are supplied. Observed inverse-log row-sum deficits are expected after fixed exclusions; the API warned but consumed the intended matrix with no `.raw` fallback.

Predictions for B were saved before independent API B inference. All label destinations matched analytic predictions; independent aggregation found no B near ties in these four queries. The subqueries are the first archive-order cells at sizes 32/128/512, specified in the pilot implementation before inference. They are exploratory implementation probes, not independent samples or the frozen main stress test.

| Query cells | A/B changed labels | Percentage | Full-pool B/subquery B changed labels | Max score residual |
|---|---:|---:|---:|---:|
| 2700 | 86 | 3.19% | 0 | 1.78e-14 |
| 32 | 3 | 9.38% | 2 | 8.88e-15 |
| 128 | 8 | 6.25% | 5 | 8.88e-15 |
| 512 | 30 | 5.86% | 13 | 1.42e-14 |

The 86/2700 result demonstrates computational sensitivity; it does not show that any new label is biologically wrong, establish calibration, or estimate population prevalence. Query-size percentages are not evidence of a monotonic law.

25 recorded numerical controls passed, including all-class reconstruction, physical C restoration, gene-order permutation, fixed-schema subquery scoring, zero non-model deletion, bounded dense/sparse agreement, fresh checkpoint/reversed order and D's per-cell shift. Six synthetic tests passed. Confidence reconstruction and threshold crossing counts are recorded. Exact changes use valid paired cells; no failed cells were counted as stable. Original pilot outputs and a fresh reproduction are retained locally. See the reproduction record for exact comparisons.

## Main-data sample and feature audit

Original Figshare file: 33,506 cells, 12,303 genes, 2,064,748,344 bytes. Its MD5 matches Figshare (`482a6a3117432d41a46effd8376b8ee0`); SHA-256 is in `manifests/artifacts.json`. Data license: CC BY 4.0. `sample_ID` denotes five merged sources; it is not the ten natural query units. `batch` plus original preparation metadata identifies those units.

Oetjen contributes 3 natural samples / 9,581 cells and Sun 4 / 8,829 cells: **18,410 eligible cells across seven samples in two studies**. All checked counts are finite, nonnegative and integer-valued; no empty count libraries were found. Available-feature library totals and fixed exclusion mass are saved. Eight recorded symbols are excluded in Oetjen and 46 in Sun using source-table ambiguity, leaving 12,295 and 12,257 scoring features respectively.

Original gene tables were recovered from the seven source GEO samples referenced in the scIB preparation notebooks. Every scIB symbol occurs in each corresponding original feature table. After removing ambiguous symbols, retained names map exactly to unique original IDs; no suffix guessing is used. This provides measured-feature support for these seven samples even without inferring it from the merge file alone. Source notebook code retains raw counts before scran/log1p; the source-count route is used. Historical AnnData 0.6.22.post1 code independently confirms an inner concatenation default; the downloaded merged artifact still supplies the restricted gene universe, not the original complete transcriptome.

10X and Freytag remain excluded from main endpoints pending their original identifier provenance; their single query units do not qualify for the practical endpoint. Villani is marked TPM and excluded from the count route. These exclusions precede main predictions. The seven eligible samples are all 10x; no Smart-seq2 or across-technology practical claim is supported. See `sample_table.csv`, original source mapping and input audit.

## Review, freeze and publication status

Private repository created at https://github.com/ashif-r/celltypist-zero-gene-audit . The user requested privacy and control over edits. Its creation is infrastructure, not a public protocol freeze. No content is made public and no message is sent to another reviewer by this report.

Before main inference: review this gate, resolve the main implementation/runtime lock, check citations/novelty for the exact final contribution, and obtain the authorized timestamped public record as required by protocol v0.2. Verify the exact files without authentication. A private GitHub commit does not satisfy the current public-freeze rule. Any requested change to that scientific transparency rule needs a disclosed prospective amendment before main predictions.

The input structure supports the planned limited two-study route. A useful preprint remains conditional on natural-sample outcomes, full controls, independent reconciliation and accurate interpretation. There is no manuscript ready for bioRxiv yet and no guaranteed posting date. No independent human scientific review has occurred.

## Reproduction

From the masterplan folder, create Python 3.12 environment and install `requirements.lock.txt`. Retrieve the files listed in artifact/source manifests and verify hashes. Run `bash run_all.sh checkpoint`, `metadata`, `test`, `pilot`; then `python scripts/04_input_audit.py`, `05_verify_gate.py`, `06_gate_report.py` in the isolated environment. The runner sets local cache and CellTypist paths. Main execution remains disabled. The gate ZIP omits model/expression downloads, credentials and conversational history; it includes code, configs, diagnostics, relevant source evidence and retrieval records.
