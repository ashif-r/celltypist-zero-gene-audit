# Filtering undetected genes changes annotations from a released CellTypist model through query-dependent score offsets

Code, analysis plan, result tables and figure data for the preprint of the same title by Ashif Rahman (2026).

CellTypist models trained with mean centring give a gene with zero counts a nonzero standardized value, so the gene contributes to every class score. Removing genes that have zero counts in every cell of a query, the cells filtered together, therefore adds a class-specific offset to every score in that query. This repository contains the analysis that measured these offsets and the resulting label changes for the released Immune_All_Low v2 model on 18,410 cells from seven donor samples in two 10x Genomics studies (GEO GSE120221 and GSE128066), together with the checks showing that the offsets account for every change. The study measured whether labels changed, not whether they were correct.

## Contents

| Path | Content |
| --- | --- |
| `reproduction/protocol.md`, `reproduction/amendments.md` | Analysis plan and its dated amendments |
| `reproduction/deliverables/` | Snapshots of the plan and implementation taken before the main queries were run |
| `reproduction/config/` | Fixed design choices: model, thresholds, seeds, partition sizes and tolerances |
| `reproduction/src/`, `reproduction/scripts/`, `reproduction/tests/` | Analysis library, stage scripts and tests |
| `reproduction/manifests/` | Input, cell, query and archive manifests with SHA-256 hashes, and the freeze record |
| `reproduction/results/` | Result tables and verification records |
| `reproduction/results/single_gene_check/` | Near-tie check of the single-gene removals in Supplementary Table S6 |
| `figures_src/` | Figure code and the data behind every figure and supplementary table |
| `SHA256SUMS.txt` | SHA-256 of every other file in the release, checked with `sha256sum -c SHA256SUMS.txt` |

The analysis plan was fixed in a private snapshot on 23 October 2026, before the main queries were run, and was not publicly registered. A later implementation revision changed only execution speed and is described in `reproduction/amendments.md`. `reproduction/manifests/freeze_record.json` and `reproduction/manifests/history/` record the SHA-256 hashes of both snapshots, and `reproduction/src/freeze.py` stops the main analysis if any frozen file has changed.

## Inputs not stored here

The original count matrices, the scIB benchmark file, the CellTypist model and the per-query score archives are not in this repository. The download scripts retrieve the public inputs and check them against the hashes in `reproduction/manifests/`, and rerunning the analysis rebuilds the score archives (`REPRODUCE.md`). Sources and their terms are listed in `THIRD_PARTY.md`.

## Checking the reported numbers

```bash
python3 reproduction/scripts/verify_paper_numbers.py reproduction/results
```

The script recomputes the counts, rates, class offsets, confidence summaries, partition results and near-tie check reported in the paper from the saved result tables and prints each comparison. It needs only Python 3 with pandas and NumPy.

## Versions

- 1.0.1 (30 October 2026): adds the near-tie check of the single-gene removals (`reproduction/results/single_gene_check/`, scripts 65 and 66) and its numbers to `verify_paper_numbers.py`, pins the packages for regenerating the figure data (`figures_src/requirements.txt`), labels the vertical axes of Figure S3c,d as cell evaluations, because each cell appears once per seed, and adds `SHA256SUMS.txt`. The analysis code, result tables and figure data of version 1.0.0 are unchanged.
- 1.0.0 (30 October 2026): first public release.

## Licence

The source code is licensed under the MIT licence (`LICENSE`). The result tables, figure data, analysis plan and documentation are licensed under CC BY 4.0 (`LICENSE-DATA.md`). Third-party data, the CellTypist model and the Inter typeface keep their own terms (`THIRD_PARTY.md`).

## Citation

Rahman A. Filtering undetected genes changes annotations from a released CellTypist model through query-dependent score offsets. Preprint, 2026. Citation metadata are in `CITATION.cff`.
