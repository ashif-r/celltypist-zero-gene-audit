# Near-tie check of the single-gene removals

Supplementary Table S6 counts the labels that removing each removed model gene on its own would change, computed by adding the gene's terms to the unfiltered scores. This folder records a check of those counts against numerical near ties.

| File | Content |
| --- | --- |
| `single_gene_screen_summary.json` | For each sample: removed model genes tested, gene-cell comparisons, comparisons within the near-tie rule, the cells involved, how many of them counted as label changes, and exact ties. Also lists any count that differs from the Table S6 data (none) |
| `single_gene_boundary_cases.csv` | The 33 gene-cell comparisons whose top-two margin was within 2(10^-6 + 10^-6 max_k \|z_k\|), the near-tie rule of the main queries, with the margin, the bound and whether the label counted as changed |
| `direct_single_gene_near_ties.json` | For each of the 33 comparisons, CellTypist's scores and label after the gene was removed from the prepared input, compared with the predicted scores and label |

The 3,529 removed model genes of the seven samples gave 8,762,606 gene-cell comparisons. Thirty-three comparisons in 17 cells fell within the near-tie rule, 18 of them counted as label changes, and none was an exact tie. CellTypist reproduced the predicted label in all 33, with a largest score difference of 5.37 x 10^-14, so no Table S6 count changes.

`scripts/65_single_gene_near_ties.py` writes the first two files from the saved score archives of the sample queries and runs no model. `scripts/66_single_gene_direct_check.py` removes each flagged gene from the prepared input and scores the cell with CellTypist. The archives and prepared inputs come from a full run (`REPRODUCE.md`). `direct_single_gene_near_ties.json` was produced by this CellTypist check in the analysis environment (macOS arm64, Python 3.12.13, CellTypist 1.7.1) on 30 October 2026.
