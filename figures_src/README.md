# Figures and supplementary table data

Source for the five main and four supplementary figures of *Filtering undetected genes changes annotations from a released CellTypist model through query-dependent score offsets*, and the data behind Supplementary Tables S3 to S6.

Each figure is drawn with D3 7.9.0 in a headless Chromium page (Playwright 1.63.0) and saved as a vector PDF, which the manuscript includes, and as a PNG preview at three times the page resolution. Figure 1 typesets its equations with KaTeX 0.16.47. Every figure is 620 CSS pixels wide, which prints at 164 mm, the text width of the manuscript, so the type sizes in the PDF are the printed sizes.

## Rebuild

```bash
python extract_figure_data.py RESULTS_DIR data
npm install
npx playwright install chromium
node build.js out
```

`RESULTS_DIR` is the analysis results folder, `reproduction/results` in the code repository. `extract_figure_data.py` reads the saved result tables and the per-query score archives and runs no model. The archives are written by the analysis and are not stored in Git because of their size. `node build.js out figure3` rebuilds a single figure. The manuscript uses the PDFs in `out/`.

## Contents

| Path | Content |
| --- | --- |
| `extract_figure_data.py` | Writes the plotted values and the supplementary table data in `data/` from the saved analysis outputs |
| `data/` | Plotted values, one file per figure panel or group of panels, and the data of Supplementary Tables S3a to S6 |
| `src/common.css`, `src/common.js` | Shared colours, type styles and drawing helpers |
| `src/figure*.js` | One drawing script per figure |
| `build.js` | Writes one HTML page per figure and renders the PDF and PNG files |
| `fonts/` | Inter 4.000 in TrueType form, its licence (`OFL.txt`) and the conversion script |
| `out/` | Rendered figures |

## Data files

| File | Figure or table | Content |
| --- | --- | --- |
| `fig2_sample_rates.json` | Figure 2 | Changed and valid cells per sample in both contrasts, with model genes removed |
| `fig3a_example_cell.json` | Figure 3a | Unfiltered and sample-filtered scores of the example cell for its three leading classes |
| `fig3b_margin_rule.csv` | Figure 3b | Certified, uncertified unchanged and changed labels in both contrasts |
| `fig4_confidence_quantiles.json` | Figure 4 | 1,001 quantiles of the confidence change per study and contrast |
| `fig5_partitions_by_seed.csv` | Figure 5 | Label-change rates and removed model genes for each seed and query size (the caption gives the largest difference between seeds) |
| `fig5_partitions_by_size.csv` | Figure 5 | Five-seed means for each study and query size |
| `figS1a_scores.csv` | Figure S1a | Predicted and observed sample-filtered scores, first 32 cells of each sample |
| `figS1b_residuals.csv` | Figure S1b | Largest score difference and model genes removed in each of the 3,794 queries |
| `figS2_label_groups.csv` | Figure S2 | Label-change rates in sample-by-label groups with at least 50 cells |
| `figS3_query_spread.csv` | Figure S3a,b | For every query: size, removed model genes, offset spread, the classes with the largest and smallest offsets, the Tem/Effector minus Tcm/Naive and Tem/Temra minus NK offset differences, and changed and certified labels in the filtering contrast |
| `figS3_spread_quintiles.csv` | Figure S3c,d | Changed and certified labels within quintiles of offset spread, for each study and partition size |
| `figS4_sensitivities.csv` | Figure S4 | Filtering-contrast rates with original counts, the scIB input and the min_cells=3 filter |
| `tableS3_rank_margin.csv` | Table S3a | Rank of the new label among the reference scores, and top-two reference margins of changed and unchanged cells, per sample and contrast |
| `tableS4_class_offsets.csv` | Table S4 | Offset spread, largest and smallest class offsets, NK offset, the two offset differences, composition offsets and NK counts per sample |
| `tableS5_offset_genes.csv` | Table S5 | Removed genes with the largest contributions to class offsets, offset differences and negative composition offsets |
| `tableS6_single_genes.csv` | Table S6 | Labels changed by removing each removed model gene on its own, for every sample |
| `tableS6_single_gene_summary.csv` | Table S6 | Per-sample summary of the single-gene removals |

In the data files, `deletion` and `mechanism` denote the filtering contrast (unfiltered versus sample-filtered), and `stress` denotes the random partitions.

## Design

Colours: Oetjen `#2a78d6`, Sun `#eb6834`, changed labels `#e34948`, neutral greys `#52514e` and `#c3c2b7`. Figure S2 uses a single-hue blue ramp. In Figure S3 the marker shape encodes partition size. Text is set in Inter at 7.1 pt for ticks and labels and 8.3 pt for panel titles. Inter was converted from its OpenType release to TrueType outlines (`fonts/otf2ttf.py`) so that the PDFs embed TrueType fonts instead of Type 3 fonts.
