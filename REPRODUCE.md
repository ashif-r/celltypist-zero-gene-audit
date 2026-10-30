# Reproducing the analysis

The analysis ran with Python 3.12.13 on macOS arm64 with the packages pinned in `reproduction/requirements.lock.txt`, including CellTypist 1.7.1, NumPy 1.26.4, SciPy 1.13.1, scikit-learn 1.6.1, Scanpy 1.11.5 and AnnData 0.11.4. Agreement on other platforms has not been tested.

## 1. Check the reported numbers

```bash
python3 reproduction/scripts/verify_paper_numbers.py reproduction/results
```

This step reads only the saved result tables and needs pandas and NumPy.

## 2. Set up the environment and run the tests

```bash
cd reproduction
uv venv --python 3.12.13 .venv
uv pip install --python .venv/bin/python -r requirements.lock.txt
export CELLTYPIST_FOLDER="$PWD/data/celltypist" MPLCONFIGDIR="$PWD/.mpl" NUMBA_CACHE_DIR="$PWD/.numba" XDG_CACHE_HOME="$PWD/.cache" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
.venv/bin/python -m pytest tests -q
```

## 3. Rebuild the analysis from the public inputs

`reproduction/REPRODUCE.md` lists the stages in order, from the downloads and hash checks, the recovery of the original counts for the scIB-selected cells and the input preparation to the main and repeated runs of the 3,794 queries and the benchmark-processed input. A wrapper runs these stages in a new folder, adds the output checks and supplementary exports (scripts 16, 17, 19, 24 to 27, 32 and 35), and never overwrites existing results:

```bash
.venv/bin/python scripts/30_reproduce_project.py --output "$PWD/runs/new_complete_run"
```

The stages were run and checked individually for the paper. The wrapper that chains them has not been run end to end. Before the main queries, `src/freeze.py` checks every frozen file against the hashes in `manifests/freeze_inputs.json` and stops if one differs.

## 4. Figures and supplementary tables

`figures_src/README.md` describes the figure build. `figures_src/extract_figure_data.py` reads the result tables and the per-query score archives written in step 3. The archives are not stored in Git because of their size, and the data it extracted for the paper are in `figures_src/data/`.

## 5. Near-tie check of the single-gene removals

After step 3, `scripts/65_single_gene_near_ties.py` recomputes the single-gene scores from the saved archives of the sample queries, checks the label-change counts against the data of Supplementary Table S6 and lists the gene-cell comparisons within the near-tie rule of the main queries. `scripts/66_single_gene_direct_check.py` then removes each listed gene from the prepared input and scores the cell with CellTypist. The outputs for the paper are in `reproduction/results/single_gene_check/`.

```bash
.venv/bin/python scripts/65_single_gene_near_ties.py --root "$PWD/runs/new_complete_run" --output "$PWD/runs/single_gene_check"
.venv/bin/python scripts/66_single_gene_direct_check.py --root "$PWD/runs/new_complete_run" --cases "$PWD/runs/single_gene_check/single_gene_boundary_cases.csv"
```
