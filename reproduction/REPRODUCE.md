# Reproduce the study

Use Python 3.12 and the exact `requirements.lock.txt`. The tested platform is macOS arm64; cross-platform agreement has not been demonstrated. Keep this directory as the project root and run commands from it. Model/data files are downloaded from official sources and checked against manifests; they are not included in the release folder.

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.lock.txt
export MPLCONFIGDIR="$PWD/.mpl"
export NUMBA_CACHE_DIR="$PWD/.numba"
export CELLTYPIST_FOLDER="$PWD/data/celltypist"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p data/geo_original data/prepared logs results/gate results/input_repair
.venv/bin/python scripts/00_download.py
bash run_all.sh checkpoint
bash run_all.sh metadata
bash run_all.sh test
bash run_all.sh pilot
.venv/bin/python scripts/07_fetch_geo.py
.venv/bin/python scripts/08_scib_masks.py
.venv/bin/python scripts/09_rebuild_inputs.py
.venv/bin/python scripts/10_public_api_control.py
.venv/bin/python scripts/12_prepare_design.py
.venv/bin/python scripts/14_scib_secondary.py --prepare
```

GEO reconstruction keeps the scIB-selected cells and matches every shared count exactly; it performs no new biological QC. Original feature tables define the common exclusion of all repeated symbols. Normalization uses all original recorded features before these exclusions. Prepared H5AD hashes are platform/library dependent; exact raw source checksums, stable cell/feature identities, normalization totals, masks and numerical reproduction checks govern scientific agreement when container bytes differ.

Main inference requires the authorized prospective private-snapshot record described in protocol v0.4. No independent external preregistration is claimed. Do not fabricate a record to bypass it. A new reproduction after publication should retrieve the published record and frozen sources, verify their integrity and clearly identify its later run time; it is not a new preregistration. The main stage checks the prepared-input and pre-analysis source hashes, then runs all natural sample/pool queries and the 32/128/512-cell partitions, seeds 0–4, with all-class analytic and restoration controls.

```bash
.venv/bin/python scripts/11_main.py --run main
.venv/bin/python scripts/13_summarize.py --run main
.venv/bin/python scripts/14_scib_secondary.py
.venv/bin/python scripts/11_main.py --run reproduction
.venv/bin/python scripts/13_summarize.py --run reproduction
.venv/bin/python scripts/14_scib_secondary.py --reproduction
```

Each inference run creates a new directory and refuses to overwrite an existing one. An unexplained primary control failure halts interpretation. Saved query archives contain all A/B/C class scores, analytic predictions, cell/class order, exact deleted genes and stability flags. D is secondary; only an identified API input-contract failure may be recorded without invalidating A/B/C. Numerical D assertion failures remain failures.

`results/main/tables/primary.csv` contains the two endpoints by sample with explicit denominators and failures. `study_summary.csv` weights cells within studies; the overall summary averages study fractions equally. `stress_summary.csv` treats artificial partitions as computational probes, not biological replicates. Complete results and failures are reported irrespective of the 1% wording threshold. Biological correctness and calibration are not established by label agreement.

This file describes the reproducibility procedure, not evidence that every listed future stage has already passed. Consult `STATUS.json` and the completed-run/reproduction records for actual execution status.
