"""Run CellTypist on every near-tie case of the single-gene removals.

For each gene-cell comparison flagged by 65_single_gene_near_ties.py, the gene is
removed from the prepared input of the cell, CellTypist scores the cell through
src.inference.api_scores (best-match mode, no majority voting), and the result is
compared with the predicted scores, the unfiltered scores plus the gene's terms.

Needs the analysis environment (requirements.lock.txt), data/Immune_All_Low.pkl,
data/prepared/{study}_original.h5ad and results/main/{study}_reference.npz.

Usage: python scripts/66_single_gene_direct_check.py [--root DIR] [--cases CSV] [--output JSON]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--root", type=Path, default=R, help="analysis folder (default: this reproduction folder)")
parser.add_argument("--cases", type=Path, help="default: ROOT/results/single_gene_check/single_gene_boundary_cases.csv")
parser.add_argument("--output", type=Path, help="default: next to the cases file")
args = parser.parse_args()
root = args.root.resolve()
cases_path = (args.cases or root / "results/single_gene_check/single_gene_boundary_cases.csv").resolve()
out = (args.output or cases_path.parent / "direct_single_gene_near_ties.json").resolve()
sys.path.insert(0, str(root))

import anndata as ad  # noqa: E402
from celltypist.models import Model  # noqa: E402
from src.inference import api_scores  # noqa: E402

with np.load(root / "results/main/tables/model_contributions.npz", allow_pickle=False) as terms:
    a_kg = terms["a_kg"]
    column = {g: j for j, g in enumerate(terms["genes"].astype(str))}
    classes = terms["classes"].astype(str)
cases = pd.read_csv(cases_path)
model = Model.load(str(root / "data/Immune_All_Low.pkl"))
report = []
for study in ["Oetjen", "Sun"]:
    selected = cases[cases["sample"].str.startswith(study)]
    if selected.empty:
        continue
    data = ad.read_h5ad(root / f"data/prepared/{study}_original.h5ad")
    genes = data.var_names.to_numpy(dtype=str)
    ids = data.obs_names.to_numpy(dtype=str)
    position = {cell: i for i, cell in enumerate(ids)}
    with np.load(root / f"results/main/{study}_reference.npz", allow_pickle=False) as ref:
        np.testing.assert_array_equal(ref["cell_ids"].astype(str), ids)
        np.testing.assert_array_equal(ref["classes"].astype(str), classes)
        unfiltered = ref["scores"]
    for case in selected.to_dict("records"):
        i = position[case["cell"]]
        (j,) = np.flatnonzero(genes == case["gene"])
        if data.X[i, j] != 0:
            raise ValueError(f"{case['gene']} has a nonzero count in {case['cell']}")
        keep = genes != case["gene"]
        observed = api_scores(data.X[i:i + 1][:, keep], genes[keep], np.array([case["cell"]]), model)[0]
        predicted = unfiltered[i] + a_kg[:, column[case["gene"]]]
        top_two = np.sort(observed)[-2:]
        report.append({**case,
                       "observed_max_residual": float(np.max(np.abs(observed - predicted))),
                       "predicted_label": str(classes[predicted.argmax()]),
                       "native_label": str(classes[observed.argmax()]),
                       "label_match": bool(predicted.argmax() == observed.argmax()),
                       "native_changed": bool(observed.argmax() != unfiltered[i].argmax()),
                       "native_margin": float(top_two[1] - top_two[0])})
        print(case["sample"], case["gene"], "label match", report[-1]["label_match"],
              "residual", report[-1]["observed_max_residual"], flush=True)
    del data

result = {"comparisons": len(report),
          "all_labels_match": all(r["label_match"] for r in report),
          "largest_score_residual": max((r["observed_max_residual"] for r in report), default=0.0),
          "cases": report}
out.write_text(json.dumps(result, indent=2))
print("comparisons", result["comparisons"], "| all labels match", result["all_labels_match"],
      "| largest score residual", result["largest_score_residual"])
if not result["all_labels_match"]:
    raise SystemExit("A CellTypist label differs from the predicted label")
