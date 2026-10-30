"""Screen the single-gene removals of the seven samples for numerical near ties.

Removing one undetected model gene g changes the score of class k by the gene's
term a_kg (Equation 2 of the paper), so the single-gene scores of a sample are
its unfiltered scores plus the terms of g. For every removed model gene of every
sample query, this script recomputes those scores, checks the number of changed
labels against the data of Supplementary Table S6 and flags each gene-cell
comparison whose top-two margin is within 2 * (1e-6 + 1e-6 * max_k |z_k|), the
near-tie rule of the main queries.

It reads saved outputs only and runs no model: results/main/queries.csv, the
observed.npz and predicted_B.npz archives of the sample queries (written by
11_main.py, not stored in Git because of their size),
results/main/tables/model_contributions.npz and the Table S6 data. The flagged
comparisons are checked with CellTypist by 66_single_gene_direct_check.py.

Usage: python scripts/65_single_gene_near_ties.py [--root DIR] [--table-s6 CSV] [--output DIR]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--root", type=Path, default=R, help="analysis folder (default: this reproduction folder)")
parser.add_argument("--table-s6", type=Path, default=R.parent / "figures_src" / "data" / "tableS6_single_genes.csv")
parser.add_argument("--output", type=Path, help="default: ROOT/results/single_gene_check")
args = parser.parse_args()
root = args.root.resolve()
out = (args.output or root / "results" / "single_gene_check").resolve()
out.mkdir(parents=True, exist_ok=True)

with np.load(root / "results/main/tables/model_contributions.npz", allow_pickle=False) as terms:
    a_kg = terms["a_kg"]
    column = {g: j for j, g in enumerate(terms["genes"].astype(str))}
table_s6 = pd.read_csv(args.table_s6)
queries = pd.read_csv(root / "results/main/queries.csv")

cases, samples, mismatches = [], [], []
for row in queries[queries.context.eq("sample")].itertuples():
    sample = row.sample.replace("Sun_sample", "Sun_")
    folder = root / "results/main/queries" / row.query_id
    with np.load(folder / "observed.npz", allow_pickle=False) as observed, \
            np.load(folder / "predicted_B.npz", allow_pickle=False) as predicted:
        scores = observed["A"]
        cells = observed["cell_ids"].astype(str)
        removed = [g for g in predicted["removed_genes"].astype(str) if g in column]
    expected = table_s6.loc[table_s6["sample"].eq(sample)].set_index("gene")["labels_changed_alone"]
    if sorted(removed) != sorted(expected.index):
        raise ValueError(f"{sample}: the removed model genes differ from the Table S6 data")
    reference = scores.argmax(axis=1)
    first = len(cases)
    exact_ties = 0
    for gene in removed:
        z = scores + a_kg[:, column[gene]]
        changed = z.argmax(axis=1) != reference
        if int(changed.sum()) != int(expected[gene]):
            mismatches.append({"sample": sample, "gene": gene, "changed": int(changed.sum()), "table_s6": int(expected[gene])})
        top_two = np.partition(z, -2, axis=1)[:, -2:]
        margin = top_two[:, 1] - top_two[:, 0]
        tolerance = 2 * (1e-6 + 1e-6 * np.abs(z).max(axis=1))
        exact_ties += int((margin == 0).sum())
        for i in np.flatnonzero(margin <= tolerance):
            cases.append({"sample": sample, "gene": gene, "cell": cells[i], "margin": float(margin[i]),
                          "twice_tolerance": float(tolerance[i]), "label_changed": bool(changed[i])})
    flagged = cases[first:]
    samples.append({"sample": sample, "query_id": row.query_id, "cells": int(len(cells)), "genes_tested": len(removed),
                    "gene_cell_comparisons": len(removed) * int(len(cells)), "near_ties": len(flagged),
                    "near_tie_cells": len({c["cell"] for c in flagged}),
                    "near_ties_counted_as_changes": sum(c["label_changed"] for c in flagged), "exact_ties": exact_ties})
    print(sample, "genes", len(removed), "near ties", len(flagged), flush=True)

frame = pd.DataFrame(cases, columns=["sample", "gene", "cell", "margin", "twice_tolerance", "label_changed"])
frame.to_csv(out / "single_gene_boundary_cases.csv", index=False)
summary = {
    "rule": "top-two margin <= 2 * (1e-6 + 1e-6 * max_k |z_k|)",
    "genes_tested": sum(s["genes_tested"] for s in samples),
    "gene_cell_comparisons": sum(s["gene_cell_comparisons"] for s in samples),
    "near_ties": len(frame),
    "near_tie_cells": int(frame.cell.nunique()),
    "near_ties_counted_as_changes": int(frame.label_changed.sum()),
    "exact_ties": sum(s["exact_ties"] for s in samples),
    "table_s6_count_mismatches": mismatches,
    "samples": samples,
}
(out / "single_gene_screen_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print("near ties", summary["near_ties"], "| cells", summary["near_tie_cells"], "| counted as changes",
      summary["near_ties_counted_as_changes"], "| exact ties", summary["exact_ties"],
      "| Table S6 mismatches", len(mismatches))
if mismatches:
    raise SystemExit("Single-gene counts differ from Table S6")
