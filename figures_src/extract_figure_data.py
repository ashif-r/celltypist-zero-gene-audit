"""Extract the data shown in the figures from the saved analysis outputs.

Usage:
    python extract_figure_data.py RESULTS_DIR OUTPUT_DIR

RESULTS_DIR is the analysis results folder that contains main/ and
scib_secondary/. No model is run: every value is read from saved tables and
score archives.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = ["Oetjen_A", "Oetjen_P", "Oetjen_U", "Sun_1_CS", "Sun_2_KC", "Sun_3_TB", "Sun_4_TC"]
CONTRASTS = [("deletion", "mechanism"), ("composition", "composition")]
PAIRS = [("Tem/Effector helper T cells", "Tcm/Naive helper T cells"), ("Tem/Temra cytotoxic T cells", "NK cells")]


def short(name):
    return str(name).replace("Sun_sample", "Sun_")


def main(results, out):
    main_dir = results / "main"
    primary = pd.read_csv(main_dir / "tables" / "primary.csv")
    cells = pd.read_csv(main_dir / "tables" / "per_cell_primary.csv.gz")
    stress = pd.read_csv(main_dir / "tables" / "stress_summary.csv")
    labels = pd.read_csv(main_dir / "tables" / "per_A_label_rates.csv")
    minimum = pd.read_csv(main_dir / "tables" / "secondary_D.csv")
    queries = pd.read_csv(main_dir / "queries.csv")
    validation = pd.read_csv(main_dir / "validation.csv")
    scib = pd.read_csv(results / "scib_secondary" / "tables" / "primary.csv")
    for frame in (primary, cells, labels, minimum, queries, scib):
        if "sample" in frame:
            frame["sample"] = frame["sample"].map(lambda s: short(s) if isinstance(s, str) else s)
    out.mkdir(parents=True, exist_ok=True)

    # Figure 2: label changes per sample in both contrasts.
    sample_queries = validation.loc[validation.context.eq("sample")].merge(
        queries[["query_id", "sample"]], on="query_id", validate="one_to_one")
    removed = dict(zip(sample_queries["sample"], sample_queries["removed_model_genes"]))
    rows = []
    for contrast, endpoint in CONTRASTS:
        for sample in SAMPLES:
            r = primary.loc[primary["sample"].eq(sample) & primary.endpoint.eq(endpoint)].iloc[0]
            rows.append({"contrast": contrast, "study": r.study, "sample": sample, "cells": int(r.valid),
                         "changed": int(r.changed), "pct": 100 * int(r.changed) / int(r.valid),
                         "deleted_model_genes": int(removed[sample])})
    (out / "fig2_sample_rates.json").write_text(json.dumps(rows, indent=1))

    # Figure 3a: the first cell, by identifier, in the most frequent filtering transition.
    changed = cells.loc[cells.A_label.ne(cells.B_sample_label)]
    transitions = (changed.groupby(["A_label", "B_sample_label"]).size().reset_index(name="n")
                   .sort_values(["n", "A_label", "B_sample_label"], ascending=[False, True, True]))
    top = transitions.iloc[0]
    cell = (changed.loc[changed.A_label.eq(top.A_label) & changed.B_sample_label.eq(top.B_sample_label)]
            .sort_values("cell_id").iloc[0])
    query_id = queries.loc[queries.context.eq("sample") & queries["sample"].eq(cell["sample"]), "query_id"].iloc[0]
    archive = np.load(main_dir / "queries" / query_id / "observed.npz", allow_pickle=True)
    i = int(np.flatnonzero(archive["cell_ids"] == cell.cell_id)[0])
    a, b, classes = archive["A"][i], archive["B"][i], archive["classes"]
    chosen = np.unique(np.r_[np.argsort(a)[-3:], np.argsort(b)[-3:]])
    chosen = chosen[np.argsort(-a[chosen])]
    example = {
        "cell_id": str(cell.cell_id), "sample": cell["sample"],
        "classes": [{"name": str(classes[j]), "unfiltered": float(a[j]), "filtered": float(b[j]),
                     "offset": float(b[j] - a[j])} for j in chosen],
        "winner_unfiltered": str(classes[int(np.argmax(a))]), "winner_filtered": str(classes[int(np.argmax(b))]),
        "conf_unfiltered": float(cell.A_confidence), "conf_filtered": float(cell.B_sample_confidence),
    }
    (out / "fig3a_example_cell.json").write_text(json.dumps(example, indent=1))

    # Figure 3b: margin-rule categories over all cells.
    rows = []
    for _, endpoint in CONTRASTS:
        z = primary.loc[primary.endpoint.eq(endpoint)]
        valid, k, safe = int(z.valid.sum()), int(z.changed.sum()), int(z.certified_stable.sum())
        rows.append({"endpoint": endpoint, "certified_stable": safe, "uncertified_unchanged": valid - safe - k,
                     "changed": k, "valid": valid})
    pd.DataFrame(rows).to_csv(out / "fig3b_margin_rule.csv", index=False)

    # Figure 4: quantiles of the confidence change in each study and contrast.
    p = np.round(np.linspace(0, 1, 1001), 3)
    series = []
    for study, g in cells.groupby("study"):
        for contrast, reference in [("deletion", "A_confidence"), ("composition", "B_pool_confidence")]:
            delta = (g.B_sample_confidence - g[reference]).to_numpy()
            series.append({"study": study, "contrast": contrast, "n": int(delta.size),
                           "quantiles": np.quantile(delta, p).tolist(), "median": float(np.median(delta)),
                           "below_zero": int((delta < 0).sum())})
    (out / "fig4_confidence_quantiles.json").write_text(json.dumps({"p": p.tolist(), "series": series}))

    # Figure 5: random partitions by seed and by target size.
    parts = validation.merge(queries[["query_id", "size", "seed"]], on="query_id", validate="one_to_one")
    parts = parts.loc[parts["size"].notna()].copy()
    parts["size"] = parts["size"].astype(int)
    parts["seed"] = parts["seed"].astype(int)
    per_seed = (parts.groupby(["study", "size", "seed"])
                .agg(queries=("query_id", "size"), mean_deleted_model_genes_per_query=("removed_model_genes", "mean"),
                     shortest_partition_cells=("cells", "min")).reset_index())
    stress = stress.rename(columns={"actual_cells": "cells"})
    stress["deletion_percent"] = 100 * stress.mechanism_changed / stress.cells
    stress["composition_percent"] = 100 * stress.composition_changed / stress.cells
    by_seed = stress[["study", "size", "seed", "cells", "deletion_percent", "composition_percent"]].merge(
        per_seed, on=["study", "size", "seed"], validate="one_to_one")
    by_seed = by_seed[["study", "size", "seed", "queries", "cells", "deletion_percent", "composition_percent",
                       "mean_deleted_model_genes_per_query", "shortest_partition_cells"]]
    by_seed = by_seed.sort_values(["study", "size", "seed"])
    by_seed.to_csv(out / "fig5_partitions_by_seed.csv", index=False)
    by_size = (by_seed.groupby(["study", "size"])
               .agg(seeds=("seed", "nunique"), queries_per_seed=("queries", "first"),
                    deletion_percent=("deletion_percent", "mean"), composition_percent=("composition_percent", "mean"),
                    shortest_partition_cells=("shortest_partition_cells", "min")).reset_index())
    genes = parts.groupby(["study", "size"]).removed_model_genes.mean().rename("mean_deleted_model_genes_per_query")
    by_size = by_size.merge(genes.reset_index(), on=["study", "size"])
    by_size = by_size[["study", "size", "seeds", "queries_per_seed", "deletion_percent", "composition_percent",
                       "mean_deleted_model_genes_per_query", "shortest_partition_cells"]]
    by_size.to_csv(out / "fig5_partitions_by_size.csv", index=False)

    # Figure S1a: predicted and observed sample-filtered scores, first 32 cells of each sample.
    rows = []
    for r in queries.loc[queries.context.eq("sample")].itertuples():
        folder = main_dir / "queries" / r.query_id
        predicted = np.load(folder / "predicted_B.npz")["scores"][:32].ravel()
        observed = np.load(folder / "observed.npz", allow_pickle=True)["B"][:32].ravel()
        rows += [(r.study, round(float(x), 6), round(float(y), 6)) for x, y in zip(predicted, observed)]
    pd.DataFrame(rows, columns=["study", "predicted", "celltypist"]).to_csv(out / "figS1a_scores.csv", index=False)

    # Figure S1b: largest score discrepancy in every query.
    validation[["query_id", "study", "context", "removed_model_genes", "B_max_score_residual"]].to_csv(
        out / "figS1b_residuals.csv", index=False)

    # Figure S2: label groups with at least 50 cells (as saved by the analysis).
    labels.to_csv(out / "figS2_label_groups.csv", index=False)

    # Figure S4: secondary comparisons per sample.
    rows = []
    for sample in SAMPLES:
        original = primary.loc[primary["sample"].eq(sample) & primary.endpoint.eq("mechanism")].iloc[0]
        reduced = scib.loc[scib["sample"].eq(sample) & scib.endpoint.eq("mechanism")].iloc[0]
        three = minimum.loc[minimum["sample"].eq(sample)].iloc[0]
        rows.append({"sample": sample, "original_pct": 100 * int(original.changed) / int(original.valid),
                     "scib_pct": 100 * float(reduced.fraction), "min3_pct": 100 * float(three.fraction)})
    pd.DataFrame(rows).to_csv(out / "figS4_sensitivities.csv", index=False)

    # Figure S3: offset spread of every query, and label changes within spread quintiles of each query size.
    stress_queries = pd.read_csv(main_dir / "tables" / "stress_queries.csv")
    spread = validation[["query_id", "study", "context", "cells", "removed_model_genes", "shift_spread",
                         "certified_stable"]].merge(queries[["query_id", "sample", "size", "seed"]], on="query_id",
                                                    validate="one_to_one")
    spread = spread.merge(stress_queries[["query_id", "mechanism_changed", "mechanism_certified"]], on="query_id",
                          how="left", validate="one_to_one")
    part = spread.context.eq("stress")
    assert (spread.loc[part, "mechanism_certified"] == spread.loc[part, "certified_stable"]).all()
    spread["changed"] = spread["mechanism_changed"]
    for r in spread.loc[spread.context.eq("sample")].itertuples():
        spread.loc[r.Index, "changed"] = int(primary.loc[primary["sample"].eq(r.sample)
                                                         & primary.endpoint.eq("mechanism"), "changed"].iloc[0])
    for r in spread.loc[spread.context.eq("pool")].itertuples():
        z = cells.loc[cells.study.eq(r.study)]
        spread.loc[r.Index, "changed"] = int(z.A_label.ne(z.B_pool_label).sum())
    spread["certified"] = spread["certified_stable"]
    for col in ["size", "seed", "changed", "certified"]:
        spread[col] = spread[col].astype("Int64")
    class_shift = pd.read_csv(main_dir / "tables" / "all_query_class_shifts.csv.gz").pivot_table(
        index="query_id", columns="class", values="score_shift")
    assert np.allclose((class_shift.max(axis=1) - class_shift.min(axis=1)).reindex(spread.query_id).to_numpy(),
                       spread.shift_spread.to_numpy(), atol=1e-9)
    extra = pd.DataFrame({"largest_class": class_shift.idxmax(axis=1), "smallest_class": class_shift.idxmin(axis=1),
                          "temeff_minus_tcmnaive": class_shift[PAIRS[0][0]] - class_shift[PAIRS[0][1]],
                          "temra_minus_nk": class_shift[PAIRS[1][0]] - class_shift[PAIRS[1][1]]})
    spread = spread.merge(extra, left_on="query_id", right_index=True, validate="one_to_one")
    spread = spread[["query_id", "study", "context", "sample", "size", "seed", "cells", "removed_model_genes",
                     "shift_spread", "largest_class", "smallest_class", "temeff_minus_tcmnaive", "temra_minus_nk",
                     "changed", "certified"]]
    spread.to_csv(out / "figS3_query_spread.csv", index=False)
    part = spread.loc[spread.context.eq("stress")].copy()
    part["size"] = part["size"].astype(int)
    part["quintile"] = part.groupby(["study", "size"]).shift_spread.transform(
        lambda v: pd.qcut(v, 5, labels=False) + 1)
    quint = (part.groupby(["study", "size", "quintile"])
             .agg(spread_lo=("shift_spread", "min"), spread_hi=("shift_spread", "max"),
                  spread_median=("shift_spread", "median"), queries=("query_id", "size"), cells=("cells", "sum"),
                  changed=("changed", "sum"), certified=("certified", "sum")).reset_index())
    quint["changed_pct"] = 100 * quint.changed / quint.cells
    quint["certified_pct"] = 100 * quint.certified / quint.cells
    quint.to_csv(out / "figS3_spread_quintiles.csv", index=False)

    # Supplementary Tables S3a to S6 use the sample and pool score archives and the per-gene terms of the offsets.
    contrib = np.load(main_dir / "tables" / "model_contributions.npz", allow_pickle=True)
    a_kg = contrib["a_kg"]
    model_genes = [str(g) for g in contrib["genes"]]
    model_classes = [str(c) for c in contrib["classes"]]
    index = {g: i for i, g in enumerate(model_genes)}
    k = {c: i for i, c in enumerate(model_classes)}
    shifts = pd.read_csv(main_dir / "tables" / "class_shifts.csv")
    shifts["sample"] = shifts["sample"].map(lambda s: short(s) if isinstance(s, str) else s)
    archives, removed_sets = {}, {}
    for r in queries.loc[queries.context.isin(["sample", "pool"])].itertuples():
        key = r.sample if r.context == "sample" else r.study
        folder = main_dir / "queries" / r.query_id
        observed = np.load(folder / "observed.npz", allow_pickle=True)
        assert [str(c) for c in observed["classes"]] == model_classes
        archives[key] = observed
        genes = [str(g) for g in np.load(folder / "predicted_B.npz", allow_pickle=True)["removed_genes"]]
        removed_sets[key] = [g for g in genes if g in index]
        offset = a_kg[:, [index[g] for g in removed_sets[key]]].sum(axis=1)
        saved = shifts.loc[shifts.query_id.eq(r.query_id)].set_index("class").score_shift.reindex(model_classes)
        assert np.abs(saved.to_numpy() - offset).max() < 1e-9, r.query_id

    def offsets(genes):
        return a_kg[:, [index[g] for g in genes]].sum(axis=1)

    # Table S3a: rank of the new label in the reference scores, and top-two reference margins.
    rows = []
    for sample in SAMPLES:
        study = "Oetjen" if sample.startswith("Oetjen") else "Sun"
        arc, pool = archives[sample], archives[study]
        ids = [str(c) for c in arc["cell_ids"]]
        pos = pd.Index([str(c) for c in pool["cell_ids"]]).get_indexer(ids)
        assert (pos >= 0).all()
        z = cells.set_index("cell_id").loc[ids]
        new = arc["B"].argmax(axis=1)
        assert (np.array(model_classes)[new] == z.B_sample_label.to_numpy()).all()
        for contrast, ref, labels in [("filtering", arc["A"], z.A_label), ("composition", pool["B"][pos], z.B_pool_label)]:
            order = np.argsort(-ref, axis=1)
            assert (np.array(model_classes)[order[:, 0]] == labels.to_numpy()).all()
            rows_ = np.arange(len(ids))
            margin = ref[rows_, order[:, 0]] - ref[rows_, order[:, 1]]
            rank = (order == new[:, None]).argmax(axis=1) + 1
            changed = rank > 1
            rows.append({"contrast": contrast, "sample": sample, "cells": len(ids), "changed": int(changed.sum()),
                         "to_second": int((rank == 2).sum()), "to_third_or_lower": int((rank >= 3).sum()),
                         "lowest_rank": int(rank.max()),
                         "median_margin_changed": float(np.median(margin[changed])),
                         "median_margin_unchanged": float(np.median(margin[~changed])),
                         "max_margin_changed": float(margin[changed].max())})
    pd.DataFrame(rows).to_csv(out / "tableS3_rank_margin.csv", index=False)

    # Table S4: class offsets and NK assignments in each sample.
    rows = []
    for sample in SAMPLES:
        study = "Oetjen" if sample.startswith("Oetjen") else "Sun"
        d = offsets(removed_sets[sample])
        eta = d - offsets(removed_sets[study])
        z = cells.loc[cells["sample"].eq(sample)]
        rows.append({"sample": sample, "removed_model_genes": len(removed_sets[sample]),
                     "classes_positive": int((d > 0).sum()), "spread": float(d.max() - d.min()),
                     "largest_class": model_classes[int(d.argmax())], "largest_offset": float(d.max()),
                     "smallest_class": model_classes[int(d.argmin())], "smallest_offset": float(d.min()),
                     "nk_offset": float(d[k["NK cells"]]),
                     "temeff_minus_tcmnaive": float(d[k[PAIRS[0][0]]] - d[k[PAIRS[0][1]]]),
                     "temra_minus_nk": float(d[k[PAIRS[1][0]]] - d[k[PAIRS[1][1]]]),
                     "composition_smallest_class": model_classes[int(eta.argmin())],
                     "composition_smallest_offset": float(eta.min()),
                     "composition_largest_offset": float(eta.max()),
                     "nk_unfiltered": int(z.A_label.eq("NK cells").sum()),
                     "nk_pool_filtered": int(z.B_pool_label.eq("NK cells").sum()),
                     "nk_sample_filtered": int(z.B_sample_label.eq("NK cells").sum()),
                     "nk_gained": int((z.B_sample_label.eq("NK cells") & z.A_label.ne("NK cells")).sum()),
                     "epithelial_unfiltered": int(z.A_label.eq("Epithelial cells").sum()),
                     "epithelial_sample_filtered": int(z.B_sample_label.eq("Epithelial cells").sum()),
                     "fibroblasts_unfiltered": int(z.A_label.eq("Fibroblasts").sum()),
                     "fibroblasts_sample_filtered": int(z.B_sample_label.eq("Fibroblasts").sum()),
                     "fibroblasts_changed": int(((z.A_label.eq("Fibroblasts") | z.B_sample_label.eq("Fibroblasts"))
                                                 & z.A_label.ne(z.B_sample_label)).sum())})
    table4 = pd.DataFrame(rows)
    table4.to_csv(out / "tableS4_class_offsets.csv", index=False)

    # Table S5: removed genes with the largest contributions to offset differences.
    rows = []
    for sample in SAMPLES:
        cols = np.array([index[g] for g in removed_sets[sample]])
        largest = table4.loc[table4["sample"].eq(sample), "largest_class"].iloc[0]
        for comparison, gain, lose in [("largest_offset", largest, None), ("temeff_minus_tcmnaive", *PAIRS[0]),
                                       ("temra_minus_nk", *PAIRS[1])]:
            diff = a_kg[k[gain], cols] - (a_kg[k[lose], cols] if lose else 0)
            order = np.argsort(-diff)[:6]
            rows.append({"comparison": comparison, "sample": sample, "gaining_class": gain, "reference_class": lose or "",
                         "total": float(diff.sum()), "top6_share": float(diff[order].sum() / diff.sum()),
                         "top_genes": ";".join(f"{removed_sets[sample][i]}:{diff[i]:.3f}" for i in order)})
    for sample in SAMPLES:
        study = "Oetjen" if sample.startswith("Oetjen") else "Sun"
        assert set(removed_sets[study]) <= set(removed_sets[sample])
        extra = sorted(set(removed_sets[sample]) - set(removed_sets[study]))
        eta = a_kg[:, [index[g] for g in extra]]
        low = int(eta.sum(axis=1).argmin())
        if eta[low].sum() >= 0:
            continue
        order = np.argsort(eta[low])[:6]
        rows.append({"comparison": "negative_composition_offset", "sample": sample, "gaining_class": model_classes[low],
                     "reference_class": "", "total": float(eta[low].sum()), "top6_share": "",
                     "top_genes": ";".join(f"{extra[i]}:{eta[low, i]:.4f}" for i in order)})
    pd.DataFrame(rows).to_csv(out / "tableS5_offset_genes.csv", index=False)

    # Table S6: labels changed by removing each removed gene on its own from the unfiltered input.
    rows, summary = [], []
    for sample in SAMPLES:
        A = archives[sample]["A"]
        reference = A.argmax(axis=1)
        genes = removed_sets[sample]
        counts = np.array([int(((A + a_kg[:, index[g]]).argmax(axis=1) != reference).sum()) for g in genes])
        rows += [{"sample": sample, "gene": g, "labels_changed_alone": int(n)} for g, n in zip(genes, counts)]
        summary.append({"sample": sample, "removed_model_genes": len(genes), "genes_changing_labels": int((counts > 0).sum()),
                        "labels_changed_all_genes": int(((A + offsets(genes)).argmax(axis=1) != reference).sum()),
                        "sum_single_gene_changes": int(counts.sum()), "max_single_gene_changes": int(counts.max()),
                        "max_effect_genes": ";".join(sorted(g for g, n in zip(genes, counts) if n == counts.max()))})
    pd.DataFrame(rows).to_csv(out / "tableS6_single_genes.csv", index=False)
    pd.DataFrame(summary).to_csv(out / "tableS6_single_gene_summary.csv", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(Path(sys.argv[1]), Path(sys.argv[2]))
