"""Recompute the numbers reported in the paper from the saved result tables.

Reads only the saved tables in RESULTS_DIR and prints one PASS or FAIL line per
reported value. Needs Python 3 with pandas and NumPy.

Usage: python3 reproduction/scripts/verify_paper_numbers.py reproduction/results
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(sys.argv[1])
T = R / "main" / "tables"
rows = []


def chk(name, paper, data, tol=0.0):
    if isinstance(paper, str) or isinstance(data, str):
        ok = str(paper) == str(data)
    else:
        ok = abs(float(paper) - float(data)) <= tol + 1e-12
    rows.append((ok, name, paper, data))


def pct(x, d=2):
    return round(100.0 * float(x), d)


def rd(path):
    """Read a CSV table and shorten Sun sample names to the form used in the paper."""
    df = pd.read_csv(path)
    if "sample" in df.columns:
        df["sample"] = df["sample"].map(lambda x: x.replace("Sun_sample", "Sun_") if isinstance(x, str) else x)
    return df


SAMPLES = ["Oetjen_A", "Oetjen_P", "Oetjen_U", "Sun_1_CS", "Sun_2_KC", "Sun_3_TB", "Sun_4_TC"]
NK, TCM, TEMEFF, TEMRA = "NK cells", "Tcm/Naive helper T cells", "Tem/Effector helper T cells", "Tem/Temra cytotoxic T cells"

pc = rd(T / "per_cell_primary.csv.gz")
pc["filt"] = pc.A_label != pc.B_sample_label
pc["comp"] = pc.B_pool_label != pc.B_sample_label
pc["pool"] = pc.A_label != pc.B_pool_label

# Cohort (Table 1)
chk("cells, total", 18410, len(pc))
chk("cells, Oetjen", 9581, (pc.study == "Oetjen").sum())
chk("cells, Sun", 8829, (pc.study == "Sun").sum())
for s, n in zip(SAMPLES, [2586, 3265, 3730, 1725, 2281, 2403, 2420]):
    chk(f"cells, {s}", n, (pc["sample"] == s).sum())

# Primary outcomes (Table 2, Results, Abstract)
F = [41, 58, 54, 36, 26, 57, 79]; Fp = [1.59, 1.78, 1.45, 2.09, 1.14, 2.37, 3.26]
C = [16, 22, 23, 21, 19, 30, 35]; Cp = [0.62, 0.67, 0.62, 1.22, 0.83, 1.25, 1.45]
fr, cr = [], []
for s, f, fp, c, cp in zip(SAMPLES, F, Fp, C, Cp):
    d = pc[pc["sample"] == s]
    chk(f"filtering changed, {s}", f, d.filt.sum()); chk(f"filtering %, {s}", fp, pct(d.filt.mean()))
    chk(f"composition changed, {s}", c, d.comp.sum()); chk(f"composition %, {s}", cp, pct(d.comp.mean()))
    fr.append(pct(d.filt.mean())); cr.append(pct(d.comp.mean()))
chk("filtering range low", 1.14, min(fr)); chk("filtering range high", 3.26, max(fr))
chk("composition range low", 0.62, min(cr)); chk("composition range high", 1.45, max(cr))
rate = {}
for st, f, fp, c, cp, pl, plp in [("Oetjen", 153, 1.60, 61, 0.64, 96, 1.00), ("Sun", 198, 2.24, 105, 1.19, 106, 1.20)]:
    d = pc[pc.study == st]
    chk(f"{st} filtering changed", f, d.filt.sum()); chk(f"{st} filtering %", fp, pct(d.filt.mean()))
    chk(f"{st} composition changed", c, d.comp.sum()); chk(f"{st} composition %", cp, pct(d.comp.mean()))
    chk(f"{st} pool vs unfiltered changed", pl, d.pool.sum()); chk(f"{st} pool vs unfiltered %", plp, pct(d.pool.mean()))
    rate[(st, "filt")] = d.filt.mean(); rate[(st, "comp")] = d.comp.mean()
chk("mean of study rates, filtering", 1.92, pct((rate[("Oetjen", "filt")] + rate[("Sun", "filt")]) / 2))
chk("mean of study rates, composition", 0.91, pct((rate[("Oetjen", "comp")] + rate[("Sun", "comp")]) / 2))
chk("all cells filtering", 351, pc.filt.sum()); chk("all cells filtering %", 1.91, pct(pc.filt.mean()))
chk("all cells composition", 166, pc.comp.sum()); chk("all cells composition %", 0.90, pct(pc.comp.mean()))


def tv(a, b):
    p, q = a.value_counts(normalize=True), b.value_counts(normalize=True)
    idx = p.index.union(q.index)
    return 0.5 * (p.reindex(idx, fill_value=0) - q.reindex(idx, fill_value=0)).abs().sum()


ratios = []
for st, tf, tc in [("Oetjen", 1.41, 0.55), ("Sun", 2.04, 1.08), ("All", 1.66, 0.77)]:
    d = pc if st == "All" else pc[pc.study == st]
    a, b = tv(d.A_label, d.B_sample_label), tv(d.B_pool_label, d.B_sample_label)
    chk(f"TV filtering, {st}", tf, pct(a)); chk(f"TV composition, {st}", tc, pct(b))
    if st != "All":
        ratios += [a / d.filt.mean(), b / d.comp.mean()]
chk("TV as % of rate, low", 87, round(100 * min(ratios))); chk("TV as % of rate, high", 91, round(100 * max(ratios)))

# Margin rule
chk("certified, filtering", 15337, pc.mechanism_certified.sum()); chk("certified %, filtering", 83.3, pct(pc.mechanism_certified.mean(), 1))
chk("certified, composition", 16974, pc.composition_certified.sum()); chk("certified %, composition", 92.2, pct(pc.composition_certified.mean(), 1))
chk("certified but changed, filtering", 0, (pc.mechanism_certified & pc.filt).sum())
chk("certified but changed, composition", 0, (pc.composition_certified & pc.comp).sum())
chk("uncertified unchanged, filtering", 2722, (~pc.mechanism_certified & ~pc.filt).sum())
chk("uncertified unchanged, composition", 1270, (~pc.composition_certified & ~pc.comp).sum())

# Confidence
dF = pc.B_sample_confidence - pc.A_confidence
chk("confidence rose, filtering", 18329, (dF > 0).sum()); chk("confidence rose %, filtering", 99.6, pct((dF > 0).mean(), 1))
chk("confidence median change, filtering", 0.0403, round(dF.median(), 4))
chk("confidence q25, filtering", 0.0081, round(dF.quantile(0.25), 4)); chk("confidence q75, filtering", 0.0917, round(dF.quantile(0.75), 4))
dec = pc[dF < 0]
chk("confidence decreases, filtering", 81, len(dec)); chk("decreases that are NK-labelled", 81, (dec.A_label == NK).sum())
chk("decreases Sun_2_KC", 11, (dec["sample"] == "Sun_2_KC").sum()); chk("decreases Sun_3_TB", 70, (dec["sample"] == "Sun_3_TB").sum())
dC = pc.B_sample_confidence - pc.B_pool_confidence
chk("confidence rose, composition", 17477, (dC > 0).sum()); chk("confidence rose %, composition", 94.9, pct((dC > 0).mean(), 1))
chk("confidence fell, composition", 933, (dC < 0).sum())
chk("confidence median change, composition", 0.0149, round(dC.median(), 4))
chk("confidence q25, composition", 0.0019, round(dC.quantile(0.25), 4)); chk("confidence q75, composition", 0.0413, round(dC.quantile(0.75), 4))
decC = pc[dC < 0].groupby(["sample", "B_pool_label"]).size()
chk("composition decreases Sun_1_CS Tcm/Naive", 851, decC.get(("Sun_1_CS", TCM), 0))
chk("composition decreases Sun_2_KC NK", 11, decC.get(("Sun_2_KC", NK), 0))
chk("composition decreases Sun_3_TB NK", 71, decC.get(("Sun_3_TB", NK), 0))
chk("composition decreases, other groups", 0, int(decC.sum()) - 851 - 11 - 71)
for name, ref, new in [("filtering", "A_confidence", "B_sample_confidence"), ("composition", "B_pool_confidence", "B_sample_confidence")]:
    expect = {"filtering": {0.5: (1185, 0), 0.9: (1445, 0), 0.95: (1316, 0)},
              "composition": {0.5: (545, 1), 0.9: (642, 1), 0.95: (548, 1)}}[name]
    for t, (up, down) in expect.items():
        chk(f"threshold {t} rose, {name}", up, ((pc[ref] <= t) & (pc[new] > t)).sum())
        chk(f"threshold {t} fell, {name}", down, ((pc[ref] > t) & (pc[new] <= t)).sum())
chk("threshold 0.50 rose %, filtering", 6.44, pct(((pc.A_confidence <= 0.5) & (pc.B_sample_confidence > 0.5)).mean()))
chk("threshold 0.90 rose %, filtering", 7.85, pct(((pc.A_confidence <= 0.9) & (pc.B_sample_confidence > 0.9)).mean()))
chk("median unfiltered confidence, changed", 0.235, round(pc[pc.filt].A_confidence.median(), 3))
chk("median unfiltered confidence, unchanged", 0.864, round(pc[~pc.filt].A_confidence.median(), 3))
chk("changed with confidence >= 0.90 before and after, filtering", 3, (pc.filt & (pc.A_confidence >= 0.9) & (pc.B_sample_confidence >= 0.9)).sum())
chk("changed with confidence >= 0.90 before and after, composition", 5, (pc.comp & (pc.B_pool_confidence >= 0.9) & (pc.B_sample_confidence >= 0.9)).sum())

# Transitions
tr = pc[pc.filt].groupby(["A_label", "B_sample_label"]).size().sort_values(ascending=False)
chk("transition Tcm/Naive to Tem/Effector", 73, tr.get((TCM, TEMEFF), 0))
chk("transition NK to Tem/Temra", 56, tr.get((NK, TEMRA), 0))
chk("two transitions combined", 129, tr.get((TCM, TEMEFF), 0) + tr.get((NK, TEMRA), 0))
chk("two transitions %", 36.8, pct((tr.get((TCM, TEMEFF), 0) + tr.get((NK, TEMRA), 0)) / pc.filt.sum(), 1))
chk("two most frequent transitions", "True", str(set(tr.index[:2]) == {(TCM, TEMEFF), (NK, TEMRA)}))
chk("transition Monocytes to Late erythroid", 14, tr.get(("Monocytes", "Late erythroid"), 0))

# Label groups (Figure S2)
grp = pc.groupby(["sample", "A_label"]).agg(n=("filt", "size"), ch=("filt", "sum"))
big = grp[grp.n >= 50].copy(); big["r"] = big.ch / big.n
chk("label groups with >= 50 cells", 75, len(big)); chk("labels spanned", 22, big.index.get_level_values(1).nunique())
chk("label-group rate low", 0.0, pct(big.r.min(), 1)); chk("label-group rate high", 32.9, pct(big.r.max(), 1))
top2 = list(big.r.sort_values(ascending=False).index[:2])
chk("two highest label-group rates are NK in Sun_4_TC and Sun_3_TB", "True", str(top2 == [("Sun_4_TC", NK), ("Sun_3_TB", NK)]))
for s, ch, n in [("Sun_4_TC", 53, 161), ("Sun_3_TB", 26, 96), ("Sun_1_CS", 15, 36), ("Sun_2_KC", 5, 15)]:
    chk(f"NK group {s}", f"{ch}/{n}", f"{int(grp.loc[(s, NK), 'ch'])}/{int(grp.loc[(s, NK), 'n'])}")
chk("NK group Sun_3_TB %", 27.1, pct(26 / 96, 1))

# NK assignments (Table S4, Results, Discussion)
sun = pc[pc.study == "Sun"]
chk("Sun NK unfiltered", 308, (sun.A_label == NK).sum()); chk("Sun NK pool-filtered", 241, (sun.B_pool_label == NK).sum())
chk("Sun NK sample-filtered", 209, (sun.B_sample_label == NK).sum())
chk("Sun NK fall %", 32.1, pct(1 - 209 / 308, 1))
chk("cells moved into NK", 0, ((pc.B_sample_label == NK) & (pc.A_label != NK)).sum())
chk("Sun NK-labelled cells that changed", 99, (sun.filt & (sun.A_label == NK)).sum())
nk = [(33, 23, 27), (2, 2, 2), (15, 15, 15), (36, 21, 25), (15, 10, 13), (96, 70, 77), (161, 108, 126)]
drops = []
for s, (a, b, p) in zip(SAMPLES, nk):
    d = pc[pc["sample"] == s]
    chk(f"NK {s} unfiltered/sample/pool", f"{a}/{b}/{p}", f"{(d.A_label == NK).sum()}/{(d.B_sample_label == NK).sum()}/{(d.B_pool_label == NK).sum()}")
    if s.startswith("Sun"):
        drops.append(pct(1 - b / a, 1))
chk("Sun per-sample NK decrease low", 27.1, min(drops)); chk("Sun per-sample NK decrease high", 41.7, max(drops))
chk("Oetjen_A NK decrease %", 30, round(100 * (1 - 23 / 33)))
chk("Epithelial cells before filtering", 0, (pc.A_label == "Epithelial cells").sum())
chk("Epithelial cells after sample filtering", 2, (pc.B_sample_label == "Epithelial cells").sum())
fib = pc[pc.A_label == "Fibroblasts"]
chk("Fibroblasts before filtering", 11, len(fib)); chk("Fibroblasts that changed", 0, fib.filt.sum())

# Class offsets (Table S4, Results)
cs = rd(T / "class_shifts.csv")
off = cs[cs.context == "sample"].pivot(index="sample", columns="class", values="score_shift").loc[SAMPLES]
pool = cs[cs.context == "pool"].pivot(index="study", columns="class", values="score_shift")
chk("classes per sample", 98, off.shape[1])
pos = (off > 0).sum(axis=1)
chk("classes raised per sample, low", 97, pos.min()); chk("classes raised per sample, high", 98, pos.max())
chk("offset minimum", -0.0182, round(off.values.min(), 4)); chk("offset maximum", 1.5512, round(off.values.max(), 4))
chk("NK smallest offset in every sample", "True", str(bool((off.idxmin(axis=1) == NK).all())))
exp_spread = [0.8098, 0.8020, 1.0002, 1.4129, 1.4378, 1.5162, 1.4747]
exp_large = [("Epithelial cells", 0.9738), ("Epithelial cells", 0.9802), ("Epithelial cells", 1.1288), ("Epithelial cells", 1.5512),
             ("Fibroblasts", 1.4196), ("Epithelial cells", 1.5161), ("Epithelial cells", 1.4970)]
exp_nk = [0.1640, 0.1783, 0.1286, 0.1383, -0.0182, -0.0001, 0.0224]
exp_d1 = [0.3424, 0.3036, 0.3387, 0.4573, 0.2150, 0.4140, 0.3359]
exp_d2 = [0.4983, 0.4167, 0.4678, 0.6874, 0.7654, 0.8899, 0.7268]
for i, s in enumerate(SAMPLES):
    o = off.loc[s]
    chk(f"offset spread {s}", exp_spread[i], round(o.max() - o.min(), 4))
    chk(f"largest offset {s}", f"{exp_large[i][0]}, {exp_large[i][1]:.4f}", f"{o.idxmax()}, {o.max():.4f}")
    chk(f"NK offset {s}", exp_nk[i], round(o[NK], 4))
    chk(f"Tem/Effector minus Tcm/Naive {s}", exp_d1[i], round(o[TEMEFF] - o[TCM], 4))
    chk(f"Tem/Temra minus NK {s}", exp_d2[i], round(o[TEMRA] - o[NK], 4))
spreads = (off.max(axis=1) - off.min(axis=1))
chk("sample spread Oetjen range", "0.80-1.00", f"{spreads[SAMPLES[:3]].min():.2f}-{spreads[SAMPLES[:3]].max():.2f}")
chk("sample spread Sun range", "1.41-1.52", f"{spreads[SAMPLES[3:]].min():.2f}-{spreads[SAMPLES[3:]].max():.2f}")
eta = off.copy()
for s in SAMPLES:
    eta.loc[s] = off.loc[s] - pool.loc[s.split("_")[0]]
chk("composition offset minimum", -0.0259, round(eta.values.min(), 4)); chk("composition offset maximum", 0.6602, round(eta.values.max(), 4))
neg = [(s, k, float(round(eta.loc[s, k], 4))) for s in SAMPLES for k in eta.columns if eta.loc[s, k] < 0]
chk("negative composition offsets", str([("Sun_1_CS", TCM, -0.0259), ("Sun_2_KC", NK, -0.0257), ("Sun_3_TB", NK, -0.0075)]), str(neg))
chk("NK composition offset Sun_4_TC", 0.0149, round(eta.loc["Sun_4_TC", NK], 4))

# Numerical checks and gene counts (validation.csv, Table 2, Table S2)
v = rd(R / "main" / "validation.csv")
chk("queries", 3794, len(v)); chk("partition queries", 3785, (v.context == "stress").sum())
chk("max predicted-minus-observed score difference (1e-14)", 2.49, round(v.B_max_score_residual.max() * 1e14, 2))
chk("max restored-minus-unfiltered difference", 0.0, v.C_max_score_residual.max())
chk("label mismatches", 0, v.B_destination_mismatches.sum()); chk("near ties", 0, v.B_near_ties.sum())
vs = v[v.context == "sample"].set_index(v[v.context == "sample"].query_id)
smp = rd(T / "model_feature_overlap.csv")
for s, n in zip(SAMPLES, [331, 313, 341, 641, 627, 651, 625]):
    chk(f"model genes removed, {s}", n, int(smp.loc[smp["sample"] == s, "B_removed_model_genes"].iloc[0]))
vp = v[v.context == "pool"].set_index("study")
chk("model genes removed, Oetjen pool", 172, vp.loc["Oetjen", "removed_model_genes"]); chk("model genes removed, Sun pool", 369, vp.loc["Sun", "removed_model_genes"])
chk("undetected genes in Oetjen pool", 10721, vp.loc["Oetjen", "removed_genes"]); chk("undetected genes in Sun pool", 12264, vp.loc["Sun", "removed_genes"])
for st, a, m in [("Oetjen", 5648, 991), ("Sun", 5256, 1383)]:
    row = smp[(smp.study == st) & (smp.context == "pool")].iloc[0]
    chk(f"model genes scored, {st}", a, row.A_matched_model_genes); chk(f"model genes absent or excluded, {st}", m, row.model_genes_absent_or_identifier_excluded_in_A)
chk("scoring columns Oetjen", 33626, smp[(smp.study == "Oetjen")].fixed_scoring_columns.iloc[0])
chk("scoring columns Sun", 32552, smp[(smp.study == "Sun")].fixed_scoring_columns.iloc[0])
chk("unmeasured model genes zero-filled in any query", "False", str(bool(smp.unmeasured_model_genes_zero_filled.any())))
chk("absent model genes % Oetjen", 14.9, pct(991 / 6639, 1)); chk("absent model genes % Sun", 20.8, pct(1383 / 6639, 1))

# Partitions (Figure 5, Results, Methods)
sq = rd(T / "stress_queries.csv")
per_seed = sq.groupby(["study", "size", "seed"]).size().unstack("seed")
for st, sizes in [("Oetjen", {32: 300, 128: 75, 512: 19}), ("Sun", {32: 276, 128: 69, 512: 18})]:
    for size, n in sizes.items():
        chk(f"partitions per seed {st} {size}", f"{n} x5", f"{int(per_seed.loc[(st, size)].iloc[0])} x{int((per_seed.loc[(st, size)] == n).sum())}")
ss = rd(T / "stress_summary.csv")
m = ss.groupby(["study", "size"])[["mechanism_fraction", "composition_fraction"]].mean()
for st, vals in [("Oetjen", [6.58, 5.42, 3.44]), ("Sun", [5.64, 4.43, 3.26])]:
    for size, val in zip([32, 128, 512], vals):
        chk(f"partition filtering mean {st} {size}", val, pct(m.loc[(st, size), "mechanism_fraction"]))
for st, size, val in [("Oetjen", 32, 5.66), ("Sun", 32, 4.70), ("Oetjen", 512, 2.49), ("Sun", 512, 2.22)]:
    chk(f"partition composition mean {st} {size}", val, pct(m.loc[(st, size), "composition_fraction"]))
chk("abstract 32-cell bone marrow", 6.6, pct(m.loc[("Oetjen", 32), "mechanism_fraction"], 1))
chk("abstract 32-cell peripheral blood", 5.6, pct(m.loc[("Sun", 32), "mechanism_fraction"], 1))
rng = ss.groupby(["study", "size"])[["mechanism_fraction", "composition_fraction"]].agg(lambda x: x.max() - x.min())
chk("largest difference between seeds (points)", 0.33, round(100 * rng.values.max(), 2))
chk("factor 512 to 32, Sun", 1.7, round(m.loc[("Sun", 32), "mechanism_fraction"] / m.loc[("Sun", 512), "mechanism_fraction"], 1))
chk("factor 512 to 32, Oetjen", 1.9, round(m.loc[("Oetjen", 32), "mechanism_fraction"] / m.loc[("Oetjen", 512), "mechanism_fraction"], 1))
vpart = v[v.context == "stress"].merge(sq[["query_id", "size"]], on="query_id")
mg = vpart.groupby(["study", "size"]).removed_model_genes.mean()
chk("mean model genes removed Oetjen 32", 2115, round(mg.loc[("Oetjen", 32)])); chk("mean model genes removed Sun 32", 2442, round(mg.loc[("Sun", 32)]))
chk("mean model genes removed Oetjen 512", 701, round(mg.loc[("Oetjen", 512)])); chk("mean model genes removed Sun 512", 1055, round(mg.loc[("Sun", 512)]))

# Partition offsets, recomputed from all_query_class_shifts (Results, Figure S3)
aq = rd(T / "all_query_class_shifts.csv.gz")
aq = aq[aq.query_id.isin(sq.query_id)]
P = aq.pivot_table(index="query_id", columns="class", values="score_shift")
info = sq.set_index("query_id").loc[P.index]
spread = P.max(axis=1) - P.min(axis=1)
d1 = P[TEMEFF] - P[TCM]; d2 = P[TEMRA] - P[NK]
med = pd.DataFrame({"spread": spread, "d1": d1, "d2": d2, "study": info.study, "size": info["size"]}).groupby(["study", "size"]).median()
cnt = info.groupby(["study", "size"]).size()
for st, vals, ns in [("Oetjen", [3.403, 2.740, 1.764], [1500, 375, 95]), ("Sun", [2.869, 2.712, 2.142], [1380, 345, 90])]:
    for size, val, n in zip([32, 128, 512], vals, ns):
        chk(f"median partition spread {st} {size}", val, round(med.loc[(st, size), "spread"], 3)); chk(f"partitions {st} {size}", n, cnt.loc[(st, size)])
chk("median Tem/Effector minus Tcm/Naive Oetjen 512", 0.575, round(med.loc[("Oetjen", 512), "d1"], 3))
chk("median Tem/Effector minus Tcm/Naive Sun 512", 0.594, round(med.loc[("Sun", 512), "d1"], 3))
chk("median Tem/Effector minus Tcm/Naive Oetjen 32", 0.982, round(med.loc[("Oetjen", 32), "d1"], 3))
chk("median Tem/Effector minus Tcm/Naive Sun 32", 1.215, round(med.loc[("Sun", 32), "d1"], 3))
chk("median Tem/Temra minus NK Oetjen 512", 1.023, round(med.loc[("Oetjen", 512), "d2"], 3))
chk("median Tem/Temra minus NK Oetjen 32", 1.858, round(med.loc[("Oetjen", 32), "d2"], 3))
chk("median Tem/Temra minus NK Sun range", "1.053-1.191", f"{med.loc['Sun', 'd2'].min():.3f}-{med.loc['Sun', 'd2'].max():.3f}")
chk("partitions with Epithelial largest and NK smallest", 3044, ((P.idxmax(axis=1) == "Epithelial cells") & (P.idxmin(axis=1) == NK)).sum())
sun32 = info[(info.study == "Sun") & (info["size"] == 32)].index
chk("Sun 32-cell partitions with Tcm/Naive smallest", 646, (P.loc[sun32].idxmin(axis=1) == TCM).sum())
part = info.assign(spread=spread).reset_index()
part["quintile"] = part.groupby(["study", "size"]).spread.transform(lambda x: pd.qcut(x, 5, labels=False) + 1)
q = part.groupby(["study", "size", "quintile"]).agg(n=("query_id", "size"), cells=("actual_cells", "sum"),
                                                     cert=("mechanism_certified", "sum"), ch=("mechanism_changed", "sum"))
q["cert_pct"] = 100 * q.cert / q.cells; q["ch_pct"] = 100 * q.ch / q.cells
chk("partitions per quintile", "18-300", f"{q.n.min()}-{q.n.max()}"); chk("cells per quintile", "7281-9728", f"{q.cells.min()}-{q.cells.max()}")
dc = [q.loc[(st, sz, 1), "cert_pct"] - q.loc[(st, sz, 5), "cert_pct"] for st in ["Oetjen", "Sun"] for sz in [32, 128, 512]]
dr = [q.loc[(st, sz, 5), "ch_pct"] - q.loc[(st, sz, 1), "ch_pct"] for st in ["Oetjen", "Sun"] for sz in [32, 128, 512]]
chk("certified fall, lowest to highest quintile", "3.6-7.0", f"{min(dc):.1f}-{max(dc):.1f}")
chk("change-rate shift, lowest to highest quintile", "-0.5-1.6", f"{min(dr):.1f}-{max(dr):.1f}")

# Benchmark-processed and min_cells=3 inputs (Figure S4, Table S7)
sb = rd(R / "scib_secondary" / "tables" / "primary.csv")
sm = sb[sb.endpoint == "mechanism"]
chk("benchmark-processed range low", 0.19, pct(sm.fraction.min())); chk("benchmark-processed range high", 0.66, pct(sm.fraction.max()))
chk("benchmark-processed contrasts coincide", "True",
    str(bool((sm.set_index("sample").changed == sb[sb.endpoint == "composition"].set_index("sample").changed).all())))
sd = rd(T / "secondary_D.csv")
chk("min_cells=3 queries passed", 9, (sd.status == "passed").sum())
exp_d = {"Oetjen": (160, 1.67), "Sun": (177, 2.00), "Oetjen_A": (76, 2.94), "Oetjen_P": (106, 3.25), "Oetjen_U": (96, 2.57),
         "Sun_1_CS": (53, 3.07), "Sun_2_KC": (49, 2.15), "Sun_3_TB": (86, 3.58), "Sun_4_TC": (122, 5.04)}
for _, r in sd.iterrows():
    key = r.study if r.context == "pool" else r["sample"]
    chk(f"min_cells=3 {key}", f"{exp_d[key][0]} ({exp_d[key][1]:.2f}%)", f"{int(r.changed)} ({pct(r.fraction):.2f}%)")
dsm = sd[sd.context == "sample"].set_index("sample").changed
chk("min_cells=3 changed more labels in all seven samples", "True", str(bool(all(dsm[s] > f for s, f in zip(SAMPLES, F)))))

# Pilot
pl = rd(R / "gate" / "pilot" / "paired_2700.csv")
chk("pilot cells", 2700, len(pl)); chk("pilot changed labels", 86, (pl.A_label != pl.B_label).sum())
chk("pilot changed %", 3.19, pct((pl.A_label != pl.B_label).mean()))

# Arithmetic stated in the text
chk("excluded scIB cells", 15096, 10727 + 3347 + 1022); chk("scIB cells minus included", 15096, 33506 - 18410)
chk("queries per seed", 757, 300 + 75 + 19 + 276 + 69 + 18); chk("partition queries total", 3785, 757 * 5)

fails = [r for r in rows if not r[0]]
for ok, name, paper, data in rows:
    print(("PASS" if ok else "FAIL"), "|", name, "| paper:", paper, "| data:", data)
print(f"\n{len(rows) - len(fails)} of {len(rows)} checks passed")
for ok, name, paper, data in fails:
    print("FAILED:", name, "| paper:", paper, "| data:", data)
sys.exit(1 if fails else 0)
