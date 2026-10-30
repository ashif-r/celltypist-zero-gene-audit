/* Figure S1. Predicted versus CellTypist scores. */
function draw(D) {
  const W = 620, H = 236;
  const svg = makeSvg(W, H);
  // ---- Panel a ----
  panelHead(svg, 0, 13, "a", "First 32 cells of each sample, all 98 classes");
  {
    const left = 44, right = 280, top = 30, bottom = H - 40;
    const dom = [-18, 13];
    const x = d3.scaleLinear().domain(dom).range([left, right]);
    const y = d3.scaleLinear().domain(dom).range([bottom, top]);
    const g = svg.append("g");
    const ticks = [-15, -10, -5, 0, 5, 10];
    const f = (d) => String(d).replace("-", "−");
    leftAxis(g, y, { x: left, ticks, format: f, grid: [left, right], title: "CellTypist score, sample-filtered", titleDx: -30 });
    bottomAxis(g, x, { y: bottom, ticks, format: f, title: "Predicted score (unfiltered + class offset)" });
    g.append("line").attr("x1", x(dom[0])).attr("y1", y(dom[0])).attr("x2", x(dom[1])).attr("y2", y(dom[1]))
      .attr("stroke", C.ink2).attr("stroke-width", 1);
    const S = D.s1a;
    const pts = S.p.map((p, i) => ({ x: x(p), y: y(S.c[i]), r: 1.1, alpha: 0.35,
      color: S.study[i] === "Oetjen" ? C.oetjen : C.sun }));
    rasterPoints(g, left - 4, top - 4, right - left + 8, bottom - top + 8, pts);
    backedNote(g, left + 8, top + 10, "224 cells × 98 classes = 21,952 scores");
  }
  // ---- Panel b ----
  panelHead(svg, 318, 13, "b", "Largest discrepancy in each of 3,794 queries");
  {
    const left = 362, right = 608, top = 30, bottom = H - 40;
    const x = d3.scaleLinear().domain([0, 3600]).range([left, right]);
    const y = d3.scaleLinear().domain([0, 2.75]).range([bottom, top]);
    const g = svg.append("g");
    leftAxis(g, y, { x: left, ticks: [0, 0.5, 1, 1.5, 2, 2.5], format: (d) => d3.format(".1f")(d), grid: [left, right],
      title: "Largest |difference|, units of 10⁻¹⁴", titleDx: -30 });
    bottomAxis(g, x, { y: bottom, ticks: [0, 1000, 2000, 3000], format: fmtInt, title: "Model genes removed in the query" });
    const pts = D.s1b.map((r) => ({ x: x(r.removed_model_genes), y: y(r.B_max_score_residual * 1e14), r: 1.6,
      alpha: 0.45, color: r.study === "Oetjen" ? C.oetjen : C.sun }));
    rasterPoints(g, left - 4, top - 4, right - left + 8, bottom - top + 8, pts);
    const lg = svg.append("g");
    legendItem(lg, right - 112, top + 6, "dot", C.oetjen, "Oetjen", { r: 3 });
    legendItem(lg, right - 52, top + 6, "dot", C.sun, "Sun", { r: 3 });
  }
}
