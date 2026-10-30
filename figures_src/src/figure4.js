/* Figure 4. Cumulative distributions of the change in CellTypist confidence. */
function draw(D) {
  const W = 620, H = 214;
  const svg = makeSvg(W, H);
  const P = D.fig4.p;
  const panels = [
    { letter: "a", title: "Filtering contrast", contrast: "deletion", x0: 0, dom: [-0.01, 0.25], ticks: [0, 0.05, 0.1, 0.15, 0.2, 0.25],
      xtitle: "Confidence change, sample-filtered − unfiltered" },
    { letter: "b", title: "Composition contrast", contrast: "composition", x0: 318, dom: [-0.008, 0.13], ticks: [0, 0.03, 0.06, 0.09, 0.12],
      xtitle: "Confidence change, sample-filtered − pool-filtered" },
  ];
  for (const p of panels) {
    panelHead(svg, p.x0, 13, p.letter, p.title);
    const left = p.x0 + 44, right = p.x0 + 300, top = 30, bottom = H - 40;
    const x = d3.scaleLinear().domain(p.dom).range([left, right]);
    const y = d3.scaleLinear().domain([0, 1]).range([bottom, top]);
    const g = svg.append("g");
    leftAxis(g, y, { x: left, ticks: [0, 0.25, 0.5, 0.75, 1], format: (d) => d3.format(".2f")(d),
      grid: [left, right], title: p.letter === "a" ? "Cumulative fraction of cells" : null, titleDx: -32 });
    bottomAxis(g, x, { y: bottom, ticks: p.ticks, format: (d) => (d === 0 ? "0" : d3.format("+.2f")(d)),
      title: p.xtitle });
    g.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", top).attr("y2", bottom)
      .attr("stroke", C.ink2).attr("stroke-width", 1).attr("shape-rendering", "crispEdges");
    const legendRows = [];
    for (const study of ["Oetjen", "Sun"]) {
      const s = D.fig4.series.find((r) => r.study === study && r.contrast === p.contrast);
      const pts = s.quantiles.map((q, i) => [x(Math.max(Math.min(q, p.dom[1]), p.dom[0])), y(P[i])]);
      g.append("path").attr("d", d3.line()(pts)).attr("fill", "none")
        .attr("stroke", study === "Oetjen" ? C.oetjen : C.sun).attr("stroke-width", 2)
        .attr("stroke-linejoin", "round").attr("stroke-linecap", "round");
      legendRows.push([study, s.median, s.below_zero, s.n]);
    }
    const lg = svg.append("g");
    legendRows.forEach(([study, med, below, n], i) => {
      const ly = bottom - 40 + i * 14;
      legendItem(lg, right - 128, ly, "line", study === "Oetjen" ? C.oetjen : C.sun,
        `${study}, median ${d3.format("+.3f")(med)}`);
    });
  }
}
