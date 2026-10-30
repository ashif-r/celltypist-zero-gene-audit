/* Figure 5. Label changes after random partitioning into smaller queries. */
function draw(D) {
  const W = 620, H = 214;
  const svg = makeSvg(W, H);
  const sizes = [32, 128, 512];
  const panels = [
    { letter: "a", title: "Filtering: unfiltered vs partition-filtered", key: "deletion_percent", x0: 0 },
    { letter: "b", title: "Composition: pool- vs partition-filtered", key: "composition_percent", x0: 318 },
  ];
  for (const p of panels) {
    panelHead(svg, p.x0, 13, p.letter, p.title);
    const left = p.x0 + 40, right = p.x0 + 290, top = 30, bottom = H - 40;
    const x = d3.scaleLog().base(2).domain([26, 630]).range([left, right]);
    const y = d3.scaleLinear().domain([0, 7]).range([bottom, top]);
    const g = svg.append("g");
    leftAxis(g, y, { x: left, ticks: [0, 1, 2, 3, 4, 5, 6, 7], format: (d) => `${d}%`, grid: [left, right],
      title: p.letter === "a" ? "Cells whose label changed" : null, titleDx: -30 });
    bottomAxis(g, x, { y: bottom, ticks: sizes, format: (d) => d, title: "Cells per query (log scale)" });
    for (const study of ["Oetjen", "Sun"]) {
      const col = study === "Oetjen" ? C.oetjen : C.sun;
      const dx = study === "Oetjen" ? -4 : 4;  // small horizontal dodge so markers do not overlap
      const mean = sizes.map((s) => D.fig5size.find((r) => r.study === study && r.size === s));
      g.append("path").attr("d", d3.line()(mean.map((r) => [x(r.size) + dx, y(r[p.key])]))).attr("fill", "none")
        .attr("stroke", col).attr("stroke-width", 2).attr("stroke-linejoin", "round");
      mean.forEach((r) => dot(g, x(r.size) + dx, y(r[p.key]), col));
    }
    if (p.letter === "a") backedNote(g, left + 6, bottom - 7, "Restored vs unfiltered: no label changed in any partition");
    const lg = svg.append("g");
    legendItem(lg, right - 120, top + 6, "line", C.oetjen, "Oetjen");
    legendItem(lg, right - 58, top + 6, "line", C.sun, "Sun");
  }
}
