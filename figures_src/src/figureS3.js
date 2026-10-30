/* Figure S3. Offsets by query size, and label outcomes within spread quintiles of each partition size. */
function draw(D) {
  const W = 620, H = 468;
  const svg = makeSvg(W, H);
  const studyCol = { Oetjen: C.oetjen, Sun: C.sun };
  const sizes = [32, 128, 512];
  const shapeOf = { 32: "circle", 128: "square", 512: "triangle" };
  const cols = [{ x0: 0, left: 44, right: 296 }, { x0: 322, left: 366, right: 618 }];
  const rows = [{ title: 13, legend: 33, top: 48, bottom: 196 }, { title: 251, legend: 271, top: 286, bottom: 428 }];

  /* Marker with a surface ring; the shape encodes partition size. */
  function mark(g, cx, cy, shape, color, r = 3) {
    if (shape === "circle") return dot(g, cx, cy, color, { r });
    if (shape === "square") {
      const s = r * 1.75;
      g.append("rect").attr("x", cx - s / 2 - 1).attr("y", cy - s / 2 - 1).attr("width", s + 2).attr("height", s + 2)
        .attr("fill", C.surface);
      g.append("rect").attr("x", cx - s / 2).attr("y", cy - s / 2).attr("width", s).attr("height", s).attr("fill", color);
      return;
    }
    const t = r * 1.25;
    const pts = [[cx, cy - t * 1.15], [cx + t, cy + t * 0.75], [cx - t, cy + t * 0.75]];
    const d = "M" + pts.map((p) => p.join(",")).join("L") + "Z";
    g.append("path").attr("d", d).attr("fill", C.surface).attr("stroke", C.surface).attr("stroke-width", 2.5)
      .attr("stroke-linejoin", "round");
    g.append("path").attr("d", d).attr("fill", color);
  }

  /* Legends: study colours and query types for the top row; study colours and partition sizes for the bottom row. */
  const lg = svg.append("g");
  let lx = legendItem(lg, 0, rows[0].legend, "line", C.oetjen, "Oetjen") + 12;
  lx = legendItem(lg, lx, rows[0].legend, "line", C.sun, "Sun") + 18;
  lx = legendItem(lg, lx, rows[0].legend, "dot", C.ink2, "Sample") + 12;
  lx = legendItem(lg, lx, rows[0].legend, "dot", C.ink2, "Study pool", { open: true }) + 14;
  [[-3, -1.5], [0, -4.5], [2.5, -0.5], [4.5, -3.5], [6.5, -1]].forEach(([dx, dy]) =>
    lg.append("circle").attr("cx", lx + 2 + dx).attr("cy", rows[0].legend - 2 + dy).attr("r", 0.9)
      .attr("fill", C.ink2).attr("fill-opacity", 0.55));
  lg.append("text").attr("class", "legend").attr("x", lx + 14).attr("y", rows[0].legend).text("Random partitions, one dot each");
  lx = legendItem(lg, 0, rows[1].legend, "line", C.oetjen, "Oetjen") + 12;
  lx = legendItem(lg, lx, rows[1].legend, "line", C.sun, "Sun") + 18;
  for (const s of sizes) {
    mark(lg, lx + 4, rows[1].legend - 3.5, shapeOf[s], C.ink2, 3);
    const t = lg.append("text").attr("class", "legend").attr("x", lx + 12).attr("y", rows[1].legend).text(`${s}-cell partitions`);
    lx += 12 + t.node().getComputedTextLength() + 12;
  }

  /* (a, b) Offset spread and the Tem/Effector minus Tcm/Naive offset difference against query size. */
  const jitter = (i) => { const v = Math.sin(i * 12.9898 + 4.1414) * 43758.5453; return v - Math.floor(v) - 0.5; };
  const scatter = [
    { letter: "a", title: "Offset spread by query size", key: "shift_spread", dom: [0, 4.3], ticks: [0, 1, 2, 3, 4],
      ytitle: "Offset spread", fmt: (d) => d },
    { letter: "b", title: "Tem/Effector minus Tcm/Naive offset", key: "temeff_minus_tcmnaive", dom: [0, 2.6],
      ticks: [0, 0.5, 1, 1.5, 2, 2.5], ytitle: "Offset difference", fmt: (d) => d.toFixed(1) },
  ];
  scatter.forEach((p, j) => {
    const c = cols[j], r = rows[0];
    panelHead(svg, c.x0, r.title, p.letter, p.title);
    const x = d3.scaleLog().domain([7, 16000]).range([c.left, c.right]);
    const y = d3.scaleLinear().domain(p.dom).range([r.bottom, r.top]);
    const g = svg.append("g");
    leftAxis(g, y, { x: c.left, ticks: p.ticks, format: p.fmt, grid: [c.left, c.right], title: p.ytitle, titleDx: -32 });
    bottomAxis(g, x, { y: r.bottom, ticks: [10, 100, 1000, 10000], format: (d) => fmtInt(d),
      title: "Cells per query (log scale)" });
    const parts = D.s3q.filter((q) => q.context === "stress");
    const points = parts.map((q, i) => {
      const dodge = q.study === "Oetjen" ? 0.8 : 1.25;
      return { x: x(q.cells * dodge * Math.exp(jitter(i) * 0.22)), y: y(q[p.key]), r: 0.9,
        color: studyCol[q.study], alpha: 0.22 };
    });
    rasterPoints(g, c.left, r.top - 4, c.right - c.left, r.bottom - r.top + 8, points);
    for (const q of D.s3q.filter((v) => v.context !== "stress")) {
      dot(g, x(q.cells), y(q[p.key]), studyCol[q.study], { r: 3.2, open: q.context === "pool" });
    }
  });

  /* (c, d) Certified and changed labels within quintiles of offset spread, by study and partition size. */
  const lines = [
    { letter: "c", title: "Certified labels by spread quintile", key: "certified_pct", dom: [36, 82], ticks: [40, 50, 60, 70, 80] },
    { letter: "d", title: "Changed labels by spread quintile", key: "changed_pct", dom: [0, 8], ticks: [0, 2, 4, 6, 8] },
  ];
  lines.forEach((p, j) => {
    const c = cols[j], r = rows[1];
    panelHead(svg, c.x0, r.title, p.letter, p.title);
    const x = d3.scaleLinear().domain([1.4, 3.8]).range([c.left, c.right]);
    const y = d3.scaleLinear().domain(p.dom).range([r.bottom, r.top]);
    const g = svg.append("g");
    leftAxis(g, y, { x: c.left, ticks: p.ticks, format: (d) => `${d}%`, grid: [c.left, c.right],
      title: "Cells in quintile", titleDx: -34 });
    bottomAxis(g, x, { y: r.bottom, ticks: [1.5, 2, 2.5, 3, 3.5], format: (d) => d.toFixed(1),
      title: "Median offset spread of quintile" });
    for (const study of ["Oetjen", "Sun"]) {
      for (const s of sizes) {
        const q = D.s3b.filter((v) => v.study === study && v.size === s).sort((a, b) => a.quintile - b.quintile);
        g.append("path").attr("d", d3.line()(q.map((v) => [x(v.spread_median), y(v[p.key])])))
          .attr("fill", "none").attr("stroke", studyCol[study]).attr("stroke-width", 1.5)
          .attr("stroke-linejoin", "round").attr("stroke-opacity", 0.85);
        q.forEach((v) => mark(g, x(v.spread_median), y(v[p.key]), shapeOf[s], studyCol[study], 2.6));
      }
    }
  });
}
