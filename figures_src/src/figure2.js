/* Figure 2. Label changes in each sample: filtering and composition contrasts. */
function draw(D) {
  const W = 620, rowH = 19, top = 40, left = 150, right = 556;
  const samples = ["Oetjen_A", "Oetjen_P", "Oetjen_U", "Sun_1_CS", "Sun_2_KC", "Sun_3_TB", "Sun_4_TC"];
  const by = (contrast, s) => D.fig2.find((r) => r.contrast === contrast && r.sample === s);

  // Row layout: a header row per study, then its samples.
  const rows = [];
  let yy = top;
  for (const study of ["Oetjen", "Sun"]) {
    rows.push({ kind: "head", study, y: yy }); yy += rowH;
    for (const s of samples.filter((x) => studyOf(x) === study)) { rows.push({ kind: "sample", s, y: yy }); yy += rowH; }
    if (study === "Oetjen") yy += 6;
  }
  const plotTop = top - 8, plotBottom = yy - 4;
  const H = plotBottom + 34;
  const svg = makeSvg(W, H);
  const x = d3.scaleLinear().domain([0, 3.5]).range([left, right]);

  // Legend (shape encodes the contrast, colour encodes the study).
  const lg = svg.append("g");
  let lx = left;
  lx = legendItem(lg, lx, 14, "dot", C.ink2, "Filtering: unfiltered vs sample-filtered") + 18;
  legendItem(lg, lx, 14, "dot", C.ink2, "Composition: pool-filtered vs sample-filtered", { open: true });

  const g = svg.append("g");
  bottomAxis(g, x, { y: plotBottom, ticks: [0, 1, 2, 3], format: (d) => `${d}%`, grid: [plotTop, plotBottom],
    title: "Cells whose predicted label changed" });
  g.append("text").attr("class", "label-2").attr("x", W - 6).attr("y", top - 6).attr("text-anchor", "end").text("Cells");

  for (const r of rows) {
    if (r.kind === "head") {
      dot(g, 8, r.y - 3.5, studyColor(r.study), { r: 3.5 });
      g.append("text").attr("class", "group").attr("x", 17).attr("y", r.y).text(STUDY_NAME[r.study]);
      continue;
    }
    const del = by("deletion", r.s), comp = by("composition", r.s);
    const col = studyColor(r.s), cy = r.y - 3.5;
    g.append("text").attr("class", "label").attr("x", 17).attr("y", r.y).text(r.s);
    g.append("line").attr("x1", x(comp.pct)).attr("x2", x(del.pct)).attr("y1", cy).attr("y2", cy)
      .attr("stroke", col).attr("stroke-opacity", 0.45).attr("stroke-width", 2).attr("stroke-linecap", "round");
    dot(g, x(comp.pct), cy, col, { open: true });
    dot(g, x(del.pct), cy, col);
    g.append("text").attr("class", "tick").attr("x", W - 6).attr("y", r.y).attr("text-anchor", "end")
      .text(fmtInt(del.cells));
  }
}
