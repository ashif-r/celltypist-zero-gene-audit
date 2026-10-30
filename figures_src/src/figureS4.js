/* Figure S4. Secondary comparisons in the seven samples. */
function draw(D) {
  const W = 620, rowH = 19, top = 52;
  const samples = ["Oetjen_A", "Oetjen_P", "Oetjen_U", "Sun_1_CS", "Sun_2_KC", "Sun_3_TB", "Sun_4_TC"];
  const yOf = (i) => top + i * rowH + (i >= 3 ? 6 : 0);
  const bottom = yOf(6) + 12;
  const H = bottom + 34;
  const svg = makeSvg(W, H);
  const panels = [
    { letter: "a", title: "Original counts vs benchmark-processed input", x0: 0, dom: [0, 3.5], ticks: [0, 1, 2, 3],
      filled: ["original_pct", "Original counts"], open: ["scib_pct", "Benchmark-processed input"] },
    { letter: "b", title: "Undetected genes vs min_cells=3 filter", x0: 318, dom: [0, 5.5], ticks: [0, 1, 2, 3, 4, 5],
      filled: ["original_pct", "Undetected genes removed"], open: ["min3_pct", "min_cells=3 filter"] },
  ];
  for (const p of panels) {
    panelHead(svg, p.x0, 13, p.letter, p.title);
    const left = p.x0 + 66, right = p.x0 + 296;
    const x = d3.scaleLinear().domain(p.dom).range([left, right]);
    const g = svg.append("g");
    bottomAxis(g, x, { y: bottom, ticks: p.ticks, format: (d) => `${d}%`, grid: [top - 12, bottom],
      title: "Cells whose label changed" });
    let lx = p.x0;
    lx = legendItem(g, lx, 31, "dot", C.ink2, p.filled[1]) + 12;
    legendItem(g, lx, 31, "dot", C.ink2, p.open[1], { open: true });
    samples.forEach((s, i) => {
      const r = D.s4.find((q) => q.sample === s);
      const cy = yOf(i), col = studyColor(s);
      g.append("text").attr("class", "label").attr("x", p.x0).attr("y", cy + 3.3).text(s);
      const a = r[p.filled[0]], b = r[p.open[0]];
      g.append("line").attr("x1", x(a)).attr("x2", x(b)).attr("y1", cy).attr("y2", cy)
        .attr("stroke", col).attr("stroke-opacity", 0.45).attr("stroke-width", 2).attr("stroke-linecap", "round");
      dot(g, x(b), cy, col, { open: true });
      dot(g, x(a), cy, col);
    });
  }
}
