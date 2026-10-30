/* Figure S2. Label changes within groups of cells sharing an unfiltered label. */
function draw(D) {
  const samples = ["Oetjen_A", "Oetjen_P", "Oetjen_U", "Sun_1_CS", "Sun_2_KC", "Sun_3_TB", "Sun_4_TC"];
  const codes = { Oetjen_A: "A", Oetjen_P: "P", Oetjen_U: "U", Sun_1_CS: "1_CS", Sun_2_KC: "2_KC", Sun_3_TB: "3_TB", Sun_4_TC: "4_TC" };
  const labels = Array.from(new Set(D.s2.map((r) => r.A_label))).sort((a, b) => a.localeCompare(b));
  const cw = 26, ch = 13, labW = 148, top = 52, gapP = 16;
  const W = 620, H = top + labels.length * ch + 30;
  const svg = makeSvg(W, H);
  const color = d3.scaleQuantize().domain([0, 35]).range(C.blue.slice(0, 12));
  const panels = [["a", "Filtering contrast", "mechanism_fraction", labW], ["b", "Composition contrast", "composition_fraction", labW + 7 * cw + gapP]];
  labels.forEach((l, i) => svg.append("text").attr("class", "label").attr("x", labW - 6).attr("y", top + i * ch + 9.3)
    .attr("text-anchor", "end").style("font-size", "9px").text(l));
  for (const [letter, title, key, x0] of panels) {
    panelHead(svg, x0, 13, letter, title);
    // column headers: study over sample code
    svg.append("text").attr("class", "note").attr("x", x0 + 1.5 * cw).attr("y", top - 18).attr("text-anchor", "middle").text("Oetjen");
    svg.append("text").attr("class", "note").attr("x", x0 + 5 * cw).attr("y", top - 18).attr("text-anchor", "middle").text("Sun");
    samples.forEach((s, j) => svg.append("text").attr("class", "tick").attr("x", x0 + j * cw + cw / 2).attr("y", top - 5)
      .attr("text-anchor", "middle").style("font-size", "8.5px").text(codes[s]));
    labels.forEach((l, i) => samples.forEach((s, j) => {
      const r = D.s2.find((q) => q.A_label === l && q.sample === s);
      const xx = x0 + j * cw, yy = top + i * ch;
      if (!r) {
        svg.append("rect").attr("x", xx + 0.75).attr("y", yy + 0.75).attr("width", cw - 1.5).attr("height", ch - 1.5)
          .attr("rx", 1.5).attr("fill", "#efeee9");
        return;
      }
      const v = 100 * r[key];
      const fill = color(v);
      svg.append("rect").attr("x", xx + 0.75).attr("y", yy + 0.75).attr("width", cw - 1.5).attr("height", ch - 1.5)
        .attr("rx", 1.5).attr("fill", fill);
      const dark = C.blue.indexOf(fill) >= 7;
      svg.append("text").attr("x", xx + cw / 2).attr("y", yy + 9.2).attr("text-anchor", "middle")
        .style("font-size", "8px").style("font-feature-settings", '"tnum" 1')
        .style("fill", dark ? "#ffffff" : C.ink).text(v === 0 ? "0" : v.toFixed(1));
    }));
  }
  // colour key
  const kx = labW + 14 * cw + gapP + 18, ky = top, kh = labels.length * ch * 0.6;
  const ks = d3.scaleLinear().domain([0, 35]).range([ky + kh, ky]);
  C.blue.slice(0, 12).forEach((c, i) => {
    const v0 = (35 / 12) * i, v1 = (35 / 12) * (i + 1);
    svg.append("rect").attr("x", kx).attr("y", ks(v1)).attr("width", 10).attr("height", ks(v0) - ks(v1) + 0.3).attr("fill", c);
  });
  [0, 10, 20, 30].forEach((t) => svg.append("text").attr("class", "tick").attr("x", kx + 14).attr("y", ks(t) + 3.3).text(`${t}%`));
  svg.append("text").attr("class", "note").attr("x", kx).attr("y", ky - 8).text("Changed");
  svg.append("rect").attr("x", kx).attr("y", ky + kh + 9).attr("width", 10).attr("height", 10).attr("rx", 1.5).attr("fill", "#efeee9");
  svg.append("text").attr("class", "note").attr("x", kx + 14).attr("y", ky + kh + 17.5).text("fewer than");
  svg.append("text").attr("class", "note").attr("x", kx + 14).attr("y", ky + kh + 28.5).text("50 cells");
}
