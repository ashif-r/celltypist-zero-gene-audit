/* Figure 1. Why undetected genes matter, how the query decides them, and the study design. */
function draw(D) {
  const W = 620, H = 360;
  const root = d3.select("#fig").style("position", "relative").style("width", W + "px").style("height", H + "px");
  const svg = makeSvg(W, H);
  const math = (tex, x, y, size = 10, anchor = "left") => {
    const div = root.append("div").style("position", "absolute").style("top", y + "px")
      .style("font-size", size + "px").style("color", C.ink).style("white-space", "nowrap").style("line-height", "1");
    div.node().innerHTML = katex.renderToString(tex, { throwOnError: true });
    const w = div.node().getBoundingClientRect().width;
    div.style("left", (anchor === "center" ? x - w / 2 : anchor === "right" ? x - w : x) + "px");
    return div;
  };
  const box = (x, y, w, h, fill = "#f4f4f1") =>
    svg.append("rect").attr("x", x).attr("y", y).attr("width", w).attr("height", h).attr("rx", 6).attr("fill", fill);
  const text = (x, y, s, cls = "label", anchor = "start") =>
    svg.append("text").attr("class", cls).attr("x", x).attr("y", y).attr("text-anchor", anchor).text(s);

  // ---- Panel a: the score contribution of a recorded zero ----
  panelHead(svg, 0, 13, "a", "A recorded zero still contributes");
  const ay = 26, aw = 136, ah = 90, gap = 10;
  box(0, ay, aw, ah);
  box(aw + gap, ay, aw, ah);
  text(10, ay + 16, "Gene kept, count 0", "group");
  text(10, ay + 33, "Standardized value", "label-2");
  math("(0-\\mu_g)/\\sigma_g", 10, ay + 38);
  text(10, ay + 66, "Contribution to class score", "label-2");
  math("-w_{kg}\\,\\mu_g/\\sigma_g", 10, ay + 71);
  const bx = aw + gap;
  text(bx + 10, ay + 16, "Gene removed", "group");
  text(bx + 10, ay + 33, "No value in the input", "label-2");
  text(bx + 10, ay + 66, "Contribution to class score", "label-2");
  math("0", bx + 10, ay + 71);
  text(0, ay + ah + 18, "Removing the gene adds the same offset to a", "label-2");
  text(0, ay + ah + 30, "class score in every cell of the query:", "label-2");
  math("\\delta_k = w_{kg}\\,\\mu_g/\\sigma_g", 196, ay + ah + 19);

  // ---- Panel b: the query decides which genes are undetected ----
  const px = 316;
  panelHead(svg, px, 13, "b", "The query decides which genes are undetected");
  const genes = ["g1", "g2", "g3", "g4", "g5", "g6"];
  const M = [
    [1, 0, 1, 0, 1, 0], [1, 1, 0, 0, 1, 0], [0, 1, 1, 0, 1, 0],
    [1, 1, 0, 1, 1, 0], [1, 0, 1, 1, 0, 0], [0, 1, 1, 0, 1, 0],
  ];
  const cs = 13, mx = px + 54, my = 40, split = 5;
  const ry = (i) => my + i * cs + (i >= 3 ? split : 0);
  genes.forEach((gname, j) => text(mx + j * cs + cs / 2, my - 4, gname, "note", "middle"));
  M.forEach((row, i) => row.forEach((v, j) => {
    svg.append("rect").attr("x", mx + j * cs + 1).attr("y", ry(i) + 1)
      .attr("width", cs - 2).attr("height", cs - 2).attr("rx", 2)
      .attr("fill", v ? C.blue[5] : "#ffffff").attr("stroke", v ? "none" : C.grid).attr("stroke-width", 1);
  }));
  text(mx - 7, ry(1) + 10, "Sample 1", "label", "end");
  text(mx - 7, ry(4) + 10, "Sample 2", "label", "end");
  const outline = (j, i0, i1, color) => svg.append("rect").attr("x", mx + j * cs - 0.5).attr("y", ry(i0) - 0.5)
    .attr("width", cs + 1).attr("height", ry(i1 - 1) + cs - ry(i0) + 1).attr("rx", 3)
    .attr("fill", "none").attr("stroke", color).attr("stroke-width", 1.5);
  outline(3, 0, 3, C.change);
  outline(5, 0, 6, C.ink);
  const kx = mx + 6 * cs + 14;
  const key = (y, color, l1, l2, l3) => {
    svg.append("line").attr("x1", kx).attr("x2", kx + 11).attr("y1", y - 3.5).attr("y2", y - 3.5)
      .attr("stroke", color).attr("stroke-width", 1.5);
    text(kx + 16, y, l1, "label");
    text(kx + 16, y + 11.5, l2, "label-2");
    if (l3) text(kx + 16, y + 23, l3, "label-2");
  };
  key(my + 12, C.change, "g4: zero in sample 1 only", "removed when sample 1", "is filtered on its own");
  key(my + 60, C.ink, "g6: zero in every cell", "removed in both cases");
  text(px, ry(5) + cs + 22, "The same sample 1 cell is scored with g4 after pool", "label-2");
  text(px, ry(5) + cs + 34, "filtering and without it after sample filtering.", "label-2");

  // ---- Panel c: inputs compared for the same cells ----
  const cy0 = 186;
  panelHead(svg, 0, cy0, "c", "Inputs compared for the same 18,410 cells");
  const by = cy0 + 14, bh = 44;
  const node = (x, w, title, lines, fill) => {
    box(x, by, w, bh, fill);
    text(x + 9, by + 15, title, "group");
    lines.forEach((l, i) => text(x + 9, by + 28 + i * 11, l, "label-2"));
  };
  node(0, 150, "Original counts", ["7 samples, 2 studies", "18,410 cells"]);
  node(170, 190, "Normalize, then exclude", ["duplicate gene symbols", "(10,000 counts per cell, log1p)"]);
  node(380, 150, "Unfiltered", ["all measured genes with", "a unique gene symbol"], "#e6effb");
  const arrowHead = (x, y, dir) => {
    const s = 4.2;
    const d = dir === "right" ? `M${x},${y} L${x - s * 1.4},${y - s} L${x - s * 1.4},${y + s} Z`
                              : `M${x},${y} L${x - s},${y - s * 1.4} L${x + s},${y - s * 1.4} Z`;
    svg.append("path").attr("d", d).attr("fill", C.muted);
  };
  const hline = (x1, x2, y) => svg.append("line").attr("x1", x1).attr("x2", x2).attr("y1", y).attr("y2", y)
    .attr("stroke", C.muted).attr("stroke-width", 1.2);
  const vline = (x, y1, y2) => svg.append("line").attr("x1", x).attr("x2", x).attr("y1", y1).attr("y2", y2)
    .attr("stroke", C.muted).attr("stroke-width", 1.2);
  hline(150, 165, by + bh / 2); arrowHead(170, by + bh / 2, "right");
  hline(360, 375, by + bh / 2); arrowHead(380, by + bh / 2, "right");
  const dy = by + bh + 30, dh = 44, dw = 140;
  // Restored is built from the sample-filtered input, so it sits to its right with its own arrow.
  const derived = [
    ["Sample-filtered", "removes genes undetected", "in the sample", 0, true],
    ["Restored", "puts the removed genes", "back as zeros", 162, false],
    ["Pool-filtered", "removes genes undetected", "in the study pool", 316, true],
    ["min_cells=3", "removes genes detected", "in fewer than 3 cells", 470, true],
  ];
  const busY = by + bh + 15;
  vline(455, by + bh, busY);
  hline(dw / 2, 470 + dw / 2, busY);
  derived.forEach(([t, l1, l2, x, fed]) => {
    if (fed) { vline(x + dw / 2, busY, dy - 5); arrowHead(x + dw / 2, dy, "down"); }
    box(x, dy, dw, dh);
    text(x + 9, dy + 15, t, "group");
    text(x + 9, dy + 28, l1, "label-2");
    text(x + 9, dy + 39, l2, "label-2");
  });
  hline(dw, 157, dy + dh / 2); arrowHead(162, dy + dh / 2, "right");
  const ly = dy + dh + 20;
  text(0, ly, "Filtering contrast: unfiltered vs sample-filtered", "label");
  text(318, ly, "Composition contrast: pool-filtered vs sample-filtered", "label");
  text(0, ly + 14, "Control: restored vs unfiltered, expected to match exactly. Secondary comparison: min_cells=3 vs unfiltered.", "label-2");
}
