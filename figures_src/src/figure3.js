/* Figure 3. (a) One cell whose label changes. (b) The margin rule across all cells. */
function draw(D) {
  const W = 620, H = 196;
  const svg = makeSvg(W, H);

  // ---- Panel a ----
  const A = D.fig3a;
  panelHead(svg, 0, 13, "a", "One cell changes label");
  svg.append("text").attr("class", "note").attr("x", 0).attr("y", 28)
    .text(`Oetjen_P cell, most frequent transition (confidence ${A.conf_unfiltered.toFixed(2)} → ${A.conf_filtered.toFixed(2)})`);
  const cls = A.classes;  // ordered by unfiltered score
  const left = 128, right = 308, top = 50, rowH = 30;
  const x = d3.scaleLinear().domain([-3.1, -1.0]).range([left, right]);
  const yOf = (i) => top + 8 + i * rowH;
  const plotBottom = top + cls.length * rowH + 2;
  const g = svg.append("g");
  bottomAxis(g, x, { y: plotBottom, ticks: [-3, -2.5, -2, -1.5, -1], format: (d) => d3.format(".1f")(d).replace("-", "−"),
    grid: [top - 6, plotBottom], title: "Decision score", titleDy: 25 });
  const short = { "Tcm/Naive helper T cells": "Tcm/Naive helper T", "Tem/Effector helper T cells": "Tem/Effector helper T", "B cells": "B cells" };
  cls.forEach((c, i) => {
    const cy = yOf(i);
    const winnerBefore = c.name === A.winner_unfiltered, winnerAfter = c.name === A.winner_filtered;
    g.append("text").attr("class", "label").attr("x", left - 8).attr("y", cy + 3.3).attr("text-anchor", "end")
      .style("font-weight", winnerBefore || winnerAfter ? 600 : 400).text(short[c.name]);
    const x0 = x(c.unfiltered), x1 = x(c.filtered);
    g.append("line").attr("x1", x0 + 5).attr("x2", x1 - 7).attr("y1", cy).attr("y2", cy)
      .attr("stroke", C.muted).attr("stroke-width", 1.2);
    g.append("path").attr("d", `M${x1 - 5},${cy} l-5,-3.2 l0,6.4 Z`).attr("fill", C.muted);
    dot(g, x0, cy, C.ink2, { open: true });
    if (winnerBefore) backedNote(g, x0 - 9, cy + 3.2, "unfiltered margin 0.002", "end");
    dot(g, x1, cy, winnerAfter ? C.change : C.ink);
    backedNote(g, (x0 + x1) / 2, cy - 6.5, `+${c.offset.toFixed(2)}`, "middle");
  });
  // legend for panel a
  const la = svg.append("g");
  let lx = 0;
  lx = legendItem(la, lx, H - 6, "dot", C.ink2, "Unfiltered", { open: true }) + 12;
  lx = legendItem(la, lx, H - 6, "dot", C.ink, "Sample-filtered") + 12;
  legendItem(la, lx, H - 6, "dot", C.change, "New top class");

  // ---- Panel b ----
  const bx = 336, bw = W - bx - 2;
  panelHead(svg, bx, 13, "b", "Margin rule across all 18,410 cells");
  const rows = [["Filtering contrast", D.fig3b.find((r) => r.endpoint === "mechanism")],
                ["Composition contrast", D.fig3b.find((r) => r.endpoint === "composition")]];
  const barH = 18;
  const xs = d3.scaleLinear().domain([0, 18410]).range([bx, bx + bw]);
  const gb = svg.append("g");
  const cats = [["certified_stable", C.keptDark, "certified"],
                ["uncertified_unchanged", C.keptLight, "unchanged"],
                ["changed", C.change, "changed"]];
  rows.forEach(([name, r], i) => {
    const y0 = 44 + i * 76;
    gb.append("text").attr("class", "group").attr("x", bx).attr("y", y0 - 6).text(name);
    let acc = 0;
    cats.forEach(([key, color], k) => {
      const v = r[key];
      const xa = xs(acc) + (k ? 1 : 0), xb = xs(acc + v) - (k < 2 ? 1 : 0);
      const w = Math.max(xb - xa, 1.5);
      if (k === 2) {
        const body = Math.max(w - 2, 0);
        gb.append("path").attr("d", `M${xa},${y0} h${body} a2,2 0 0 1 2,2 v${barH - 4} a2,2 0 0 1 -2,2 h${-body} Z`).attr("fill", color);
      } else {
        gb.append("rect").attr("x", xa).attr("y", y0).attr("width", w).attr("height", barH).attr("fill", color);
      }
      if (k === 0) gb.append("text").attr("class", "label").attr("x", xa + 7).attr("y", y0 + barH / 2 + 3.3)
        .style("fill", "#ffffff").text(`${(100 * v / r.valid).toFixed(1)}% certified by the rule`);
      acc += v;
    });
    // Two legend lines in bar order: certified; then not certified but unchanged, and changed.
    legendItem(gb, bx, y0 + barH + 15, "swatch", cats[0][1], `${fmtInt(r[cats[0][0]])} certified`);
    const lx = legendItem(gb, bx, y0 + barH + 29, "swatch", cats[1][1], `${fmtInt(r[cats[1][0]])} not certified, unchanged`) + 12;
    legendItem(gb, lx, y0 + barH + 29, "swatch", cats[2][1], `${fmtInt(r[cats[2][0]])} changed`);
  });
}
