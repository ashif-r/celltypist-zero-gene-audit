/* Shared tokens and helpers for all figures. */
const C = {
  surface: "#ffffff",
  ink: "#0b0b0b",
  ink2: "#52514e",
  muted: "#898781",
  grid: "#e1e0d9",
  axis: "#c3c2b7",
  oetjen: "#2a78d6",   // categorical slot 1
  sun: "#eb6834",      // categorical slot 2
  change: "#e34948",   // accent for changed labels
  keptDark: "#52514e",
  keptLight: "#c3c2b7",
  blue: ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
         "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"],
};
const STUDY_NAME = { Oetjen: "Oetjen (bone marrow)", Sun: "Sun (peripheral blood)" };
const studyOf = (s) => (s.startsWith("Oetjen") ? "Oetjen" : "Sun");
const studyColor = (s) => (studyOf(s) === "Oetjen" ? C.oetjen : C.sun);

function makeSvg(w, h) {
  return d3.select("#fig").append("svg")
    .attr("xmlns", "http://www.w3.org/2000/svg")
    .attr("width", w).attr("height", h).attr("viewBox", `0 0 ${w} ${h}`);
}

function panelHead(svg, x, y, letter, title) {
  svg.append("text").attr("class", "panel-letter").attr("x", x).attr("y", y).text(letter);
  if (title) svg.append("text").attr("class", "panel-title").attr("x", x + 14).attr("y", y).text(title);
}

/* Horizontal axis at the bottom of a plot area with optional vertical gridlines. */
function bottomAxis(g, x, { y, ticks, format = (d) => d, grid = null, title = null, titleDy = 26 }) {
  if (grid) {
    g.append("g").selectAll("line").data(ticks).join("line")
      .attr("x1", (d) => x(d)).attr("x2", (d) => x(d))
      .attr("y1", grid[0]).attr("y2", grid[1])
      .attr("stroke", C.grid).attr("stroke-width", 1).attr("shape-rendering", "crispEdges");
  }
  g.append("line").attr("x1", x.range()[0]).attr("x2", x.range()[1]).attr("y1", y).attr("y2", y)
    .attr("stroke", C.axis).attr("stroke-width", 1).attr("shape-rendering", "crispEdges");
  g.append("g").selectAll("text").data(ticks).join("text")
    .attr("class", "tick").attr("x", (d) => x(d)).attr("y", y + 12).attr("text-anchor", "middle")
    .text(format);
  if (title) {
    g.append("text").attr("class", "axis-title")
      .attr("x", (x.range()[0] + x.range()[1]) / 2).attr("y", y + titleDy)
      .attr("text-anchor", "middle").text(title);
  }
}

/* Vertical axis on the left with optional horizontal gridlines. */
function leftAxis(g, y, { x, ticks, format = (d) => d, grid = null, title = null, titleDx = -34 }) {
  if (grid) {
    g.append("g").selectAll("line").data(ticks).join("line")
      .attr("x1", grid[0]).attr("x2", grid[1])
      .attr("y1", (d) => y(d)).attr("y2", (d) => y(d))
      .attr("stroke", C.grid).attr("stroke-width", 1).attr("shape-rendering", "crispEdges");
  }
  g.append("line").attr("x1", x).attr("x2", x).attr("y1", y.range()[0]).attr("y2", y.range()[1])
    .attr("stroke", C.axis).attr("stroke-width", 1).attr("shape-rendering", "crispEdges");
  g.append("g").selectAll("text").data(ticks).join("text")
    .attr("class", "tick").attr("x", x - 5).attr("y", (d) => y(d) + 3.3).attr("text-anchor", "end")
    .text(format);
  if (title) {
    const cy = (y.range()[0] + y.range()[1]) / 2;
    g.append("text").attr("class", "axis-title")
      .attr("transform", `translate(${x + titleDx},${cy}) rotate(-90)`)
      .attr("text-anchor", "middle").text(title);
  }
}

/* Filled marker with a 2px surface ring. */
function dot(g, cx, cy, color, { r = 4, open = false } = {}) {
  g.append("circle").attr("cx", cx).attr("cy", cy).attr("r", r + 1)
    .attr("fill", C.surface);
  g.append("circle").attr("cx", cx).attr("cy", cy).attr("r", open ? r - 0.75 : r)
    .attr("fill", open ? C.surface : color)
    .attr("stroke", color).attr("stroke-width", open ? 1.5 : 0);
}

function legendItem(g, x, y, kind, color, text, opts = {}) {
  if (kind === "dot") dot(g, x + 4, y - 3.5, color, opts);
  else if (kind === "line") {
    g.append("line").attr("x1", x).attr("x2", x + 14).attr("y1", y - 3.5).attr("y2", y - 3.5)
      .attr("stroke", color).attr("stroke-width", opts.width || 2).attr("stroke-linecap", "round")
      .attr("opacity", opts.opacity || 1);
  } else if (kind === "swatch") {
    g.append("rect").attr("x", x).attr("y", y - 8.5).attr("width", 10).attr("height", 10)
      .attr("rx", 2).attr("fill", color);
  }
  const t = g.append("text").attr("class", "legend").attr("x", x + (kind === "line" ? 19 : 14)).attr("y", y).text(text);
  return x + (kind === "line" ? 19 : 14) + t.node().getComputedTextLength();
}

/* A muted note on a white backing box, so it stays legible over gridlines and marks
   (a box instead of a stroked text halo keeps the PDF fonts as TrueType). */
function backedNote(g, x, y, s, anchor = "start") {
  const t = g.append("text").attr("class", "note").attr("x", x).attr("y", y).attr("text-anchor", anchor).text(s);
  const b = t.node().getBBox();
  g.insert("rect", () => t.node()).attr("x", b.x - 2).attr("y", b.y - 0.5)
    .attr("width", b.width + 4).attr("height", b.height + 1).attr("fill", C.surface);
  return t;
}

/* Dense point layers are drawn on a high-resolution canvas and placed in the SVG
   as one image, which keeps the PDF small; axes and text stay vector. */
function rasterPoints(g, x0, y0, w, h, points, scale = 6) {
  const cv = document.createElement("canvas");
  cv.width = Math.round(w * scale);
  cv.height = Math.round(h * scale);
  const ctx = cv.getContext("2d");
  ctx.scale(scale, scale);
  for (const p of points) {
    ctx.globalAlpha = p.alpha;
    ctx.fillStyle = p.color;
    ctx.beginPath();
    ctx.arc(p.x - x0, p.y - y0, p.r, 0, 2 * Math.PI);
    ctx.fill();
  }
  g.append("image").attr("x", x0).attr("y", y0).attr("width", w).attr("height", h)
    .attr("preserveAspectRatio", "none").attr("href", cv.toDataURL("image/png"));
}

const pct = (v, digits = 0) => `${v.toFixed(digits)}%`;
const fmtInt = d3.format(",");
