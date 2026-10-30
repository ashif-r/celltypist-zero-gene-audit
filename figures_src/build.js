/* Build data.js from the saved figure data, write one HTML page per figure,
   then render each page to a vector PDF (for LaTeX) and a 3x PNG preview.
   Usage: node build.js [OUTPUT_DIR] [figure names...] */
const fs = require("fs");
const path = require("path");

const ROOT = __dirname;
const DATA_DIR = path.join(ROOT, "data");
const OUT_DIR = process.argv[2] || path.join(ROOT, "out");
const ONLY = process.argv.slice(3);

/* Locate an installed package folder (local node_modules first). */
function modDir(name) {
  const roots = [path.join(ROOT, "node_modules"), ...module.paths, "/opt/npm-tools/node_modules"];
  for (const r of roots) {
    const dir = path.join(r, name);
    if (fs.existsSync(path.join(dir, "package.json"))) return dir;
  }
  throw new Error(`Package ${name} not found; run npm install`);
}
const { chromium } = require(modDir("playwright"));
const url = (name, file) => "file://" + path.join(modDir(name), file);

function csv(name) {
  const [head, ...rows] = fs.readFileSync(path.join(DATA_DIR, name), "utf8").trim().split(/\r?\n/);
  const keys = head.split(",");
  return rows.map((r) => {
    const v = r.split(",");
    const o = {};
    keys.forEach((k, i) => {
      const x = v[i];
      o[k] = x !== "" && !isNaN(Number(x)) ? Number(x) : x;
    });
    return o;
  });
}
const json = (name) => JSON.parse(fs.readFileSync(path.join(DATA_DIR, name), "utf8"));

const s1aRows = csv("figS1a_scores.csv");
const DATA = {
  fig2: json("fig2_sample_rates.json"),
  fig3a: json("fig3a_example_cell.json"),
  fig3b: csv("fig3b_margin_rule.csv"),
  fig4: json("fig4_confidence_quantiles.json"),
  fig5seed: csv("fig5_partitions_by_seed.csv"),
  fig5size: csv("fig5_partitions_by_size.csv"),
  s1a: { study: s1aRows.map((r) => r.study), p: s1aRows.map((r) => r.predicted), c: s1aRows.map((r) => r.celltypist) },
  s1b: csv("figS1b_residuals.csv"),
  s2: csv("figS2_label_groups.csv"),
  s3q: csv("figS3_query_spread.csv"),
  s3b: csv("figS3_spread_quintiles.csv"),
  s4: csv("figS4_sensitivities.csv"),
};
fs.mkdirSync(path.join(ROOT, "html"), { recursive: true });
fs.writeFileSync(path.join(ROOT, "html", "data.js"), "window.DATA = " + JSON.stringify(DATA) + ";\n");

const FIGS = ["figure1", "figure2", "figure3", "figure4", "figure5", "figureS1", "figureS2", "figureS3", "figureS4"];
const page = (f) => `<!doctype html>
<html><head><meta charset="utf-8">
<title>${f.startsWith("figureS") ? "Supplementary Figure S" + f.slice(7) : "Figure " + f.slice(6)}</title>
<link rel="stylesheet" href="${url("katex", "dist/katex.min.css")}">
<link rel="stylesheet" href="../src/common.css">
<script src="${url("d3", "dist/d3.min.js")}"></script>
<script src="${url("katex", "dist/katex.min.js")}"></script>
<script src="../src/common.js"></script>
<script src="data.js"></script>
<script src="../src/${f}.js"></script>
</head><body><div id="fig"></div>
<script>
(async () => {
  await Promise.all(["400", "600", "700"].map((w) => document.fonts.load(w + " 10px InterFig")));
  await Promise.all(["10px KaTeX_Main", "italic 10px KaTeX_Math"].map((f) => document.fonts.load(f)));
  await document.fonts.ready;
  if (!document.fonts.check("400 10px InterFig")) console.error("InterFig font did not load");
  draw(window.DATA);
  await document.fonts.ready;
  window.__done = true;
})();
</script>
</body></html>`;
fs.mkdirSync(path.join(ROOT, "html"), { recursive: true });
for (const f of FIGS) fs.writeFileSync(path.join(ROOT, "html", `${f}.html`), page(f));

(async () => {
  fs.mkdirSync(OUT_DIR, { recursive: true });
  const browser = await chromium.launch();
  const ctx = await browser.newContext({ deviceScaleFactor: 3, viewport: { width: 800, height: 800 } });
  const pg = await ctx.newPage();
  pg.on("console", (m) => console.log(`[${m.type()}] ${m.text()}`));
  pg.on("pageerror", (e) => console.log("PAGE ERROR:", e.message));
  for (const f of FIGS) {
    if (ONLY.length && !ONLY.includes(f)) continue;
    await pg.goto("file://" + path.join(ROOT, "html", `${f}.html`));
    await pg.waitForFunction("window.__done === true", null, { timeout: 30000 });
    const el = await pg.$("#fig");
    const box = await el.boundingBox();
    const w = Math.ceil(box.width), h = Math.ceil(box.height);
    await pg.setViewportSize({ width: w, height: h });
    await pg.pdf({ path: path.join(OUT_DIR, `${f}.pdf`), width: `${w}px`, height: `${h}px`,
      printBackground: true, margin: { top: "0", right: "0", bottom: "0", left: "0" }, pageRanges: "1" });
    await el.screenshot({ path: path.join(OUT_DIR, `${f}.png`) });
    console.log(`${f}: ${w} x ${h}px`);
  }
  await browser.close();
})();
