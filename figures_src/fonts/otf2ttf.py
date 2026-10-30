"""Convert Inter OpenType (CFF) fonts to TrueType outlines with fontTools.

Usage:
    python otf2ttf.py SRC_DIR DST_DIR Regular SemiBold Bold

Chromium embeds TrueType-outline fonts in PDF output as TrueType, so the
figure PDFs contain no Type 3 fonts. Glyphs are converted with cu2qu
(maximum error 1 font unit); no glyph is otherwise changed.
"""
import sys
from fontTools.ttLib import TTFont, newTable
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.ttGlyphPen import TTGlyphPen

def convert(src, dst):
    f = TTFont(src)
    assert f.sfntVersion == "OTTO" and "CFF " in f
    order = f.getGlyphOrder()
    gs = f.getGlyphSet()
    quad = {}
    for name in gs.keys():
        tt = TTGlyphPen(gs)
        gs[name].draw(Cu2QuPen(tt, 1.0, reverse_direction=True))
        quad[name] = tt.glyph()
    f["loca"] = newTable("loca")
    f["glyf"] = glyf = newTable("glyf")
    glyf.glyphOrder = order
    glyf.glyphs = quad
    del f["CFF "]
    if "VORG" in f:
        del f["VORG"]
    glyf.compile(f)
    hmtx = f["hmtx"]
    for n, g in glyf.glyphs.items():
        if hasattr(g, "xMin"):
            hmtx[n] = (hmtx[n][0], g.xMin)
    f["maxp"] = maxp = newTable("maxp")
    maxp.tableVersion = 0x00010000
    for a in ("maxZones",): setattr(maxp, a, 1)
    for a in ("maxTwilightPoints", "maxStorage", "maxFunctionDefs", "maxInstructionDefs",
              "maxStackElements", "maxSizeOfInstructions"): setattr(maxp, a, 0)
    maxp.maxComponentElements = max(len(getattr(g, "components", []) or []) for g in glyf.glyphs.values())
    maxp.compile(f)
    post = f["post"]
    post.formatType = 2.0
    post.extraNames = []
    post.mapping = {}
    post.glyphOrder = order
    try:
        post.compile(f)
    except OverflowError:
        post.formatType = 3.0
    f.sfntVersion = "\000\001\000\000"
    f.save(dst)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    for weight in sys.argv[3:]:
        convert(f"{src}/Inter-{weight}.otf", f"{dst}/Inter-{weight}.ttf")
        print("converted", weight)
