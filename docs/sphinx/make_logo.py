# /// script
# requires-python = ">=3.12"
# dependencies = ["fonttools==4.55.0"]
# ///
"""Writes the Numerixx logo files into _static/: run `uv run docs/sphinx/make_logo.py` from the project root.

The logo is "nxx" in IBM Plex Mono SemiBold, four square dots in Lumina's accent, a rule and "numerixx" in Plex Mono
Regular, drawn from the font's outlines (SIL Open Font License) so that no file needs the font. Outputs:
numerixx-logo-{light,dark}.svg (site header, README), numerixx-mark-{light,dark}.svg (the mark alone), favicon.svg,
and numerixx-logo.pdf (the PDF manual's title page; pdflatex includes PDF, not SVG). The fonts are downloaded from the
pinned npm release and checked against their SHA-256.
"""

import hashlib
import io
import sys
import urllib.request
from pathlib import Path

from fontTools.pens.basePen import BasePen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

FONT_URL = "https://cdn.jsdelivr.net/npm/@ibm/plex-mono@2.5.0/fonts/complete/woff/IBMPlexMono-{}.woff"
FONT_SHA256 = {
    "Regular": "3a0d6e1587e8dc784be5a75e1c55b1a5533919b04072a29c750bb29cfe4a4b74",
    "SemiBold": "7fb3daa74f03d2c9f4e092b9f1b8f2cdb301eafeba0d0af375d856d44dce3327",
}
ACCENT = "#10b981"  # Lumina's --lumina-accent
THEMES = {
    "light": {"fg": "#18181b", "rule": "#d4d4d8", "accent": ACCENT},
    "dark": {"fg": "#e6ebe8", "rule": "#3a443f", "accent": ACCENT},
}


def num(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s == "-0" else s


class Outline(BasePen):
    """A drawing recorded as closed contours of lines and cubic curves, y down, in px.

    BasePen turns the font's quadratic curves into cubics, which both SVG and PDF draw.
    """

    def __init__(self, glyph_set=None):
        super().__init__(glyph_set)
        self.contours = []

    def _moveTo(self, p):
        self.contours.append([("M", p)])

    def _lineTo(self, p):
        self.contours[-1].append(("L", p))

    def _curveToOne(self, a, b, c):
        self.contours[-1].append(("C", a, b, c))

    def _closePath(self):
        pass

    def rect(self, x, y, w, h):
        self.contours.append([("M", (x, y)), ("L", (x + w, y)), ("L", (x + w, y + h)), ("L", (x, y + h))])


def bounds(outlines):
    pen = BoundsPen(None)
    for outline in outlines:
        for contour in outline.contours:
            for op, *points in contour:
                {"M": pen.moveTo, "L": pen.lineTo, "C": pen.curveTo}[op](*points)
            pen.closePath()
    return pen.bounds


def load_fonts():
    fonts = {}
    for weight, digest in FONT_SHA256.items():
        with urllib.request.urlopen(FONT_URL.format(weight)) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != digest:
            sys.exit(f"IBMPlexMono-{weight}.woff: unexpected SHA-256, so not the pinned release")
        fonts[weight] = TTFont(io.BytesIO(data))
    return fonts


def text(font, s, size, x0, baseline, tracking_em):
    """The outline of s set at (x0, baseline) and the end of its advance."""
    glyphs, cmap = font.getGlyphSet(), font.getBestCmap()
    scale = size / font["head"].unitsPerEm
    outline, x = Outline(glyphs), x0
    for ch in s:
        glyph = glyphs[cmap[ord(ch)]]
        glyph.draw(TransformPen(outline, (scale, 0, 0, -scale, x, baseline)))
        x += glyph.width * scale + tracking_em * size
    return outline, x


def x_height(font, size):
    # The flat top of "x", which has no overshoot, so the dots and the centre line match the letters exactly.
    pen = BoundsPen(font.getGlyphSet())
    font.getGlyphSet()[font.getBestCmap()[ord("x")]].draw(pen)
    return pen.bounds[3] * size / font["head"].unitsPerEm


def lockup(fonts):
    """The mark's parts and the rest of the lockup, as (role, outline).

    Laid out by ink, not by text boxes: both words' x-height bands and the rule share the centre line y = 0, and the
    rule is the same distance from the last dot and from the first stroke of "numerixx". Sizes as on the approved
    design: "nxx" 56 px with tracking -0.02em; dots 0.16em square, centred 0.2em and 0.6em after a 0.04em gap, one
    pair from the x-height down and one up from the baseline; rule 2 x 48 px, 20 px from each side; "numerixx" 32 px
    with tracking -0.01em.
    """
    semi, reg = fonts["SemiBold"], fonts["Regular"]
    size, name_size, gap = 56, 32, 20
    xh = x_height(semi, size)
    baseline = xh / 2
    nxx, x_end = text(semi, "nxx", size, 0, baseline, -0.02)
    dots, u, side = Outline(), size / 100, 0.16 * size
    for centre in (20, 60):
        for top in (baseline - xh, baseline - side):
            dots.rect(x_end + 0.04 * size + (centre - 8) * u, top, side, side)
    rule_x = bounds([dots])[2] + gap
    rule = Outline()
    rule.rect(rule_x, -24, 2, 48)
    name_baseline = x_height(reg, name_size) / 2
    probe, _ = text(reg, "numerixx", name_size, 0, name_baseline, -0.01)
    name, _ = text(reg, "numerixx", name_size, rule_x + 2 + gap - bounds([probe])[0], name_baseline, -0.01)
    return [("fg", nxx), ("accent", dots)], [("rule", rule), ("fg", name)]


def svg_document(shapes, colours):
    x0, y0, x1, y1 = bounds([outline for _, outline in shapes])
    w, h = x1 - x0, y1 - y0

    def d(outline):
        return "".join(op + " ".join(f"{num(x)} {num(y)}" for x, y in points) for contour in outline.contours
                       for op, *points in [*contour, ("Z",)])

    paths = "".join(f'<path fill="{colours[role]}" d="{d(outline)}"/>\n' for role, outline in shapes)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{num(x0)} {num(y0)} {num(w)} {num(h)}" '
            f'width="{num(w)}" height="{num(h)}" role="img" aria-label="Numerixx">\n<title>Numerixx</title>\n'
            f'{paths}</svg>\n')


def pdf_document(shapes, colours):
    """A one-page PDF of the filled outlines, 1 px = 1 pt, written by hand: no timestamps, so a rerun is identical."""
    x0, y0, x1, y1 = bounds([outline for _, outline in shapes])
    operators = {"M": "m", "L": "l", "C": "c"}
    content = []
    for role, outline in shapes:
        rgb = (int(colours[role][i:i + 2], 16) / 255 for i in (1, 3, 5))
        content.append(" ".join(f"{c:.3f}" for c in rgb) + " rg")
        for contour in outline.contours:
            content += [" ".join(f"{num(x - x0)} {num(y1 - y)}" for x, y in points) + " " + operators[op]
                        for op, *points in contour]
            content.append("h")
        content.append("f")
    stream = "\n".join(content).encode("ascii")
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {num(x1 - x0)} {num(y1 - y0)}] /Resources << >> "
        f"/Contents 4 0 R >>".encode("ascii"),
        b"<< /Length %d >>\nstream\n%s\nendstream" % (len(stream), stream),
    ]
    pdf, offsets = bytearray(b"%PDF-1.4\n"), []
    for number, body in enumerate(objects, 1):
        offsets.append(len(pdf))
        pdf += b"%d 0 obj\n%s\nendobj\n" % (number, body)
    xref = len(pdf)
    pdf += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objects) + 1)
    pdf += b"".join(b"%010d 00000 n \n" % offset for offset in offsets)
    pdf += b"trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objects) + 1, xref)
    return bytes(pdf)


def main():
    out = Path(__file__).resolve().parent / "_static"
    mark, rest = lockup(load_fonts())
    for theme, colours in THEMES.items():
        (out / f"numerixx-logo-{theme}.svg").write_text(svg_document(mark + rest, colours), "ascii")
        (out / f"numerixx-mark-{theme}.svg").write_text(svg_document(mark, colours), "ascii")
    # Printed on white paper, so the light colours.
    (out / "numerixx-logo.pdf").write_bytes(pdf_document(mark + rest, THEMES["light"]))
    # The favicon: the four dots on Lumina's dark green (--lumina-accent-light in dark mode), in a 64-unit tile.
    squares = "".join(f"M{x} {y}h15v15h-15z" for y in (14, 35) for x in (14, 35))
    (out / "favicon.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="64" height="64" role="img" '
        'aria-label="Numerixx">\n<title>Numerixx</title>\n<rect width="64" height="64" rx="12" fill="#1c3329"/>\n'
        f'<path fill="{ACCENT}" d="{squares}"/>\n</svg>\n', "ascii")


if __name__ == "__main__":
    main()
