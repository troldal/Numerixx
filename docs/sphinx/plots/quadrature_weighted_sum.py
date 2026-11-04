"""The figure of quadrature/index.md: a quadrature rule as a weighted sum. GaussLegendre<5> on 4 / (1 + x^2) over
[0, 1], drawn as five bars side by side, the i-th of width w_i and height f(x_i), against the area under the curve
(PHASE8_PLAN.md).

The bars and the curve come from examples/quadrature_weighted_sum.cpp, a program built and run as a test, which
compares them with the table in its raw string; this script reads that table, so the figure shows what the test
checked. Sphinx runs it through matplotlib's plot directive (conf.py), from this folder. The colours are those of
docs/CLAUDE.md, which pass the dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_weighted_sum.cpp"
text = source.read_text(encoding="utf-8").split('weighted_sum_table = R"(\n', 1)[1].split(')";', 1)[0]
bars_part, curve_part = text.strip().split("\n\n")
bars = [[float(v) for v in line.split(",")] for line in bars_part.splitlines()[1:]]  # node, weight, left, right, f
curve = [[float(v) for v in line.split(",")] for line in curve_part.splitlines()[1:]]

grey = "#757470"  # 4.7:1 on white paper, 3.7:1 on the dark site
blue = "#3987e5"
orange = "#d95926"
plt.rcParams.update({
    "font.size": 9,
    "text.color": grey,
    "axes.labelcolor": grey,
    "axes.edgecolor": grey,
    "xtick.color": grey,
    "ytick.color": grey,
    "savefig.transparent": True,
})

# A little lower than the manual's 3.2 inches, so that the figure fits on the chapter's first page of the PDF.
fig, ax = plt.subplots(figsize=(7.0, 2.6), layout="constrained")
xs = [x for x, _ in curve]
fs = [f for _, f in curve]
area = ax.fill_between(xs, fs, color=grey, alpha=0.4, linewidth=0, label="area under $f(x)$: the integral")
rects = ax.bar([b[2] for b in bars], [b[4] for b in bars], width=[b[1] for b in bars], align="edge", color=blue,
              alpha=0.3, edgecolor=blue, linewidth=1.0,
              label=r"bars of width $w_i$ and height $f(x_i)$: the rule's sum")
ax.plot(xs, fs, color=grey, linewidth=1.2)
(dots,) = ax.plot([b[0] for b in bars], [b[4] for b in bars], linestyle="none", marker="o", markersize=5,
                  color=orange, label=r"nodes $x_i$, where $f$ is evaluated")

ax.set_xlim(curve[0][0], curve[-1][0])
ax.set_ylim(0, 1.25 * max(fs))
ax.set_xlabel("x")
ax.set_ylabel("f(x)")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(handles=[area, rects, dots], loc="upper right", frameon=False)
