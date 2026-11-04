"""The figure of quadrature/infinite-intervals.md: where the sinh-sinh map looks. The points TanhSinh places on the
whole line in the first five levels, which every run takes, and the peak 1 / (1 + (x - x0)^2) at x0 = 0, 10 and 40
(PHASE8_PLAN.md).

The points and the peaks come from examples/quadrature_line.cpp, a program built and run as a test, which compares
them with the tables in its raw string; this script reads those tables, so the figure shows what the test checked.
Sphinx runs it through matplotlib's plot directive (conf.py), from this folder. The colours are those of
docs/CLAUDE.md, which pass the dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_line.cpp"
text = source.read_text(encoding="utf-8").split('line_table = R"(\n', 1)[1].split(')";', 1)[0]
points_part, peaks_part = text.strip().split("\n\n")
points = [float(v) for v in points_part.splitlines()[1:]]
peaks = {}
for line in peaks_part.splitlines()[1:]:
    x0, x, f = (float(v) for v in line.split(","))
    peaks.setdefault(x0, []).append((x, f))

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

fig, ax = plt.subplots(figsize=(7.0, 3.0), layout="constrained")
# Each point a line up to just above the peaks, leaving room for the labels and the legend.
reach = 1.15
ax.vlines(points, 0, reach, color=blue, linewidth=0.8, label="the points of the first five levels")
for i, (x0, curve) in enumerate(peaks.items()):
    ax.plot([x for x, _ in curve], [f for _, f in curve], color=orange, linewidth=1.5,
            label=r"the peak $1/(1 + (x - x_0)^2)$" if i == 0 else None)
    ax.annotate(f"$x_0 = {x0:g}$", xy=(x0, reach), xytext=(0, 3), textcoords="offset points", ha="center")

ax.set_xlim(min(points) * 1.03, max(points) * 1.03)  # the outermost lines unclipped
ax.set_ylim(0, 1.55)
ax.set_xlabel("x")
ax.set_yticks([])
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper left", frameon=False, ncols=2)
