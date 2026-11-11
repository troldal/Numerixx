"""The first figure of roots/index.md: how the methods close in. The width of the bracket after each evaluation of fn,
for bisection and Brent's method on x^2 - 2 over [1, 2] (PHASE10_PLAN.md).

The widths come from examples/roots_closing_in.cpp, a program built and run as a test, which compares them with the
table in its raw string; this script reads that table, so the figure shows what the test checked. Sphinx runs it
through matplotlib's plot directive (conf.py), from this folder. The colours are those of docs/CLAUDE.md, which pass
the dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "roots_closing_in.cpp"
table = source.read_text(encoding="utf-8").split('closing_in_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [line.split(",") for line in table.strip().splitlines()[1:]]
series = {}
for method, evaluations, width in rows:
    series.setdefault(method, []).append((int(evaluations), float(width)))

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

fig, ax = plt.subplots(figsize=(7.0, 3.2), layout="constrained")
for method, colour, name in (("bisection", blue, "bisection"), ("brent", orange, "Brent's method")):
    points = series[method]
    ax.plot([n for n, _ in points], [w for _, w in points], color=colour, marker="o", markersize=3, linewidth=1.2,
            label=f"{name}: {points[-1][0]} evaluations")

widths = [w for points in series.values() for _, w in points]
ax.set_yscale("log")
ax.set_ylim(min(widths) / 4, max(widths) * 4)
ax.set_xlim(0, max(n for points in series.values() for n, _ in points) + 1)
ax.set_xlabel("evaluations of fn")
ax.set_ylabel("width of the bracket")
ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
ax.minorticks_off()
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper right", frameon=False)
