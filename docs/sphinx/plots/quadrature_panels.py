"""The figure of quadrature/index.md: the panels the default method ends with on the peak 1 / ((x - 0.3)^2 + 1e-4) over
[0, 1], at the default request and at a relative tolerance of 1e-12, below the curve (PHASE8_PLAN.md).

The panels and the curve come from examples/quadrature_panels.cpp, a program built and run as a test, which compares
them with the tables in its raw string; this script reads those tables, so the figure shows what the test checked.
Sphinx runs it through matplotlib's plot directive (conf.py), from this folder. The colours are those of
docs/CLAUDE.md, which pass the dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_panels.cpp"
text = source.read_text(encoding="utf-8").split('panels_table = R"(\n', 1)[1].split(')";', 1)[0]
panels_part, curve_part = text.strip().split("\n\n")
panels = [line.split(",") for line in panels_part.splitlines()[1:]]
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

fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(7.0, 3.2), layout="constrained",
                                  gridspec_kw={"height_ratios": [2.2, 1]})

top.plot([x for x, _ in curve], [f for _, f in curve], color=grey, linewidth=1.2)
top.set_ylabel("f(x)")
top.set_ylim(bottom=0)
top.grid(True, color=grey, alpha=0.2, linewidth=0.5)

# One row per request, a tick at each panel's end: where the ticks crowd, the panels are narrow.
rows = []
requests = list(dict.fromkeys(p[0] for p in panels))  # "default", then the other, in the table's order
for request, colour in zip(requests, (blue, orange)):
    label = "default request" if request == "default" else f"rel_tol = {request}"
    ends = sorted({float(p[1]) for p in panels if p[0] == request} | {float(p[2]) for p in panels if p[0] == request})
    count = sum(1 for p in panels if p[0] == request)
    rows.append((ends, colour, f"{label}: {count} panels"))
bottom.eventplot([ends for ends, _, _ in rows], colors=[c for _, c, _ in rows], lineoffsets=[1, 0],
                 linelengths=0.7, linewidths=1.0)
bottom.set_yticks([1, 0], [label for _, _, label in rows])
bottom.set_ylim(-0.6, 1.6)
bottom.tick_params(axis="y", length=0)
bottom.set_xlim(min(float(p[1]) for p in panels), max(float(p[2]) for p in panels))
bottom.set_xlabel("x")

for ax in (top, bottom):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
bottom.spines["left"].set_visible(False)
