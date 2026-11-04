"""The figure of quadrature/tolerances.md: the rounding floor. The default method on the peak 1 / ((x - 0.3)^2 + 1e-4)
over [0, 1] at every relative tolerance from 1e-2 to 1e-16: its estimate and its actual error, relative to the
integral, against the request and the floor (PHASE8_PLAN.md).

The numbers come from examples/quadrature_floor.cpp, a program built and run as a test, which compares them with the
table in its raw string; this script reads that table, so the figure shows what the test checked. Sphinx runs it
through matplotlib's plot directive (conf.py), from this folder. The colours are those of docs/CLAUDE.md, which pass
the dataviz palette checks on the white page and on the dark site; the request and the floor, references rather than
data, are grey.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_floor.cpp"
table = source.read_text(encoding="utf-8").split('floor_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [[float(v) for v in line.split(",")] for line in table.strip().splitlines()[1:]]
requests = [r[0] for r in rows]
estimates = [r[2] for r in rows]
errors = [r[3] for r in rows]
floor = rows[0][4]

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
ax.plot(requests, requests, color=grey, linestyle="--", linewidth=1.0, label="the request, rel_tol")
ax.axhline(floor, color=grey, linestyle=":", linewidth=1.2, label=r"the rounding floor / integral, $50\,\varepsilon$")
ax.plot(requests, estimates, color=blue, marker="o", markersize=4, linewidth=1.2,
        label="the estimate, abs_error / integral")
# The actual error; where it is 0, the result is the double nearest the integral, drawn on the bottom edge, and the
# line breaks there (NaN) rather than passing over it.
bottom = min(e for e in errors if e > 0) / 30
ax.plot(requests, [e if e > 0 else float("nan") for e in errors], color=orange, marker="o", markersize=4,
        linewidth=1.2, label="the actual error / integral")
zeros = [q for q, e in zip(requests, errors) if e == 0]
ax.plot(zeros, [bottom] * len(zeros), linestyle="none", marker="v", markersize=6, color=orange, clip_on=False)
for q in zeros:
    ax.annotate("error 0", xy=(q, bottom), xytext=(6, 2), textcoords="offset points")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(max(requests) * 2, min(requests) / 2)  # tighter requests to the right
ax.set_ylim(bottom, max(estimates) * 30)
ax.set_xlabel("requested relative error, rel_tol")
ax.set_ylabel("relative to the integral")
ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper right", frameon=False)
