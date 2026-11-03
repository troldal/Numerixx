"""The figure of quadrature/choose-a-method.md: where the default and TanhSinh evaluate 1/sqrt(x) on [0, 1], counted
per decade of x (PHASE7_PLAN.md).

The counts come from examples/quadrature_sampling.cpp, a program built and run as a test, which compares them with the
table in its raw string; this script reads that table, so the figure shows what the test checked. Sphinx runs it
through matplotlib's plot directive (conf.py), from this folder.

One figure serves the dark site and the white page: a transparent background, grey axes and text, and two series
colours that pass the dataviz palette checks on both surfaces; the grey reads on both too.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_sampling.cpp"
table = source.read_text(encoding="utf-8").split('sampling_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [line.split(",") for line in table.strip().splitlines()[1:]]

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

# Each decade [10^d, 10^(d + 1)) is drawn at its geometric middle; only the decades a method evaluated in.
totals = [sum(int(r[column]) for r in rows) for column in (1, 2)]
series = [(1, f"default, GaussKronrod<15>: {totals[0]} evaluations", blue),
          (2, f"TanhSinh: {totals[1]} evaluations", orange)]
left = -45  # the axis ends here; TanhSinh's two points below it are named in a note

fig, ax = plt.subplots(figsize=(7.0, 3.2), layout="constrained")
for column, label, colour in series:
    points = [(10.0 ** (int(r[0]) + 0.5), int(r[column])) for r in rows if int(r[column]) > 0 and int(r[0]) >= left]
    ax.plot([x for x, _ in points], [n for _, n in points], linestyle="none", marker="o", markersize=5,
            color=colour, label=label)
beyond = [int(r[0]) for r in rows if int(r[0]) < left and int(r[2]) > 0]
ax.annotate(f"TanhSinh: {len(beyond)} more, near $10^{{{beyond[-1]}}}$ and $10^{{{beyond[0]}}}$",
            xy=(10.0 ** (left + 0.5), 1.0), xytext=(10.0 ** (left + 3), 2.2),
            arrowprops={"arrowstyle": "->", "color": grey, "linewidth": 1})

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10.0 ** left, 1.0)
ax.set_ylim(0.7, 300)
ax.set_xticks([10.0 ** d for d in range(left + 5, 1, 5)])
ax.set_xlabel("x, by decade")
ax.set_ylabel("evaluations in the decade")
ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper left", frameon=False)
