"""The second figure of roots/index.md: a function that cancels near its root. The computed values of (x - 1)^3
written out, x^3 - 3x^2 + 3x - 1 by Horner's rule, and as a product, at 201 points across [1 - 1e-5, 1 + 1e-5]
(PHASE10_PLAN.md).

The values come from examples/roots_cancellation.cpp, a program built and run as a test, which compares them with the
table in its raw string; this script reads that table, so the figure shows what the test checked. Sphinx runs it
through matplotlib's plot directive (conf.py), from this folder. The colours are those of docs/CLAUDE.md, which pass
the dataviz palette checks on the white page and on the dark site; the zero line, a reference rather than data, is
grey.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "roots_cancellation.cpp"
table = source.read_text(encoding="utf-8").split('cancellation_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [[float(v) for v in line.split(",")] for line in table.strip().splitlines()[1:]]
# The axes' units: x - 1 in millionths, the values in units of 1e-16.
offsets = [offset * 1e6 for offset, _, _ in rows]
written_out = [y * 1e16 for _, y, _ in rows]
product = [y * 1e16 for _, _, y in rows]

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
ax.axhline(0.0, color=grey, linewidth=0.8)
ax.plot(offsets, written_out, color=orange, marker="o", markersize=2.5, linestyle="none",
        label=r"written out, $x^3 - 3x^2 + 3x - 1$")
ax.plot(offsets, product, color=blue, linewidth=1.4, label=r"product, $(x - 1)^3$")

ax.set_xlim(min(offsets), max(offsets))
ax.set_xlabel(r"$x - 1$, in units of $10^{-6}$")
ax.set_ylabel(r"computed value, in units of $10^{-16}$")
ax.grid(True, color=grey, alpha=0.2, linewidth=0.5)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper left", frameon=False)
