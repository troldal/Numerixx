"""The figure of quadrature/high-precision.md: evaluations against digits. The default method and TanhSinh on
1 / (1 + x^2) over [0, 1] in cpp_bin_float_50, at requests of 4 to 48 digits (PHASE8_PLAN.md).

The counts come from examples/quadrature_digits.cpp, a program built and run as a test, which compares them with the
table in its raw string and checks each result against pi / 4; this script reads that table, so the figure shows what
the test checked. Sphinx runs it through matplotlib's plot directive (conf.py), from this folder. The colours are those
of docs/CLAUDE.md, which pass the dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_digits.cpp"
table = source.read_text(encoding="utf-8").split('digits_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [[int(v) for v in line.split(",")] for line in table.strip().splitlines()[1:]]  # digits, default, tanh_sinh
digits = [r[0] for r in rows]

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
ax.plot(digits, [r[1] for r in rows], color=blue, marker="o", markersize=3.5, linewidth=1.2,
        label="default, GaussKronrod<15>")
ax.plot(digits, [r[2] for r in rows], color=orange, marker="o", markersize=3.5, linewidth=1.2, label="TanhSinh")

ax.set_yscale("log")
# Plain ticks at 1, 2 and 5 times the powers of ten the counts span, so that a count can be read off.
ticks = [m * 10 ** e for e in range(2, 6) for m in (1, 2, 5)]
counts = [c for r in rows for c in r[1:]]
# The axis starts at the tick below the lowest count, so that the plateaus of the first evaluations can be read.
low = max(t for t in ticks if t <= min(counts) / 1.3)
high = max(counts) * 1.3
ax.set_ylim(low, high)
ax.set_yticks([t for t in ticks if low <= t <= high], [f"{t:,}" for t in ticks if low <= t <= high])
ax.minorticks_off()
ax.set_xlim(digits[0] - 1, digits[-1] + 1)
ax.set_xlabel("digits requested, d: rel_tol = $10^{-d}$")
ax.set_ylabel("evaluations")
ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper left", frameon=False)
