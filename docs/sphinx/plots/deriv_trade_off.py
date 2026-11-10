"""The figure of deriv/index.md: the error of the central difference against its step. The relative error of
Central<2> on sqrt x at 2 for steps from 1 to 1e-13, the default step, and the error of Ridders' method
(PHASE9_PLAN.md).

The numbers come from examples/deriv_trade_off.cpp, a program built and run as a test, which compares them with the
table in its raw string; this script reads that table, so the figure shows what the test checked. Sphinx runs it
through matplotlib's plot directive (conf.py), from this folder. The colours are those of docs/CLAUDE.md, which pass
the dataviz palette checks on the white page and on the dark site; the default step, a reference rather than data, is
grey.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "deriv_trade_off.cpp"
text = source.read_text(encoding="utf-8").split('trade_off_table = R"(\n', 1)[1].split(')";', 1)[0]
curve_part, summary_part = text.strip().split("\n\n")
curve = [[float(v) for v in line.split(",")] for line in curve_part.splitlines()[1:]]
step, _, ridders_error, _ = (float(v) for v in summary_part.splitlines()[1].split(","))
steps = [h for h, _ in curve]
errors = [e for _, e in curve]

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
ax.plot(steps, errors, color=blue, marker="o", markersize=3, linewidth=1.2, label="central difference, Central<2>")
ax.axvline(step, color=grey, linestyle="--", linewidth=1.0, label="the default step")
ax.axhline(ridders_error, color=orange, linewidth=1.2, label="Ridders' method")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(max(steps) * 2, min(steps) / 2)  # smaller steps to the right
ax.set_ylim(ridders_error / 10, max(errors) * 10)
ax.set_xlabel("step h")
ax.set_ylabel("relative error")
ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
ax.minorticks_off()
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
# Right of the default step's line, above the rounding branch, so that the line does not cross the legend.
ax.legend(loc="upper left", bbox_to_anchor=(0.43, 1.0), frameon=False)
