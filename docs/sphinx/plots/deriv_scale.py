"""The figure of deriv/step.md: the valley moves with the function's scale. The relative error of the central
difference against its step on sqrt x at 2 and on the peak x / (x^2 + a^2), a = 1e-4, at 2e-4, for steps from 1 to
1e-14, each with its default step marked (PHASE9_PLAN.md).

The numbers come from examples/deriv_scale.cpp, a program built and run as a test, which compares them with the table
in its raw string; this script reads that table, so the figure shows what the test checked. Sphinx runs it through
matplotlib's plot directive (conf.py), from this folder. The colours are those of docs/CLAUDE.md, which pass the
dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "deriv_scale.cpp"
text = source.read_text(encoding="utf-8").split('scale_table = R"(\n', 1)[1].split(')";', 1)[0]
curve_part, default_part = text.strip().split("\n\n")
curve = [[float(v) for v in line.split(",")] for line in curve_part.splitlines()[1:]]  # h, sqrt, peak
defaults = {name: (float(h), float(e)) for name, h, e in (line.split(",") for line in default_part.splitlines()[1:])}
steps = [r[0] for r in curve]

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
series = (("sqrt", 1, blue, r"$\sqrt{x}$ at 2"), ("peak", 2, orange, r"$x / (x^2 + a^2)$, $a = 10^{-4}$, at $2a$"))
for name, column, colour, label in series:
    ax.plot(steps, [r[column] for r in curve], color=colour, marker="o", markersize=3, linewidth=1.2, label=label)
    h, e = defaults[name]
    ax.plot([h], [e], linestyle="none", marker="o", markersize=9, markerfacecolor="none", markeredgecolor=colour,
            markeredgewidth=1.5)
# One legend entry for both marks of the default step, in the text grey.
handles, labels = ax.get_legend_handles_labels()
handles.append(Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor="none",
                      markeredgecolor=grey, markeredgewidth=1.5))
labels.append("the default step")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(max(steps) * 2, min(steps) / 2)  # smaller steps to the right
errors = [e for r in curve for e in r[1:]]
ax.set_ylim(min(errors) / 10, max(errors) * 10)
ax.set_xlabel("step h")
ax.set_ylabel("relative error")
ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
ax.minorticks_off()
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.45, 1.0), frameon=False)
