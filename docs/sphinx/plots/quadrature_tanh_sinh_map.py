"""The figure of quadrature/rough-integrands.md: why TanhSinh copes with a singularity at an end. Left, the substitution
x = (1 + tanh(sinh t)) / 2 onto [0, 1], with points at equal steps of t and where they land in x; right, the integrand
1 / sqrt x after it, f(x(t)) dx/dt, on a logarithmic scale (PHASE8_PLAN.md).

The curves come from examples/quadrature_tanh_sinh_map.cpp, a program built and run as a test, which compares them with
the table in its raw string and checks the map against the points TanhSinh evaluates; this script reads that table, so
the figure shows what the test checked. Sphinx runs it through matplotlib's plot directive (conf.py), from this folder.
The colours are those of docs/CLAUDE.md, which pass the dataviz palette checks on the white page and on the dark site.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# The plot directive runs the script from its folder, without __file__; a run by hand has __file__.
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "quadrature_tanh_sinh_map.cpp"
table = source.read_text(encoding="utf-8").split('tanh_sinh_map_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [line.split(",") for line in table.strip().splitlines()[1:]]  # t, x, g, marked
ts = [float(r[0]) for r in rows]
xs = [float(r[1]) for r in rows]
gs = [float(r[2]) for r in rows]
marked = [(float(r[0]), float(r[1]), float(r[2])) for r in rows if r[3] == "1"]

grey = "#757470"  # 4.7:1 on white paper, 3.7:1 on the dark site
blue = "#3987e5"
plt.rcParams.update({
    "font.size": 9,
    "text.color": grey,
    "axes.labelcolor": grey,
    "axes.edgecolor": grey,
    "xtick.color": grey,
    "ytick.color": grey,
    "savefig.transparent": True,
})

fig, (left, right) = plt.subplots(1, 2, figsize=(7.0, 3.0), layout="constrained")

# Left: the map, and each marked point carried across to the vertical x axis, where the lines crowd towards 0 and 1.
for t, x, _ in marked:
    left.plot([ts[0], t], [x, x], color=blue, linewidth=0.6, alpha=0.8)
left.plot(ts, xs, color=grey, linewidth=1.2)
left.plot([t for t, _, _ in marked], [x for _, x, _ in marked], linestyle="none", marker="o", markersize=4,
          color=blue, clip_on=False)
left.set_xlim(ts[0], ts[-1])
left.set_ylim(min(xs), max(xs))
left.set_xlabel("t")
left.set_ylabel("x")

# Right: the integrand in t, finite everywhere and falling off double-exponentially at both ends.
right.plot(ts, gs, color=grey, linewidth=1.2)
right.plot([t for t, _, _ in marked], [g for _, _, g in marked], linestyle="none", marker="o", markersize=4,
           color=blue, clip_on=False)
right.set_yscale("log")
right.set_xlim(ts[0], ts[-1])
right.set_xlabel("t")
right.set_ylabel(r"$f(x(t))\;dx/dt$")

for ax in (left, right):
    ax.grid(True, which="major", color=grey, alpha=0.2, linewidth=0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
