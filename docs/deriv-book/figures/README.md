# figures/

One file per figure, `<chapter-id>-<name>.tex`, holding only the `tikzpicture` (the float and its caption stay in the
chapter, so that the author places them). Every plot uses the `eqtufte` style (thin grey axes, no frame or grid,
curves labelled directly rather than with a legend) and reads its data from `data/*.csv` with
`\addplot table[col sep=comma, comment chars={\#}]`: no hand-typed coordinates for computed results. Diagrams that
are not data plots (a stencil's nodes, the Ridders tableau) are TikZ in the same visual language: thin lines, the
`accent` colour for the one thing to notice, `muted` for the rest.
