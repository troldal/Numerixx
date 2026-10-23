#!/usr/bin/env bash
# Compiles one chapter on its own, with the book's preamble and macros, for an author working on that chapter.
# Usage: bash docs/deriv-book/tools/check_chapter.sh chapters/<file>.tex
# Prints COMPILED OK when the chapter compiles. References to other chapters and citations are undefined on their
# own and are not errors here; LaTeX errors and overfull boxes wider than 1 pt are.
set -euo pipefail
book="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
chapter="${1:?usage: check_chapter.sh chapters/<file>.tex}"
name="$(basename "$chapter" .tex)"
out="$book/../../build/deriv-book-chapters/$name"
mkdir -p "$out"
cd "$book"
cat > "$out/wrap.tex" <<EOF
\\documentclass[10pt,a4paper,oneside]{report}
\\input{preamble}
\\input{macros}
\\begin{document}
\\input{${chapter%.tex}}
\\end{document}
EOF
if pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$out" "$out/wrap.tex" > /dev/null 2>&1 &&
   pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$out" "$out/wrap.tex" > /dev/null 2>&1; then
    over=$(grep -E 'Overfull \\[hv]box \(([1-9][0-9]*|[1-9])\.?[0-9]*pt' "$out/wrap.log" || true)
    if [ -n "$over" ]; then echo "$over"; echo "overfull boxes: fix them"; exit 1; fi
    echo "COMPILED OK ($out/wrap.pdf)"
else
    grep -E '^!' -A4 "$out/wrap.log" | head -40
    echo "FAILED: see $out/wrap.log"
    exit 1
fi
