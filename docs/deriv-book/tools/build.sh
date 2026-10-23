#!/usr/bin/env bash
# Builds the book's PDF and checks it. Usage, from anywhere:  bash docs/deriv-book/tools/build.sh [outdir]
# The default output directory is build/deriv-book in the repository. The CMake target numerixx_book does the same
# after regenerating data/ (cmake --build --preset <preset> --target numerixx_book); this script does not regenerate
# data/, so run that target (or numerixx_book_data) after a change to the library.
set -euo pipefail
book="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out="${1:-$book/../../build/deriv-book}"
mkdir -p "$out"
cd "$book"
latexmk -interaction=nonstopmode -halt-on-error -outdir="$out" main.tex > "$out/latexmk.txt" 2>&1 || {
    grep -E '^!' -A4 "$out/main.log" | head -40
    echo "build failed: see $out/main.log"
    exit 1
}
python tools/check.py "$out/main.log"
echo "PDF: $out/main.pdf"
