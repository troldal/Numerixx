# data/

Written only by `tools/generate_data.cpp` (CMake target `numerixx_book_data`, also run by `numerixx_book`); every file
carries the generator's header, and `tools/check.py` rejects one that does not. Never edit these files by hand: change
the generator and rerun it. The generated files are committed, so that a diff shows when a number in the book changes.

- `values.tex`: one `\val<Name>` macro per number quoted in the text, read by `macros.tex`.
- `<figure>.csv`: one file per plotted data set; the first line is a `#` comment, the second the column names.

Which figure or chapter uses which file: keep this list current.

| File | Used by |
|---|---|
| `template-step-error.csv`, `\valTemplateDefaultError` | the layout test (remove both with it) |
