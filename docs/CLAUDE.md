# Documentation

The root `CLAUDE.md` describes the documentation's layout and its build targets; the deriv book
(`docs/deriv-book/`) has its own `CLAUDE.md`. This file holds how the user manual (`docs/sphinx/`) makes its figures.

## Plots in the manual

A figure is drawn by matplotlib when the manual is built, on the site as SVG and in the PDF as vector graphics, from
data a tested program checked. A model to follow is the figure of the default method's panels in `quadrature/index.md`,
with `examples/quadrature_panels.cpp` and `plots/quadrature_panels.py` (`PHASE8_PLAN.md`); `PHASE7_PLAN.md` records why
figures are built this way. The manual's rule holds for figures too: every number comes from a program built and run as
a test.

### 1. The data: an example program

`docs/sphinx/examples/<name>.cpp`, built and run as the test `numerixx_docs_example_<name>` (the glob in
`tests/CMakeLists.txt`; no CMake change). It computes the data with the library, prints it as CSV and compares it
with the table it holds in a raw string, as the programs compare the output their pages quote (`examples/printed.hpp`):

```cpp
if (!docs::capture(argv[0])) return 100;
// ... compute, then print the table: a header line, then one row per point
std::println("request,left,right");
constexpr std::string_view panels_table = R"(
request,left,right
...
)";
if (!(docs::printed() == docs::quoted(panels_table))) return 1;
```

- No file I/O beyond `printed.hpp`: the examples also run under node (Emscripten), where a file of the source tree
  cannot be read. The table lives in the program, and the script reads it from there.
- The table must be the same on every toolchain, since CI compares it on all of them: compute with `+ - * /` and
  `std::sqrt` (correctly rounded), not a libm function's last bits (`std::log10`, `std::exp`), and keep data that do
  depend on them out of the table, or plot the bounds the program asserts.
- To fill the table: build the program and run it once; it prints the table to its `.out` file (`printed.hpp`) and to
  stderr, and fails until the raw string matches. Paste the printed lines between `R"(` and `)";`.
- `main` also asserts every claim the page's text, alt text and caption make about the figure, each with its own
  return code: the table pins the numbers, the asserts pin the sentences.

### 2. The drawing: a script

`docs/sphinx/plots/<name>.py`, run by `matplotlib.sphinxext.plot_directive` (`conf.py`). It reads the table from the
program's raw string, splitting on the marker, and draws only that:

```python
here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
source = here.parent / "examples" / "<name>.cpp"
table = source.read_text(encoding="utf-8").split('<name>_table = R"(\n', 1)[1].split(')";', 1)[0]
rows = [line.split(",") for line in table.strip().splitlines()[1:]]
```

- The plot directive runs the script from its own folder, without `__file__`; a run by hand has it.
- No number typed into the script that the table holds (a legend's total is summed from the rows).
- One figure serves the dark site and the white page: `savefig.transparent`, text, axes and ticks in grey `#757470`
  (4.7:1 on white, 3.7:1 on the dark site), font size 9, about 7 x 3.2 inches. Series colours must pass the dataviz
  skill's validator (`validate_palette.js`) with `--mode light` and `--mode dark`: blue `#3987e5` and orange `#d95926`
  do; a third colour is validated with them first. Text in the text grey, never in a series colour; a legend for two
  series or more.
- Preview before building: run the script through `uv run --no-project --managed-python --python 3.14
  --with-requirements docs/sphinx/requirements.txt` and save it on `#fcfcfb` and on `#1a1a19` to see both surfaces.

### 3. The page: the directive

Inside an `{eval-rst}` block, the path from `docs/sphinx/`:

````md
```{eval-rst}
.. plot:: plots/<name>.py
   :alt: Plain text: what the figure shows, with its numbers, true of the table.
   :caption: The caption may use :math:`1/\sqrt{x}` and ``code``.
```
````

- Not MyST's `{plot}`: matplotlib's directive inserts reStructuredText into the document, which MyST's directives do
  not support (`MockStateMachine has not yet implemented attribute 'insert_input'`).
- The path is from `docs/sphinx/` (`plot_basedir` in `conf.py`): under MyST the directive does not resolve one from
  the page's folder, and a leading `/` is the drive's root on Windows.
- The directive takes no `:name:`, so a figure cannot be cited by number: the text says "the figure", in the
  paragraph after it. The PDF keeps each figure where the text puts it (`figure_align` H in `conf.py`); the site's CSS
  (`_static/numerixx.css`) gives it the content column's full width.
- Words for both media (voice rule 7 in `PHASE4_PLAN.md`): "the figure", never "below" for something not adjacent.

### 4. Checks before a commit

- The program on the five toolchains (`ctest --preset <preset> -R <name>`; WSL from a copy, as the root file says).
- `numerixx_docs` and `numerixx_docs_pdf`, with no warning; then look at the figure on the site
  (`build/<preset>/docs/public/html`) and in the PDF (`build/<preset>/docs/public/pdf/numerixx.pdf`).
- `docs-reviewer` on the page, the program and the script.

### Dependencies

matplotlib and its dependencies are pinned in `docs/sphinx/requirements.txt` with the rest, each release at least two
weeks old. To update them: pick the versions, then check that `uv pip compile --python 3.14` on the pinned list
resolves to exactly that list (nothing added, nothing changed) before the docs build.
