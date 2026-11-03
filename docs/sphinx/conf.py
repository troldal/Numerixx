"""Sphinx settings for the Numerixx documentation (docs/sphinx/PHASE2_PLAN.md).

The CMake target numerixx_docs writes Doxygen's XML and runs Sphinx through uv with requirements.txt, passing the XML
folder in NUMERIXX_DOXYGEN_XML. Breathe renders the reference from that XML, with the patch and the checks of the
local extension _ext/numerixx_breathe.py.
"""

import os
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

project = "Numerixx"
html_title = "Numerixx"
# Lumina (sphinx-lumina-theme): dark by default; Sphinx's own search, since the theme's default, Pagefind, needs a
# tool of its own; no SEO metadata or sitemap until the site has a public URL.
html_theme = "lumina"
html_theme_options = {
    "dark_mode_default": "dark",
    "search_backend": "sphinx",
    "disable_seo": "true",
    "wide_layout": "always",
    # The logo already spells the name, so the header shows it alone; Lumina swaps the files with its theme toggle.
    "light_logo": "numerixx-logo-light.svg",
    "dark_logo": "numerixx-logo-dark.svg",
    "logo_only": "true",
}
html_favicon = "_static/favicon.svg"
# Wider still, so that the examples' code blocks need no horizontal scroll bar (_static/numerixx.css).
html_static_path = ["_static"]
html_css_files = ["numerixx.css"]
# The theme loads sphinx-sitemap, which warns without html_baseurl, the site's public URL.
suppress_warnings = ["sitemap.configuration"]
extensions = ["breathe", "myst_parser", "numerixx_breathe", "numerixx_latex"]
exclude_patterns = ["PHASE*.md", "examples", "_ext"]

# One source, two manuals (PHASE6_PLAN.md). The PDF build passes -t pdf: its root is manual.md, whose toctree makes
# each page a chapter, with introduction.md (index.md's text without its toctrees) first. index.md cannot be the PDF's
# root, since the LaTeX builder would put every page inside its last section, and each root is excluded from the
# other build.
if tags.has("pdf"):
    root_doc = "manual"
    exclude_patterns += ["index.md"]
else:
    exclude_patterns += ["manual.md", "introduction.md"]
latex_documents = [("manual", "numerixx.tex", "Numerixx", "Kenneth Troldal Balslev", "manual", True)]
# A4, no blank pages before chapters, and code small enough that the examples' 120-column lines fit the text width,
# as on the site. A link to a chapter or a section prints its page, since a printed reader cannot follow it
# (_ext/numerixx_latex.py keeps the page off the C++ names). The CMake target passes the project's version as
# release, for the title page.
latex_elements = {
    "papersize": "a4paper",
    "extraclassoptions": "openany,oneside",
    "sphinxsetup": "hmargin=2cm, vmargin=2.5cm",
    "fvset": r"\fvset{fontsize=\footnotesize}",
    # The default chapter style (fncychap's Bjarne) sets titles in capitals, which would misspell nxx::real.
    "fncychap": "",
    # The index's two columns are narrower than names such as std::formatter<nxx::quadrature::Failure<T>>::format.
    "printindex": r"\footnotesize\raggedright\printindex",
}
latex_show_pagerefs = True
# Formulas in $...$ and $$...$$, rendered by MathJax from its CDN (a site read offline shows the TeX); anchors on
# headings down to level 3, so that links such as basics.md#results resolve.
myst_enable_extensions = ["dollarmath"]
myst_heading_anchors = 3

breathe_projects = {"numerixx": os.environ["NUMERIXX_DOXYGEN_XML"]}
breathe_default_project = "numerixx"
breathe_show_include = False
primary_domain = "cpp"
highlight_language = "cpp"

# Every C++ name must resolve, except those with nothing to link to: the namespace nxx, the standard library, the
# detail concepts the requires clauses name (Doxygen's HTML shows them as text too), and the formatters' base class,
# which Doxygen writes without its std::.
nitpicky = True
nitpick_ignore_regex = [(r"cpp:identifier", r"nxx|std(::.*)?|(nxx::)?detail(::.*)?|formatter<string_view>|string_view")]

# Doxygen 1.18 writes a <conceptparts> element that Breathe 5.1 does not know; it reads the rest of the concept.
warnings.filterwarnings("ignore", message=r'.*unexpected element "conceptparts"')
