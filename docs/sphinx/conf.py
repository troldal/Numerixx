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
}
# Wider still, so that the examples' code blocks need no horizontal scroll bar (_static/numerixx.css).
html_static_path = ["_static"]
html_css_files = ["numerixx.css"]
# The theme loads sphinx-sitemap, which warns without html_baseurl, the site's public URL.
suppress_warnings = ["sitemap.configuration"]
extensions = ["breathe", "myst_parser", "numerixx_breathe"]
exclude_patterns = ["PHASE*.md", "examples", "_ext"]

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
