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
html_theme = "furo"
extensions = ["breathe", "myst_parser", "numerixx_breathe"]
exclude_patterns = ["PHASE*.md", "examples", "_ext"]

breathe_projects = {"numerixx": os.environ["NUMERIXX_DOXYGEN_XML"]}
breathe_default_project = "numerixx"
breathe_show_include = False
primary_domain = "cpp"
highlight_language = "cpp"

# Every C++ name must resolve, except those with nothing to link to: the namespace nxx, the detail concepts the
# requires clauses name (Doxygen's HTML shows them as text too), the standard library, and the formatters' base
# class, which Doxygen writes without its std::.
nitpicky = True
nitpick_ignore_regex = [(r"cpp:identifier", r"nxx|(nxx::)?detail(::.*)?|std::.*|formatter<string_view>|string_view")]

# Doxygen 1.18 writes a <conceptparts> element that Breathe 5.1 does not know; it reads the rest of the concept.
warnings.filterwarnings("ignore", message=r'.*unexpected element "conceptparts"')
