"""Sphinx settings for the Numerixx documentation (docs/sphinx/PHASE2_PLAN.md).

The CMake target numerixx_docs writes Doxygen's XML and runs Sphinx through uv with requirements.txt, passing the XML
folder in NUMERIXX_DOXYGEN_XML. Breathe renders the reference from that XML, with the patch below.
"""

import os
import re
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import breathe
from breathe import _parser
from breathe.renderer import sphinxrenderer
from docutils import nodes

# The patch below replaces Breathe internals; another version may have moved them. The checks at the end would catch
# that, but a wrong environment is better reported for what it is.
if breathe.__version__ != "5.1.0":
    raise RuntimeError(f"Breathe {breathe.__version__} found, 5.1.0 needed: run Sphinx with docs/sphinx/requirements.txt")

project = "Numerixx"
html_title = "Numerixx"
html_theme = "furo"
extensions = ["breathe", "myst_parser"]
exclude_patterns = ["PHASE*.md", "examples"]

xml_dir = Path(os.environ["NUMERIXX_DOXYGEN_XML"])
breathe_projects = {"numerixx": str(xml_dir)}
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

# The patch (breathe#482, open since 2020): Breathe parses a requires clause and a concept's definition from the XML
# but renders neither. Each part counts what it emits, for the checks at the end.
emitted = {"clauses": 0, "concepts": 0}

# 1. The requires clause after the template parameter list, where Sphinx's C++ domain accepts it since Sphinx 5.2.
#    Relies on SphinxRenderer.create_template_prefix, which every templated declaration goes through.
_template_prefix = sphinxrenderer.SphinxRenderer.create_template_prefix


def _template_prefix_with_clause(self, decl):
    prefix = _template_prefix(self, decl)
    clause = getattr(decl, "requiresclause", None)
    if prefix and clause is not None:
        prefix += " requires " + " ".join("".join(n.astext() for n in self.render(clause)).split())
        emitted["clauses"] += 1
    return prefix


sphinxrenderer.SphinxRenderer.create_template_prefix = _template_prefix_with_clause

# 2. A concept's definition, Doxygen's <initializer>, as a code block before its description; Sphinx's cpp:concept
#    shows only the template head and the name. Relies on the handler table SphinxRenderer.node_handlers, which Breathe
#    fills from its visit_ methods when the class is created, so replacing the method itself would have no effect. The
#    text comes from the XML file itself: Breathe's rendering of the initializer drops the spaces and line breaks
#    between its parts.
_compounddef = sphinxrenderer.SphinxRenderer.node_handlers[_parser.Node_compounddefType]


def _compounddef_with_definition(self, node):
    result = _compounddef(self, node)
    if node.kind == _parser.DoxCompoundKind.concept and node.initializer is not None:
        initializer = ET.parse(xml_dir / f"{node.id}.xml").getroot().find("compounddef/initializer")
        text = "".join(initializer.itertext())
        result.insert(0, nodes.literal_block(text, text, language="cpp"))
        emitted["concepts"] += 1
    return result


sphinxrenderer.SphinxRenderer.node_handlers[_parser.Node_compounddefType] = _compounddef_with_definition


# The checks, after a successful build; a mismatch fails it.
def _scope(name):
    """The enclosing scope of a qualified name, ignoring the :: inside template arguments."""
    depth = 0
    for i in range(len(name) - 1, 0, -1):
        depth += {">": 1, "<": -1}.get(name[i], 0)
        if depth == 0 and name[i - 1 : i + 1] == "::":
            return name[: i - 1]
    return ""


def _directives(srcdir):
    """{directive: {normalised argument}} over every page."""
    found = {}
    for page in Path(srcdir).rglob("*.md"):
        for kind, arg in re.findall(r"^```\{doxygen(\w+)\}\s*(.+)$", page.read_text(encoding="utf-8"), re.M):
            found.setdefault(kind, set()).add(re.sub(r"\s+", "", arg))
    return found


def _check(app, exception):
    if exception is not None:
        return
    problems = []

    # The patch ran once for every clause and concept in the XML: none dropped, none rendered twice.
    clauses = concepts = 0
    for f in xml_dir.glob("*.xml"):
        if f.name != "index.xml":
            root = ET.parse(f).getroot()
            clauses += len(root.findall(".//requiresclause"))
            concepts += len(root.findall(".//compounddef[@kind='concept']"))
    if emitted["clauses"] != clauses:
        problems.append(f"{emitted['clauses']} requires clauses rendered, {clauses} in the XML")
    if emitted["concepts"] != concepts:
        problems.append(f"{emitted['concepts']} concept definitions rendered, {concepts} in the XML")

    # Coverage: every symbol in the XML is reached by a directive. A doxygennamespace renders its namespace's members
    # and classes, not its concepts; a doxygenstruct or doxygenclass renders its members.
    directives = _directives(app.srcdir)
    namespaces = directives.get("namespace", set())
    named = directives.get("struct", set()) | directives.get("class", set())

    def reached(name, kind):
        key = re.sub(r"\s+", "", name)
        if kind == "concept":
            return key in directives.get("concept", set())
        scope = key
        while scope:
            if scope in named or (scope != key and scope in namespaces):
                return True
            scope = _scope(scope)
        return False

    index = ET.parse(xml_dir / "index.xml").getroot()
    for compound in index.findall("compound"):
        kind, name = compound.get("kind"), compound.findtext("name")
        if kind in ("struct", "class", "union", "concept"):
            if not reached(name, kind):
                problems.append(f"{kind} {name} is on no page")
        elif kind == "namespace":
            for member in compound.findall("member"):
                if member.get("kind") != "enumvalue" and name not in namespaces:
                    problems.append(f"{member.get('kind')} {name}::{member.findtext('name')} is on no page")

    if problems:
        raise RuntimeError("documentation checks failed:\n  " + "\n  ".join(problems))


def setup(app):
    app.connect("build-finished", _check)
