"""A patch to Breathe 5.1.0 for the Numerixx reference, and the checks that the reference is complete.

Breathe parses a C++ requires clause and a concept's definition from Doxygen's XML but renders neither (breathe#482,
open since 2020). This extension adds both, and after the build fails it unless every class, concept and declaration
in the XML was rendered exactly once, with its clause and its definition. Loaded by conf.py after Breathe; the XML is
Breathe's project breathe_default_project.
"""

import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import breathe
from breathe import _parser
from breathe.renderer import sphinxrenderer
from docutils import nodes

# The patch replaces Breathe internals; another version may have moved them. The checks would catch that, but a wrong
# environment is better reported for what it is.
if breathe.__version__ != "5.1.0":
    raise RuntimeError(f"Breathe {breathe.__version__} found, 5.1.0 needed: run Sphinx with docs/sphinx/requirements.txt")

xml_dir = Path()  # set from the configuration, before any page is read

# The Doxygen IDs of what reached the site: every declaration and class Breathe rendered, and those the two parts of
# the patch completed.
rendered = Counter()
emitted = {"clauses": Counter(), "concepts": Counter()}

# 1. The requires clause after the template parameter list, where Sphinx's C++ domain accepts it since Sphinx 5.2.
#    Relies on SphinxRenderer.create_template_prefix, which every templated declaration goes through.
_template_prefix = sphinxrenderer.SphinxRenderer.create_template_prefix


def _template_prefix_with_clause(self, decl):
    prefix = _template_prefix(self, decl)
    clause = getattr(decl, "requiresclause", None)
    if prefix and clause is not None:
        prefix += " requires " + " ".join("".join(n.astext() for n in self.render(clause)).split())
        emitted["clauses"][decl.id] += 1
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
        emitted["concepts"][node.id] += 1
    return result


sphinxrenderer.SphinxRenderer.node_handlers[_parser.Node_compounddefType] = _compounddef_with_definition

# 3. What Breathe rendered, recorded in SphinxRenderer.render after the directive's filters (such as :desc-only:) have
#    applied, so that the checks see what was rendered, not what a page asked for.
_render = sphinxrenderer.SphinxRenderer.render


def _render_recorded(self, node, context=None, tag=None):
    result = _render(self, node, context, tag)
    if result and isinstance(node, (_parser.Node_memberdefType, _parser.Node_compounddefType)):
        rendered[node.id] += 1
    return result


sphinxrenderer.SphinxRenderer.render = _render_recorded


def _find_xml(app, config):
    global xml_dir
    xml_dir = Path(config.breathe_projects[config.breathe_default_project])


def _check(app, exception):
    """After a successful build, fail it unless every declaration, class and concept in the XML was rendered exactly
    once, and the patch completed every requires clause and concept definition exactly once."""
    if exception is not None:
        return
    declarations, clauses, concepts, names = set(), set(), set(), {}
    for f in xml_dir.glob("*.xml"):
        compound = ET.parse(f).getroot().find("compounddef")
        if compound is None or compound.get("kind") not in ("namespace", "struct", "class", "union", "concept"):
            continue
        if compound.get("kind") != "namespace":
            declarations.add(compound.get("id"))
        if compound.get("kind") == "concept":
            concepts.add(compound.get("id"))
        names[compound.get("id")] = compound.findtext("compoundname")
        for member in compound.iter("memberdef"):
            declarations.add(member.get("id"))
            names[member.get("id")] = f"{compound.findtext('compoundname')}::{member.findtext('name')}"
        clauses |= {e.get("id") for e in compound.iter() if e.find("requiresclause") is not None}

    problems = []
    for what, expected, counted in (("rendered", declarations, rendered),
                                    ("given its requires clause", clauses, emitted["clauses"]),
                                    ("given its definition", concepts, emitted["concepts"])):
        for doxygen_id in sorted(expected):
            if counted[doxygen_id] != 1:
                problems.append(f"{names[doxygen_id]} ({doxygen_id}) {what} {counted[doxygen_id]} times, not once")
    if problems:
        raise RuntimeError("documentation checks failed:\n  " + "\n  ".join(problems))


def setup(app):
    app.connect("config-inited", _find_xml)
    app.connect("build-finished", _check)
    return {"parallel_read_safe": False}
