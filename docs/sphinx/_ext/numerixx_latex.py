"""The PDF manual's LaTeX writer and its check of code width (PHASE6_PLAN.md, step 3).

latex_show_pagerefs prints "(page n)" after every internal link. A printed reader needs it after a link in the text to
a chapter or a section ("Working with Numerixx (page 3)"), but Sphinx also adds it to every C++ name the reference
and the text link, template parameters included, which crowds the signatures past the margin. Those names keep their
link, without the page; the index gives their pages.

Code inside an API description, which the reference indents and whose lines come from the headers, is set a size
smaller than the guides' code. Two closing quotes in text, as Sphinx makes f'' into, are kept apart, since TeX would
join them into a closing double quote.

The guides' code fits 115 columns (8 pt on A4 with 2 cm margins, conf.py), the reference's 122 (7 pt, at its
indent), both measured in the PDF; a longer line would wrap there, so every build warns, which -W makes an error, on a
code block with a longer line.
"""

import re

from docutils import nodes
from sphinx import addnodes
from sphinx.util import logging
from sphinx.writers.latex import LaTeXTranslator

logger = logging.getLogger(__name__)
_code_columns = 115
_description_columns = 122

# The anchor of a link to a C++ name: Sphinx's C++ domain ("_CPPv4...") or, for a name Doxygen links in the
# reference's text, a Doxygen ID ("structnxx_1_1...", "namespacenxx_1_1..._1a..."). A section's anchor is a slug of
# words and hyphens, and a link to a whole document has none. The target decides, not the link's text, which may be
# code for a section too.
_cpp_anchor = re.compile(r"_CPPv|(namespace|struct|class|union|concept)\w*_1")


def _links_to_text(node):
    anchor = (node.get("refuri") or node.get("refid", "")).rpartition("#")[2]
    return not _cpp_anchor.match(anchor)


def _inside_description(node):
    while node is not None:
        if isinstance(node, addnodes.desc_content):
            return True
        node = node.parent
    return False


class Translator(LaTeXTranslator):
    def visit_reference(self, node):
        # LaTeXTranslator closes an internal link with "}}} (\autopageref...)" when latex_show_pagerefs is on, as the
        # last context entry, which depart_reference writes; drop the page from that entry for a C++ name.
        super().visit_reference(node)
        if self.context and self.context[-1].startswith("}}} (") and not _links_to_text(node):
            self.context[-1] = "}}}"

    def visit_Text(self, node):
        # A group between two closing quotes stops TeX's ligature; Text writes its whole node, so only the body's
        # new part is changed.
        start = len(self.body)
        super().visit_Text(node)
        self.body[start:] = [re.sub("\u2019(?=\u2019)", "\u2019{}", part) for part in self.body[start:]]

    def visit_literal_block(self, node):
        # LaTeXTranslator writes a code block whole and skips its children; the group keeps the smaller size to it.
        if not _inside_description(node):
            return super().visit_literal_block(node)
        self.body.append("\n\\begingroup\\fvset{fontsize=\\scriptsize}")
        try:
            super().visit_literal_block(node)
        except nodes.SkipNode:
            self.body.append("\\endgroup\n")
            raise
        node["numerixx_group"] = True

    def depart_literal_block(self, node):
        super().depart_literal_block(node)
        if node.get("numerixx_group"):
            self.body.append("\\endgroup\n")


def _check_code_width(app, doctree, docname):
    for block in doctree.findall(nodes.literal_block):
        limit = _description_columns if _inside_description(block) else _code_columns
        # LaTeX's code sets a tab to the next multiple of 8 columns.
        widest = max((len(line.expandtabs(8)) for line in block.astext().split("\n")), default=0)
        if widest > limit:
            logger.warning(f"a code line of {widest} columns wraps in the PDF, which fits {limit} here",
                           location=block)


def setup(app):
    app.connect("doctree-resolved", _check_code_width)
    app.set_translator("latex", Translator)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
