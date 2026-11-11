"""Prose words per manual page, by the budget rule of docs/sphinx/CLAUDE.md, so that every agent counts the same way.

Fenced blocks (code, output, eval-rst, toctree, raw, tabularcolumns) are excluded; notes, headings and table cells
count; an inline code span or a formula is one word; a link's text counts, its target does not.

Usage, from the project root: python docs/sphinx/count_words.py docs/sphinx/roots [more pages or folders]
"""
import glob
import os
import re
import sys

ADMONITIONS = ("{note}", "{warning}", "{tip}", "{important}", "{attention}")


def count(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    out = []
    skip = False
    in_admonition = False
    for ln in lines:
        m = re.match(r"^\s*(`{3,})\s*(\{?[\w-]*\}?)", ln)
        if m:
            kind = m.group(2)
            if skip:
                skip = False
                continue
            if in_admonition:
                in_admonition = False
                continue
            if kind in ADMONITIONS:
                in_admonition = True
                continue
            skip = True
            continue
        if skip or ln.startswith("%"):
            continue
        out.append(ln)
    t = "\n".join(out)
    t = re.sub(r"\$\$.*?\$\$", " FORMULA ", t, flags=re.S)
    t = re.sub(r"\$[^$\n]+\$", " FORMULA ", t)
    t = re.sub(r"`[^`\n]+`", " CODE ", t)
    t = re.sub(r"\{cpp:[a-z]+\}", "", t)
    t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", t)
    t = re.sub(r"^\s*#+\s*", "", t, flags=re.M)
    t = re.sub(r"^\s*\|?[-|: ]+\|?\s*$", "", t, flags=re.M)
    t = t.replace("|", " ")
    return len([w for w in re.split(r"\s+", t) if re.search(r"[A-Za-z0-9]", w)])


def main(args):
    for arg in args:
        pages = sorted(glob.glob(os.path.join(arg, "*.md"))) if os.path.isdir(arg) else [arg]
        total = 0
        for page in pages:
            n = count(page)
            total += n
            print(f"{page:50} {n}")
        if len(pages) > 1:
            print(f"{arg} total {total}\n")


if __name__ == "__main__":
    main(sys.argv[1:] or ["docs/sphinx/deriv", "docs/sphinx/roots", "docs/sphinx/quadrature"])
