"""Check that a language pass (PLAN.md section 12.9) changed only the wording of a chapter: compare each chapter
file at a git revision with the working file.

Strict, must be identical and in the same order (exit status 1 otherwise):
- the headings (\\chapter, \\section, \\subsection, \\paragraph) with their titles;
- the environments, \\begin and \\end in order (figures, tables, sidebars, lists, displayed mathematics, ...);
- the text of every displayed equation, punctuation included;
- the \\val macros, \\snippet and \\snippetheader, \\cite keys, \\label (and sidebar labels), \\cref and \\Cref
  targets, \\input, \\tablebody and \\tablerows.
Reported for review (a language pass may change them, the editor reads every difference):
- the number of margin notes (\\mnote);
- inline mathematics and \\code{}, as multisets: what was dropped and what was added;
- plain numbers typed in the prose, as a multiset (layout arguments such as \\needspace are ignored).
Comments are ignored throughout.

Usage, from docs/deriv-book:  python tools/check_language_pass.py <revision> chapters/01-why.tex ...
"""
import collections
import difflib
import pathlib
import re
import subprocess
import sys

BOOK = pathlib.Path(__file__).resolve().parent.parent
REPO = BOOK.parent.parent
DISPLAY_ENVS = ('equation', 'equation*', 'align', 'align*', 'gather', 'gather*', 'multline', 'multline*')
DISPLAY = re.compile(r'\\begin\{(' + '|'.join(re.escape(e) for e in DISPLAY_ENVS) + r')\}(.*?)\\end\{\1\}', re.S)
BRACKETS = re.compile(r'\\\[(.*?)\\\]', re.S)
INLINE = re.compile(r'(?<!\\)\$(.+?)(?<!\\)\$', re.S)
CODE = re.compile(r'\\code\{((?:[^{}]|\{[^{}]*\})*)\}')


def strip_comments(text: str) -> str:
    return '\n'.join(re.sub(r'(?<!\\)%.*', '', line) for line in text.splitlines())


def norm(s: str) -> str:
    return re.sub(r'\s+', ' ', s).strip()


def without_display(text: str) -> str:
    return BRACKETS.sub(' ', DISPLAY.sub(' ', text))


def strict_sequences(text: str) -> dict:
    heading = r'\\(chapter|section|subsection|paragraph)\*?(?:\[[^\]]*\])?\{((?:[^{}]|\{[^{}]*\})*)\}'
    return {
        'headings': [f'{m.group(1)}: {norm(m.group(2))}' for m in re.finditer(heading, text)],
        'environments': [m.group(1) + ' ' + m.group(2) for m in re.finditer(r'\\(begin|end)\{([^}]*)\}', text)],
        'display maths': [f'{m.group(1)}: {norm(m.group(2))}' for m in DISPLAY.finditer(text)]
                         + [f'\\[: {norm(m.group(1))}' for m in BRACKETS.finditer(text)],
        'val macros': re.findall(r'\\val[A-Za-z]+', text),
        'snippets': [f'{a}{{{b}}}{{{c}}}' for a, b, c in
                     re.findall(r'\\(snippet|snippetheader)\{([^}]*)\}\{([^}]*)\}', text)],
        'cite keys': [k.strip() for g in re.findall(r'\\cite\{([^}]*)\}', text) for k in g.split(',')],
        'labels': re.findall(r'\\label\{([^}]*)\}', text) + re.findall(r'label=([a-z]+:[A-Za-z0-9-]+)', text),
        'cref targets': [k.strip() for g in re.findall(r'\\[cC]ref\{([^}]*)\}', text) for k in g.split(',')],
        'inputs': [f'{a}{{{b}}}' for a, b in re.findall(r'\\(input|tablebody|tablerows)\{([^}]*)\}', text)],
    }


def inline_maths(text: str) -> collections.Counter:
    return collections.Counter(norm(m) for m in INLINE.findall(without_display(text)))


def code_names(text: str) -> collections.Counter:
    return collections.Counter(CODE.findall(text))


def prose_numbers(text: str) -> collections.Counter:
    text = INLINE.sub(' ', without_display(text))
    text = CODE.sub(' ', text)
    text = re.sub(r'\\(label|cref|Cref|cite|snippet|snippetheader|input|tablebody|tablerows|begin|end|needspace)'
                  r'(\{[^}]*\})+', ' ', text)
    text = re.sub(r'\\[A-Za-z]+', ' ', text)
    text = re.sub(r'\[[^\]]*\]', ' ', text)
    return collections.Counter(re.findall(r'(?<![\w.])\d+(?:\.\d+)?', text))


def compare(revision: str, chapter: str) -> bool:
    path = (BOOK / chapter).resolve()
    relative = path.relative_to(REPO).as_posix()
    old = strip_comments(subprocess.run(['git', '-C', str(REPO), 'show', f'{revision}:{relative}'],
                                        capture_output=True, text=True, check=True).stdout)
    new = strip_comments(path.read_text(encoding='ascii'))
    print(f'== {relative} ({revision} -> working tree)')
    ok = True
    so, sn = strict_sequences(old), strict_sequences(new)
    for key in so:
        if so[key] == sn[key]:
            print(f'  same    {key}: {len(so[key])}, in the same order')
        else:
            ok = False
            print(f'  DIFFER  {key}: {len(so[key])} -> {len(sn[key])}')
            for line in difflib.unified_diff(so[key], sn[key], lineterm='', n=1):
                print('          ' + line)
    print(f'  review  margin notes: {len(re.findall(r"\\mnote\{", old))} -> {len(re.findall(r"\\mnote\{", new))}')
    for name, count in (('inline maths', inline_maths), ('code', code_names), ('prose numbers', prose_numbers)):
        dropped = sorted((count(old) - count(new)).elements())
        added = sorted((count(new) - count(old)).elements())
        print(f'  review  {name}: dropped {dropped or "none"}; added {added or "none"}')
    words = [len(re.findall(r'[A-Za-z]+', re.sub(r'\\[A-Za-z]+', ' ', t))) for t in (old, new)]
    print(f'  info    words, roughly: {words[0]} -> {words[1]}')
    return ok


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    results = [compare(sys.argv[1], chapter) for chapter in sys.argv[2:]]
    print('RESULT:', 'identical content sequences' if all(results) else 'CONTENT SEQUENCES DIFFER')
    return 0 if all(results) else 1


if __name__ == '__main__':
    sys.exit(main())
