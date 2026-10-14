"""Everything the editor can emit must be a command LaTeX actually knows.

The saved-file rule was "ASCII TeX commands, never a bare Unicode glyph", and
test_symbols.py checks exactly that.  Being ASCII is not the same as being
real: \\nsubset round-tripped perfectly through the parser for as long as this
program has existed and is not a command in LaTeX or in any package, so every
equation containing one failed to typeset.  Nothing noticed, because nothing
had ever put the output in front of pdflatex.

This does.  Every symbol in the table and every template goes through
pdflatex, and anything undefined is a failure -- except the handful listed in
PACKAGES, which are real commands that need a package the saved fragment does
not carry.  Those are recorded with the package that supplies them, and the
list is checked both ways: a new dependency fails here, and one that quietly
stops needing its package fails too, so the list cannot drift out of date.

Run:  python tests\\test_tex_compiles.py
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

import eqnedit_core as E  # noqa: E402

# Real commands that need a package beyond amsmath.  A saved equation is a
# fragment, so the host document supplies these; docs/GUI_SPEC.md says so.
PACKAGES = {
    r"\therefore": "amssymb",
    r"\because": "amssymb",
    r"\nexists": "amssymb",
    r"\varkappa": "amssymb",
    r"\mathbb": "amsfonts",
    r"\mathfrak": "amsfonts",
    r"\bm": "bm",
    r"\cancel": "cancel",
    # \oiiint is deliberately absent: the closed volume integral exists only
    # in packages that replace the whole math font, so no template offers it
    # and the editor never originates one.  Pasted \oiiint still displays.
}

BASE = ["amsmath"]
EXTRA = sorted(set(PACKAGES.values()))


def pieces():
    out = [("symbol %s" % s, s) for s in E.symbol_commands()]
    out.extend([
        ("alphabet mathrm", r"\mathrm{x}"),
        ("alphabet mathit", r"\mathit{x}"),
        ("alphabet mathbf", r"\mathbf{x}"),
        ("alphabet mathsf", r"\mathsf{x}"),
        ("alphabet mathtt", r"\mathtt{x}"),
        ("alphabet mathcal", r"\mathcal{X}"),
        ("alphabet mathbb", r"\mathbb{R}"),
        ("alphabet mathfrak", r"\mathfrak{F}"),
        ("alphabet bm", r"\bm{\alpha}"),
        ("alphabet mathnormal", r"\mathnormal{x}"),
    ])
    for kind in E.Equation.templates():
        equation = E.Equation()
        equation.load_latex("")
        if equation.insert_template(kind):
            equation.insert_text("x")
            out.append(("template %s" % kind, equation.latex()))
    return out


def undefined(items, packages, workdir):
    """Which items pdflatex cannot typeset, and the token it choked on."""
    head = ["\\documentclass[12pt]{article}"]
    head += ["\\usepackage{%s}" % p for p in packages]
    head.append("\\begin{document}")
    body = []
    for label, tex in items:
        body.append("%% CHECK %s" % label)
        body.append("$%s$" % tex)
        body.append("")
    text = "\n".join(head + body + ["\\end{document}", ""])

    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "check.tex"), "w", encoding="utf-8") as fh:
        fh.write(text)
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "check.tex"],
                   capture_output=True, cwd=workdir)
    log = os.path.join(workdir, "check.log")
    if not os.path.exists(log):
        raise RuntimeError("pdflatex produced no log")
    with open(log, encoding="utf-8", errors="replace") as fh:
        contents = fh.read()

    lines = text.splitlines()
    bad = {}
    for match in re.finditer(r"! Undefined control sequence\.(.*?)l\.(\d+)",
                             contents, re.S):
        line_no = int(match.group(2))
        label = None
        for back in range(min(line_no, len(lines)) - 1, -1, -1):
            if lines[back].startswith("% CHECK "):
                label = lines[back][len("% CHECK "):]
                break
        token = re.findall(r"(\\[A-Za-z]+)\s*$", match.group(1).strip())
        bad[label or ("line %d" % line_no)] = token[-1] if token else "?"
    return bad


def main() -> int:
    if not shutil.which("pdflatex"):
        print("skip  pdflatex is not installed; cannot check that the output "
              "compiles")
        return 0

    work = os.path.join(os.environ.get("TEMP", "."), "eqncompiletest")
    items = pieces()
    failures = []

    # With every package loaded, nothing at all may be undefined.
    with_all = undefined(items, BASE + EXTRA, work)
    for label, token in sorted(with_all.items()):
        failures.append("%s emits %s, which no package defines" % (label, token))

    # With amsmath alone, exactly the recorded dependencies may fail.
    bare = undefined(items, BASE, work)
    seen = set()
    for label, token in sorted(bare.items()):
        if token in PACKAGES:
            seen.add(token)
            continue
        if label not in with_all:      # already reported above
            failures.append(
                "%s emits %s, which needs a package that is not recorded "
                "in PACKAGES" % (label, token))
    for token, package in sorted(PACKAGES.items()):
        if token not in seen and token not in [t for t in with_all.values()]:
            failures.append(
                "%s is recorded as needing %s but compiles without it -- "
                "remove it from PACKAGES" % (token, package))

    total = len(items)
    if failures:
        print("FAIL  %d of %d" % (len(failures), total))
        for line in failures:
            print("  " + line)
        return 1
    print("ok    %d symbols and templates compile; %d need a package (%s)"
          % (total, len(PACKAGES), ", ".join(EXTRA)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
