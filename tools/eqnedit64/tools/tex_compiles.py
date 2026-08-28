"""Does every command we can emit actually compile?

The saved-file rule is "ASCII TeX commands, never a bare Unicode glyph", and
tests/test_symbols.py checks exactly that.  Being ASCII is not the same as
being real: a command can round-trip perfectly through the parser and still
be undefined in LaTeX, and nothing here has ever put the output in front of
pdflatex to find out.

This runs every symbol in the table, and every template, through pdflatex in
one document and reports which ones LaTeX does not know -- and with which
package they would need, so the choice is between fixing the command and
documenting a dependency rather than being surprised by it later.

    python tools\\tex_compiles.py           # amsmath only, as a document gets
    python tools\\tex_compiles.py --amssymb # also load amssymb
"""
from __future__ import annotations

import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

import eqnedit_core as E  # noqa: E402


def build_document(pieces, extra_packages):
    head = ["\\documentclass[12pt]{article}", "\\usepackage{amsmath}"]
    head += ["\\usepackage{%s}" % p for p in extra_packages]
    head.append("\\begin{document}")
    body = []
    for label, tex in pieces:
        # A marker line so a failure can be tied back to its command: TeX
        # reports the line number, not what we were trying.
        body.append("%% CHECK %s" % label)
        body.append("$%s$" % tex)
        body.append("")
    return "\n".join(head + body + ["\\end{document}", ""])


def run(pieces, extra_packages, workdir):
    os.makedirs(workdir, exist_ok=True)
    name = "compiles"
    path = os.path.join(workdir, name + ".tex")
    text = build_document(pieces, extra_packages)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    subprocess.run(["pdflatex", "-interaction=nonstopmode", name + ".tex"],
                   capture_output=True, cwd=workdir)
    log = os.path.join(workdir, name + ".log")
    if not os.path.exists(log):
        raise RuntimeError("pdflatex produced no log")
    with open(log, encoding="utf-8", errors="replace") as fh:
        contents = fh.read()

    lines = text.splitlines()
    bad = {}
    # "! Undefined control sequence." then a line whose tail is the token,
    # then "l.<n> ..."
    for match in re.finditer(r"! Undefined control sequence\.(.*?)l\.(\d+)",
                             contents, re.S):
        line_no = int(match.group(2))
        label = None
        for back in range(line_no - 1, -1, -1):
            if back < len(lines) and lines[back].startswith("% CHECK "):
                label = lines[back][len("% CHECK "):]
                break
        token = re.findall(r"(\\[A-Za-z]+)\s*$", match.group(1).strip())
        bad[label or ("line %d" % line_no)] = token[-1] if token else "?"
    return bad


def main(with_amssymb):
    pieces = [("symbol %s" % s, s) for s in E.symbol_commands()]
    for kind in E.Equation.templates():
        equation = E.Equation()
        equation.load_latex("")
        if equation.insert_template(kind):
            equation.insert_text("x")
            pieces.append(("template %s" % kind, equation.latex()))

    work = os.path.join(os.environ.get("TEMP", "."), "eqncompile")
    packages = ["amssymb"] if with_amssymb else []
    bad = run(pieces, packages, work)

    print("checked %d symbols and %d templates with %s"
          % (len(E.symbol_commands()), len(E.Equation.templates()),
             ", ".join(["amsmath"] + packages)))
    if not bad:
        print("ok    every one of them compiles")
        return 0
    print("\n%d do not compile:" % len(bad))
    for label in sorted(bad):
        print("  %-28s undefined: %s" % (label, bad[label]))
    return 1


if __name__ == "__main__":
    raise SystemExit(main("--amssymb" in sys.argv))
