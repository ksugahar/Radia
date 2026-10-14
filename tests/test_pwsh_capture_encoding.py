"""A pwsh or git capture must name its encoding.

pwsh 7 and git write UTF-8; a `text=True` capture without `encoding=` decodes
with the machine code page instead, which is cp932 on the lab's Windows hosts.
On 2026-09-20 that failed the main CI job on a single byte of a localised git
error message inside a subprocess reader thread -- the step under test was
fine, only the reading of its diagnostic was not.  Three separate tests have
now been repaired for this same reason, so the rule is checked rather than
remembered.

The scope is deliberately narrow: only captures of pwsh, PowerShell or git,
because those are the programs that emit localised, non-ASCII text on these
hosts.  This is not a general style rule about subprocess.
"""

from __future__ import annotations

import ast
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCANNED = ("tests", "tools")
NAMES = re.compile(r"\b(pwsh|powershell|git)\b", re.I)


def _unnamed_captures():
    bad = []
    for area in SCANNED:
        for path in sorted((ROOT / area).rglob("*.py")):
            if "fixtures" in path.parts:
                continue
            source = path.read_text(encoding="utf-8", errors="replace")
            try:
                tree = ast.parse(source)
            except SyntaxError:
                continue
            lines = source.splitlines(keepends=True)
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                if not ast.unparse(node.func).startswith("subprocess"):
                    continue
                keywords = {k.arg for k in node.keywords if k.arg}
                decodes = "text" in keywords or "universal_newlines" in keywords
                if not decodes or "encoding" in keywords:
                    continue
                segment = "".join(
                    lines[node.lineno - 1:(node.end_lineno or node.lineno)])
                if NAMES.search(segment):
                    bad.append("%s:%d" % (path.relative_to(ROOT).as_posix(),
                                          node.lineno))
    return bad


def test_every_pwsh_or_git_capture_names_its_encoding():
    unnamed = _unnamed_captures()
    assert not unnamed, (
        "these captures decode pwsh/git output with the machine code page, "
        "which raises UnicodeDecodeError on a Japanese Windows host; pass "
        'encoding="utf-8": ' + ", ".join(unnamed))


def test_the_scan_can_actually_see_such_a_capture(tmp_path):
    """Guard the guard: a checker that finds nothing proves nothing."""
    global ROOT
    sample = tmp_path / "tools"
    sample.mkdir()
    (sample / "probe.py").write_text(
        "import subprocess" + chr(10) +
        'subprocess.run(["git", "status"], capture_output=True, text=True)' + chr(10),
        encoding="utf-8")
    original = ROOT
    try:
        ROOT = tmp_path
        assert _unnamed_captures() == ["tools/probe.py:2"]
    finally:
        ROOT = original
