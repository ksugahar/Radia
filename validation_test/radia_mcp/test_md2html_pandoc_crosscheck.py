"""Cross-check md2html's math handling against pandoc as an independent parser.

Why this exists (2026-08-19, commit 7f571e83d): md2html was originally written
because "pandoc had bugs", then md2html itself accumulated seven defect
classes, and the migration argument had become circular.  Measurement broke
the circle: pandoc 3.9's own markdown reader handles every math invariant the
md2html audit defined -- it parses math as grammar, so the regex-layer defect
classes cannot exist in it -- while the LAB requirements ([N] citation links,
||x||->\\Vert, cp932 input) are genuinely absent from stock pandoc.

This file keeps that verdict executable instead of anecdotal:

* the INVARIANT tests assert that md2html and pandoc, fed the same source,
  both satisfy the math invariants -- a true independent-implementation
  cross-check of the arithmatex-based converter;
* the GAP test asserts the LAB features are still md2html-only.  If a future
  pandoc grows them, the assert flips and the migration calculus has changed.

Lane placement: validation_test (optional external binary, environment
specific).  Package tests in packages/radia-mcp/tests/test_md2html.py stay
the fast, dependency-light goldens.
"""
from __future__ import annotations

import html as html_lib
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

markdown = pytest.importorskip("markdown")
pytest.importorskip("pymdownx")

from radia_mcp.md2html import converter as md2html_converter  # noqa: E402
from radia_mcp.md2html import md_to_html  # noqa: E402

if not hasattr(md2html_converter, "_normalize_math_fences"):
    pytest.skip(
        "installed radia-mcp predates the arithmatex delegation (7f571e83d); "
        "point the environment at a tree that contains it",
        allow_module_level=True,
    )


def _find_pandoc() -> str | None:
    exe = shutil.which("pandoc")
    if exe:
        return exe
    local = os.path.expandvars(r"%LOCALAPPDATA%\Pandoc\pandoc.exe")
    if os.name == "nt" and os.path.isfile(local):
        return local
    return None


PANDOC = _find_pandoc()
if PANDOC is None:
    pytest.skip("pandoc not installed (PATH and %LOCALAPPDATA%\\Pandoc checked)",
                allow_module_level=True)


def _pandoc_html(src: Path) -> str:
    """pandoc's own markdown reader: raw_tex / tex_math_dollars /
    pipe_tables are all native grammar, which is the whole point of the
    comparison."""
    out = src.with_suffix(".pandoc.html")
    r = subprocess.run(
        [PANDOC, "-f", "markdown", "-t", "html5", "-s", "--mathjax",
         str(src), "-o", str(out)],
        capture_output=True, text=True, encoding="utf-8", timeout=120)
    assert r.returncode == 0, f"pandoc failed: {r.stderr[:200]}"
    return out.read_text(encoding="utf-8")


def _md2html_html(src: Path) -> str:
    out = src.with_suffix(".md2html.html")
    md_to_html(str(src), str(out), "crosscheck")
    return out.read_text(encoding="utf-8-sig")


def _both(tmp_path: Path, name: str, text: str) -> tuple[str, str]:
    src = tmp_path / name
    src.write_text(text, encoding="utf-8")
    return _md2html_html(src), _pandoc_html(src)


def _visible(h: str) -> str:
    return re.sub(r"\s+", " ", html_lib.unescape(re.sub(r"<[^>]+>", "", h)))


# ---------------------------------------------------------------------------
# Invariants: both implementations must satisfy them on identical input.
# ---------------------------------------------------------------------------

def test_bare_align_keeps_row_separators(tmp_path):
    ours, ref = _both(tmp_path, "align.md",
                      "\\begin{align}\na &= b \\\\\nc &= d\n\\end{align}\n")
    for label, h in (("md2html", ours), ("pandoc", ref)):
        body = h[h.find("\\begin{align}"):h.find("\\end{align}")]
        assert "\\\\" in body, f"{label}: align row separator collapsed"


def test_dollar_in_inline_code_stays_code(tmp_path):
    ours, ref = _both(tmp_path, "code.md",
                      "use `$PATH` here, and math $E=mc^2$ there\n")
    for label, h in (("md2html", ours), ("pandoc", ref)):
        assert re.search(r"<code[^>]*>\$PATH</code>", h), \
            f"{label}: `$PATH` left the code span"
        assert "E=mc^2" in h, f"{label}: sibling math lost"


def test_single_pipe_math_in_table_keeps_row(tmp_path):
    ours, ref = _both(tmp_path, "pipe.md",
                      "| a | b |\n|---|---|\n| x | $|H_t|$ |\n")
    for label, h in (("md2html", ours), ("pandoc", ref)):
        n_td = len(re.findall(r"<td[^>]*>", h))
        assert n_td == 2, f"{label}: table row split on math pipe (td={n_td})"
        assert "|H_t|" in _visible(h), f"{label}: math content lost"


def test_less_than_in_math_survives(tmp_path):
    ours, ref = _both(tmp_path, "lt.md", "For $\\xi < 0.5$: expand.\n")
    for label, h in (("md2html", ours), ("pandoc", ref)):
        assert "< 0.5" in _visible(h), \
            f"{label}: '< 0.5' swallowed as a tag opening"
        assert ": expand." in _visible(h), f"{label}: trailing prose lost"


# ---------------------------------------------------------------------------
# The gap that keeps md2html alive: LAB features stock pandoc does not have.
# If pandoc ever grows one of these, this test flips and the migration
# calculus should be revisited.
# ---------------------------------------------------------------------------

def test_lab_features_remain_md2html_only(tmp_path):
    text = ("cite [1] and norm $||v||$\n\n"
            "## References\n\n1. Author, 2026.\n")
    ours, ref = _both(tmp_path, "lab.md", text)

    # md2html: both features present.
    assert 'href="#ref1"' in ours and 'id="ref1"' in ours
    assert "\\Vert" in ours

    # pandoc: both absent (documented gap, not a defect).
    assert 'href="#ref' not in ref, \
        "pandoc now emits [N] citation anchors -- revisit the migration calculus"
    assert "\\Vert" not in ref, \
        "pandoc now rewrites ||x|| -- revisit the migration calculus"
