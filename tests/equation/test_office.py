"""Markdown -> Word, checked at the boundary that matters.

What is verified is that the equations arrive as *equations*: the .docx carries
one ``<m:oMath>`` per equation in the source, not a picture and not literal
``$x$`` text.  Whether Word lays them out beautifully is Word's business.
"""

from __future__ import annotations

import zipfile

import pytest

equation = pytest.importorskip("radia.equation")
pytest.importorskip("docx")
pytest.importorskip("lxml")

SOURCE = """# Title

Text with $a_{j}$ inline.

$$\\sigma_{f}(x) = \\sum_{j} B_{fj}(x) a_{j}$$

- bullet with $\\frac{a}{b}$
- bullet with `$HOME` which is code, not math

Prices are $5 and $6.

```bash
echo "$PATH"
```

End.
"""

EXPECTED_MATH = [
    "a_{j}",
    "\\sigma_{f}(x) = \\sum_{j} B_{fj}(x) a_{j}",
    "\\frac{a}{b}",
]


def test_split_math_finds_the_equations_and_nothing_else():
    pieces = equation.split_math(SOURCE)
    assert [p.latex for p in pieces if p.is_math] == EXPECTED_MATH


def test_display_and_inline_are_distinguished():
    pieces = equation.split_math(SOURCE)
    assert [p.display for p in pieces if p.is_math] == [False, True, False]


def test_docx_carries_native_equations(tmp_path):
    out = tmp_path / "probe.docx"
    equation.markdown_to_docx(SOURCE, str(out))
    assert equation.office.count_equations(str(out)) == len(EXPECTED_MATH)


@pytest.mark.parametrize("literal", ["$5", "$6", "$HOME", "$PATH"])
def test_dollars_that_are_not_math_survive_as_text(tmp_path, literal):
    out = tmp_path / "probe.docx"
    equation.markdown_to_docx(SOURCE, str(out))
    with zipfile.ZipFile(out) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    assert literal in xml


def test_a_code_fence_keeps_its_content_and_drops_its_markers(tmp_path):
    out = tmp_path / "probe.docx"
    equation.markdown_to_docx(SOURCE, str(out))
    with zipfile.ZipFile(out) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    assert "```" not in xml
