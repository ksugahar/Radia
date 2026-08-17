"""RTF: the form an equation takes on the clipboard.

Word offers "Rich Text Format" when an equation is copied, and the maths inside
it is OMML spelled as control words.  These checks assert that spelling
directly -- it was measured by copying each construct out of Word, so a
deviation here means the transcription drifted, not that Word changed.

The end-to-end proof (put it on the clipboard, paste into Word, read back
<m:oMath> from the saved .docx) needs Word itself and lives in
validation_test; what is checked here is everything that does not.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

# LaTeX in, the control words Word writes for it.  Only the distinguishing
# fragments are asserted, not the whole payload: the run formatting Word emits
# around each piece is its business, not ours.
SPELLING = [
    (r"\frac{a}{b}",        [r"{\mf", r"{\mfPr", r"{\mctrlPr}", r"{\mnum", r"{\mden"]),
    (r"a^{2}",              [r"{\msSup", r"{\msSupPr", r"{\me", r"{\msup"]),
    (r"a_{i}",              [r"{\msSub", r"{\msSubPr", r"{\msub"]),
    (r"a_{i}^{2}",          [r"{\msSubSup", r"{\msub", r"{\msup"]),
    (r"\sqrt{x}",           [r"{\mrad", r"{\mradPr", r"{\mdegHide on}", r"{\mdeg}"]),
    (r"\sqrt[3]{x}",        [r"{\mrad", r"{\mdeg"]),
    (r"\left(x\right)",     [r"{\md", r"{\mdPr"]),
    (r"\left[x\right]",     [r"{\mbegChr [}", r"{\mendChr ]}"]),
    (r"\sum_{i}^{n} a",     [r"{\mnary", r"{\mnaryPr", r"{\mchr", r"{\mlimLoc undOvr}"]),
    (r"\int_{0}^{1} a",     [r"{\mlimLoc subSup}"]),
    (r"\oint_{C} a",        [r"{\msupHide on}"]),
    (r"\begin{pmatrix} a & b \\ c & d \end{pmatrix}",
                            [r"{\mm", r"{\mmPr", r"{\mmcs", r"{\mmc", r"{\mcount 2}",
                             r"{\mmcJc center}", r"{\mmr"]),
    (r"\hat{x}",            [r"{\macc", r"{\maccPr"]),
    (r"\overline{x}",       [r"{\mbar", r"{\mpos top}"]),
]


@pytest.mark.parametrize("latex,expected", SPELLING, ids=[c[0] for c in SPELLING])
def test_control_words(latex, expected):
    rtf = equation.tex_to_rtf(latex)
    for word in expected:
        assert word in rtf, f"{word!r} missing from the RTF for {latex}"


def test_is_a_complete_rtf_document():
    """A bare fragment is not something Word will paste."""
    rtf = equation.tex_to_rtf(r"a+b")
    assert rtf.startswith(r"{\rtf1")
    assert rtf.endswith("}")
    assert r"Cambria Math;" in rtf          # the maths font must be declared
    assert r"{\*\moMath" in rtf
    assert rtf.count("{") == rtf.count("}"), "unbalanced braces"


def test_display_uses_its_own_line():
    opts = equation.RtfOptions()
    opts.display = True
    assert r"{\*\moMathPara" in equation.tex_to_rtf(r"a+b", opts)
    assert r"{\*\moMathPara" not in equation.tex_to_rtf(r"a+b")


def test_runs_carry_their_style():
    """Word writes \\msty2 for an italic variable and \\msty0 for an upright
    operator; losing that makes every symbol look like a variable."""
    rtf = equation.tex_to_rtf("x+1")
    assert r"\msty2 x" in rtf
    assert r"\msty0 +" in rtf
    assert r"\msty0 1" in rtf


def test_non_ascii_is_escaped_not_dropped():
    """RTF is a byte format: a Greek letter has to become a \\uN escape."""
    rtf = equation.tex_to_rtf(r"\sigma")
    assert r"\u963?" in rtf                 # U+03C3
    assert "σ" not in rtf              # and not the raw character


def test_braces_in_text_are_escaped():
    rtf = equation.tex_to_rtf(r"\left\{x\right\}")
    # The delimiters reach RTF as escaped literals, not as group syntax.
    assert r"\{" in rtf and r"\}" in rtf
    assert rtf.count("{") == rtf.count("}")


@pytest.mark.parametrize("latex", [c[0] for c in SPELLING])
def test_rtf_and_omml_describe_the_same_equation(latex):
    """Both come from one walk; a divergence means the walk was duplicated."""
    import re

    rtf = equation.tex_to_rtf(latex)
    omml = equation.tex_to_omml(latex)
    for elem in ["f", "nary", "rad", "d", "m", "acc", "bar", "sSub", "sSup"]:
        # A control word ends at a non-letter, or "{\md" also matches "{\mden".
        in_rtf = re.search(r"\{\\m" + elem + r"(?![A-Za-z])", rtf) is not None
        in_omml = f"<m:{elem}>" in omml
        assert in_rtf == in_omml, \
            f"{latex}: m:{elem} appears in only one of the two outputs"
