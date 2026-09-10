"""presentation_check_text_box_overflow / _check_quoted_figure_credit.

Both come from building the IGTE'26 deck (2026-09-06).  The overflow check
is the one that would have caught the defect that recurred five times there:
the builder shrinks the title band and the takeaway banner to fit, but not
the body, so one line too many slides under the banner and is only visible
once the slide is rendered.  The credit check is the rule that a figure
reproduced from a paper must name its source on the same slide.
"""
import json

import pytest

pptx = pytest.importorskip("pptx")

from pptx.util import Inches, Pt  # noqa: E402

from radia_mcp.presentation.tools import (  # noqa: E402
    presentation_check_quoted_figure_credit,
    presentation_check_text_box_overflow,
)


def _deck(tmp_path, paragraphs, *, box=(0.6, 1.05, 12.13, 5.35), size_pt=24,
          name="deck.pptx"):
    prs = pptx.Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    slide = prs.slides.add_slide(prs.slide_layouts[6])   # blank
    left, top, width, height = (Inches(v) for v in box)
    shape = slide.shapes.add_textbox(left, top, width, height)
    frame = shape.text_frame
    frame.word_wrap = True
    for i, text in enumerate(paragraphs):
        para = frame.paragraphs[0] if i == 0 else frame.add_paragraph()
        run = para.add_run()
        run.text = text
        run.font.size = Pt(size_pt)
    out = tmp_path / name
    prs.save(str(out))
    return out


def test_a_body_that_fits_is_not_flagged(tmp_path):
    deck = _deck(tmp_path, ["A short line.", "Another short line."])
    r = presentation_check_text_box_overflow(str(deck))
    assert r["overflow_count"] == 0
    assert r["checked_text_frames"] == 1
    assert r["off_slide_count"] == 0


def test_a_body_one_line_too_long_is_flagged(tmp_path):
    # 5.35 in = 385 pt of box; 24 pt at line spacing 1.2 is 28.8 pt a line,
    # so 14 lines fit and 20 paragraphs of two lines each cannot.
    deck = _deck(tmp_path, ["This paragraph is long enough to wrap onto a "
                            "second line of the body box." ] * 20)
    r = presentation_check_text_box_overflow(str(deck))
    assert r["overflow_count"] == 1
    bad = r["overflows"][0]
    assert bad["slide"] == 1
    assert bad["needed_pt"] > bad["available_pt"]
    assert bad["overflow_pt"] > 0
    assert bad["estimated_lines"] >= 20


def test_a_two_line_overflow_is_caught(tmp_path):
    """The real defect scale: a body two lines past its box, not twenty.

    14 lines of 24 pt fill the 5.35 in body box; 16 single-line paragraphs
    are two lines over, which is more than the one-line allowance.
    """
    deck = _deck(tmp_path, ["A single short line of body text."] * 16)
    r = presentation_check_text_box_overflow(str(deck))
    assert r["overflow_count"] == 1
    assert 1.0 < r["overflows"][0]["overflow_lines"] < 4.0


def test_a_one_line_strip_is_left_to_its_author(tmp_path):
    """A banner only ever holds one line; the estimate has nothing to add."""
    deck = _deck(tmp_path, ["Bulk and surface are different function spaces. Use both."],
                 box=(0.39, 6.445, 12.50, 0.54), size_pt=28)
    r = presentation_check_text_box_overflow(str(deck))
    assert r["overflow_count"] == 0
    assert r["skipped"] and "one-line strip" in r["skipped"][0]["why"]


def test_wrapping_uses_the_box_width(tmp_path):
    """The same text in a narrow box needs more lines than in a wide one."""
    text = ["Bulk and surface are different function spaces, so use both of them."] * 6
    wide = presentation_check_text_box_overflow(
        str(_deck(tmp_path, text, box=(0.6, 1.05, 12.13, 5.35), name="wide.pptx")))
    narrow = presentation_check_text_box_overflow(
        str(_deck(tmp_path, text, box=(0.6, 1.05, 6.0, 5.35), name="narrow.pptx")))
    assert wide["overflow_count"] == 0
    assert narrow["overflow_count"] == 0          # still fits, but with more lines
    assert narrow["checked_text_frames"] == wide["checked_text_frames"] == 1


def test_a_box_hanging_off_the_slide_is_reported(tmp_path):
    deck = _deck(tmp_path, ["One line."], box=(0.6, 7.2, 6.0, 1.5))
    r = presentation_check_text_box_overflow(str(deck))
    assert r["off_slide_count"] == 1
    assert r["off_slide"][0]["bottom_pt"] > r["off_slide"][0]["slide_height_pt"]


def test_japanese_width_is_the_callers_choice(tmp_path):
    """A CJK body needs about one em a character, so the caller widens it.

    Thirteen 27-character lines fit a 6 in x 5.35 in box at the Latin
    calibration and cannot fit it at one em a character.
    """
    deck = _deck(tmp_path, ["これは日本語の本文で、一文字が約一エムの幅を占めます。"] * 13,
                 box=(0.6, 1.05, 6.0, 5.35))
    latin = presentation_check_text_box_overflow(str(deck))
    cjk = presentation_check_text_box_overflow(str(deck), char_width_em=1.0)
    assert latin["overflow_count"] == 0
    assert cjk["overflow_count"] == 1
    assert cjk["overflows"][0]["estimated_lines"] == 26


def test_missing_file_and_bad_arguments(tmp_path):
    assert "error" in presentation_check_text_box_overflow(str(tmp_path / "nope.pptx"))
    deck = _deck(tmp_path, ["x"])
    assert "error" in presentation_check_text_box_overflow(str(deck), char_width_em=0)


# --------------------------------------------------------------------------
# quoted-figure credit
# --------------------------------------------------------------------------
def _deck_with_picture(tmp_path, body_text, name="q.pptx"):
    png = tmp_path / "fig.png"
    # a 2x2 white PNG, written by hand so the test needs no image library
    png.write_bytes(bytes.fromhex(
        "89504e470d0a1a0a0000000d4948445200000002000000020802000000fdd49a"
        "73000000164944415408d763fcffff3f0323030303130000f3ff01fd3d0a0f00"
        "00000049454e44ae426082"))
    prs = pptx.Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    shape = slide.shapes.add_textbox(Inches(0.6), Inches(1.0), Inches(6.0), Inches(3.0))
    shape.text_frame.text = body_text
    pic = slide.shapes.add_picture(str(png), Inches(7.0), Inches(1.0), Inches(4.0))
    out = tmp_path / name
    prs.save(str(out))
    return out, pic.name


def _evidence(tmp_path, shape_name, quoted: bool):
    entry = {"slide": 1, "shape": shape_name,
             "source_evidence": {"minimum_source_font_pt": 7.0, "source_width_cm": 4.4}}
    if quoted:
        entry["source_evidence"]["quoted_from"] = "K. Kuriyama et al., TMag 2019, Fig. 9(a)"
    path = tmp_path / "figure_text_source_evidence.json"
    path.write_text(json.dumps({"pictures": [entry]}), encoding="utf-8")
    return path


def test_a_quoted_figure_without_a_credit_is_a_violation(tmp_path):
    deck, shape = _deck_with_picture(tmp_path, "The band follows the second point.")
    _evidence(tmp_path, shape, quoted=True)
    r = presentation_check_quoted_figure_credit(str(deck))
    assert r["violation_count"] == 1
    assert r["violations"][0]["slide"] == 1
    assert "Kuriyama" in r["violations"][0]["quoted_from"][0]


def test_a_quoted_figure_with_a_credit_passes(tmp_path):
    deck, shape = _deck_with_picture(
        tmp_path, "The band follows the second point.\nFigure: [2] Fig. 9(a), © 2019 IEEE.")
    _evidence(tmp_path, shape, quoted=True)
    r = presentation_check_quoted_figure_credit(str(deck))
    assert r["violation_count"] == 0
    assert r["credited"] and r["credited"][0]["slide"] == 1


def test_a_credit_on_an_own_figure_is_reported_as_unverified(tmp_path):
    deck, shape = _deck_with_picture(tmp_path, "Our own plot.\n© 2019 IEEE.")
    _evidence(tmp_path, shape, quoted=False)
    r = presentation_check_quoted_figure_credit(str(deck))
    assert r["violation_count"] == 0
    assert r["unverified_credit"] and r["unverified_credit"][0]["slide"] == 1


def test_without_evidence_nothing_can_be_claimed(tmp_path):
    deck, _ = _deck_with_picture(tmp_path, "No evidence file beside this deck.")
    r = presentation_check_quoted_figure_credit(str(deck))
    assert r["evidence"] is None
    assert r["violation_count"] == 0
    assert r["quoted_pictures"] == 0


# --------------------------------------------------------------------------
# mathematics (2026-09-06): python-pptx reads only a:r runs, so an equation
# was zero characters to the estimate.  Two IGTE'26 slides whose bodies are
# half mathematics were reported clean and rendered under the takeaway band.
# --------------------------------------------------------------------------
_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
_M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
_A14 = "http://schemas.microsoft.com/office/drawing/2010/main"


def _omml_paragraph(text, size_pt=24, nary=False, prose=""):
    """One a:p carrying OMML, optionally with prose beside it."""
    from lxml import etree
    runs = "".join(f"<m:r><m:t>{ch}</m:t></m:r>" for ch in text)
    if nary:
        runs = ("<m:nary><m:naryPr><m:chr m:val=\"\u2211\"/></m:naryPr>"
                f"<m:e>{runs}</m:e></m:nary>")
    lead = (f"<a:r><a:rPr sz=\"{int(size_pt * 100)}\"/>"
            f"<a:t>{prose}</a:t></a:r>") if prose else ""
    xml = (f'<a:p xmlns:a="{_A}" xmlns:m="{_M}" xmlns:a14="{_A14}">'
           f'<a:pPr><a:defRPr sz="{int(size_pt * 100)}"/></a:pPr>'
           f'{lead}<a14:m><m:oMathPara><m:oMath>{runs}'
           f'</m:oMath></m:oMathPara></a14:m></a:p>')
    return etree.fromstring(xml)


def _math_deck(tmp_path, n, *, chars=70, nary=False, prose="", name="m.pptx"):
    prs = pptx.Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    shape = slide.shapes.add_textbox(Inches(0.6), Inches(1.05),
                                     Inches(6.0), Inches(5.35))
    frame = shape.text_frame
    frame.word_wrap = True
    body = frame.paragraphs[0]._p
    for _ in range(n):
        body.addprevious(_omml_paragraph("x" * chars, nary=nary, prose=prose))
    body.getparent().remove(body)
    out = tmp_path / name
    prs.save(str(out))
    return out


def test_mathematics_counts_towards_the_box(tmp_path):
    """Six 70-character equations cannot fit where the estimate saw nothing."""
    deck = _math_deck(tmp_path, 6)
    r = presentation_check_text_box_overflow(str(deck))
    assert r["checked_text_frames"] == 1          # not "an empty frame"
    assert r["overflow_count"] == 1
    assert r["overflows"][0]["estimated_lines"] >= 6


def test_a_short_equation_still_fits(tmp_path):
    """The fix must not turn every equation into an overflow."""
    deck = _math_deck(tmp_path, 2, chars=20, name="short.pptx")
    assert presentation_check_text_box_overflow(str(deck))["overflow_count"] == 0


def test_a_display_sum_is_taller_than_a_line_of_prose(tmp_path):
    """A built-up operator raises the line box; tall_math_scale says by how much."""
    plain = _math_deck(tmp_path, 5, chars=30, name="plain.pptx")
    tall = _math_deck(tmp_path, 5, chars=30, nary=True, name="tall.pptx")
    assert presentation_check_text_box_overflow(str(plain))["overflow_count"] == 0
    assert presentation_check_text_box_overflow(str(tall))["overflow_count"] == 1


def test_inline_mathematics_does_not_stretch_a_prose_paragraph(tmp_path):
    """A sum inside a sentence must not scale that paragraph's whole height.

    Scaling every paragraph that merely CONTAINS a tall structure reported
    five clean IGTE'26 slides as overflowing; the factor belongs to display
    equations, which stand on a line of their own.
    """
    deck = _math_deck(tmp_path, 5, chars=10, nary=True,
                      prose="Inline, in a sentence: ", name="inline.pptx")
    assert presentation_check_text_box_overflow(str(deck))["overflow_count"] == 0


def test_math_char_em_must_be_positive(tmp_path):
    deck = _deck(tmp_path, ["x"], name="args.pptx")
    assert "error" in presentation_check_text_box_overflow(str(deck), math_char_em=0)
    assert "error" in presentation_check_text_box_overflow(str(deck),
                                                           tall_math_scale=0)
