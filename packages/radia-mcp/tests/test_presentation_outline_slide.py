"""Major sections must be announced where each section begins."""
import pytest

pptx = pytest.importorskip("pptx")

from pptx.util import Inches  # noqa: E402

from radia_mcp.presentation._outline import presentation_check_outline_slide  # noqa: E402


ARC = [
    ("The ladder", ["Kameari 2018 builds it."], "A circuit, exact at DC."),
    ("Its far end", ["Ten rungs."], "A finite ladder cannot make that slope."),
    ("One space", ["Inspired by XFEM [10]."], "Keep the ladder, add one surface mode."),
    ("A cylinder", ["Copper, 5 mm."], "Two unknowns hold the band to 0.06 %."),
    ("Summary", ["Bulk plus surface."], "Keep the ladder, add the surface mode."),
]

MOTIVATION = ("1 Motivation", ["Why the finite ladder misses the tail."], "")
METHOD = ("2 Proposed method", ["Two spaces in one Galerkin system."], "")
RESULTS = ("3 Results", ["Accuracy for four geometries."], "")
AGENDA = ("Outline", ["Motivation", "Proposed method", "Results"], "")


def _deck(tmp_path, rows, name="outline.pptx"):
    prs = pptx.Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    for title, body, takeaway in [("Cover", ["Talk"], "")] + rows:
        slide = prs.slides.add_slide(prs.slide_layouts[6])
        head = slide.shapes.add_textbox(Inches(0.6), Inches(0.1),
                                        Inches(12.0), Inches(0.7))
        head.text_frame.text = title
        for line in body:
            box = slide.shapes.add_textbox(Inches(0.6), Inches(1.5),
                                           Inches(12.0), Inches(3.0))
            box.text_frame.text = line
        if takeaway:
            band = slide.shapes.add_textbox(Inches(0.4), Inches(6.4),
                                            Inches(12.5), Inches(0.6))
            band.text_frame.paragraphs[0].add_run().text = takeaway
    out = tmp_path / name
    prs.save(str(out))
    return out


def test_deck_without_section_dividers_fails(tmp_path):
    result = presentation_check_outline_slide(str(_deck(tmp_path, ARC)))
    assert result["score"] == 0.0
    assert result["section_dividers"] == []


def test_single_front_agenda_is_not_recurring_navigation(tmp_path):
    result = presentation_check_outline_slide(str(_deck(tmp_path, [AGENDA] + ARC)))
    assert result["section_coverage"]["motivation"] == [2]
    assert result["section_coverage"]["method"] == []
    assert result["section_coverage"]["results"] == []
    assert result["score"] == 3.3


def test_section_dividers_at_each_transition_pass(tmp_path):
    rows = [MOTIVATION] + ARC[:2] + [METHOD] + ARC[2:3] + [RESULTS] + ARC[3:]
    result = presentation_check_outline_slide(str(_deck(tmp_path, rows)))
    assert result["score"] == 10.0
    assert result["section_coverage"] == {
        "motivation": [2], "method": [5], "results": [7]
    }


def test_section_names_are_examples_and_japanese_labels_work(tmp_path):
    rows = [
        ("背景と課題", ["有限梯子の限界"], ""), *ARC[:2],
        ("提案手法", ["二つの空間を結合"], ""), *ARC[2:3],
        ("検証結果", ["四形状の精度"], ""), *ARC[3:],
    ]
    result = presentation_check_outline_slide(str(_deck(tmp_path, rows)))
    assert result["score"] == 10.0


def test_deck_too_short_to_need_sections_says_so(tmp_path):
    result = presentation_check_outline_slide(str(_deck(tmp_path, ARC[:2])))
    assert "too few" in result["error"]
