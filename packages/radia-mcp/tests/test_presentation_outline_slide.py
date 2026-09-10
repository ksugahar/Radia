"""途中の目次: 転の直前に境目の枚があるか。

松尾先生 2026-09-10:「途中に目次/アウトラインを挟むと分かりやすい。どこまでが
現状の問題で、どこからが今回の解決法かを明示する方が聴衆には分かりやすい」。

だからこの検査の要は「目次があるか」ではなく「**境目に**あるか」で、冒頭に
一度だけ置いた目次は不合格になる ―― それがまさに、聞き手が現在地を見失う
場面で役に立たない置き方だからである。
"""
import pytest

pptx = pytest.importorskip("pptx")

from pptx.util import Inches  # noqa: E402

from radia_mcp.presentation._outline import (  # noqa: E402
    presentation_check_outline_slide,
)

# 起, 承, 転, 結, 結 -- 起承転結 の検査と同じ形。転は 4 枚目に落ちる。
ARC = [
    ("The ladder", ["Kameari 2018 builds it."], "A circuit, exact at DC."),
    ("Its far end", ["Ten rungs."], "A finite ladder cannot make that slope."),
    ("One space", ["Inspired by XFEM [10]: basis and enrichment in one matrix."],
     "Keep the ladder, add one surface mode."),
    ("A cylinder", ["Copper, 5 mm."], "Two unknowns hold the band to 0.06 %."),
    ("Summary", ["Bulk plus surface."], "Keep the ladder, add the surface mode."),
]
# 区切りの枚。下端の主張文は持たない（持たせると筋の判定が動いてしまう）。
DIVIDER = ("Outline", ["Background and the problem: slides 2-3",
                       "Our method and its evidence: slides 5-7"], "")
# 章の名前だけを並べ、どちらが問題でどちらが解決法かを言わない目次。
BARE = ("Outline", ["The ladder", "One space", "A cylinder"], "")


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


def test_a_deck_with_no_outline_fails_and_is_told_where_to_put_one(tmp_path):
    r = presentation_check_outline_slide(str(_deck(tmp_path, ARC)))
    assert r["outline_slides"] == []
    assert r["score"] == 0.0
    # 転は 4 枚目。助言はその直前を指す。
    assert r["turn_slide"] == 4
    assert r["suggested_outline"]["insert_before_slide"] == 4
    assert [s["slide"] for s in r["suggested_outline"]["problem"]] == [2, 3]
    assert 4 in [s["slide"] for s in r["suggested_outline"]["solution"]]


def test_an_outline_only_at_the_front_is_not_at_the_boundary(tmp_path):
    rows = [DIVIDER] + ARC          # 2 枚目が目次、転は 5 枚目
    r = presentation_check_outline_slide(str(_deck(tmp_path, rows)))
    assert [f["slide"] for f in r["outline_slides"]] == [2]
    assert r["outline_at_boundary"] == []
    assert r["checks"]["目次・区切りの枚がある"] is True
    assert r["checks"]["その一枚が転の直前にある（冒頭だけではない）"] is False
    assert any("境目" in c for c in r["comments"])


def test_a_divider_just_before_the_turn_passes(tmp_path):
    rows = ARC[:2] + [DIVIDER] + ARC[2:]     # 4 枚目が区切り、転は 5 枚目
    r = presentation_check_outline_slide(str(_deck(tmp_path, rows)))
    assert r["turn_slide"] == 5
    assert r["outline_at_boundary"] == [4]
    assert r["score"] == 10.0


def test_an_outline_that_does_not_name_the_two_sides_is_flagged(tmp_path):
    rows = ARC[:2] + [BARE] + ARC[2:]
    r = presentation_check_outline_slide(str(_deck(tmp_path, rows)))
    assert r["outline_at_boundary"] == [4]
    assert r["checks"]["問題の側と解決法の側を両方名指ししている"] is False
    assert any("どこまでが問題" in c for c in r["comments"])


def test_a_divider_is_found_by_its_content_when_the_title_is_its_own(tmp_path):
    # 題が目次語でなくても、デッキ自身の題を並べていれば区切りと分かる。
    own = ("From problem to proposal", ["The ladder and its far end",
                                       "One space, two bases",
                                       "A cylinder"], "")
    rows = ARC[:2] + [own] + ARC[2:]
    r = presentation_check_outline_slide(str(_deck(tmp_path, rows)))
    assert [f["slide"] for f in r["outline_slides"]] == [4]
    assert "lists" in r["outline_slides"][0]["detected_by"]


def test_a_deck_too_short_to_have_sections_says_so(tmp_path):
    r = presentation_check_outline_slide(str(_deck(tmp_path, ARC[:2])))
    assert "too few" in r["error"]
