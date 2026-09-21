"""presentation_check_negated_premise: a negation needs an introduced premise."""
import pytest

from radia_mcp.presentation._claim_premises import (
    _stem,
    check_negated_premises,
    presentation_check_negated_premise,
)


def test_stem_keeps_the_inflections_of_one_verb_together():
    assert {_stem(w) for w in ("fit", "fits", "fitted", "fitting")} == {"fit"}
    assert {_stem(w) for w in ("tune", "tuned", "tuning")} == {"tun"}
    assert _stem("crossover") != _stem("cross")


def test_premise_first_seen_inside_the_negation_is_flagged():
    slides = [
        {"slide": 1, "title": "The ladder", "text": "A Krylov basis at s = 0."},
        {"slide": 2, "title": "Schur complement, not a fit",
         "text": "The crossover comes out of the projection. Nothing is fitted."},
    ]
    f = check_negated_premises(slides)
    assert [x["premise"] for x in f] == ["fitted"]
    assert f[0]["slide"] == 2 and f[0]["negation"].lower() == "nothing"


def test_premise_introduced_earlier_is_not_flagged():
    slides = [
        {"slide": 1, "title": "The known repair",
         "text": "The usual repair fits a crossover frequency by hand."},
        {"slide": 2, "title": "This talk", "text": "Nothing is fitted: the projection decides."},
    ]
    assert check_negated_premises(slides) == []


def test_premise_earlier_on_the_same_slide_counts():
    slides = [{"slide": 1, "title": "Time domain",
               "text": "Earlier models tuned a knob. Here: no tuning."}]
    assert check_negated_premises(slides) == []


def test_japanese_negation():
    slides = [
        {"slide": 1, "title": "手法", "text": "バルクモードと表面モードを同じ空間に置く。"},
        {"slide": 2, "title": "円柱", "text": "未知数二つ、当てはめなし。"},
    ]
    f = check_negated_premises(slides)
    assert [x["premise"] for x in f] == ["当てはめ"]
    slides.insert(1, {"slide": 1, "title": "従来", "text": "従来は接続を当てはめで決めていた。"})
    assert check_negated_premises(slides) == []


def test_qualifier_between_negation_and_premise():
    """'no empirical crossover parameter' denies the crossover, not empiricism."""
    slides = [
        {"slide": 1, "title": "Coupling",
         "text": "The crossover is a result of the projection. No empirical crossover parameter enters."},
        {"slide": 2, "title": "Time domain", "text": "Twenty-two poles, no extra tuning."},
    ]
    f = check_negated_premises(slides)
    assert [x["premise"] for x in f] == ["tuning"]


def test_stop_words_and_free_suffix():
    slides = [{"slide": 1, "title": "Cost",
               "text": "No other body was tried. The model is mesh-free."}]
    f = check_negated_premises(slides)
    assert [x["premise"] for x in f] == ["mesh"]


def test_on_a_pptx(tmp_path):
    pptx = pytest.importorskip("pptx")
    prs = pptx.Presentation()
    layout = prs.slide_layouts[1]
    for title, body in (("The ladder", "Krylov basis at s = 0."),
                        ("Result", "Two unknowns. Nothing fitted.")):
        s = prs.slides.add_slide(layout)
        s.shapes.title.text = title
        s.placeholders[1].text = body
    out = tmp_path / "deck.pptx"
    prs.save(str(out))
    r = presentation_check_negated_premise(str(out))
    assert r["ok"] is False and r["n_findings"] == 1
    assert r["findings"][0]["premise"] == "fitted" and r["findings"][0]["slide"] == 2
    assert "hint" in r
