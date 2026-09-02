"""直訳調・AI調の検査。2026-09-02 の基盤C通読で残った実例が種になっている。"""
from radia_mcp._shared.translationese import check_translationese
from radia_mcp.grant_writing import tools as gw
from radia_mcp.paper_writing import cross_lint as paper
from radia_mcp.presentation import cross_lint as pres


def test_flags_intransitive_verb_used_transitively_as_grammar():
    result = check_translationese(
        "四名の資産を誘導加熱へ展開し、その判定則を加速器電磁石設計へ発展する。"
    )

    grammar = [r for r in result["risks"] if r["type"] == "intransitive_verb_used_transitively"]
    assert grammar and grammar[0]["severity"] == "HIGH"
    assert grammar[0]["examples"][0]["suggestion"] == "を発展させる"


def test_accepts_the_causative_form():
    result = check_translationese(
        "四名の手法を誘導加熱で結合し、その判定則を加速器電磁石設計へ発展させる。"
    )

    assert all(r["type"] != "intransitive_verb_used_transitively" for r in result["risks"])


def test_flags_lowercase_english_gloss_but_not_acronym_definition():
    result = check_translationese(
        "AIが利用する知識・実行インターフェース（MCP）は、各機関の技術を結ぶ接着層（glue）である。"
    )

    gloss = [r for r in result["risks"] if r["type"] == "english_gloss_after_japanese_term"]
    assert gloss and gloss[0]["examples"][0]["gloss"] == "glue"
    assert gloss[0]["count"] == 1


def test_flags_calques_with_a_replacement():
    result = check_translationese(
        "特定ソフトウェアの提供や劇的な差の発生に成否を依存させない。"
        "本手法は高速な設計探索を可能にする。"
    )

    calque = next(r for r in result["risks"] if r["type"] == "calque_phrase")
    labels = {hit["why"] for hit in calque["examples"]}
    assert "dramatic の直訳" in labels
    assert "enable の直訳" in labels
    assert all(hit["better"] for hit in calque["examples"])


def test_clean_japanese_scores_ten_and_english_is_inapplicable():
    clean = check_translationese(
        "手法差による設計量の変動を求め、その変動を考慮しても順位を確定できる条件を示す。"
    )
    assert clean["applicable"] and clean["score"] == 10.0 and clean["risks"] == []

    english = check_translationese("This proposal develops a coupling method.")
    assert english["applicable"] is False


def test_wrappers_share_one_implementation_and_strip_latex():
    tex = "\\section{研究目的}\n% 劇的な comment must not count\n判定則を設計へ発展する。"

    grant = gw.grant_writing_translationese_check(tex)
    assert grant["applicable"]
    assert [r["type"] for r in grant["risks"]] == ["intransitive_verb_used_transitively"]

    assert paper.paper_writing_translationese_check("判定則を設計へ発展する。")["risk_count"] == 1
    assert pres.presentation_translationese_check("判定則を設計へ発展する。")["risk_count"] == 1


def test_health_report_reports_translationese_as_a_finding():
    report = gw.grant_writing_health_report(
        "本研究の学術的問いはXである。判定則を設計へ発展する。接着層（glue）である。",
        program="kaken_generic",
    )

    finding = next(f for f in report["findings"] if f["name"] == "translationese_check")
    assert finding["severity"] == "HIGH"
    assert "translationese" in report["detailed_scores"]

    skipped = gw.grant_writing_health_report("判定則を設計へ発展する。", skip="translationese")
    assert "translationese" not in skipped["detailed_results"]
