"""Review regressions: lexical signals must not imply verified capability."""

import pytest

from radia_mcp.grant_writing import tools as t


@pytest.mark.parametrize("body", [
    "H行列加速、NGSolve連携、HDiv-MMMを統合する公開電磁界解析基盤。",
    "Radiaは磁場を計算し、磁極形状を評価する。",
    "Radiaに統合し、磁場を計算する。",
])
def test_nonpast_operation_is_ambiguous_not_future(body):
    result = t.grant_writing_capability_status_map(body, capability_names="Radia,HDiv-MMM")
    statements = [s for c in result["capabilities"] for s in c["statements"]]
    assert statements and all(s["status"] == "ambiguous" for s in statements)


@pytest.mark.parametrize("body", [
    "Radiaはまだ実装済みではない。", "Radia is not yet implemented.",
    "Radiaは実装していない。", "Radiaは検証済みではなかった。",
])
def test_negated_completion_never_becomes_current(body):
    c = t.grant_writing_capability_status_map(body, capability_names="Radia")["capabilities"][0]
    assert not c["has_current_statement"]
    assert c["statements"][0]["status"] == "negated_or_mixed"


@pytest.mark.parametrize("body", ["Radiaを実装し、公開した。", "Radiaを開発してきた。"])
def test_past_completion_is_a_current_statement(body):
    c = t.grant_writing_capability_status_map(body, capability_names="Radia")["capabilities"][0]
    assert c["has_current_statement"] and not c["has_future_statement"]


@pytest.mark.parametrize("name", ["NGSolve", "Gmsh"])
def test_product_name_is_not_an_operation(name):
    entry = t.grant_writing_named_software_first_use_check(
        f"解析ソフトウェア{name}を用いる。"
    )["entries"][0]
    assert not entry["functional_identity_present"]
    assert not entry["operation_hits"]


def test_prestige_label_is_not_a_category():
    entry = t.grant_writing_named_software_first_use_check(
        "公開研究基盤Radiaを用い、FFAG磁石の形状を最適化する。"
    )["entries"][0]
    assert not entry["functional_identity_present"]


@pytest.mark.parametrize("prefix", [
    "# 公開研究基盤Radiaによる磁場設計\n",
    "% Radiaを用いる。\n",
    r"\section{Radiaによる磁場設計}" + "\n",
    "| Radia | 開発予定 |\n| --- | --- |\n",
    "```text\nRadiaを用いる。\n```\n",
])
def test_nonprose_is_not_first_use(prefix):
    body = "Radiaは磁場を計算するソフトウェアである。"
    entry = t.grant_writing_named_software_first_use_check(prefix + body)["entries"][0]
    assert entry["functional_identity_present"]
    assert entry["first_use"] == body


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_soft_wrap_and_version_do_not_truncate_definition(newline):
    entry = t.grant_writing_named_software_first_use_check(
        "Radia v4.95は磁石の形状から" + newline + "磁場を計算するソフトウェアである。"
    )["entries"][0]
    assert entry["functional_identity_present"]
    assert "v4.95" in entry["first_use"]
    assert entry["first_use"].endswith("ソフトウェアである。")


def test_ascii_periods_in_url_and_abbreviation_do_not_split():
    from radia_mcp.grant_writing._draft_checks import _named_diagnostic_sentences
    text = "Radia (https://example.org/v4.95; IEEE Trans.) は磁場を計算するソフトウェアである。"
    assert _named_diagnostic_sentences(text) == [{"line": 1, "text": text}]


def test_english_whole_word_operation_still_matches():
    entry = t.grant_writing_named_software_first_use_check(
        "Radia is software that computes magnetic fields."
    )["entries"][0]
    assert entry["functional_identity_present"]


def test_latex_table_is_not_a_capability_statement():
    text = r"\begin{tabular}{ll}Radia & planned \\ \end{tabular}" + "\nRadiaを公開した。"
    item = t.grant_writing_capability_status_map(text, capability_names="Radia")["capabilities"][0]
    assert item["mention_count"] == 1
    assert item["has_current_statement"] and not item["has_future_statement"]


def test_plain_percentage_does_not_start_a_tex_comment():
    entry = t.grant_writing_named_software_first_use_check(
        "Radiaは誤差2%以下で磁場を計算するソフトウェアである。"
    )["entries"][0]
    assert entry["functional_identity_present"]
