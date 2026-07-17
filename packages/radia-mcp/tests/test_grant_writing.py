from radia_mcp.document_meta.tools import document_meta_lint_all
from radia_mcp.grant_writing import tools as gw
from radia_mcp.meta.catalog import CATALOG


KDDI_SAMPLE = (
    "本研究の主題はパワーエレクトロニクス基板CAE-AI環境の社会実装である。"
    "社会的課題は1000人以下の地域製造業のパワーエレクトロニクス設計である。"
    "本提案は生成AIとMCPを用いてLTspice、SPICE、Radia、NGSolve、PEEC、CAEを接続し、"
    "回路・寄生・EMC・熱・インダクタンスを協調して評価する。"
    "商用CAEを直ちに置換するものではなく、Python-nativeで現代的なAPIを持つ入口を作る。"
    "三菱電機でのEMC経験とIH熱解析、RadiaとLTspiceの実績を基盤に、"
    "厚銅基板のPoC試作、計測評価、OSSレポジトリ公開、MotorAIを含む導入候補への技術プレゼンを行う。"
    "1年目、2年目、3年目の年度スケジュールを定め、"
    "Claude、Codex、Fable、MDXの計算資源と基板評価費を予算化する。"
    "予算は助成上限額に近い申請額とし、単価、数量、月数、年度配分、見積根拠を積算する。"
)


def test_grant_writing_kddi_health_report_runs():
    report = gw.grant_writing_health_report(KDDI_SAMPLE, program="kddi_digital")

    assert report["program"] == "kddi_digital"
    assert report["overall_score"] > 0
    assert "kddi_digital" in report["detailed_results"]
    assert "power_electronics_focus" in report["detailed_results"]
    assert "budget" in report["detailed_results"]


def test_grant_writing_kddi_power_electronics_focus_check():
    result = gw.grant_writing_kddi_power_electronics_focus_check(KDDI_SAMPLE)

    assert result["score"] >= 8
    assert result["axis_results"]["main_theme_specificity"]["ok"]
    assert result["axis_results"]["commercial_positioning"]["ok"]


def test_grant_writing_kddi_power_electronics_focus_warns_on_generic_cae():
    result = gw.grant_writing_kddi_power_electronics_focus_check(
        "本研究の主題は1000人以下の会社に一般的なCAE導入を行うことである。"
        "商用CAEの代替を作る。"
    )

    assert result["score"] < 8
    assert result["comments"]


def test_grant_writing_budget_alignment_accepts_near_ceiling_itemization():
    result = gw.grant_writing_budget_alignment_check(KDDI_SAMPLE)

    assert result["score"] >= 9
    assert result["axis_results"]["near_ceiling_strategy"]["ok"]
    assert result["axis_results"]["itemized_calculation"]["ok"]
    assert "上限" in result["budget_policy"]


def test_grant_writing_budget_alignment_requires_ceiling_and_calculation():
    result = gw.grant_writing_budget_alignment_check(
        "Claude、Codex、Fable、MDXの計算資源と基板試作、計測評価、発表旅費を予算化する。"
    )

    assert result["score"] < 8
    assert "near_ceiling_strategy" in result["missing_axes"]
    assert "itemized_calculation" in result["missing_axes"]
    assert any("上限" in comment for comment in result["comments"])
    assert any("単価" in comment for comment in result["comments"])


def test_grant_writing_reexports_ja_lint_helpers():
    result = gw.grant_writing_lint_bedrock("これは重要であると考えられる。")

    assert isinstance(result, dict)
    assert "issue_count" in result


def test_document_meta_grant_domain_uses_public_grant_writing(tmp_path):
    draft = tmp_path / "grant.md"
    draft.write_text(KDDI_SAMPLE, encoding="utf-8")

    result = document_meta_lint_all(str(draft), domain="grant")

    assert "grant_writing_health_report" in result["lints_run"]
    assert result["lints_unavailable"] == []
    assert "grant lint stays" not in str(result)


def test_meta_catalog_lists_document_writing_trio_and_merged_presentation():
    for key in ("paper-writing", "figure", "grant-writing"):
        assert key in CATALOG
        assert CATALOG[key]["entry_point"].startswith("mcp-server-")
    # 2026-07-17: presentation merged into paper-writing (standalone
    # server retired); the old name must still resolve for discovery
    # and the tools must be advertised on the paper-writing entry.
    assert "presentation" not in CATALOG
    assert "presentation_usage" in CATALOG["paper-writing"]["primary_tools"]
    from radia_mcp.meta.catalog import _resolve
    assert _resolve("presentation") == "paper-writing"


def test_paper_writing_server_serves_merged_presentation_tools():
    from radia_mcp.paper_writing import server as pw_server
    assert pw_server._N_PRESENTATION_TOOLS > 60
