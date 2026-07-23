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

KAKEN_OSS_SAMPLE = (
    "本研究の学術的な問いは、日本語の技術報告をAI時代の実行可能な共同研究資産へ"
    "変換する条件は何かである。研究室単位の属人化を属研究室化と定義し、"
    "AIによる開発加速が重複実装を増幅することを防ぐ。"
    "JP-MARsを研究基盤としてGitHub上で運営し、"
    "issue、pull request、CI、ライセンス、ドキュメント、貢献者表記を整備する。"
    "問題仕様、参照実装、ベンチマーク、software.html、elemag/index.php、"
    "geometry.php、mcp-server教材を公開する。AIとPython/MCPによる再実行で、"
    "境界条件、メッシュ、収束条件、失敗例として暗黙知を顕在化する。"
    "新規実装前にmeshio等の既存OSSとの機能重複を調べ、再利用または上流貢献する。"
    "試験データの由来、期待値、許容誤差、既知制約をCIで確認し、別機関が再実行した"
    "コードだけを参照実装としてreleaseする。"
    "伊田氏とのHACApK・CEFCの国内共同実装、ADVENTUREとの接続、Hollaus氏・"
    "TU WienとのIGTE共著、TU GrazのopenCFSとNGSolveの知見を起点とする。"
    "国内外の学生と若手がコードレビューを通じて実装力と検証力を高める。"
    "企業のWindows設計環境とLinux計算環境、mdx、スパコンを同一仕様で接続し、"
    "GPU版ICCGをA100とH200で比較する。ARMを含む将来のアーキテクチャにも"
    "移植できる可搬な環境構築手順を整備する。Radiaは異なる手法を融合した先行実証とする。"
)


def test_grant_writing_kddi_health_report_runs():
    report = gw.grant_writing_health_report(KDDI_SAMPLE, program="kddi_digital")

    assert report["program"] == "kddi_digital"
    assert report["overall_score"] > 0
    assert "kddi_digital" in report["detailed_results"]
    assert "power_electronics_focus" in report["detailed_results"]
    assert "budget" in report["detailed_results"]


def test_grant_writing_kaken_oss_health_report_runs():
    report = gw.grant_writing_health_report(KAKEN_OSS_SAMPLE, program="kaken_oss")

    assert report["program"] == "kaken_oss"
    assert report["overall_score"] > 0
    assert "kaken_oss_platform" in report["detailed_results"]
    assert "budget" in report["detailed_results"]


def test_grant_writing_kaken_oss_platform_check():
    result = gw.grant_writing_kaken_oss_platform_check(KAKEN_OSS_SAMPLE)

    assert result["score"] >= 8
    assert result["axis_results"]["technical_reports_as_source"]["ok"]
    assert result["axis_results"]["jpmars_github_governance"]["ok"]
    assert result["axis_results"]["environment_portability"]["ok"]
    assert result["axis_results"]["lab_silo_and_ai_urgency"]["ok"]
    assert result["axis_results"]["reuse_and_upstream_first"]["ok"]
    assert result["axis_results"]["scientific_quality_gate"]["ok"]
    assert result["domestic_evidence_hits"]
    assert result["overseas_evidence_hits"]
    assert all(result["environment_results"].values())
    assert all(result["quality_results"].values())
    assert result["reuse_hits"]
    assert result["radia_integration"]["positioned_as_integration_evidence"]
    assert result["platform_focus"]["ok"]


def test_grant_writing_kaken_oss_accepts_ai_acceleration_variant():
    text = KAKEN_OSS_SAMPLE.replace(
        "AIによる開発加速が重複実装を増幅することを防ぐ。",
        "AIにより開発が加速する現在、重複実装を防ぐ。",
    )

    result = gw.grant_writing_kaken_oss_platform_check(text)

    assert "AIにより開発が加速" in result["ai_urgency_hits"]
    assert not any("why-now" in comment for comment in result["comments"])


def test_grant_writing_kaken_oss_accepts_negated_catch_up_frame():
    result = gw.grant_writing_kaken_oss_platform_check(
        KAKEN_OSS_SAMPLE + "本研究は欧州に追いつく計画ではなく、国内外の相互交流を作る。"
    )

    assert not any("catch-up" in comment for comment in result["comments"])


def test_grant_writing_kaken_oss_warns_on_hardware_only_plan():
    result = gw.grant_writing_kaken_oss_platform_check(
        "本研究ではGPUを購入し、計算機を購入する。"
    )

    assert result["score"] < 8
    assert any("hardware acquisition" in comment for comment in result["comments"])


def test_grant_writing_kaken_oss_warns_when_radia_hides_jpmars():
    result = gw.grant_writing_kaken_oss_platform_check(
        KAKEN_OSS_SAMPLE
        + "Radiaで実装し、Radiaを公開し、Radiaを研究基盤として運営する。"
    )

    assert not result["platform_focus"]["ok"]
    assert any("JP-MARs" in comment for comment in result["comments"])


def test_grant_writing_kaken_oss_rejects_repository_dump_mindset():
    result = gw.grant_writing_kaken_oss_platform_check(
        "新しい変換ソフトをOSSとしてGitHubに公開する。検証は利用者の自己責任とする。"
    )

    assert result["score"] < 8
    assert any("repository dump" in comment for comment in result["comments"])
    assert any("no-warranty" in comment for comment in result["comments"])
    assert any("upstream-first" in comment for comment in result["comments"])


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


def test_internal_evidence_to_external_scale_accepts_transfer_and_validation():
    result = gw.grant_writing_internal_evidence_to_external_scale_check(
        "予備成果として、研究室内で教員と学生が知識基盤を利用している。"
        "問題仕様、参照実装、入力データ、検証手順を公開し、"
        "他機関の第三者が再実行して結果を比較・検証する。"
    )

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["missing_axes"] == []


def test_internal_evidence_to_external_scale_warns_on_internal_only_success():
    result = gw.grant_writing_internal_evidence_to_external_scale_check(
        "予備成果として研究室内で運用している。"
    )

    assert result["applicable"]
    assert result["score"] < 8
    assert "transferable_unit" in result["missing_axes"]
    assert "external_actor" in result["missing_axes"]
    assert result["comments"]


def test_internal_evidence_to_external_scale_is_optional_without_prior_result():
    result = gw.grant_writing_internal_evidence_to_external_scale_check(
        "本研究では基礎理論を構築し、新しい数値解法の性質を明らかにする。"
    )

    assert not result["applicable"]
    assert result["score"] is None
    assert result["comments"] == []


def test_grant_writing_reexports_ja_lint_helpers():
    result = gw.grant_writing_lint_bedrock("これは重要であると考えられる。")

    assert isinstance(result, dict)
    assert "issue_count" in result


def test_sentence_analysis_ignores_latex_scaffolding():
    tex = r"""
% This deliberately long template comment must not be treated as proposal prose even when it exceeds the sentence threshold by a wide margin.
\section{研究目的}
\input{pieces/template_header}
本研究は技術報告と参照実装を接続する。第三者が結果を検証する。
"""

    result = gw.grant_writing_analyze_sentences(tex, max_len=50)

    assert result["total_sentences"] == 2
    assert result["over_threshold_count"] == 0
    assert all("template" not in item["head"] for item in result["over_threshold_examples"])


def test_document_meta_grant_domain_uses_public_grant_writing(tmp_path):
    draft = tmp_path / "grant.md"
    draft.write_text(KDDI_SAMPLE, encoding="utf-8")

    result = document_meta_lint_all(str(draft), domain="grant")

    assert "grant_writing_health_report" in result["lints_run"]
    assert result["lints_unavailable"] == []
    assert "grant lint stays" not in str(result)


def test_meta_catalog_lists_document_writing_pair_and_merged_presentation_figure():
    for key in ("paper-writing", "grant-writing"):
        assert key in CATALOG
        assert CATALOG[key]["entry_point"].startswith("mcp-server-")
    from radia_mcp.meta.catalog import _resolve
    # 2026-07-17: presentation merged into paper-writing (standalone
    # server retired); the old name must still resolve for discovery
    # and the tools must be advertised on the paper-writing entry.
    assert "presentation" not in CATALOG
    assert "presentation_usage" in CATALOG["paper-writing"]["primary_tools"]
    assert _resolve("presentation") == "paper-writing"
    # 2026-07-18: figure merged into paper-writing the same way.
    assert "figure" not in CATALOG
    assert "figure_style_guide" in CATALOG["paper-writing"]["primary_tools"]
    assert _resolve("figure") == "paper-writing"
    assert _resolve("mcp-server-figure") == "paper-writing"


def test_paper_writing_server_serves_merged_presentation_and_figure_tools():
    from radia_mcp.paper_writing import server as pw_server
    assert pw_server._N_PRESENTATION_TOOLS > 60
    assert pw_server._N_FIGURE_TOOLS > 5
