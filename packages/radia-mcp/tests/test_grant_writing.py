import asyncio
import json
import os
import sys
import zipfile
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
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
    "成果発表の旅費は国内2回、国際1回の単価と泊数から積算し、外国旅費へ計上する。"
    "外部費用は公式料金表URL、料金年度、参照日、税込区分、最低購入単位、有効期限、"
    "為替換算、端数処理を記録する。"
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

COLLABORATIVE_INTEGRATION_SAMPLE = (
    "本研究の学術的問いは、異種解析を結合して新しい設計知を導く条件は何かである。"
    "MCPは仮説を検証する実装手段とする。所有者側の自己記述・試験の初期整備と保守、"
    "利用側の作業を含む総負担を、中核ペアと別問題への再利用回数に対して評価する。"
    "既存基盤を置き換えず再利用し、その境界にある科学的自己記述を研究対象とする。"
    "中核実証の成立条件と、独立課題として試す発展候補を分ける。結合不能や反例も、"
    "原因と適用境界を同定できれば成果とする。各担当者には共著の既往成果、役割、"
    "直ちに着手できる準備がある。分析単位は結合課題であり、個人の能力を評価しない。"
    "手作業時間を記録する前に倫理該当性を確認する。資産の権利と保守主体を確認し、"
    "利用不能時は公開ベンチマークと参照実装を用いる。"
)

READABLE_JAPANESE_GRANT = (
    "磁気機器の設計では、解析条件の選択に長い時間を要する。"
    "本研究では、設計候補の順位が一致する条件を明らかにする。"
    "二つの解析法を同じ指標で比較し、適用範囲を判定する。"
    "予備試験では解析解との一致を確認した。"
    "これにより、設計者は必要な解析法を選択できる。"
)

DENSE_JAPANESE_GRANT = (
    "本研究ではHDiv-MMM、HCurl eddy-bubble、AGE、CLN、MCP及びAIを"
    "統合的に高度化することによって磁性導電体を含む複雑な電磁機器に"
    "対する高精度かつ高速で汎用的な解析基盤の構築及び社会実装を実現し"
    "さらに国内外の研究者や企業との連携を促進することでこれまで困難で"
    "あった設計最適化を可能にすることを目指す。"
)

RESEARCH_MEETING_MANUSCRIPT = (
    "本稿では、軸対称磁界解析に高次要素を適用した。"
    "提案法を解析解と比較した結果、磁束密度の相対誤差は一パーセント未満であった。"
    "図二にメッシュ収束を示す。"
    "以上から、曲面要素が境界近傍の誤差を低減することを確認した。"
)


def test_grant_writing_kddi_health_report_runs():
    report = gw.grant_writing_health_report(KDDI_SAMPLE, program="kddi_digital")

    assert report["program"] == "kddi_digital"
    assert 0 <= report["defect_score"] <= report["score_max"]
    assert report["tools_run"]
    assert "kddi_digital" in report["detailed_results"]
    assert "power_electronics_focus" in report["detailed_results"]
    assert "budget" in report["detailed_results"]


def test_grant_writing_kaken_oss_health_report_runs():
    report = gw.grant_writing_health_report(KAKEN_OSS_SAMPLE, program="kaken_oss")

    assert report["program"] == "kaken_oss"
    assert 0 <= report["defect_score"] <= report["score_max"]
    assert report["tools_run"]
    assert "kaken_oss_platform" in report["detailed_results"]
    assert "kaken_basic_research_positioning" in report["detailed_results"]
    assert "named_software_abstraction" in report["detailed_results"]
    assert "reviewer_vocabulary" in report["detailed_results"]
    assert "persuasion_quality" in report["detailed_results"]
    assert "adjacent_reviewer_readability" in report["detailed_results"]
    assert "reviewer_momentum" in report["detailed_results"]
    assert "literature_gap_evidence" in report["detailed_results"]
    assert "budget" in report["detailed_results"]


def test_grant_writing_server_exposes_adjacent_reviewer_readability():
    from radia_mcp.grant_writing.server import mcp

    tools = asyncio.run(mcp.list_tools())
    by_name = {tool.name: tool for tool in tools}
    names = set(by_name)

    assert "grant_writing_adjacent_reviewer_readability_check" in names
    assert "grant_writing_reviewer_momentum_check" in names
    assert "grant_writing_japanese_genre_contract" in names
    assert "grant_writing_japanese_readability_score" in names
    assert "grant_writing_kaken_basic_research_positioning_check" in names
    assert "grant_writing_budget_source_consistency_check" in names
    assert "bib_path" in by_name["grant_writing_publication_list"].inputSchema[
        "required"
    ]
    assert "document_type" in by_name[
        "grant_writing_japanese_readability_score"
    ].inputSchema["required"]
    assert mcp._mcp_server.instructions
    domain_tools = [
        tool
        for tool in tools
        if tool.name.startswith("grant_writing_")
        and not tool.name.endswith("_reload_code")
    ]
    assert all(tool.title for tool in tools)
    assert all(tool.annotations is not None for tool in tools)
    assert all(tool.annotations.readOnlyHint for tool in domain_tools)
    assert all(not tool.annotations.destructiveHint for tool in domain_tools)


def test_japanese_readability_scores_clear_grant_prose():
    result = gw.grant_writing_japanese_readability_score(
        READABLE_JAPANESE_GRANT,
        document_type="grant_proposal",
    )

    assert result["applicable"]
    assert result["status"] == "pass"
    assert result["score"] >= 90
    assert sum(
        axis["score_max"] for axis in result["scoring_axes"].values()
    ) == 100
    assert "English is neither scored nor averaged" in result["scoring_policy"]


def test_japanese_readability_rejects_long_method_inventory():
    result = gw.grant_writing_japanese_readability_score(
        DENSE_JAPANESE_GRANT,
        document_type="grant_proposal",
    )

    assert result["status"] == "fail"
    assert result["score"] < 70
    assert result["revision_priorities"][0]["axis"] == (
        "one_claim_sentence_rhythm"
    )
    lexical = result["scoring_axes"]["lexical_and_concept_load"]
    assert lexical["score"] < lexical["score_max"]


def test_japanese_readability_does_not_score_english_grant():
    result = gw.grant_writing_japanese_readability_score(
        "This proposal develops a readable method for magnetic design. "
        "It validates the method against measurements.",
        document_type="grant_proposal",
    )

    assert not result["applicable"]
    assert result["status"] == "not_applicable"
    assert result["score"] is None
    assert "Japanese grant prose only" in result["scoring_policy"]


def test_japanese_genre_contract_separates_grant_and_manuscript_review():
    grant = gw.grant_writing_japanese_genre_contract("助成金申請")
    manuscript = gw.grant_writing_japanese_genre_contract(
        "research_meeting_manuscript"
    )

    assert grant["status"] == "supported"
    assert grant["review_owner"] == "grant-writing"
    assert "reviewer_visible_problem_and_why_now" in (
        grant["grant_proposal_criteria"]
    )
    assert manuscript["status"] == "wrong_genre"
    assert manuscript["review_owner"] == "paper-writing"
    assert manuscript["route_to"]["server"] == "mcp-server-paper-writing"
    assert "result_figure_table_and_citation_traceability" in (
        manuscript["research_manuscript_criteria"]
    )
    assert set(grant["grant_proposal_criteria"]).isdisjoint(
        manuscript["research_manuscript_criteria"]
    )


def test_japanese_grant_score_rejects_research_meeting_manuscript_genre():
    result = gw.grant_writing_japanese_readability_score(
        RESEARCH_MEETING_MANUSCRIPT,
        document_type="research_meeting_manuscript",
    )

    assert not result["applicable"]
    assert result["status"] == "wrong_genre"
    assert result["score"] is None
    assert result["genre_contract"]["review_owner"] == "paper-writing"
    assert "paper_writing_bilingual_readability_check" in (
        result["genre_contract"]["route_to"]["tools"]
    )


def test_health_report_exposes_japanese_readability_without_merit_claim():
    result = gw.grant_writing_health_report(READABLE_JAPANESE_GRANT)

    assert result["japanese_readability_score"] >= 90
    assert result["japanese_readability_status"] == "pass"
    readability = result["detailed_results"]["japanese_readability"]
    assert "does not assess novelty" in readability["interpretation"]


async def _probe_japanese_readability_stdio() -> dict:
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.grant_writing.server"],
        cwd=str(package_root),
        env=env,
    )
    async with (
        stdio_client(params) as (read_stream, write_stream),
        ClientSession(read_stream, write_stream) as session,
    ):
        initialized = await session.initialize()
        listed = await session.list_tools()
        called = await session.call_tool(
            "grant_writing_japanese_readability_score",
            {
                "text": READABLE_JAPANESE_GRANT,
                "document_type": "grant_proposal",
            },
        )
        payload = json.loads(called.content[0].text)
        rejected = await session.call_tool(
            "grant_writing_japanese_readability_score",
            {
                "text": RESEARCH_MEETING_MANUSCRIPT,
                "document_type": "research_meeting_manuscript",
            },
        )
        rejected_payload = json.loads(rejected.content[0].text)
        return {
            "server_name": initialized.serverInfo.name,
            "listed": any(
                tool.name == "grant_writing_japanese_readability_score"
                for tool in listed.tools
            ),
            "is_error": bool(called.isError),
            "status": payload["status"],
            "score": payload["score"],
            "rejected_is_error": bool(rejected.isError),
            "rejected_status": rejected_payload["status"],
            "rejected_score": rejected_payload["score"],
        }


def test_japanese_readability_passes_real_stdio_protocol():
    result = asyncio.run(asyncio.wait_for(
        _probe_japanese_readability_stdio(),
        timeout=45,
    ))
    assert result["server_name"] == "mcp-server-grant-writing"
    assert result["listed"]
    assert not result["is_error"]
    assert result["status"] == "pass"
    assert result["score"] >= 90
    assert not result["rejected_is_error"]
    assert result["rejected_status"] == "wrong_genre"
    assert result["rejected_score"] is None


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


def test_kaken_basic_research_positioning_keeps_tools_below_question():
    result = gw.grant_writing_kaken_basic_research_positioning_check(
        "本研究の学術的問いは、異なる解析手法を結合しても設計判断を保存する"
        "成立条件と適用境界は何かである。MCPは各機関の技術を接続するglueという"
        "検証手段とし、GitHubは版と試験履歴を保持する。JMAGと実機測定による"
        "独立検証から選択則を得る。これにより工学へ波及し、国際競争力に貢献する。"
    )

    assert result["applicable"]
    assert result["risks"] == []
    assert result["evidence"]["academic_question"]
    assert result["evidence"]["tool_as_means"]


def test_kaken_basic_research_positioning_warns_on_tool_or_industry_as_goal():
    result = gw.grant_writing_kaken_basic_research_positioning_check(
        "本研究の目的はMCP付きOSS基盤をGitHubに構築し、自動車産業を強化することである。"
    )

    assert result["applicable"]
    assert {risk["type"] for risk in result["risks"]} >= {
        "tool_without_academic_question",
        "tool_role_not_subordinated",
        "engineering_context_without_hierarchy",
    }


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
    assert result["axis_results"]["pricing_provenance"]["ok"]
    assert "上限" in result["budget_policy"]
    assert "常時稼働" in result["budget_policy"]
    assert "短期集中" in result["budget_policy"]
    assert "ホストCPU" in result["budget_policy"]
    assert "設計判断" in result["budget_policy"]
    assert "参照日" in result["budget_policy"]
    assert "最低購入単位" in result["budget_policy"]
    assert "間接経費" in result["budget_policy"]
    assert "泊数" in result["budget_policy"]


def test_grant_writing_budget_alignment_requires_ceiling_and_calculation():
    result = gw.grant_writing_budget_alignment_check(
        "Claude、Codex、Fable、MDXの計算資源と基板試作、計測評価、発表旅費を予算化する。"
    )

    assert result["score"] < 8
    assert "near_ceiling_strategy" in result["missing_axes"]
    assert "itemized_calculation" in result["missing_axes"]
    assert "pricing_provenance" in result["missing_axes"]
    assert any("上限" in comment for comment in result["comments"])
    assert any("単価" in comment for comment in result["comments"])
    assert any("公式料金表" in comment for comment in result["comments"])


def test_grant_writing_budget_alignment_requires_price_provenance():
    result = gw.grant_writing_budget_alignment_check(
        "生成AIとGPUを予算化する。助成上限まで、単価、数量、月数、年度配分を積算する。"
    )

    assert "pricing_provenance" in result["missing_axes"]
    assert not result["axis_results"]["pricing_provenance"]["ok"]
    assert any("参照日" in comment for comment in result["comments"])


def test_budget_source_consistency_reconciles_rows_and_exact_totals(tmp_path):
    header = (
        "費目区分/Expenditure Categories,年度/FY,品名・仕様/Item (Specification),"
        "設置機関/Place,品目/Item,数量/Qty,単価/Unit Price,金額/Amount\n"
    )
    body = (
        "B,2027,,,SSD・HDD,,,50\n"
        "G,2027,,,JMAG Plan 3,,,1023\n"
        "G,2028,,,JMAG maintenance,,,165\n"
        "G,2029,,,JMAG maintenance,,,165\n"
        "G,2027,,,other,,,1221\n"
        "G,2028,,,other,,,1195\n"
        "G,2029,,,other,,,1151\n"
    )
    source = tmp_path / "budget.csv"
    comparison = tmp_path / "budget-copy.csv"
    source.write_text(header + body, encoding="utf-8-sig")
    comparison.write_text(header + body, encoding="utf-8-sig")

    result = gw.grant_writing_budget_source_consistency_check(
        str(source),
        comparison_source=str(comparison),
        expected_total_thousand_yen=4970,
        expected_year_totals_json='{"2027": 2294, "2028": 1360, "2029": 1316}',
        expected_category_totals_json='{"B": 50, "G": 4920}',
    )

    assert result["consistent"]
    assert result["canonical"]["row_count"] == 7
    assert result["comparison"]["missing_from_comparison"] == []


def test_budget_source_consistency_reports_exact_delta(tmp_path):
    source = tmp_path / "budget.csv"
    source.write_text(
        "category,年度/FY,x,x,品目/Item,x,x,金額/Amount\n"
        "G,2027,,,AI subscription,,,605\n",
        encoding="utf-8",
    )

    result = gw.grant_writing_budget_source_consistency_check(
        str(source), expected_total_thousand_yen=600
    )

    assert not result["consistent"]
    assert result["differences"][0]["delta"] == 5


def test_budget_source_consistency_ignores_trailing_amount_only_summary_rows(tmp_path):
    source = tmp_path / "budget.csv"
    source.write_text(
        "category,年度/FY,x,x,品目/Item,x,x,金額/Amount\n"
        "F,2027,,,AI subscription,,,202\n"
        "\n"
        ",,,,,,,202\n"
        ",,,,,,,202\n",
        encoding="utf-8",
    )
    result = gw.grant_writing_budget_source_consistency_check(
        str(source), expected_total_thousand_yen=202
    )
    assert result["consistent"]
    assert result["canonical"]["row_count"] == 1


def test_budget_source_consistency_keeps_amount_only_rows_inside_ledger_strict(tmp_path):
    source = tmp_path / "budget.csv"
    source.write_text(
        "category,年度/FY,x,x,品目/Item,x,x,金額/Amount\n"
        "F,2027,,,AI subscription,,,202\n"
        ",,,,,,,100\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="invalid budget row 3"):
        gw.grant_writing_budget_source_consistency_check(str(source))


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


def test_domain_outcome_chain_accepts_field_decision_rule():
    text = (
        "磁気浮上機器の平衡高さ、復元剛性、損失を設計量とする。"
        "候補間の性能差から、低費用解析と高忠実度解析の選択則と適用境界を求める。"
        "区間が分離すれば設計順位を確定し、包含条件を満たさなければ反証する。"
        "MCPとGitHubは問いではなく、この仮説を検証する手段である。"
    )

    result = gw.grant_writing_domain_outcome_chain_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["missing_axes"] == []


def test_domain_outcome_chain_rejects_platform_as_outcome():
    result = gw.grant_writing_domain_outcome_chain_check(
        "本研究ではMCPとGitHubを用いたOSS研究基盤を構築し、公開する。"
    )

    assert result["applicable"]
    assert result["score"] < 5
    assert "measurable_domain_quantity" in result["missing_axes"]
    assert "conditional_knowledge_product" in result["missing_axes"]
    assert "falsifiable_gate" in result["missing_axes"]


ARGUMENT_TRACEABLE = (
    "中心の問いは、異なる解析手法でも設計候補の順位を確定できる条件は何か、である。"
    "既往研究では個別手法が高度化した一方、手法間の比較条件は体系化されていない。"
    "二つの解析結果を同じ設計量へ射影し、候補順位を比較する。"
    "順位の一致率と許容差で成立条件を判定し、反例から適用境界を定める。"
    "成果として解析経路の選択則と設計指針を示す。"
    "予備実証では結合系を実装し、解析解との一致を確認した。"
    "この実績により、研究項目2の結合試験から着手できる。"
    "佐藤は最適化、比留間は行列解法、長嶺は数理理論、菅原は統合を担当する。"
    "結合が成立しない場合も、不能理由と適用境界を知識成果とする。"
)


def test_argument_evidence_map_indexes_each_role_without_scoring_it():
    result = gw.grant_writing_argument_evidence_map(ARGUMENT_TRACEABLE)

    assert result["applicable"]
    assert result["untraced_roles"] == []
    assert all(
        item["candidate_count"] > 0 for item in result["evidence_map"].values()
    )
    assert "score" not in result
    assert "risks" not in result


def test_argument_evidence_map_treats_absence_as_a_review_prompt_not_a_defect():
    result = gw.grant_writing_argument_evidence_map(
        "本研究では公開リポジトリと実行環境を整備する。"
        "複数機関がソースコードを共同編集する。"
    )

    assert "central_question" in result["untraced_roles"]
    assert "decision_rule" in result["untraced_roles"]
    assert "未検出は欠陥を意味しない" in result["manual_review_prompts"][0]


def test_argument_evidence_map_prompts_for_preparation_to_plan_link():
    result = gw.grant_writing_argument_evidence_map(
        "誘導加熱解析を実装し、解析解との一致を確認した。"
        "大学間予備実証を完了した。"
    )

    assert result["evidence_map"]["preliminary_evidence"]["candidate_count"] > 0
    assert result["evidence_map"]["preparation_plan_link"]["candidate_count"] == 0
    assert any(
        "研究項目への橋渡し文" in prompt
        for prompt in result["manual_review_prompts"]
    )


def test_health_report_exposes_argument_map_only_for_manual_review():
    result = gw.grant_writing_health_report(ARGUMENT_TRACEABLE)

    assert "argument_evidence_map" in result["detailed_results"]
    assert result["manual_review_prompts"]
    assert all(item["name"] != "argument_evidence_map" for item in result["findings"])


def test_health_report_rejects_unknown_program_instead_of_underchecking():
    with pytest.raises(ValueError, match="unknown grant program"):
        gw.grant_writing_health_report("本研究の目的を述べる。", program="kaken-oss")


def test_kaken_generic_uses_general_review_axes_without_oss_theme_check():
    text = (
        "本研究の学術的重要性と独創性を示す。研究方法と検証手順を定める。"
        "研究実績、実施体制および研究環境を整えている。"
        "国際共同研究を通じて世界の研究の発展に貢献する。"
    )
    report = gw.grant_writing_health_report(text, program="kaken_generic")

    assert report["program"] == "kaken_generic"
    assert report["detailed_results"]["sections"]["missing_axes"] == []
    assert "kaken_oss_platform" not in report["detailed_results"]


def test_kaken_review_axes_separates_main_axes_internationality_and_budget():
    result = gw.grant_writing_kaken_review_axes()

    assert result["scheme"] == "科研費 基盤研究(B・C)(一般)"
    assert len(result["research_plan_axes"]) == 3
    assert result["separate_scored_axis"]["id"] == "internationality"
    assert "同格" in result["budget_validity"]["role"]
    assert "複数" in result["budget_validity"]["consequence"]
    assert result["review_process"]["reviewers_for_scientific_research_c"] == 3
    assert len(result["sources"]) >= 4


def test_kaken_review_axes_records_budget_itemization_requirements():
    result = gw.grant_writing_kaken_review_axes()
    requirements = "".join(result["budget_entry_requirements"])

    assert "一式" in requirements
    assert "積算根拠" in requirements
    assert "90%" in requirements
    assert "研究代表者・研究分担者本人" in requirements


def test_health_report_rejects_unknown_skip_id():
    with pytest.raises(ValueError, match="unknown grant-writing skip"):
        gw.grant_writing_health_report("本研究の目的を述べる。", skip="sentences")


def test_derived_metric_validation_accepts_frozen_holdout_design():
    text = (
        "求解差、離散化差、連成差の合成則から設計判定区間を定義する。"
        "9校正ケースを校正集合とし、最終候補を得る前に式と安全係数を凍結する。"
        "校正に用いない保留データ15点で検証し、全点包含を合格条件とする。"
        "1点でも満たさない場合は再校正せず、適用境界として主張しない。"
    )

    result = gw.grant_writing_derived_metric_validation_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["missing_axes"] == []


def test_derived_metric_validation_flags_same_data_tuning():
    result = gw.grant_writing_derived_metric_validation_check(
        "残差から新しい設計判定区間を構成し、結果が合うよう安全係数を調整する。"
    )

    assert result["applicable"]
    assert result["score"] < 6
    assert "pretest_freeze" in result["missing_axes"]
    assert "heldout_validation" in result["missing_axes"]
    assert "failure_consequence" in result["missing_axes"]


def test_cross_organization_pilot_accepts_artifact_result_and_limit():
    text = (
        "大学間予備実証として、A大学提供の行列ソルバー実装をB大学が同一問題へ統合した。"
        "C大学提供のモデルで再実行し、168反復、真の残差8e-11で収束した。"
        "別機関レビューによる設計採否は未検証であり、本研究で実施する。"
    )

    result = gw.grant_writing_cross_organization_pilot_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["missing_axes"] == []
    assert result["evidence_level"].startswith("L2")
    assert not result["independent_review_or_adoption"]


def test_cross_organization_pilot_rejects_link_and_user_count_only():
    result = gw.grant_writing_cross_organization_pilot_check(
        "予備成果としてリポジトリを公開し、公式ページからリンクされた。"
        "研究室内の利用者が増えた。"
    )

    assert result["applicable"]
    assert result["score"] < 4
    assert "cross_organization_actor" in result["missing_axes"]
    assert "bounded_task" in result["missing_axes"]
    assert "observed_outcome" in result["missing_axes"]


def test_literature_gap_evidence_flags_bounded_corpus_generalization():
    text = (
        "電気機器の電磁界解析は細分化・高度化し、全てを一研究室で再実装できない。"
        "有限要素・数値解法・実機解析を継続的に扱う電気学会技術報告"
        "No. 1043・1317・1391・1471・1543・1547（2006–2023年）を調査対象とした。"
        "同一語彙で全文検索し、検出頁を目視確認した。ICCGは全6冊に現れた。"
        "一方、Maxwell問題向け補助空間前処理hypre AMSを確認できなかった。"
        "回転境界ではスライドメッシュ系が継承される一方、海外で発展した"
        "Air-Gap Elementは同一問題の比較候補にない。高次要素は以前から論じられる。"
        "しかしNo. 1543・1547の実機例は一次辺要素である。"
        "国内全実装の未利用や海外手法の優位を意味せず、実装言語、自由度対応、"
        "補助作用素、検証情報の不足により、科学的適合性を調べる前に候補から外れる"
        "未解決の研究障壁を示す。"
    )

    result = gw.grant_writing_literature_gap_evidence_check(text)

    assert result["applicable"]
    assert result["score"] < 5
    assert {risk["type"] for risk in result["risks"]} == {
        "field_generalization",
        "unsupported_causal_inference",
        "absence_as_academic_gap",
    }
    assert result["rewrite_strategy"]


def test_literature_gap_evidence_disclaimer_does_not_cure_same_inference():
    text = (
        "調査対象は6報告に限り、国内全体を断定しない。"
        "しかし全文検索でXを確認できなかったことは、国内では言語障壁により"
        "Xが普及していないことを示す。"
    )

    result = gw.grant_writing_literature_gap_evidence_check(text)

    assert result["applicable"]
    assert result["risk_count"] >= 2
    assert any(risk["type"] == "field_generalization" for risk in result["risks"])


def test_literature_gap_evidence_accepts_case_selection_only():
    text = (
        "対象10件ではXを確認できなかった。"
        "本調査は網羅調査ではなく、実証対象の選定にのみ用いる。"
        "学術的問いは、異種定式化の結合時に設計判断を保証する条件は何かである。"
    )

    result = gw.grant_writing_literature_gap_evidence_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["risks"] == []


def test_literature_gap_evidence_ignores_independent_literature_gap():
    text = (
        "既往研究は個々の手法の高速化を進めてきた。"
        "一方、異なる定式化間で設計判断の再現性を保証する条件は明らかでない。"
    )

    result = gw.grant_writing_literature_gap_evidence_check(text)

    assert not result["applicable"]
    assert result["score"] is None
    assert result["risks"] == []


def test_named_software_abstraction_flags_core_research_framing():
    text = r"""
\subsection{研究背景と学術的問い}
海外にはNGSolveの高次要素実装がある。
本研究の中心の問いは、ONELABで異種コードを結合できるかである。
"""

    result = gw.grant_writing_named_software_abstraction_check(text)

    assert result["applicable"]
    assert result["risk_count"] == 2
    assert {risk["software"] for risk in result["risks"]} == {
        "NGSolve",
        "ONELAB",
    }
    assert result["recommendations"]


def test_named_software_abstraction_accepts_category_level_framing():
    text = (
        "海外には高次要素を備えたOSS実装がある。"
        "中心の問いは、既存OSS連成基盤で異種モジュールを結合しても、"
        "設計判断を再現可能に保つ条件は何かである。"
    )

    result = gw.grant_writing_named_software_abstraction_check(text)

    assert not result["applicable"]
    assert result["score"] is None
    assert result["risks"] == []


def test_named_software_abstraction_allows_methods_and_prior_evidence():
    text = r"""
\subsection{研究方法と達成指標}
実証ではNGSolveの高次要素を用いて比較する。
\subsection{着想を支える準備状況}
代表者はONELABとの接続試験を実施し、国際会議で発表した。
"""

    result = gw.grant_writing_named_software_abstraction_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["risks"] == []
    assert {item["software"] for item in result["allowed_mentions"]} == {
        "NGSolve",
        "ONELAB",
    }


def test_reviewer_vocabulary_flags_unexplained_and_internal_terms():
    text = (
        "TU Grazとの共同研究でOSS、LLM、MCPを用いる。"
        "MMMとH(curl)を接続して解析する。"
    )

    result = gw.grant_writing_reviewer_vocabulary_check(text)

    assert result["applicable"]
    risk_types = {risk["type"] for risk in result["risks"]}
    assert "foreign_institution_alias" in risk_types
    assert "unexplained_acronym" in risk_types
    assert "domain_shorthand" in risk_types
    assert not result["term_results"]["OSS"]["explained"]
    assert not result["term_results"]["MCP"]["explained"]


def test_reviewer_vocabulary_accepts_japanese_role_first():
    text = (
        "代表者はグラーツ工科大学で共同研究を行った。"
        "オープンソースソフトウェア（OSS）を版管理し、"
        "AIが利用する知識・実行インターフェース（Model Context Protocol; MCP）"
        "へ写像する。大規模言語モデルは候補探索を支援する。"
        "磁気モーメント法と高次辺要素を比較する。"
    )

    result = gw.grant_writing_reviewer_vocabulary_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["risks"] == []
    assert result["term_results"]["OSS"]["explained"]
    assert result["term_results"]["MCP"]["explained"]


def test_reviewer_vocabulary_rejects_mcp_as_a_storage_location():
    result = gw.grant_writing_reviewer_vocabulary_check(
        "AIが利用する知識・実行インターフェース（MCP）を整備する。"
        "MCPには技術報告の実装判断・検証手順を蓄積している。"
    )

    assert any(
        risk["type"] == "mcp_described_as_storage"
        for risk in result["risks"]
    )


def test_reviewer_vocabulary_accepts_mcp_as_an_access_interface():
    result = gw.grant_writing_reviewer_vocabulary_check(
        "AIが利用する知識・実行インターフェース（MCP）を整備する。"
        "技術報告に基づく実装判断・検証手順をMCPサーバーから利用可能にする。"
    )

    assert not any(
        risk["type"] == "mcp_described_as_storage"
        for risk in result["risks"]
    )


def test_reviewer_vocabulary_rejects_benchmark_as_engineering_significance():
    text = "本研究の独創性はTEAM Problem 28を高精度に解くことである。"

    result = gw.grant_writing_reviewer_vocabulary_check(text)

    risk_types = {risk["type"] for risk in result["risks"]}
    assert "benchmark_unframed" in risk_types
    assert "benchmark_without_limit" in risk_types
    assert "benchmark_as_significance" in risk_types
    assert result["benchmark_result"]["used_as_significance"]


def test_reviewer_vocabulary_accepts_benchmark_as_initial_check_only():
    text = (
        "単純化した磁気浮上公開基準問題（TEAM Problem 28）は、"
        "解析経路の初期検証に限って用いる。"
        "結合前後の誤差を4%に収め、整合確認基準を得た。"
        "主実証では制約付き機器設計の設計判断を評価する。"
    )

    result = gw.grant_writing_reviewer_vocabulary_check(text)

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["risks"] == []
    assert result["benchmark_result"]["framed_as_verification"]
    assert result["benchmark_result"]["limitation_stated"]
    assert result["benchmark_result"]["engineering_value_distinguished"]
    assert not result["benchmark_result"]["negative_disclaimer_present"]


def test_reviewer_vocabulary_rejects_self_negating_benchmark_disclaimer():
    text = (
        "単純化した公開基準問題（TEAM Problem 28）は初期検証に用いる。"
        "この一致は工学的有用性の根拠にしない。"
    )

    result = gw.grant_writing_reviewer_vocabulary_check(text)

    risk_types = {risk["type"] for risk in result["risks"]}
    assert "benchmark_self_negating_disclaimer" in risk_types
    assert "benchmark_without_limit" in risk_types
    assert result["benchmark_result"]["negative_disclaimer_present"]


def test_persuasion_quality_flags_self_negating_evidence():
    result = gw.grant_writing_persuasion_quality_check(
        "公開基準問題で誤差を4%に収めた。"
        "この一致は工学的有用性の根拠にしない。"
    )

    assert result["score"] < 10
    assert any(
        risk["type"] == "self_negating_evidence" for risk in result["risks"]
    )


def test_persuasion_quality_accepts_positive_evidence_bridge():
    result = gw.grant_writing_persuasion_quality_check(
        "公開基準問題で誤差を4%に収め、結合前後の整合確認基準を得た。"
        "主実証では制約付き機器設計に移し、設計順位が保存される条件を検証する。"
    )

    assert result["score"] == 10.0
    assert result["risks"] == []


def test_persuasion_quality_flags_equation_with_unintroduced_symbols():
    text = r"""
設計判定区間を次式で定義する。
\[
\Delta_q(d)=\gamma_s|q_s-q_s^+|+\gamma_h|q_{h,p}-q_{h/2,p+1}|
+\gamma_c|q_{1way}-q_{iter}|
\]
この区間で候補順位を判断する。
"""

    result = gw.grant_writing_persuasion_quality_check(text)

    symbol_risks = [
        risk
        for risk in result["risks"]
        if risk["type"] == "equation_symbols_not_introduced"
    ]
    assert symbol_risks
    assert "q_s" in symbol_risks[0]["missing_symbols"]
    assert "q_{1way}" in symbol_risks[0]["missing_symbols"]


def test_persuasion_quality_accepts_equation_on_ramp_and_interpretation():
    text = r"""
候補$d$の低費用解析値を$\hat q_d$とする。停止条件を強めた変化を
求解差$\delta_s$、高次化による変化を離散化差$\delta_h$、反復連成による
変化を連成差$\delta_c$とする。値が動き得る幅$\Delta_q$と判定区間$I_q$を
\[
\Delta_q(d)=\gamma_s\delta_s(d)+\gamma_h\delta_h(d)+\gamma_c\delta_c(d),
\qquad I_q(d)=[\hat q_d-\Delta_q(d),\hat q_d+\Delta_q(d)]
\]
と表す。区間が分離すれば順位を確定し、重なれば高忠実度解析へ進む。
$\gamma_s,\gamma_h,\gamma_c$は校正で定める安全係数である。
"""

    result = gw.grant_writing_persuasion_quality_check(text)

    assert result["score"] == 10.0
    assert result["counts"]["display_equations"] == 1
    assert result["risks"] == []


def test_persuasion_quality_flags_inline_condition_before_physical_meaning():
    result = gw.grant_writing_persuasion_quality_check(
        r"加速器電磁石では、$\eta_{\rm out}=0$を満たす形状を設計する。"
    )

    risk = next(
        risk
        for risk in result["risks"]
        if risk["type"] == "inline_condition_before_physical_meaning"
    )
    assert risk["symbol"] == r"\eta_{\rm out}"
    assert result["counts"]["unexplained_inline_conditions"] == 1


def test_persuasion_quality_accepts_explained_inline_condition():
    text = (
        "粒子の運動量が異なっても出口で同じ位置へ戻す。"
        r"運動量差に対する軌道位置の変化率を分散関数$\eta$と呼ぶ。"
        r"目標は出口分散$\eta_{\rm out}=0$である。"
    )

    result = gw.grant_writing_persuasion_quality_check(text)

    assert result["score"] == 10.0
    assert result["counts"]["unexplained_inline_conditions"] == 0
    assert result["risks"] == []


def test_persuasion_quality_ignores_latex_linebreak_spacing():
    text = r"""
\begin{center}
図を配置する。
\\[-0.5zw]{\small 図1　研究構想}
\end{center}
"""

    result = gw.grant_writing_persuasion_quality_check(text)

    assert result["score"] == 10.0
    assert result["counts"]["display_equations"] == 0
    assert result["risks"] == []


def test_sentence_analysis_does_not_count_display_equation_source():
    text = r"""
設計判定区間を定義する。
\[
I_q(d)=[\hat q_d-\gamma_s\delta_s(d)-\gamma_h\delta_h(d)-
\gamma_c\delta_c(d),\hat q_d+\gamma_s\delta_s(d)+
\gamma_h\delta_h(d)+\gamma_c\delta_c(d)]
\]
区間が分離すれば順位を確定する。
"""

    result = gw.grant_writing_analyze_sentences(text)

    assert result["over_threshold_count"] == 0


def test_adjacent_reviewer_readability_flags_short_dense_layer_mix():
    text = (
        "異種解析モジュールの機能、入出力物理量、実行条件、判定区間を"
        "MCP、GitHub、CIで版管理し、設計候補の順位確定と第三者反証へ対応付ける。"
        "XFEM、HACApK、ESIMをCauer縮約へ接続し、損失と効率を判定する。"
        "再現、反証、判定、採否、凍結、版管理、許容差、適用限界を記録する。"
        "別機関が同じ結論を再構成し、反証結果を判定する。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)

    assert result["applicable"]
    assert result["score"] is None
    assert result["risk_count"] >= 2
    types = {risk["type"] for risk in result["risks"]}
    assert "compressed_concept_density" in types
    assert "notation_or_method_pile" in types
    assert "three_layer_paragraph" in types


def test_adjacent_reviewer_readability_accepts_plain_causal_sequence():
    text = (
        "誘導加熱では、コイルに流す電流が品物を発熱させる。"
        "しかし、解析手法を替えると発熱量の予測が変わる。"
        "そこで二つの予測を実測値と比べ、設計順位を保てる条件を明らかにする。"
        "実行履歴は最後に保存し、別機関が結果を確認できるようにする。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)

    assert result["applicable"]
    assert result["risk_count"] == 0


def test_adjacent_reviewer_readability_flags_section_claim_compression():
    text = (
        "独創性は、機器設計則を、校正・保留検証済み区間として与える点にある。"
        "さらに異なるコード系譜間で判断を再現する条件を示す。"
        "必達範囲は誘導加熱と加速器電磁石の二課題とする。"
        "本研究は一方へ統一せず、同一設計量で採否する。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "result_representation_type_mismatch" in types
    assert "ambiguous_relation_or_decision_object" in types
    assert "required_scope_without_deliverable" in types


def test_adjacent_reviewer_readability_accepts_explicit_section_claims():
    text = (
        "必要な解析忠実度を選び、候補順位を確定できる条件を示す。"
        "開発母体と内部形式が異なる解析コードを結合する。"
        "二課題で結合条件を実証することを必達範囲とする。"
        "各手法を単一の内部形式へ統一せず、適用の可否を判断する。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "result_representation_type_mismatch" not in types
    assert "ambiguous_relation_or_decision_object" not in types
    assert "required_scope_without_deliverable" not in types


def test_adjacent_reviewer_readability_flags_applicant_internal_metaphors():
    text = (
        "比較対象が自研究室のコード系譜に偏る。"
        "本研究では四名の資産を誘導加熱へ結び、判定則を加速器設計へ移す。"
        "助成期間前に共同実装を進める準備が整っている。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "applicant_internal_abstraction" in types
    assert "vague_readiness_status" in types
    assert result["metrics"]["applicant_internal_abstraction_count"] == 2
    assert result["metrics"]["vague_readiness_status_count"] == 1


def test_adjacent_reviewer_readability_accepts_concrete_keiko_revision():
    text = (
        "比較できる対象が自研究室の既存コードで扱える手法に限られやすい。"
        "本研究では、誘導加熱と加速器電磁石の二課題について検証する。"
        "四名の資産を誘導加熱へ展開し、判定則を加速器設計へ発展する。"
        "以降は、3人の分担者の役割を説明する。"
        "本助成期間前に共同実装を進める準備を整備済みである。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "applicant_internal_abstraction" not in types
    assert "vague_readiness_status" not in types


def test_adjacent_reviewer_readability_flags_implementation_lineage_jargon():
    text = (
        "二つのプログラムの系譜を比較する。"
        "ソースコード譜系の違いを検証する。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    risks = {risk["type"]: risk for risk in result["risks"]}

    assert "applicant_internal_abstraction" in risks
    assert result["metrics"]["applicant_internal_abstraction_count"] == 2
    recommendation = risks["applicant_internal_abstraction"]["recommendation"]
    assert "既存プログラム" in recommendation
    assert "ソースコード" in recommendation


def test_adjacent_reviewer_readability_accepts_concrete_implementation_terms():
    text = (
        "比較できる対象は既存プログラムで扱える手法に限られる。"
        "独立に開発された二つのソースコードで同じ問題を解析する。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "applicant_internal_abstraction" not in types


def test_adjacent_reviewer_flags_abstract_feasibility_evidence_labels():
    text = (
        "遂行可能性は、共著実績、複数機関資産の再実行、担当者の公開変更、"
        "研究室内MCP運用に基づく。"
    )
    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    risks = {risk["type"]: risk for risk in result["risks"]}
    assert "feasibility_evidence_not_observable" in risks
    assert result["metrics"]["feasibility_evidence_not_observable_count"] == 1
    assert "実施者、対象、完了した操作" in risks[
        "feasibility_evidence_not_observable"
    ]["recommendation"]


def test_adjacent_reviewer_accepts_observable_feasibility_evidence():
    text = (
        "遂行可能性は、共著実績、複数機関資産の再実行、担当者のビルド・CI"
        "改善4件を公開版へ取り込んだ実績、研究室内MCP接続・実行実績に基づく。"
    )
    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}
    assert "feasibility_evidence_not_observable" not in types


def test_adjacent_reviewer_readability_flags_takeaway_after_evidence():
    text = (
        "菅原・長嶺らは、Cauer縮約をGalerkin系へ実装した。"
        "解析解に対する誤差が1%未満であることを確認し、IGTE 2026で発表した。"
        "これらの誘導加熱実績により、研究項目2は実装済みの結合系から開始できる。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    risks = {risk["type"]: risk for risk in result["risks"]}

    assert "takeaway_after_evidence" in risks
    assert risks["takeaway_after_evidence"]["severity"] == "HIGH"
    assert risks["takeaway_after_evidence"]["rewrite_order"][0] == (
        "reviewer_takeaway"
    )
    assert result["metrics"]["takeaway_after_evidence_count"] == 1
    assert result["revision_protocol"]["sequence"] == [
        "reviewer_takeaway",
        "plain_language_role",
        "specific_method_or_evidence",
        "limit_or_remaining_question",
    ]


def test_adjacent_reviewer_readability_accepts_takeaway_first_evidence():
    text = (
        "誘導加熱の結合系は実装・検証済みであり、研究項目2を直ちに開始できる。"
        "具体的には、Cauer縮約をGalerkin系へ実装した。"
        "解析解に対する誤差が1%未満であることを確認し、IGTE 2026で発表した。"
    )

    result = gw.grant_writing_adjacent_reviewer_readability_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "takeaway_after_evidence" not in types
    assert result["metrics"]["takeaway_after_evidence_count"] == 0


def test_reviewer_momentum_accepts_concrete_problem_to_payoff():
    text = (
        "粒子線がん治療では、患者のがん形状に合わせた高精度な照射が求められる。"
        "しかし、磁場をそろえる初期化運転には半年を要し、治療装置の運用を"
        "圧迫している。そこで本研究では、磁気ヒステリシスを定量化して磁場を"
        "制御する。初期化運転を行わず、設定電流を即座に算出できる技術を確立する。"
    )

    result = gw.grant_writing_reviewer_momentum_check(text)

    assert result["applicable"]
    assert result["score"] is None
    assert result["risk_count"] == 0
    assert result["metrics"]["arc_complete"]


def test_reviewer_momentum_accepts_unused_opportunity_arc():
    text = (
        "金属3Dプリンタの実用化により、加熱コイルの設計自由度が急速に拡大した。"
        "しかし、任意形状の自動メッシュ生成は60年間実現されておらず、設計は"
        "熟練者の試行錯誤に依存する。そこで本研究では、AIエージェントと"
        "電磁界・熱解析を結ぶ。製造自由度を使い切るコイル設計を自動化し、"
        "設計期間の短縮と加熱効率の向上を実現する。"
    )

    result = gw.grant_writing_reviewer_momentum_check(text)

    assert result["risk_count"] == 0
    assert result["metrics"]["arc_complete"]
    assert result["metrics"]["productive_tension"]


def test_reviewer_momentum_flags_method_first_inventory():
    text = (
        "FEM、RNA、CLN、MCP、OSS、APIを統合するモデルを構築する。"
        "本研究では解析モジュールを接続する。"
        "これにより設計条件を示す。"
    )

    result = gw.grant_writing_reviewer_momentum_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "reviewer_stakes_delayed" in types
    assert "bottleneck_not_concrete" in types
    assert "method_before_problem" in types
    assert "lead_method_inventory" in types


def test_reviewer_momentum_flags_unsupported_hype():
    text = (
        "本研究は世界初の画期的で革新的な解析基盤を構築する。"
        "新しいモデルを開発する。"
        "設計へ応用する。"
    )

    result = gw.grant_writing_reviewer_momentum_check(text)
    types = {risk["type"] for risk in result["risks"]}

    assert "unsupported_excitement_language" in types


def test_subject_predicate_distance_reads_fullwidth_japanese_comma():
    text = (
        "本研究は，異なる手法で得た多数の候補について設計量を比較し，"
        "第三者が同じ順位を再現できる条件を明らかにする。"
    )

    result = gw.grant_writing_check_subject_predicate_distance(text, max_chars=20)

    assert result["analyzed_sentences"] == 1
    assert result["violation_count"] == 1


def test_sentence_analysis_does_not_join_figure_caption_to_aims():
    text = r"""
\begin{center}
\includegraphics[width=0.8\linewidth]{concept.png}
\\[-0.5zw]{\small\textbf{図1　研究機関別の解析資産を結合する構想}}
\end{center}
\subsection{研究目的}
目的は、磁気浮上機構の設計判断を再現する条件を示すことである。
"""

    result = gw.grant_writing_analyze_sentences(text, max_len=50)

    assert result["over_threshold_count"] == 0


def test_persuasion_quality_reads_cp932_tex(tmp_path):
    path = tmp_path / "proposal.tex"
    path.write_bytes(
        (
            "公開基準問題で誤差を4%に収めた。"
            "この一致は工学的有用性の根拠にしない。"
        ).encode("cp932")
    )

    result = gw.grant_writing_persuasion_quality_check(str(path))

    assert any(
        risk["type"] == "self_negating_evidence" for risk in result["risks"]
    )


def test_persuasion_quality_flags_equation_without_decision_on_ramp():
    text = r"""
候補について計算する。
\[
x_a=x_b+x_c
\]
以上である。
"""

    result = gw.grant_writing_persuasion_quality_check(text)

    risk_types = {risk["type"] for risk in result["risks"]}
    assert "equation_without_on_ramp" in risk_types
    assert "equation_without_interpretation" in risk_types


def test_persuasion_quality_flags_defensive_core_and_optional_branch():
    text = (
        "本研究は統一ソルバーを目的としない。規格策定を前提にしない。"
        "商用コードは含めない。単なる基盤整備ではない。"
        "条件付き追加検証は必達成果には含めない。"
    )

    result = gw.grant_writing_persuasion_quality_check(text)

    risk_types = {risk["type"] for risk in result["risks"]}
    assert "defensive_paragraph" in risk_types
    assert "optional_branch_in_core_plan" in risk_types


def test_persuasion_quality_flags_internal_memo_shorthand():
    result = gw.grant_writing_persuasion_quality_check(
        "実証Aと実証BをL2比較する。A・Bの仕様を固定する。"
        "国際会議で報告する（年度末: 結合条件を確定）。"
    )

    memo_risks = [
        risk
        for risk in result["risks"]
        if risk["type"] == "internal_memo_shorthand"
    ]
    assert {risk["memo_type"] for risk in memo_risks} == {
        "lettered_experiment",
        "letter_pair",
        "coded_stage",
        "parenthetical_milestone",
    }
    assert result["counts"]["memo_shorthand"] == 5


def test_persuasion_quality_allows_formal_numbered_research_items():
    result = gw.grant_writing_persuasion_quality_check(
        "研究項目1では磁気浮上設計を検証する。"
        "研究項目2では行列解法の選択条件を求める。"
    )

    assert result["score"] == 10.0
    assert result["risks"] == []


def test_persuasion_quality_flags_acronym_pile():
    result = gw.grant_writing_persuasion_quality_check(
        "OSS、MCP、LLM、GPU、CPU、XFEMをGitHubで接続して検証する。"
    )

    acronym_risk = next(
        risk for risk in result["risks"] if risk["type"] == "acronym_pile"
    )
    assert set(acronym_risk["acronyms"]) >= {
        "OSS",
        "MCP",
        "LLM",
        "GPU",
        "CPU",
        "XFEM",
    }


def test_collaborative_integration_risk_check_accepts_complete_plan():
    result = gw.grant_writing_collaborative_integration_risk_check(
        COLLABORATIVE_INTEGRATION_SAMPLE
    )

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["missing_axes"] == []
    assert all(axis["ok"] for axis in result["axis_results"].values())


def test_collaborative_integration_risk_check_exposes_hidden_costs_and_boundaries():
    result = gw.grant_writing_collaborative_integration_risk_check(
        "MCPで複数研究室のソフトウェアを統合して効率化する。"
    )

    assert result["applicable"]
    assert result["score"] < 5
    assert "provider_and_reuse_cost" in result["missing_axes"]
    assert "existing_ecosystem_boundary" in result["missing_axes"]
    assert "negative_result_value" in result["missing_axes"]
    assert result["comments"]


def test_collaborative_integration_risk_check_is_optional_for_unrelated_plan():
    result = gw.grant_writing_collaborative_integration_risk_check(
        "本研究では新しい磁性材料の物性を測定する。"
    )

    assert not result["applicable"]
    assert result["score"] is None


def test_kaken_review_format_flags_color_only_figure():
    result = gw.grant_writing_kaken_review_format_check(
        "図1に赤線で提案法、青線で従来法の損失を示す。"
    )

    assert result["applicable"]
    risk = next(
        r for r in result["risks"] if r["type"] == "color_dependent_figure"
    )
    assert risk["severity"] == "HIGH"
    assert result["score"] < 10


def test_kaken_review_format_accepts_monochrome_safe_figure():
    result = gw.grant_writing_kaken_review_format_check(
        "図1では提案法を実線、従来法を破線の線種で区別し、"
        "白黒印刷でも判別できるようにする。"
    )

    assert not any(
        r["type"] == "color_dependent_figure" for r in result["risks"]
    )


def test_kaken_review_format_requires_safeguards_for_human_subjects():
    result = gw.grant_writing_kaken_review_format_check(
        "設計者へのアンケート調査で使い勝手を評価する。"
    )

    risk = next(
        r
        for r in result["risks"]
        if r["type"] == "human_subjects_without_safeguard"
    )
    assert risk["severity"] == "HIGH"


def test_kaken_review_format_accepts_safeguarded_survey():
    result = gw.grant_writing_kaken_review_format_check(
        "アンケート調査は倫理委員会の承認を得て、同意を取得し匿名化して実施する。"
    )

    assert not any(
        r["type"] == "human_subjects_without_safeguard" for r in result["risks"]
    )


def test_kaken_review_format_ignores_wrapped_ethics_form_example():
    # Extracted old PDF forms wrap one instruction across lines, removing the
    # polite ending from the fragments that carry the trigger words.
    result = gw.grant_writing_kaken_review_format_check(
        "例えば、個人情報を伴う、アンケート調査・インタビュー調査、\n"
        "提供を受けた試料の使用、ヒト遺伝子解析研究など、\n"
        "承認手続が必要となる調査・研究・実験などが対象となります。\n"
        "本研究は公開済みの数値データだけを解析する。"
    )

    assert not any(
        r["type"] == "human_subjects_without_safeguard" for r in result["risks"]
    )


def test_kaken_review_format_ignores_old_form_and_admin_survey_labels():
    result = gw.grant_writing_kaken_review_format_check(
        "個人情報を伴うアンケート調査・インタビュー調査、提供を受けた試料の使用\n"
        "パワーアカデミー研究助成に関するアンケート\n"
        "本研究は公開済みの数値データだけを解析する。"
    )

    assert not any(
        r["type"] == "human_subjects_without_safeguard" for r in result["risks"]
    )


def test_kaken_review_format_box_heading_is_not_a_safeguard():
    # The box heading itself contains 「遵守」; quoting it must not
    # suppress the missing-safeguard check.
    result = gw.grant_writing_kaken_review_format_check(
        "人権の保護及び法令等の遵守への対応: 設計者へのアンケート調査を行う。"
    )

    assert any(
        r["type"] == "human_subjects_without_safeguard" for r in result["risks"]
    )


def test_kaken_review_format_wants_rationale_next_to_not_applicable():
    result = gw.grant_writing_kaken_review_format_check(
        "人権の保護及び法令等の遵守への対応: 該当なし。"
    )

    assert any(
        r["type"] == "not_applicable_without_rationale" for r in result["risks"]
    )


def test_kaken_review_format_accepts_rationale_with_not_applicable():
    result = gw.grant_writing_kaken_review_format_check(
        "本研究は数値解析のみで人や動物を対象としないため、該当なし。"
    )

    assert not any(
        r["type"] == "not_applicable_without_rationale" for r in result["risks"]
    )


def test_kaken_review_format_ignores_not_applicable_in_an_unrelated_box():
    result = gw.grant_writing_kaken_review_format_check(
        "研究計画と進捗評価を受けた研究課題の関連性: 該当なし。"
    )

    assert not any(
        r["type"] == "not_applicable_without_rationale" for r in result["risks"]
    )


def test_kaken_review_format_flags_final_year_not_applicable_text():
    result = gw.grant_writing_kaken_review_format_check(
        r"""
        \section{研究計画最終年度前年度応募を行う場合の記述事項}
        \newcommand{\最終年度研究種目名}{該当なし}
        \newcommand{\最終年度研究課題番号}{該当なし}
        本応募は研究計画最終年度前年度応募には該当しない。
        """
    )

    risks = {risk["type"]: risk for risk in result["risks"]}

    assert "final_year_non_applicant_field_not_blank" in risks
    assert risks["final_year_non_applicant_field_not_blank"]["severity"] == "HIGH"
    assert result["final_year_blank_rule_checked"]
    assert len(result["final_year_blank_violations"]) == 3


def test_kaken_review_format_accepts_blank_final_year_fields():
    result = gw.grant_writing_kaken_review_format_check(
        r"""
        \section{研究計画最終年度前年度応募を行う場合の記述事項}
        \newcommand{\最終年度研究種目名}{}
        \newcommand{\最終年度研究課題番号}{}
        \newcommand{\最終年度研究課題名}{}
        \newcommand{\最終年度研究期間}{}
        \textbf{当初研究計画及び研究成果}
        \textbf{前年度応募する理由}
        """
    )

    assert not any(
        risk["type"] == "final_year_non_applicant_field_not_blank"
        for risk in result["risks"]
    )
    assert result["final_year_blank_rule_checked"]
    assert result["final_year_blank_violations"] == []


def test_kaken_review_format_ignores_final_year_form_instruction():
    result = gw.grant_writing_kaken_review_format_check(
        "研究計画最終年度前年度応募を行う場合の記述事項。"
        "該当しない場合は記述欄を削除することなく、空欄のまま提出すること。"
    )

    assert not any(
        risk["type"] == "final_year_non_applicant_field_not_blank"
        for risk in result["risks"]
    )


def test_kaken_review_format_flags_unidentifiable_publications():
    result = gw.grant_writing_kaken_review_format_check(
        "研究遂行能力: 関連する主要論文があり、実行可能性は高い。"
    )

    assert any(
        r["type"] == "publication_not_identifiable" for r in result["risks"]
    )


def test_kaken_review_format_accepts_identifiable_publications():
    result = gw.grant_writing_kaken_review_format_check(
        "研究遂行能力: 代表論文はIEEE Trans. Magn., vol. 54, 2018に掲載された。"
    )

    assert not any(
        r["type"] == "publication_not_identifiable" for r in result["risks"]
    )


def test_kaken_review_format_checks_funding_overlap_box():
    result = gw.grant_writing_kaken_review_format_check(
        "応募中の研究費: 基盤研究C(代表)。"
    )

    risk = next(
        r for r in result["risks"] if r["type"] == "funding_overlap_format"
    )
    assert any("相違点" in part for part in risk["missing_parts"])
    assert any("役職" in part for part in risk["missing_parts"])


def test_kaken_review_format_accepts_complete_funding_overlap_box():
    result = gw.grant_writing_kaken_review_format_check(
        "応募中の研究費: 挑戦的研究(萌芽)。本応募課題との相違点は解析対象で、"
        "応募する理由は検証装置の整備である。研究代表者は近畿大学教授。"
    )

    assert not any(
        r["type"] == "funding_overlap_format" for r in result["risks"]
    )


def test_kaken_review_format_carries_briefing_notes():
    result = gw.grant_writing_kaken_review_format_check("本研究の目的を述べる。")

    assert result["score"] == 10.0
    assert any("100件" in note for note in result["briefing_notes"])
    assert any("充足率" in note for note in result["briefing_notes"])
    assert any("白黒" in note for note in result["briefing_notes"])


def test_kaken_review_format_skips_full_draft_heuristics_on_fragments():
    result = gw.grant_writing_kaken_review_format_check(
        "本研究の目的を述べる。"
    )

    assert not result["full_draft_heuristics_applied"]
    assert result["criteria_axis_results"] == {}


def test_kaken_review_format_runs_in_health_report():
    report = gw.grant_writing_health_report(KDDI_SAMPLE, program="kddi_digital")

    assert "kaken_review_format" in report["detailed_results"]
    assert "kaken_review_format" in report["detailed_scores"]


ORIGINALITY_QUESTION_DRAFT = (
    "本研究の中心の問いは、異種の電磁界解析モジュールを内部形式のまま結合したとき、"
    "設計候補の順位を確定できる条件をどのように定量化できるか、である。"
    "既往研究では、個々の手法の高精度化・高速化が進められてきた。"
    "一方、手法間の差が設計量へ及ぼす影響を定量化する枠組みは体系化されていない。"
    "本研究の独自性は、精度を一律に高めず性能差に応じて解析経路を選ぶ点にある。"
)


def test_question_originality_accepts_a_contrastive_pair():
    # Real proposals split the contrast across two sentences joined by 一方.
    # Requiring the single-sentence form would fail correct Japanese.
    result = gw.grant_writing_question_originality_check(
        ORIGINALITY_QUESTION_DRAFT
    )

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["gap_statements"][0]["form"] == "contrastive_pair"
    assert result["originality_markers"]


def test_question_originality_accepts_a_single_sentence_contrast():
    result = gw.grant_writing_question_originality_check(
        "本研究の中心の問いは、異種解析を結合したとき設計候補の順位を確定できる"
        "条件をどのように定量化できるか、である。"
        "既往研究は個々の手法の高速化を進めてきたが、その差を設計量へ伝播させる"
        "方法は確立していない。本研究の新規性はこの伝播則にある。"
    )

    assert result["score"] == 10.0
    assert result["gap_statements"][0]["form"] == "single_sentence"


def test_question_originality_flags_a_question_with_no_position():
    # The review item three of five reviewers marked down: the question is
    # stated, but nothing says what is new or what prior work leaves open.
    result = gw.grant_writing_question_originality_check(
        "本研究の中心の問いは、異種の電磁界解析モジュールを内部形式のまま結合"
        "したとき、設計候補の順位を確定できる条件をどのように定量化できるか、"
        "である。誘導加熱と加速器電磁石の二課題で検証し、結果を公開する。"
    )

    assert result["applicable"]
    types = {r["type"] for r in result["risks"]}
    assert "no_originality_claim" in types
    assert "no_gap_against_prior_work" in types
    assert result["score"] < 6


def test_question_originality_is_not_applicable_without_a_question():
    result = gw.grant_writing_question_originality_check(
        "誘導加熱の発熱量を評価し、損失分布を求める。"
    )

    assert not result["applicable"]
    assert result["score"] is None


def test_template_residue_check_flags_unfilled_placeholders():
    # Measured on a real 2026 draft: the money boxes still read ○○○○千円 when
    # it reached the co-investigator. A form office sends this back before any
    # reviewer sees it.
    result = gw.grant_writing_template_residue_check(
        "(A) 設備備品費　…　小計：○○○○千円（税込）\n"
        "(1) 氏名：○○ ○○\n"
        "研究の目的を平易に述べる。"
    )

    assert result["applicable"]
    assert result["risk_count"] == 2
    assert all(r["type"] == "unfilled_placeholder" for r in result["risks"])
    assert result["score"] < 10


def test_template_residue_check_counts_instructions_without_calling_them_defects():
    # The rule that failed the detector test: an experienced PI deleted 6
    # instruction sentences and deliberately kept 13, and flat text cannot
    # tell the two groups apart. So they are counted and handed back as a
    # question, never reported as located defects.
    result = gw.grant_writing_template_residue_check(
        "パワーアカデミー研究マップをご覧いただき、課題を選択してください。\n"
        "性別を回答したくない方は「該当しない」を選択してください。\n"
        "本研究の目的は損失低減である。"
    )

    assert result["risk_count"] == 0
    assert result["instruction_sentence_count"] == 2
    assert result["questions"]
    assert result["score"] == 10.0


def test_template_residue_check_is_clean_on_a_finished_draft():
    result = gw.grant_writing_template_residue_check(
        "設備備品費は1,200千円である。氏名は菅原賢悟である。"
    )

    assert result["risk_count"] == 0
    assert result["score"] == 10.0


def test_health_report_separates_findings_from_questions():
    report = gw.grant_writing_health_report(KDDI_SAMPLE, program="kddi_digital")

    assert "findings" in report and "questions" in report
    assert "defect_counts" in report
    assert report["defect_counts"]["total"] == len(report["findings"])
    # Questions carry no severity and no score: they are prompts, not defects.
    for q in report["questions"]:
        assert q["kind"] == "question"
        assert "severity" not in q
        assert "score" not in q
    # The score reflects located defects only.
    assert "defect_score" in report


def test_vague_claim_verb_check_flags_integration_without_a_mechanism():
    # The wording an experienced PI removes first: it promises a result and
    # names no operation. Verified against a real 2026 proposal rewrite where
    # 「統合し」 went 3 -> 0 and 「双方向に連成」 went 0 -> 5.
    result = gw.grant_writing_vague_claim_verb_check(
        "本研究では、研究者三者が有する非線形磁気モデリング、物理ベース等価回路抽出、"
        "高周波損失測定の技術を統合し、回路シミュレータ上で利用可能なモデルを構築する。"
    )

    assert result["applicable"]
    risk = next(
        r for r in result["risks"]
        if r["type"] == "claim_verb_without_mechanism"
    )
    assert risk["verb"] == "統合し"
    assert result["score"] < 10


def test_vague_claim_verb_check_accepts_a_named_operation():
    result = gw.grant_writing_vague_claim_verb_check(
        "両モデルを巻線電流と誘起電圧を介して双方向に連成し、"
        "PoL変換回路の損失と動作波形を統一的に解析する。"
    )

    assert result["risks"] == []
    assert result["concrete_uses"] == [] or result["score"] == 10.0


def test_vague_claim_verb_check_credits_a_mechanism_in_the_same_sentence():
    result = gw.grant_writing_vague_claim_verb_check(
        "二つの解析を、電圧と電流を介して双方向に統合する。"
    )

    assert result["risks"] == []
    assert result["concrete_uses"]
    assert result["concrete_uses"][0]["mechanism_markers"]


def test_vague_claim_verb_check_is_not_applicable_without_such_verbs():
    result = gw.grant_writing_vague_claim_verb_check(
        "誘導加熱の発熱量を評価し、損失分布を求める。"
    )

    assert not result["applicable"]
    assert result["score"] is None


DIVERGENT_CLAIM_DRAFT = (
    "本研究は「異なる研究室の解析手法を、内部形式を統一せずに連携・差し替えたとき、"
    "解析手法の違いによって設計候補の順位が変わる境界を、どう定量化し、"
    "第三者が検証可能な形で示せるか」を問う。"
    "電気機器の電磁界解析では手法が高度化されてきた。"
    "中心の問いは、異なる研究室が独自形式で実装した解析を、内部形式を統一せずに結合し、"
    "その差が設計候補の優劣を覆さない条件を、どのように記述・検証できるか、である。"
)

UNIFIED_CLAIM_DRAFT = (
    "本研究は次を問う。異なる研究室の解析手法を、内部形式を統一せずに連携・差し替えたとき、"
    "解析手法の違いを考慮しても設計候補の順位を確定できる条件を、どのように定量化し、"
    "第三者が検証可能な形で示せるか。"
    "電気機器の電磁界解析では手法が高度化されてきた。"
    "中心の問いは次である。こうして定義した解析モジュールを、内部形式のまま連携・差し替える。"
    "このとき解析手法が異なっても設計候補の順位を確定できる条件を、どのように定量化できるか。"
)

DISTINCT_ROLE_CLAIM_DRAFT = (
    "本研究の目的は、解析手法が異なっても設計候補の順位を確定できる条件を定量化し、"
    "その条件を高忠実度解析へ進む判断基準へ展開することである。"
    "中心の問いは、設計量の変動幅から順位を確定できる条件をどのように定量化できるか、である。"
)


def test_central_claim_check_flags_a_question_restated_with_other_nouns():
    # The defect a keyword-coverage checker cannot see: every required word is
    # present, but the summary promises a 境界 while the body promises a 条件,
    # so a reviewer cannot tell whether there is one question or two.
    result = gw.grant_writing_central_claim_consistency_check(
        DIVERGENT_CLAIM_DRAFT
    )

    assert result["applicable"]
    assert result["statement_count"] == 2
    risk = next(
        r for r in result["risks"] if r["type"] == "outcome_noun_divergence"
    )
    assert risk["severity"] == "HIGH"
    assert "境界" in risk["comment"] and "条件" in risk["comment"]
    assert result["score"] < 10


def test_central_claim_check_accepts_a_question_restated_with_the_same_nouns():
    result = gw.grant_writing_central_claim_consistency_check(
        UNIFIED_CLAIM_DRAFT
    )

    assert result["applicable"]
    assert result["statement_count"] == 2
    assert result["risks"] == []
    assert result["score"] == 10.0


def test_central_claim_check_allows_distinct_terms_in_distinct_roles():
    result = gw.grant_writing_central_claim_consistency_check(
        DISTINCT_ROLE_CLAIM_DRAFT
    )

    assert result["applicable"]
    assert not any(
        risk["type"] == "outcome_noun_divergence"
        for risk in result["risks"]
    )


def test_central_claim_check_does_not_swallow_the_second_statement():
    # A greedy multi-sentence window merged both claims into one and reported
    # a single statement, which silently disabled the whole check.
    result = gw.grant_writing_central_claim_consistency_check(
        DIVERGENT_CLAIM_DRAFT
    )

    markers = [s["marker"] for s in result["statements"]]
    assert "を問う" in markers
    assert "中心の問い" in markers


def test_central_claim_check_ignores_headings_and_passing_mentions():
    result = gw.grant_writing_central_claim_consistency_check(
        "研究背景と学術的問い。手法を評価する。"
    )

    assert not result["applicable"]
    assert result["statements"] == []


def test_central_claim_check_runs_in_health_report():
    report = gw.grant_writing_health_report(DIVERGENT_CLAIM_DRAFT)

    assert "central_claim_consistency" in report["detailed_results"]
    issue = next(
        i for i in report["findings"]
        if i["name"] == "central_claim_consistency_check"
    )
    assert issue["severity"] == "HIGH"


RESEARCH_PLAN_SECTION = (
    "本研究の学術的問いは、異種解析を結合しても設計判断を保つ条件は何かである。"
    "低費用解析で候補順位を確定できる領域を評価し、順位が定まらない候補だけを"
    "高忠実度解析で再評価する。AIが開発を加速する現在、共有仕様と試験がなければ"
    "重複実装が増える。mdxで最終評価を行い、成果を国際会議で発表する。"
    "MCPには検証手順を蓄積し、教員・学生が利用している。"
)


def test_budget_check_is_not_applicable_to_a_research_plan_section():
    # A methods/plan section carries no budget by design. Reporting it as a
    # thin budget was a false HIGH that dragged the whole health report down.
    result = gw.grant_writing_budget_alignment_check(RESEARCH_PLAN_SECTION)

    assert result["applicable"] is False
    assert result["score"] is None
    assert result["comments"] == []


def test_budget_axes_require_a_cost_token_near_the_keyword():
    # 評価 and AI appear in every methods section; alone they are not evidence
    # that PoC work or AI usage was costed.
    result = gw.grant_writing_budget_alignment_check(
        "予算として旅費を計上する。評価を行い、AIが開発を加速する。"
    )

    assert result["applicable"] is True
    assert not result["axis_results"]["poc_experiment"]["ok"]
    assert not result["axis_results"]["ai_agent_costs"]["ok"]
    assert result["axis_results"]["dissemination"]["ok"]


def test_budget_axes_accept_keywords_stated_as_costs():
    result = gw.grant_writing_budget_alignment_check(
        "その他2,400千円の内訳は、生成AI費1,404千円とmdx計算資源664千円である。"
        "基板試作費と計測評価費を消耗品費へ計上する。"
    )

    assert result["applicable"] is True
    assert result["axis_results"]["ai_agent_costs"]["ok"]
    assert result["axis_results"]["compute_resources"]["ok"]
    assert result["axis_results"]["poc_experiment"]["ok"]


def test_health_report_omits_an_inapplicable_budget_from_scoring():
    report = gw.grant_writing_health_report(RESEARCH_PLAN_SECTION)

    assert report["detailed_results"]["budget"]["applicable"] is False
    assert "budget" not in report["detailed_scores"]
    assert not any(
        issue["name"] == "budget_alignment_check"
        for issue in report["findings"]
    )


def test_ethics_axis_ignores_merely_naming_students():
    result = gw.grant_writing_collaborative_integration_risk_check(
        "異種解析を結合する。MCPには検証手順を蓄積し、教員・学生が利用している。"
    )

    assert result["axis_results"]["evaluation_unit_and_ethics"]["not_applicable"]


def test_ethics_axis_still_fires_when_people_are_measured():
    result = gw.grant_writing_collaborative_integration_risk_check(
        "異種解析を結合する。学生の作業を記録し、手順ごとに比較する。"
    )

    assert not result["axis_results"]["evaluation_unit_and_ethics"].get(
        "not_applicable"
    )


def test_core_vs_optional_scope_accepts_ripple_effect_wording():
    result = gw.grant_writing_collaborative_integration_risk_check(
        "異種解析を結合する。必達範囲は二課題とし、NVH等は波及効果とする。"
    )

    assert result["axis_results"]["core_vs_optional_scope"]["ok"]


def test_budget_policy_mentions_fill_rate_strategy():
    result = gw.grant_writing_budget_alignment_check(KDDI_SAMPLE)

    assert "充足率" in result["budget_policy"]
    assert "挑戦的研究" in result["budget_policy"]


def test_grant_writing_reexports_ja_lint_helpers():
    result = gw.grant_writing_lint_bedrock("これは重要であると考えられる。")

    assert isinstance(result, dict)
    assert "issue_count" in result


def test_acronym_audit_accepts_gloss_and_ignores_hyphenated_project_name():
    result = gw.grant_writing_acronym_usage_audit(
        "知識・実行インターフェース（Model Context Protocol; MCP）を用いる。"
        "MCPで機能を提示し、MCPで実行する。JP-MARsで版管理する。"
        "READMEに再現手順を記す。"
    )

    by_name = {item["acronym"]: item for item in result["findings"]}
    assert by_name["MCP"]["verdict"] == "ok"
    assert "JP" not in by_name
    assert "README" not in by_name


def test_undefined_acronyms_accept_gloss_and_hyphenated_project_name():
    result = gw.grant_writing_find_undefined_acronyms(
        "知識・実行インターフェース（Model Context Protocol; MCP）を用いる。"
        "JP-MARsで版管理する。"
    )

    assert result["undefined"] == []


def test_sentence_analysis_ignores_latex_scaffolding():
    tex = r"""
% This deliberately long template comment must not be treated as proposal prose even when it exceeds the sentence threshold by a wide margin.
\section{研究目的}
\input{pieces/template_header}
本研究は技術報告と参照実装を接続する。第三者が結果を検証する。
"""

    result = gw.grant_writing_analyze_sentences(tex, max_len=50)

    # Fixed form headings are not applicant prose; the two body sentences are.
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


def test_international_standing_accepts_named_two_way_evidence():
    result = gw.grant_writing_international_standing_check(
        "日本発のCauer縮約と欧州発の高次要素を相互検証し、双方へ還流する。"
        "グラーツ工科大学での共同研究とIGTEでの共著がその起点である。"
    )

    assert result["applicable"]
    assert result["score"] == 10.0
    assert result["named_counterparts"]
    assert result["national_value_markers"]


def test_international_standing_flags_a_declaration_with_no_partner():
    # The axis a 2025 disclosure scored 1.60 against an adopted 2.70.
    result = gw.grant_writing_international_standing_check(
        "本研究は将来的に国際的な展開を目指す。海外の研究者とも連携したい。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "no_named_counterpart" in types
    assert result["score"] < 6


def test_international_standing_flags_a_catch_up_frame():
    result = gw.grant_writing_international_standing_check(
        "海外の研究者と共著を進め、国際会議で発表する。"
        "グラーツ工科大学と交流し、世界水準の技術に追いつくことを目指す。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "one_way_catch_up_frame" in types


def test_international_standing_is_not_applicable_without_the_subject():
    result = gw.grant_writing_international_standing_check(
        "誘導加熱の発熱量を評価し、損失分布を求める。"
    )

    assert not result["applicable"]
    assert result["score"] is None


def test_international_standing_flags_outputs_that_are_all_planned():
    # A network being formed and one already established must not read alike.
    # Stating an accepted paper as accepted is what separates them.
    result = gw.grant_writing_international_standing_check(
        "グラーツ工科大学と相互に共同研究を進める。"
        "国際会議での共著発表を目指す。日本発の手法を還流する予定である。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "international_output_all_planned" in types
    assert result["planned_output_sentences"]


def test_international_standing_credits_an_accepted_output():
    result = gw.grant_writing_international_standing_check(
        "ウィーン工科大学のHollaus氏と相互に共同研究を進め、"
        "日本発の手法との相互検証をIGTE Symposium 2026採択共著論文として得た。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "international_output_all_planned" not in types
    assert result["achieved_output_sentences"]


def test_international_standing_records_adopted_jsps_visit_as_preparation():
    result = gw.grant_writing_international_standing_check(
        "ウィーン工科大学のHollaus氏との共同研究に向け、日本学術振興会の"
        "外国人研究者招へい事業に採択され、31日間招へいすることが決定した。"
        "日本発の手法と欧州の手法を相互検証し、双方へ還流する。"
    )

    types = {risk["type"] for risk in result["risks"]}
    assert result["preparation_evidence_sentences"]
    assert "採択" in result["preparation_status_markers"]
    assert "31日間" in result["preparation_scale_markers"]
    assert "international_preparation_status_unclear" not in types
    assert "international_preparation_scale_missing" not in types


def test_irreplaceability_flags_intent_without_an_asset_or_demand():
    # The question a funder on either side asks: why this partner rather than
    # someone closer. Shared enthusiasm answers neither direction.
    result = gw.grant_writing_collaboration_irreplaceability_check(
        "ウィーン工科大学と国際共同研究を進め、相互に交流する。"
        "国際会議での共著発表を目指す。"
    )

    assert result["applicable"]
    types = {r["type"] for r in result["risks"]}
    assert "no_asset_this_side_holds" in types
    assert "no_evidence_partner_wants_it" in types
    assert result["score"] < 6


def test_irreplaceability_accepts_a_named_asset_the_partner_asked_for():
    result = gw.grant_writing_collaboration_irreplaceability_check(
        "日本発の階層行列ライブラリを用いた積分方程式解法について、"
        "ミラノ工科大学より議論の招請を受けている。"
        "独立に発展した別系統との相互検証は、国内の近隣機関では代替できない。"
    )

    assert result["risks"] == []
    assert result["score"] == 10.0
    assert result["asset_markers"] and result["demand_markers"]


def test_irreplaceability_reports_the_missing_half_only():
    # The current 基盤(C) draft states the asset and why no substitute works,
    # but never shows the counterpart asking. That single gap is the finding.
    result = gw.grant_writing_collaboration_irreplaceability_check(
        "日本発のCauer縮約と欧州発の高次要素は異なる系譜であり、"
        "相互検証によってのみ妥当性を確認できる。グラーツ工科大学と行う。"
    )

    types = {r["type"] for r in result["risks"]}
    assert types == {"no_evidence_partner_wants_it"}


def test_irreplaceability_is_not_applicable_without_a_named_partner():
    result = gw.grant_writing_collaboration_irreplaceability_check(
        "誘導加熱の発熱量を評価し、損失分布を求める。"
    )

    assert not result["applicable"]
    assert result["score"] is None


def test_budget_narrative_flags_amounts_without_a_calculation_basis():
    # A bare total duplicates the table without satisfying the current JSPS
    # request for a calculation basis.
    result = gw.grant_writing_budget_narrative_check(
        "鉄心を改造する必要がある。そのためコイルを巻き直す費用が60万円、"
        "またXYステージを70万円として見積もっている。"
    )

    assert result["applicable"]
    risk = next(
        r for r in result["risks"]
        if r["type"] == "amount_without_calculation_basis"
    )
    assert "60万円" in risk["amounts"]
    assert result["score"] < 10


def test_budget_narrative_accepts_the_editor_rewrite():
    result = gw.grant_writing_budget_narrative_check(
        "加速器を想定してギャップ部の分布を詳細に計測する必要があるため、"
        "ギャップ付き鉄心に改造する費用を計上している。"
        "成果は電気学会で発表する。"
    )

    assert result["risks"] == []
    assert result["score"] == 10.0


def test_budget_narrative_accepts_an_itemized_calculation_basis():
    result = gw.grant_writing_budget_narrative_check(
        "国内旅費は、学会発表のため1件15万円×3回=45万円として計上する。"
        "会場未定の国際会議は1件60万円として暫定計上する。"
    )

    assert result["risks"] == []
    assert result["score"] == 10.0


def test_budget_narrative_flags_travel_with_no_venue():
    result = gw.grant_writing_budget_narrative_check(
        "共同研究先へ出張する必要があるため、旅費を計上した。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "travel_without_dissemination_plan" in types


def test_budget_narrative_is_not_applicable_without_a_necessity_section():
    result = gw.grant_writing_budget_narrative_check(
        "誘導加熱の発熱量を評価し、損失分布を求める。"
    )

    assert not result["applicable"]
    assert result["score"] is None


def _write_form_pdf(path, fields, overflow_notice=""):
    """Build a stand-in compiled form: one page per entry of ``fields``."""
    import fitz

    doc = fitz.open()
    for title, pages in fields:
        for offset in range(pages):
            page = doc.new_page()
            head = title if offset == 0 else f"[{title} (tsuzuki)]"
            page.insert_text((72, 72), head, fontname="helv", fontsize=11)
            page.insert_text((72, 700), "body", fontname="helv", fontsize=11)
    if overflow_notice:
        page = doc.new_page()
        # Helvetica cannot encode the Japanese notice; use the built-in CJK font.
        page.insert_text((72, 72), overflow_notice, fontname="japan", fontsize=11)
    doc.save(str(path))
    doc.close()


def _write_form_tex(path, title, max_pages):
    path.write_text(
        "\\section{" + title + "}\n"
        f"%    <<最大 {max_pages}ページ>>\n\n本文。\n",
        encoding="utf-8",
    )


def test_page_limit_flags_a_field_past_its_allowance(tmp_path):
    _write_form_tex(tmp_path / "a.tex", "PURPOSE", 2)
    _write_form_pdf(tmp_path / "form.pdf", [("PURPOSE", 3)])

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    assert result["applicable"]
    over = next(r for r in result["risks"] if r["severity"] == "CRITICAL")
    assert "3ページ占めている" in over["comment"]
    assert result["fields"][0]["used_pages"] == 3


def test_page_limit_reads_the_notice_the_form_prints_itself(tmp_path):
    # No .tex declaration at all: the template's own notice must still fire.
    _write_form_pdf(
        tmp_path / "form.pdf",
        [("PURPOSE", 1)],
        overflow_notice="「PURPOSE」は4ページ以内で書いてください。",
    )

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    assert result["applicable"]
    assert any(r["severity"] == "CRITICAL" for r in result["risks"])
    assert "様式が超過を印字している" in result["comments"][0]


def test_page_limit_reports_an_allowance_left_unused(tmp_path):
    _write_form_tex(tmp_path / "a.tex", "PURPOSE", 4)
    _write_form_pdf(tmp_path / "form.pdf", [("PURPOSE", 2)])

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    risk = next(r for r in result["risks"] if r["severity"] == "MEDIUM")
    assert "4ページ許容のうち2ページ" in risk["comment"]


def test_page_limit_keeps_quiet_on_a_short_single_page_field(tmp_path):
    # A one-page compliance field is often short because the honest answer is.
    _write_form_tex(tmp_path / "a.tex", "RIGHTS", 1)
    _write_form_pdf(tmp_path / "form.pdf", [("RIGHTS", 1)])

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    assert result["risks"] == []
    assert result["score"] == 10.0


def test_health_report_finds_a_single_sibling_pdf(tmp_path):
    _write_form_tex(tmp_path / "a.tex", "PURPOSE", 1)
    _write_form_pdf(tmp_path / "form.pdf", [("PURPOSE", 2)])

    result = gw.grant_writing_health_report(str(tmp_path / "a.tex"))

    assert "page_limit" in result["detailed_results"]
    assert any(f["name"] == "page_limit_check" for f in result["findings"])


def test_health_report_leaves_page_limits_alone_when_the_pdf_is_ambiguous(tmp_path):
    _write_form_tex(tmp_path / "a.tex", "PURPOSE", 1)
    _write_form_pdf(tmp_path / "one.pdf", [("PURPOSE", 2)])
    _write_form_pdf(tmp_path / "two.pdf", [("PURPOSE", 2)])

    result = gw.grant_writing_health_report(str(tmp_path / "a.tex"))

    page_result = result["detailed_results"]["page_limit"]
    assert not page_result["applicable"]
    assert "ambiguous" in page_result["reason"]
    assert len(page_result["pdf_candidates"]) == 2


def test_publication_list_is_not_linted_as_one_long_sentence():
    proposal = (
        "本研究は結合条件を明らかにする。\n"
        "\\begin{enumerate}\n"
        "\\item K. Sugahara, ``Electromagnetic Analysis of Eddy Current Testing"
        " With Kelvin Transformation,'' IEEE Trans. Magn., 58(9), 1--6 (2022).\n"
        "\\item H. Nagamine, T. Yamaguchi, and K. Sugahara, ``A Pullback-Based"
        " Formulation of Kelvin Transformation,'' CEFC 2026.\n"
        "\\end{enumerate}\n"
    )

    result = gw.grant_writing_analyze_sentences(gw._prose_for_lint(proposal))

    assert result["max_length"] < 60
    assert result["over_threshold_count"] == 0


def test_main_tex_assembles_sibling_inputs_but_not_template_pieces(tmp_path):
    # The JSPS template keeps one file per field beside the main file and its
    # own scaffolding in pieces/. Run on a field file, every whole-proposal
    # check reported the other fields as missing (2026-09-02).
    (tmp_path / "pieces").mkdir()
    (tmp_path / "pieces" / "form_header.tex").write_text(
        "本欄には研究目的を記述すること。TEMPLATE_PIECE\n", encoding="utf-8"
    )
    (tmp_path / "kiban_c_01_purpose_plan.tex").write_text(
        "\\input{pieces/form_header}\n\\section{研究目的}\n本研究の学術的問いはXである。\n",
        encoding="utf-8",
    )
    (tmp_path / "kiban_c_03_abilities.tex").write_text(
        "\\section{遂行能力}\n菅原は開境界解析を研究してきた。\n", encoding="utf-8"
    )
    (tmp_path / "comment_only.tex").write_text(
        "COMMENTED_INPUT_MUST_NOT_APPEAR\n", encoding="utf-8"
    )
    (tmp_path / "kiban_c.tex").write_text(
        "\\documentclass{jarticle}\n"
        "%\\input{kiban_c_03_abilities} commented out, must not be read twice\n"
        "\\begin{document}\n"
        "\\input{kiban_c_01_purpose_plan}\n"
        "\\input{kiban_c_03_abilities.tex} % \\input{comment_only} trailing comment\n"
        "\\input{kiban_c}\n"
        "\\end{document}\n",
        encoding="utf-8",
    )

    text = gw._read_text_if_path(str(tmp_path / "kiban_c.tex"))

    assert "学術的問いはX" in text
    assert text.count("開境界解析") == 1
    assert "TEMPLATE_PIECE" not in text
    assert "COMMENTED_INPUT_MUST_NOT_APPEAR" not in text
    assert "\\input{pieces/form_header}" in text

    # A field file on its own still reads as before.
    alone = gw._read_text_if_path(str(tmp_path / "kiban_c_03_abilities.tex"))
    assert "学術的問い" not in alone


def test_cross_organization_pilot_reads_the_limit_from_the_whole_paragraph():
    # The honest closing sentence of a real 準備状況 paragraph sat four
    # sentences after the last trigger word and fell outside the window.
    text = (
        "大学間予備実証では、A大学提供のソルバー実装をB大学が同一問題へ統合した。"
        "C大学提供のモデルで再実行した。"
        "197,395自由度で168反復、真の相対残差8e-11で収束した。"
        "D大学は微分形式から材料則の再構成を導いた。"
        "これをCEFC 2026で発表した。"
        "ただし、別機関が共通の設計量から候補順位を判定できるかは未検証である。"
    )

    result = gw.grant_writing_cross_organization_pilot_check(text)

    assert result["applicable"]
    assert "remaining_gap" not in result["missing_axes"]


def test_cross_organization_pilot_treats_wrapped_lines_as_one_paragraph():
    text = (
        "大学間予備実証では、A大学提供の実装をB大学が統合した。\n"
        "C大学のモデルで再実行した。\n"
        "197,395自由度で168反復、相対残差8e-11で収束した。\n"
        "微分形式から材料則の再構成を導いた。\n"
        "CEFCで成果を発表した。\n"
        "別機関が候補順位を判定できるかは未検証である。\n\n"
        "次の独立した段落で研究計画を述べる。"
    )

    result = gw.grant_writing_cross_organization_pilot_check(text)

    assert "remaining_gap" not in result["missing_axes"]


def test_budget_source_consistency_accepts_japanese_category_headings(tmp_path):
    header = (
        "費目区分/Expenditure Categories,年度/FY,品名・仕様/Item (Specification),"
        "設置機関/Place,品目/Item,数量/Qty,単価/Unit Price,金額/Amount\n"
    )
    source = tmp_path / "budget.csv"
    source.write_text(
        header
        + "B,2027,,,SSD,,,50\n"
        + "E,2027,,,研究協力者謝金,,,150\n"
        + "F,2027,,,JMAG,,,1023\n",
        encoding="utf-8-sig",
    )

    result = gw.grant_writing_budget_source_consistency_check(
        str(source),
        expected_category_totals_json=(
            '{"消耗品費": 50, "人件費・謝金": 150, "その他": 1023, "設備備品費": 0}'
        ),
    )

    # 設備備品費 declared as zero with no ledger row is agreement, not a gap.
    assert result["consistent"]
    assert result["canonical"]["category_labels"] == {
        "B": "消耗品費", "E": "人件費・謝金", "F": "その他",
    }

    result = gw.grant_writing_budget_source_consistency_check(
        str(source),
        expected_category_totals_json='{"消耗品費": 50, "国内旅費": 300}',
    )

    # A non-zero total for a category the ledger never mentions is reported
    # as absent, not as a mismatch against null.
    assert [d["type"] for d in result["differences"]] == ["category_not_in_ledger"]
    assert result["differences"][0]["category"] == "C"

    result = gw.grant_writing_budget_source_consistency_check(
        str(source),
        expected_category_totals_json='{"消耗品費": 50, "人件費・謝金": 150, "F": 1000}',
    )

    assert [d["type"] for d in result["differences"]] == ["category_total_mismatch"]
    assert result["differences"][0]["delta"] == 23


def test_proper_noun_load_lists_unplaced_singletons_but_keeps_venues_and_roles():
    # The 2026-09-02 foreign-matter case: FEMM and Meeker appear once, with no
    # role in the plan. Hollaus is named once too, but with 招へい beside the
    # name; COMPUMAG is a venue in the yearly plan.
    text = (
        "比留間の拡張有限要素法（XFEM）を混合Galerkin縮約へ接続する。"
        "Hollaus氏の31日間招へいが採択された。"
        "COMPUMAG 2027で結合手法を発表する。"
        "軸対称解析でも、公開ソフトFEMMの作者Meeker氏の断片コードを2次要素へ拡張実装した。"
    )

    result = gw.grant_writing_proper_noun_load_check(text)

    names = {s["name"]: s for s in result["singletons"]}
    # 「作者Meeker氏」 attributes a program; it does not place him in the plan.
    assert "Meeker" in names and names["Meeker"]["kind"] == "person"
    assert not names["Meeker"]["role_stated"]
    assert names["Hollaus"]["role_stated"]
    # All-caps tokens (FEMM, XFEM, COMPUMAG) belong to the acronym audit.
    assert "FEMM" not in names and "COMPUMAG" not in names
    assert "Galerkin" not in names
    assert result["unplaced_singleton_count"] == 1
    assert result["comments"] and "Meeker" in result["comments"][0]

    dated_nonvenue = gw.grant_writing_proper_noun_load_check(
        "2027年度にMeeker氏の実装断片を解析へ加える。"
    )
    assert dated_nonvenue["unplaced_singleton_count"] == 1

    report = gw.grant_writing_health_report(text, program="kaken_generic")
    assert any(q["name"] == "proper_noun_load_check" for q in report["questions"])
    assert all(f["name"] != "proper_noun_load_check" for f in report["findings"])

    clean = gw.grant_writing_proper_noun_load_check(
        "手法差による設計量の変動を求め、順位を確定できる条件を示す。"
    )
    assert clean["applicable"] is False


def test_prose_list_items_survive_but_stay_separate():
    proposal = (
        "\\begin{itemize}\n"
        "\\item 誘導加熱の発熱量を評価する\n"
        "\\item 加速器電磁石の出口位置ずれを抑える\n"
        "\\end{itemize}\n"
    )

    prose = gw._prose_for_lint(proposal)

    assert "誘導加熱の発熱量を評価する。" in prose
    assert "加速器電磁石の出口位置ずれを抑える。" in prose


def test_a_parenthesised_gloss_is_not_an_unfilled_placeholder():
    # From a submitted proposal: 入力 opens an ordinary technical gloss.
    result = gw.grant_writing_template_residue_check(
        "設計技術を高度化し，高いビーム効率（入力したエネルギーに対する"
        "ビーム強度）を実現する。"
    )

    assert [r for r in result["risks"] if r["type"] == "unfilled_placeholder"] == []


def test_a_real_placeholder_parenthetical_still_fires():
    result = gw.grant_writing_template_residue_check(
        "申請金額（記入してください）を確定する。研究期間は（未定）である。"
    )

    matches = [
        r["match"] for r in result["risks"] if r["type"] == "unfilled_placeholder"
    ]
    assert matches


def test_an_undecided_future_conference_venue_is_not_a_placeholder():
    result = gw.grant_writing_template_residue_check(
        "電気学会（開催地未定）へ4名が参加する。"
        "開催地未定のため、学内旅費規程に基づき1名7万円で暫定積算する。"
    )

    assert [r for r in result["risks"] if r["type"] == "unfilled_placeholder"] == []


def test_a_conference_venue_is_not_an_international_claim():
    # A country name inside a travel line names where a meeting is held.
    result = gw.grant_writing_international_standing_check(
        "本研究は負ミュオン核変換の磁場設計技術を確立する。"
        "国際会議 Conference（2026/5/17~22,フランス）：50万円を計上する。"
    )

    assert not result["applicable"]


def test_a_worldwide_problem_does_not_need_a_named_counterpart():
    result = gw.grant_writing_international_standing_check(
        "放射性廃棄物の処理は世界的な社会課題であり，長寿命核分裂生成物の"
        "低減が求められている。加速器電磁石の磁場精度を高める。"
    )

    assert not result["applicable"]


def test_a_domestic_partnership_is_not_a_foreign_one():
    result = gw.grant_writing_international_standing_check(
        "本事業は株式会社MotorAIと近畿大学の共同開発として実施する。"
        "想定する国内、海外市場の規模を調査する。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "no_named_counterpart" not in types


def test_alternating_current_is_not_academic_exchange():
    # 交流 in an electrical proposal is AC. A glossary row tripped this.
    result = gw.grant_writing_international_standing_check(
        "誘導加熱は交流電流によって導体を加熱する技術である。"
        "海外市場は年間4%で成長すると予測されている。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "no_named_counterpart" not in types


def test_a_conference_acronym_is_not_a_collaboration_partner():
    result = gw.grant_writing_collaboration_irreplaceability_check(
        "成果はCOMPUMAG 2027およびCEFC 2028で発表する。"
        "誘導加熱コイルの設計最適化を行う。"
    )

    assert not result["applicable"]


def test_a_named_institution_still_starts_the_partner_question():
    result = gw.grant_writing_collaboration_irreplaceability_check(
        "ミラノ工科大学と階層行列圧縮の適用について共同研究を行う。"
    )

    assert result["applicable"]


def test_a_compact_form_is_not_judged_on_the_three_review_criteria():
    # One 要旨 box plus keywords and amounts: the form offers nowhere to
    # write 研究遂行能力, so demanding it is a finding its author would argue.
    form = (
        "研究テーマ 負ミュオン核変換実現のための加速器用電磁石の磁場計算。"
        "研究テーマの要旨 " + "加速器電磁石の磁場分布を数値的に求める。" * 20
        + "キーワード 加速器、電磁石。申請金額 230万円。"
    )

    result = gw.grant_writing_kaken_review_format_check(form)

    types = {r["type"] for r in result["risks"]}
    assert "review_criteria_axis_missing" not in types


def test_a_proposal_body_missing_one_axis_is_still_reported():
    body = (
        "本研究の独自性は異種解析の結合条件を定量化する点にある。" * 25
        + "妥当性は解析解との比較で検証し、実証する。" * 25
        + "波及効果として設計手順が再利用できる。" * 25
    )

    result = gw.grant_writing_kaken_review_format_check(body)

    types = {r["type"] for r in result["risks"]}
    assert "review_criteria_axis_missing" in types


def test_a_year_by_task_matrix_is_not_one_long_sentence():
    # An adopted proposal's 年度計画 is a matrix of short cells with no full
    # stop. Joined, it was reported as a single 455-character sentence.
    plan = (
        "［研究計画］\n令和2年度\n令和3年度\n令和4年度\n"
        "マルチスケールモデル縮約\n定式化・実装\n（汎用シミュレータ実装）\n"
        "マルチフィジクスモデル縮約\n実現方法の検討・定式化\n実装\n"
        "モータモデル縮約\n回転機への応用\nマルチスケール化検討\n"
        "非線形化\nマルチフィジクス化\nシミュレータ化\n実証用モータ実験"
    )

    result = gw.grant_writing_analyze_sentences(plan)

    assert result["max_length"] < 30
    assert result["over_threshold_count"] == 0


def test_a_capability_handed_to_a_non_member_is_reported():
    # The line is from a rejected 基盤C whose novelty was machine learning.
    proposal = (
        "本研究では機械学習による電気機器設計の自動化を行う。"
        "トポロジー最適化に機械学習を用いる点が新しい。"
        "機械学習の観点からは最適化計算の高速化を検討する。\n"
        "研究代表者\t菅原賢悟：研究統括，電磁界解析技術全般担当\n"
        "連携研究者\t浅川伸一：機械学習に関する専門知識の供与\n"
    )

    result = gw.grant_writing_capability_responsibility_check(proposal)

    assert result["applicable"]
    risk = result["risks"][0]
    assert risk["role"] == "連携研究者"
    assert "機械学習" in risk["terms"]
    assert "浅川伸一" not in risk["terms"]
    assert "連携研究者" not in risk["terms"]


def test_a_team_of_members_only_is_not_judged():
    # Every name in an adopted proposal's 役割分担 carried a 分担 role.
    proposal = (
        "本研究はモデル縮約法を開発する。モデル縮約の要素技術は多い。"
        "モデル縮約は物理現象の本質を抽出する技術である。\n"
        "松尾哲司：とりまとめ，モデル縮約定式化\n"
        "高橋康人：モータモデル縮約法の実装・実証\n"
    )

    result = gw.grant_writing_capability_responsibility_check(proposal)

    assert not result["applicable"]


def test_prose_mentioning_a_collaborator_is_not_an_assignment():
    # 「有能な研究協力者を有する」 describes; it does not hand anyone a job.
    proposal = (
        "計算電磁気学の理論研究のため計算機環境は整備されており，"
        "大学院生など有能な研究協力者を有する。"
        "計算機環境は理論研究に用いる。計算機環境の増強も進めている。"
    )

    result = gw.grant_writing_capability_responsibility_check(proposal)

    assert not result["applicable"]


def test_role_description_filler_is_not_a_capability():
    proposal = (
        "本研究は誘導加熱コイルを設計する。誘導加熱の発熱量を評価する。"
        "誘導加熱の損失分布を求める。\n"
        "アドバイザー\t伊藤英昭：全般に関する助言の提供\n"
    )

    result = gw.grant_writing_capability_responsibility_check(proposal)

    assert result["risks"] == []
    assert result["score"] == 10.0


def test_form_instruction_text_is_not_the_applicants_prose():
    # Verbatim from an adopted 令和2年度 S-14 form.
    form = (
        "本欄には、本研究の目的と方法などについて、３頁以内で記述すること。\n"
        "冒頭にその概要を簡潔にまとめて記述し、本文には、(1)本研究の学術的背景、"
        "研究課題の核心をなす学術的「問い」、(2)本研究の目的および学術的独自性と"
        "創造性について具体的かつ明確に記述すること。\n"
        "本研究はモデル縮約法を開発する。\n"
    )

    prose = gw._prose_for_lint(form)

    assert prose == "本研究はモデル縮約法を開発する。"


def test_polite_form_sentences_are_read_as_the_forms_voice():
    # A proposal body is written in である調; the form speaks in ですます調.
    mixed = (
        "承認手続が必要となる調査・研究・実験などが対象となります。\n"
        "本研究は誘導加熱の発熱量を評価する。\n"
        "損失分布を求める。解析経路を選ぶ。候補順位を確定する。\n"
    )

    prose = gw._prose_for_lint(mixed)

    assert "対象となります" not in prose
    assert "本研究は誘導加熱の発熱量を評価する。" in prose


def test_a_document_written_in_polite_form_keeps_its_text():
    polite = (
        "本研究では誘導加熱の発熱量を評価します。\n"
        "損失分布を求めます。\n"
        "コイルのインピーダンスを測定します。\n"
    )

    prose = gw._prose_for_lint(polite)

    assert "発熱量を評価します" in prose


def test_a_furigana_field_is_not_a_foreign_counterpart():
    # A form puts スガハラ ケンゴ on one line and 氏名 on the next.
    form = (
        "（フリガナ）スガハラ ケンゴ\n氏名 菅原賢悟\n"
        "本事業は誘導加熱コイルの設計最適化を行う。海外市場の規模を調査する。"
    )

    assert gw._NAMED_PARTNER.search(form) is None
    assert not gw.grant_writing_international_standing_check(form)["applicable"]


def test_citing_foreign_prior_work_is_not_a_standing_claim():
    survey = (
        "電磁界のモデル縮約については有力な数学的手法がいくつか存在する。"
        "フランスの研究グループによる直交分解法が知られており、"
        "その適用範囲を検証した報告がある。"
    )

    assert not gw.grant_writing_international_standing_check(survey)["applicable"]


def test_an_ieee_publication_counts_as_international_output():
    record = (
        "ウィーン工科大学と共同研究を行う。"
        "成果は IEEE Transactions on Magnetics に発表した。"
    )

    result = gw.grant_writing_international_standing_check(record)

    types = {r["type"] for r in result["risks"]}
    assert "no_international_output" not in types


def test_a_name_on_its_own_line_does_not_join_the_paragraph_below():
    doc = (
        "菅原賢悟（研究分担者）\n"
        "(1)これまでの研究活動 三菱電機在職中に実務経験を積んだ後、"
        "豊富な知識を活用して先端技術開発を行っている（S1,2）\n"
    )

    result = gw.grant_writing_vague_claim_verb_check(doc)

    assert result["risks"] == []


def test_an_adnominal_claim_verb_is_not_a_promise():
    # 活用する modifies 産業分野; it says who uses the technology.
    sentence = "本研究開発の成果は、誘導加熱技術を活用する幅広い産業分野に波及する。"

    result = gw.grant_writing_vague_claim_verb_check(sentence)

    assert result["risks"] == []


def test_a_forward_claim_with_a_vague_verb_still_fires():
    # No mechanism named: a sentence that says 縮約 or 出力 would be concrete.
    sentence = "異なる研究室の資産を活用する。"

    result = gw.grant_writing_vague_claim_verb_check(sentence)

    assert any(r["type"] == "claim_verb_without_mechanism" for r in result["risks"])


def test_a_software_inventory_is_not_an_acronym_pile():
    inventory = (
        "本研究の遂行に必要な計算機資源を有する。"
        "Adventure, CST Studio, ELF, Elmer, EMCoS, EMSolution, FastCap, "
        "JMAG, COMSOL を保有している。"
    )

    result = gw.grant_writing_persuasion_quality_check(inventory)

    assert not [r for r in result["risks"] if r.get("type") == "acronym_pile"]


def test_a_price_charged_is_not_a_cost_incurred():
    plan = (
        "ライセンスビジネスの権利付与型は、顧客の採算が取れる1件あたり"
        "5,000千円に設定し、共同型は当社の人件費も計上するため1件あたり"
        "10,000千円に設定する。"
    )

    result = gw.grant_writing_budget_narrative_check(plan)

    assert not [
        r for r in result["risks"]
        if r["type"] == "amount_without_calculation_basis"
    ]


# Everything a real application form contains that the applicant did not
# argue: the form's instructions, a furigana field, a publication list, a
# year-by-task matrix, a software inventory, a price table, and headings.
# Eight separate false-positive families were traced to one of these being
# read as prose, so the suite must have nothing to say about a document made
# of nothing else.
NON_PROSE_ONLY = """１　研究目的、研究方法など
本欄には、本研究の目的と方法などについて、４頁以内で記述すること。
冒頭にその概要を簡潔にまとめて記述し、本文には、(1)本研究の学術的背景、研究課題の核心をなす学術的「問い」、(2)本研究の目的および学術的独自性と創造性について具体的かつ明確に記述すること。
本研究計画調書は「小区分」の審査区分で審査されます。
承認手続が必要となる調査・研究・実験などが対象となります。
なお、該当しない場合には、その旨記述してください。
記入に当たっては、基盤研究（Ｃ）（一般）研究計画調書作成・記入要領を参照してください。
（フリガナ）スガハラ ケンゴ
氏名 菅原賢悟
研究業績
\\begin{enumerate}
\\item K. Sugahara, ``Electromagnetic Analysis of Eddy Current Testing With Kelvin Transformation,'' IEEE Transactions on Magnetics, 58(9), 1--6 (2022).
\\item S. Hiruma and H. Igarashi, ``Fast 3-D Analysis of Eddy Current in Litz Wire Using Integral Equation,'' IEEE Transactions on Magnetics, 53(6), 1--4 (2017).
\\item H. Nagamine, T. Yamaguchi, and K. Sugahara, ``A Pullback-Based Formulation of Kelvin Transformation,'' CEFC 2026.
\\end{enumerate}
［研究計画］
令和9年度
令和10年度
令和11年度
モジュール結合
定式化・実装
検証
二課題移転
実装
再検証
研究環境
Adventure, CST Studio, ELF, Elmer, EMCoS, EMSolution, FastCap, JMAG, COMSOL
経費明細
設備備品費 1,000千円
消耗品費 500千円
旅費 800千円
"""


def test_the_suite_has_nothing_to_say_about_a_document_with_no_prose():
    report = gw.grant_writing_health_report(NON_PROSE_ONLY, program="kaken_oss")

    # section_presence also strips the form and correctly reports that a form
    # skeleton has no applicant argument in it.
    findings = [f for f in report["findings"] if f["name"] != "section_presence"]
    assert findings == [], [
        (f["name"], f["comments"][:1]) for f in findings
    ]


def test_kaken_form_instructions_do_not_satisfy_review_axis_presence():
    form = (
        "本欄には、研究課題の学術的重要性、研究方法の妥当性、研究遂行能力及び"
        "研究環境、国際性について記述してください。"
    )

    result = gw.grant_writing_section_presence(form, program="kaken_generic")

    assert set(result["missing_axes"]) == {
        "academic_importance",
        "method_validity",
        "feasibility_environment",
        "internationality",
    }


def test_direct_prose_tools_do_not_lint_the_form_or_publication_list():
    results = [
        gw.grant_writing_lint_bedrock(NON_PROSE_ONLY),
        gw.grant_writing_check_misuse_japanese(NON_PROSE_ONLY),
        gw.grant_writing_check_subject_predicate_distance(NON_PROSE_ONLY),
        gw.grant_writing_suggest_redundancy_fixes(NON_PROSE_ONLY),
    ]

    assert results[0]["issue_count"] == 0
    assert results[1]["total_matches"] == 0
    assert results[2]["violation_count"] == 0
    assert results[3]["total_matches"] == 0


def test_a_document_with_no_prose_leaves_no_sentences_to_measure():
    prose = gw._prose_for_lint(NON_PROSE_ONLY)
    result = gw.grant_writing_analyze_sentences(prose)

    assert result.get("over_threshold_count", 0) == 0
    # What survives is the software inventory line, ~75 characters. The
    # citation list, the year matrix and the instruction block are gone; left
    # in, any one of them alone exceeds the 90-character threshold.
    assert result.get("max_length", 0) < 90


def test_misuse_check_is_alive_even_though_proposals_never_trip_it():
    # Inherited from the shared Japanese table, which targets speech and
    # email: よろしかったでしょうか, のほう, こんにちわ. Eight real proposals
    # score zero, and that is the genre, not a broken check.
    result = gw.grant_writing_check_misuse_japanese(
        "資料のほうを送付いたします。こんにちわ。お連絡ありがとうございます。"
    )

    assert result["total_matches"] >= 2


def test_misuse_check_rejects_novelty_placement_meta_phrase():
    bad = gw.grant_writing_check_misuse_japanese(
        "本研究では物理量制約と高速求解の統合に新規性を置く。"
    )
    clean = gw.grant_writing_check_misuse_japanese(
        "本研究では、先行法にない物理量制約を低ランク求解へ組み込む。"
    )

    assert bad["total_matches"] == 1
    assert clean["total_matches"] == 0


def test_an_asserted_absence_needs_an_account_of_the_search():
    # Present in four of eight real proposals, three of them adopted.
    result = gw.grant_writing_literature_gap_evidence_check(
        "本提案事業に関して、類似する計画は存在しない。"
        "誘導加熱コイルの設計分野では、直接的な競合製品は存在しない。"
    )

    assert result["applicable"]
    assert [r["type"] for r in result["risks"]] == [
        "absence_claimed_without_search"
    ]
    assert result["unbacked_absence_claims"]


def test_an_absence_with_a_stated_search_is_not_flagged_as_unbacked():
    result = gw.grant_writing_literature_gap_evidence_check(
        "IEEE Xplore を対象文献として「stream function coil」で検索した結果、"
        "鉄心を含む設計例は存在しない。"
    )

    types = {r["type"] for r in result["risks"]}
    assert "absence_claimed_without_search" not in types


def test_a_proposal_that_claims_no_absence_is_not_judged_on_it():
    result = gw.grant_writing_literature_gap_evidence_check(
        "既往研究では、個々の手法の高精度化・高速化が進められてきた。"
        "本研究は手法間の差を設計量へ写す。"
    )

    assert not result["applicable"]


def test_inline_latex_comments_and_commented_inputs_are_not_linted(tmp_path):
    (tmp_path / "a.tex").write_text("A本文。\n", encoding="utf-8")
    (tmp_path / "b.tex").write_text("B本文。\n", encoding="utf-8")
    main = tmp_path / "main.tex"
    main.write_text(
        "\\input{a} % \\input{b} TODO 赤線 \\[x=y\\]\n",
        encoding="utf-8",
    )

    assembled = gw._read_text_if_path(str(main))
    assert "A本文" in assembled
    assert "B本文" not in assembled
    assert "TODO" not in gw._prose_for_lint(assembled)

    persuasion = gw.grant_writing_persuasion_quality_check(str(main))
    assert persuasion["counts"]["display_equations"] == 0
    residue = gw.grant_writing_template_residue_check(str(main))
    assert residue["risks"] == []
    review_format = gw.grant_writing_kaken_review_format_check(str(main))
    assert all(r["type"] != "color_dependent_figure" for r in review_format["risks"])


def test_fixed_form_headings_are_not_applicant_prose():
    prose = gw._prose_for_lint(
        "\\section{人権の保護及び法令等の遵守への対応}\n"
        "本研究では個人情報を扱わない。\n"
    )

    assert "人権の保護及び法令等の遵守への対応" not in prose
    assert prose == "本研究では個人情報を扱わない。"


def test_raw_comment_only_triggers_do_not_activate_prose_checks():
    commented = (
        "% 連携 統合 Radia 予算 内訳 単価 数量 見積書 100千円\n"
        "本研究では磁気浮上の成立条件を明らかにする。\n"
    )

    assert not gw.grant_writing_collaborative_integration_risk_check(commented)["applicable"]
    assert not gw.grant_writing_budget_alignment_check(commented)["applicable"]
    software = gw.grant_writing_named_software_abstraction_check(commented)
    assert not software["applicable"]


def test_short_keyword_matching_respects_word_and_character_boundaries():
    assert gw._contains_any("domain decomposition", ["ai", "oss", "api", "ci", "arm"]) == []
    assert gw._contains_any("core loss", ["oss"]) == []
    assert gw._contains_any("磁力と熱伝導", ["力", "熱"]) == []
    assert gw._contains_any("AI と OSS を用いる", ["ai", "oss"]) == ["ai", "oss"]

    basic = gw.grant_writing_kaken_basic_research_positioning_check(
        "domain decomposition improves a core loss model."
    )
    assert not basic["applicable"]
    domain = gw.grant_writing_domain_outcome_chain_check("core loss is measured.")
    assert not domain["applicable"]
    assert not gw._mentions_cost_nearby("server maintenance 費を計上する。", "ai")


def test_partner_matching_has_boundaries_and_excludes_the_home_institution():
    assert gw._NAMED_PARTNER.search("METHOD is compared with another method.") is None
    assert gw._NAMED_PARTNER.search("Kindai University develops the model.") is None
    assert gw._NAMED_PARTNER.search("TU Wien develops the independent solver.")


def test_electromagnetic_terms_do_not_satisfy_international_relationship_words():
    result = gw.grant_writing_international_standing_check(
        "ウィーン工科大学と海外共同研究を行う。"
        "円筒導体の交流磁界を研究分担者が解析する。"
    )

    assert result["applicable"]
    assert "交流" not in result["reciprocity_markers"]
    assert "分担" not in result["reciprocity_markers"]


def test_radia_name_does_not_match_radial_or_radiation():
    result = gw.grant_writing_kaken_oss_platform_check(
        "JP-MARsを研究基盤とする。radial field and radiation are evaluated."
    )

    assert not result["radia_integration"]["mentioned"]


def test_overlapping_self_negation_patterns_report_one_finding():
    result = gw.grant_writing_persuasion_quality_check(
        "この成果は工学的有用性を示すものではなく、根拠にはしない。"
    )

    assert result["counts"]["self_negating_evidence"] == 1


def test_form_role_words_are_not_parsed_as_people():
    result = gw.grant_writing_proper_noun_load_check(
        "同氏は研究代表者である。氏名 菅原賢悟。"
    )

    names = {item["name"] for item in result.get("singletons", [])}
    assert "同" not in names
    assert "究代表者" not in names


def test_two_item_enumeration_before_nado_is_not_a_hedge():
    result = gw.grant_writing_count_weak_expressions("高次要素、補助空間前処理などを比較する。")

    assert "など" not in result["by_pattern"]


def test_past_modifier_does_not_turn_a_future_claim_into_a_record():
    result = gw.grant_writing_vague_claim_verb_check(
        "統合したシステムを活用する。"
    )

    assert any(r["type"] == "claim_verb_without_mechanism" for r in result["risks"])


def test_bare_item_number_is_not_concrete_momentum_evidence():
    result = gw.grant_writing_reviewer_momentum_check(
        "設計判断が困難である。研究項目1を設定する。"
        "そこで本研究では成立条件を明らかにする。候補順位を確定する。"
    )

    assert not result["metrics"]["bottleneck_is_specific"]


def test_calendar_year_alone_does_not_identify_a_publication():
    result = gw.grant_writing_kaken_review_format_check(
        "2026年度の研究業績として主要論文がある。"
    )

    assert any(r["type"] == "publication_not_identifiable" for r in result["risks"])


def test_page_limit_reports_a_declared_field_missing_from_the_pdf(tmp_path):
    _write_form_tex(tmp_path / "a.tex", "PURPOSE", 2)
    _write_form_pdf(tmp_path / "form.pdf", [("OTHER", 4)])

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    assert result["applicable"]
    assert result["score"] < 10
    assert any(r["type"] == "declared_field_heading_not_found" for r in result["risks"])


def test_page_limit_does_not_consume_the_next_section_declaration(tmp_path):
    (tmp_path / "a.tex").write_text(
        "\\section{FIRST}\n本文。\n\\section{SECOND}\n% <<最大 1ページ>>\n本文。\n",
        encoding="utf-8",
    )
    _write_form_pdf(tmp_path / "form.pdf", [("FIRST", 1), ("SECOND", 1)])

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    assert [field["field"] for field in result["fields"]] == ["SECOND"]


def test_page_limit_handles_two_field_headings_on_the_same_page(tmp_path):
    import fitz

    (tmp_path / "a.tex").write_text(
        "\\section{FIRST}\n% <<最大 1ページ>>\n"
        "\\section{SECOND}\n% <<最大 1ページ>>\n",
        encoding="utf-8",
    )
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 72), "FIRST SECOND", fontname="helv", fontsize=11)
    doc.save(str(tmp_path / "form.pdf"))
    doc.close()

    result = gw.grant_writing_page_limit_check(str(tmp_path / "form.pdf"))

    assert all(field["used_pages"] >= 1 for field in result["fields"])
    assert any(r["type"] == "field_headings_share_page" for r in result["risks"])


def test_health_report_uses_tex_directory_for_an_explicit_pdf_elsewhere(tmp_path):
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "output"
    source_dir.mkdir()
    output_dir.mkdir()
    tex = source_dir / "proposal.tex"
    _write_form_tex(tex, "PURPOSE", 1)
    pdf = output_dir / "compiled.pdf"
    _write_form_pdf(pdf, [("PURPOSE", 2)])

    result = gw.grant_writing_health_report(str(tex), pdf=str(pdf))

    assert result["detailed_results"]["page_limit"]["applicable"]
    assert any(f["name"] == "page_limit_check" for f in result["findings"])


def test_health_report_prefers_a_same_stem_pdf_among_siblings(tmp_path):
    tex = tmp_path / "proposal.tex"
    _write_form_tex(tex, "PURPOSE", 1)
    _write_form_pdf(tmp_path / "proposal.pdf", [("PURPOSE", 2)])
    _write_form_pdf(tmp_path / "template.pdf", [("OTHER", 1)])

    result = gw.grant_writing_health_report(str(tex))

    assert result["detailed_results"]["page_limit"]["applicable"]
    assert any(f["name"] == "page_limit_check" for f in result["findings"])


def test_budget_parser_raises_on_a_non_numeric_candidate_row(tmp_path):
    source = tmp_path / "budget.csv"
    source.write_text(
        "expenditure categories/費目,FY/年度,item/品目,,amount/金額\n"
        "A,2026,PC,,300\n"
        "A,2027,server,,=100*2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="invalid budget row 3"):
        gw.grant_writing_budget_source_consistency_check(str(source))


def test_budget_csv_preserves_a_quoted_newline_and_decimal_total(tmp_path):
    source = tmp_path / "budget.csv"
    source.write_text(
        "expenditure categories,FY,item,amount\n"
        '人件費,2026,"RA\nwork",0.1\n'
        "謝金,2026,lecture,0.2\n",
        encoding="utf-8",
    )

    result = gw.grant_writing_budget_source_consistency_check(
        str(source),
        expected_total_thousand_yen=0.3,
        expected_category_totals_json='{"E": 0.3}',
    )

    assert result["consistent"]
    assert result["canonical"]["rows"][0]["item"] == "RA work"


def test_xlsx_reader_ignores_phonetics_and_honours_row_numbers_and_merges(tmp_path):
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    doc_rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    pkg_rel = "http://schemas.openxmlformats.org/package/2006/relationships"
    shared = (
        f'<sst xmlns="{main_ns}" count="1" uniqueCount="1">'
        '<si><t>設備備品費</t><rPh sb="0" eb="5"><t>セツビビヒンヒ</t></rPh></si>'
        "</sst>"
    )
    sheet = (
        f'<worksheet xmlns="{main_ns}"><sheetData>'
        '<row r="1"><c r="A1" t="inlineStr"><is><t>expenditure categories</t></is></c>'
        '<c r="B1" t="inlineStr"><is><t>FY</t></is></c>'
        '<c r="C1" t="inlineStr"><is><t>item</t></is></c>'
        '<c r="D1" t="inlineStr"><is><t>amount</t></is></c></row>'
        '<row r="3"><c r="A3" t="s"><v>0</v></c><c r="B3"><v>2026</v></c>'
        '<c r="C3" t="inlineStr"><is><t>PC</t></is></c><c r="D3"><v>300</v></c></row>'
        '<row r="4"><c r="B4"><v>2026</v></c>'
        '<c r="C4" t="inlineStr"><is><t>GPU</t></is></c><c r="D4"><v>200</v></c></row>'
        '</sheetData><mergeCells count="1"><mergeCell ref="A3:A4"/></mergeCells></worksheet>'
    )
    workbook = (
        f'<workbook xmlns="{main_ns}" xmlns:r="{doc_rel}"><sheets>'
        '<sheet name="s" sheetId="1" r:id="rId1"/></sheets></workbook>'
    )
    rels = (
        f'<Relationships xmlns="{pkg_rel}">'
        '<Relationship Id="rId1" Type="x" Target="/xl/worksheets/sheet1.xml"/>'
        "</Relationships>"
    )
    path = tmp_path / "budget.xlsx"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", rels)
        archive.writestr("xl/sharedStrings.xml", shared)
        archive.writestr("xl/worksheets/sheet1.xml", sheet)

    rows = gw._budget_ledger_rows(path, "s")

    assert [row["category"] for row in rows] == ["設備備品費", "設備備品費"]
    assert [row["source_row"] for row in rows] == [3, 4]


def test_health_report_scores_translationese_and_surfaces_non_scoring_diagnostics():
    text = (
        "設計判断には長い時間が必要である。"
        "この劇的な方法は解析条件を設計へ発展する。"
        "異種解析モジュールの機能、入出力物理量、実行条件、判定区間を"
        "MCP、GitHub、CIで版管理し、設計候補の順位確定と第三者反証へ対応付ける。"
        "XFEM、HACApK、ESIMをCauer縮約へ接続し、損失と効率を判定する。"
    )
    result = gw.grant_writing_health_report(text)

    assert "translationese" in result["detailed_scores"]
    assert result["defect_score"] < 10
    question_names = {question["name"] for question in result["questions"]}
    assert "adjacent_reviewer_readability_check" in question_names
    assert "reviewer_momentum_check" in question_names


def test_japanese_readability_requires_the_grant_genre_and_excludes_english():
    wrong = gw.grant_writing_japanese_readability_score(
        "解析結果を比較し、適用範囲を明らかにした。",
        document_type="research_meeting_manuscript",
    )
    english = gw.grant_writing_japanese_readability_score(
        "This proposal is written in English and must not receive a Japanese score.",
        document_type="grant_proposal",
    )
    grant = gw.grant_writing_japanese_readability_score(
        "設計条件の選択には時間を要する。"
        "本研究では候補順位が一致する条件を明らかにする。"
        "二つの解析法を比較し、適用範囲を判定する。",
        document_type="grant_proposal",
    )

    assert wrong["status"] == "wrong_genre" and wrong["score"] is None
    assert english["status"] == "not_applicable" and english["score"] is None
    assert grant["applicable"] and 0 <= grant["score"] <= 100


def _write_publication_bib(path):
    path.write_text(
        """@article{old,
  author = {Sugahara, Kengo}, title = {Older}, journal = {Journal},
  year = {2022}, volume = {1}, pages = {1--2}, doi = {10.1/old}
}
@article{future,
  author = {Sugahara, Kengo}, title = {Future}, journal = {Journal},
  year = {2027}, volume = {2}, pages = {3--4}, doi = {10.1/future}
}
@inproceedings{conf,
  author = {Sugahara, Kengo}, title = {Conference}, booktitle = {Domestic Meeting},
  year = {2025}, pages = {1--2}
}
@conference{conf2,
  author = {Sugahara, Kengo}, title = {Conference Two}, booktitle = {Meeting},
  year = {2024}, pages = {3--4}, number = {SA-27-0xx}
}
@article{edited,
  editor = {Sugahara, Kengo}, author = {Someone, Else}, title = {Edited},
  journal = {Journal}, year = {2025}
}
""",
        encoding="utf-8",
    )


def test_publication_claims_do_not_duplicate_or_infer_peer_review_or_internationality(tmp_path):
    bib = tmp_path / "papers.bib"
    _write_publication_bib(bib)

    peer = gw.grant_writing_achievement_count_check(
        "査読付き学術論文を12件発表した。",
        bib_path=str(bib),
    )
    international = gw.grant_writing_achievement_count_check(
        "国際会議発表を1件行った。",
        bib_path=str(bib),
    )

    assert peer.count("[要確認]") == 1
    assert "[不一致]" not in peer
    assert international.count("[要確認]") == 1
    assert "BibTeX種別だけでは国際・国内を判定できない" in international


def test_publication_relative_scope_uses_the_current_year_and_lists_each_group_once(tmp_path):
    bib = tmp_path / "papers.bib"
    _write_publication_bib(bib)

    count = gw.grant_writing_achievement_count_check(
        "過去5年の学術論文を2件発表した。",
        bib_path=str(bib),
    )
    listing = gw.grant_writing_publication_list(bib_path=str(bib))

    assert "学術論文 2" in count
    assert listing.count("## 国際会議・研究発表") == 1
    assert "Edited" not in listing
    assert "プレースホルダを含む書誌項目" in listing
    assert "将来年の業績" in listing
