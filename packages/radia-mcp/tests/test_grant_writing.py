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
    assert "named_software_abstraction" in report["detailed_results"]
    assert "reviewer_vocabulary" in report["detailed_results"]
    assert "persuasion_quality" in report["detailed_results"]
    assert "literature_gap_evidence" in report["detailed_results"]
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
