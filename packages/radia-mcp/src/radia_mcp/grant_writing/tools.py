"""Tool functions for grant proposal writing.

All functions are plain callables; ``server.py`` wraps them as MCP tools.
The Japanese technical-prose lint helpers are re-exported from the
grant-writing implementation that was already preserved inside
``paper_writing._ja_lint`` during the public radia-mcp migration.
"""
from __future__ import annotations

import pathlib
import re

from radia_mcp.paper_writing._ja_lint import (  # noqa: F401
    grant_writing_acronym_usage_audit,
    grant_writing_check_kanji_ratio,
    grant_writing_check_misuse_japanese,
    grant_writing_check_notation_variants,
    grant_writing_check_subject_predicate_distance,
    grant_writing_find_undefined_acronyms,
    grant_writing_lint_bedrock,
    grant_writing_suggest_redundancy_fixes,
)

from .._shared.hedges import HEDGE_PATTERNS, scan_hedges


_HERE = pathlib.Path(__file__).resolve().parent


def _load_skill() -> str:
    return (_HERE / "skill.md").read_text(encoding="utf-8")


def _read_text_if_path(text_or_path: str) -> str:
    """Treat a short existing .md/.tex/.txt path as a file, else as text."""
    if text_or_path is None:
        return ""
    s = str(text_or_path)
    if "\n" in s or "\r" in s:
        return s
    if len(s) > 320:
        return s
    try:
        p = pathlib.Path(s)
        if p.exists() and p.suffix.lower() in {".md", ".tex", ".txt"}:
            return p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        pass
    return s


def _prose_for_lint(text: str) -> str:
    """Remove common LaTeX scaffolding before prose-oriented checks.

    Section-presence and program checks still inspect the original source. This
    normalization is only for sentence length, hedge, and Japanese prose lint;
    otherwise template comments and commands are reported as applicant prose.
    """
    if "\\" not in text and not re.search(r"(?m)^\s*%", text):
        return text

    text = re.sub(r"(?m)^\s*%.*$", " ", text)
    text = re.sub(r"\$[^$]*\$", " 数式 ", text)
    text = re.sub(
        r"\\(?:textbf|textit|emph|underline|section|subsection|subsubsection)"
        r"\*?\{([^{}]*)\}",
        r" \1 ",
        text,
    )
    text = re.sub(
        r"\\(?:begin|end|input|include|includegraphics|bibliography|bibliographystyle)"
        r"\*?(?:\[[^\]]*\])?\{[^{}]*\}",
        " ",
        text,
    )
    text = re.sub(r"\\[A-Za-z@]+\*?(?:\[[^\]]*\])?", " ", text)
    text = text.replace("\\\\", " ").replace("{", " ").replace("}", " ")
    return re.sub(r"\s+", " ", text).strip()


def _contains_any(text_lower: str, keywords: list[str]) -> list[str]:
    return [kw for kw in keywords if kw.lower() in text_lower]


def _severity_from_score(score: float | None) -> str:
    if score is None:
        return "UNKNOWN"
    if score < 4:
        return "CRITICAL"
    if score < 6:
        return "HIGH"
    if score < 8:
        return "MEDIUM"
    return "LOW"


def grant_writing_usage() -> str:
    """Return the grant-writing guide."""
    return _load_skill()


def grant_writing_analyze_sentences(text: str, max_len: int = 90) -> dict:
    """Analyze Japanese sentence length for grant proposals.

    Grant drafts can tolerate denser prose than slides, but application
    reviewers still need a clear one-claim-per-sentence rhythm.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    sentences = [s.strip() for s in re.split(r"[。．!?！？]", text) if s.strip()]
    if not sentences:
        return {"error": "no sentences found"}
    lengths = [len(s) for s in sentences]
    long_ones = [
        {"index": i, "length": n, "head": s[:50] + ("..." if len(s) > 50 else "")}
        for i, (s, n) in enumerate(zip(sentences, lengths))
        if n > max_len
    ]
    return {
        "total_sentences": len(sentences),
        "avg_length": round(sum(lengths) / len(lengths), 1),
        "max_length": max(lengths),
        "threshold": max_len,
        "over_threshold_count": len(long_ones),
        "over_threshold_examples": long_ones[:8],
        "target": "avg <= 70-80 chars; one reviewer-relevant claim per sentence",
        "warning": "sentences above 90 chars usually hide two claims",
    }


def grant_writing_count_weak_expressions(text: str) -> dict:
    """Count hedges and grant-specific non-commitment phrases."""
    text = _prose_for_lint(_read_text_if_path(text))
    total, by_pattern = scan_hedges(text)
    patterns = dict(HEDGE_PATTERNS)
    patterns.update({
        "検討する": r"検討(?:する|します|を行う|を進める)",
        "目指す": r"目指(?:す|します)",
        "努める": r"努め(?:る|ます)",
        "など": r"など",
        "今後": r"今後",
    })
    # Recompute with the grant-specific additions included.
    total = 0
    by_pattern = {}
    for name, pat in patterns.items():
        n = len(re.findall(pat, text))
        if n:
            by_pattern[name] = n
            total += n
    return {
        "total_weak_expressions": total,
        "by_pattern": by_pattern,
        "target": "0 in aims / outcomes; replace with deliverables and dates",
        "warning": ">=5 suggests the plan reads as intent rather than execution",
    }


_GENERIC_AXES = {
    "social_or_scientific_problem": ["課題", "問題", "ニーズ", "必要", "社会", "産業"],
    "objective": ["目的", "本研究", "本提案", "実現", "構築", "開発"],
    "novelty": ["新規", "独自", "先行", "特徴", "差別化", "初"],
    "method": ["手法", "方法", "解析", "実験", "検証", "評価"],
    "feasibility": ["実績", "経験", "環境", "体制", "遂行", "可能"],
    "schedule": ["年度", "計画", "スケジュール", "1年目", "2年目", "3年目"],
    "outcomes": ["成果", "公開", "発表", "論文", "レポート", "デモ"],
    "budget": ["予算", "経費", "費用", "助成", "使途", "計算資源"],
}

_KDDI_DIGITAL_AXES = {
    "social_issue": ["社会的課題", "地域", "中小企業", "製造業", "人材", "現場"],
    "social_implementation": ["社会実装", "実装", "PoC", "実証", "試作", "技術プレゼン"],
    "digital_use": ["デジタル", "AI", "生成AI", "MCP", "LTspice", "SPICE", "CAE"],
    "field_validation": ["基板", "計測", "評価", "実験", "熱", "EMC", "厚銅"],
    "open_outputs": ["公開", "OSS", "レポジトリ", "デモ", "教材", "報告書"],
    "schedule": ["年度", "1年目", "2年目", "3年目", "スケジュール", "マイルストーン"],
    "budget_alignment": ["予算", "Claude", "Codex", "Fable", "MDX", "計算資源"],
    "applicant_fit": ["三菱電機", "EMC", "IH", "Radia", "NGSolve", "LTspice", "実績"],
}

_KAKEN_OSS_PLATFORM_AXES = {
    "academic_question": [
        "学術的な問い",
        "学術的問い",
        "研究課題",
        "問い",
        "共同研究サイクル",
    ],
    "technical_reports_as_source": [
        "技術報告",
        "報告書",
        "日本語",
        "知識源",
        "既存知",
    ],
    "lab_silo_and_ai_urgency": [
        "属研究室化",
        "研究室単位",
        "重複実装",
        "重複開発",
        "AIが開発を加速",
        "AIによる開発加速",
    ],
    "jpmars_github_governance": [
        "JP-MARs",
        "GitHub",
        "issue",
        "pull request",
        "CI",
        "ライセンス",
    ],
    "executable_public_outputs": [
        "問題仕様",
        "参照実装",
        "ベンチマーク",
        "教材",
        "mcp-server",
        "ワークショップ",
    ],
    "reuse_and_upstream_first": [
        "既存OSS",
        "再利用",
        "上流貢献",
        "上流で改良",
        "meshio",
        "機能重複",
    ],
    "scientific_quality_gate": [
        "試験データ",
        "期待値",
        "許容誤差",
        "既知制約",
        "参照実装",
        "別機関",
        "CI",
    ],
    "ai_machine_reexecution": [
        "AI",
        "機械実行",
        "再実行",
        "Python/MCP",
        "環境構築支援",
        "依存関係",
    ],
    "tacit_knowledge_conversion": [
        "暗黙知",
        "解析条件",
        "境界条件",
        "メッシュ",
        "収束",
        "失敗例",
    ],
    "domestic_international_collaboration": [
        "国内外",
        "国際",
        "共同実装",
        "コードレビュー",
        "TU Wien",
        "TU Graz",
        "IGTE",
    ],
    "education_and_capacity": [
        "学生",
        "若手",
        "博士",
        "教育",
        "技術力",
        "実装力",
        "検証力",
    ],
    "environment_portability": [
        "mdx",
        "Windows",
        "Linux",
        "スパコン",
        "GPU",
        "ARM",
        "アーキテクチャ",
        "可搬",
        "移植",
    ],
    "maintenance_governance": [
        "運営",
        "保守",
        "CI",
        "リリース",
        "ドキュメント",
        "再現性",
        "貢献者",
    ],
}

_BUDGET_POLICY = (
    "申請予算は遠慮して小さく見せず、助成上限に近い額まで必要な計画として組む。"
    "ただし、上限近くでも不自然に見えないよう、単価・数量・月数/回数・年度配分・"
    "見積根拠を具体的に積算し、検証ループと社会実装に直結する経費として説明する。"
)

_BUDGET_AXIS_COMMENTS = {
    "ai_agent_costs": (
        "AI agent / LLM costs are thin. Tie Claude, Codex, Fable, or related tools "
        "to proposal drafting, design automation, validation, or public deliverables."
    ),
    "compute_resources": (
        "Compute costs are thin. Tie servers, cloud, GPU, or MDX-like resources to "
        "the verification loop and expected run volume."
    ),
    "poc_experiment": (
        "PoC / experiment costs are thin. Tie boards, parts, consumables, measurement, "
        "or evaluation work to concrete implementation evidence."
    ),
    "dissemination": (
        "Dissemination costs are thin. Tie travel, presentations, workshops, or reports "
        "to handoff and social implementation."
    ),
    "near_ceiling_strategy": (
        "予算を遠慮して小さく見せるのではなく、助成上限額に近い申請額が必要である方針を明記する。"
    ),
    "itemized_calculation": (
        "上限近くでも不自然に見えないよう、単価、数量、月数/回数、年度配分、見積根拠を具体的に積算する。"
    ),
}

_POWER_ELECTRONICS_FOCUS_TRIGGERS = [
    "パワーエレクトロニクス",
    "パワエレ",
    "厚銅",
    "LTspice",
    "Radia",
    "NGSolve",
    "CAE-AI",
]

_POWER_ELECTRONICS_FOCUS_AXES = {
    "main_theme_specificity": [
        "パワーエレクトロニクス基板",
        "パワエレ基板",
        "電源基板",
        "厚銅基板",
        "CAE-AI",
    ],
    "board_physics": [
        "寄生",
        "EMC",
        "熱",
        "電流集中",
        "放熱",
        "浮遊容量",
        "インダクタンス",
    ],
    "ai_mcp_loop": ["MCP", "AI", "LLM", "生成AI", "自律駆動", "自然言語"],
    "tool_chain": ["LTspice", "SPICE", "Radia", "NGSolve", "PEEC"],
    "social_target": ["1000人以下", "中小企業", "地域製造業", "製造業"],
    "poc_handoff": ["PoC", "試作", "技術プレゼン", "試験導入", "導入候補", "MotorAI"],
    "commercial_positioning": ["商用CAE", "置換", "入口", "届かない", "習得コスト"],
    "llm_native_advantage": ["Python-native", "現代的", "API", "ツール呼び出し"],
}


def _section_axes_for_program(program: str) -> dict[str, list[str]]:
    if program == "kddi_digital":
        return _KDDI_DIGITAL_AXES
    if program in {"kaken_oss", "kaken_oss_platform"}:
        return _KAKEN_OSS_PLATFORM_AXES
    return _GENERIC_AXES


def grant_writing_section_presence(text: str, program: str = "generic") -> dict:
    """Check whether a proposal draft contains the expected review axes."""
    text = _read_text_if_path(text)
    low = text.lower()
    axes = _section_axes_for_program(program)
    axis_results = {}
    missing = []
    for axis, keywords in axes.items():
        matches = _contains_any(low, keywords)
        ok = bool(matches)
        axis_results[axis] = {
            "ok": ok,
            "matches": matches[:8],
            "keywords": keywords,
        }
        if not ok:
            missing.append(axis)
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "program": program,
        "score": score,
        "axes_total": len(axes),
        "axes_present": len(axes) - len(missing),
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "target": "all core review axes present at least once",
    }


def grant_writing_kddi_digital_check(text: str) -> dict:
    """KDDI Foundation Digital Innovation / social implementation check."""
    text = _read_text_if_path(text)
    presence = grant_writing_section_presence(text, program="kddi_digital")
    comments = []
    for axis in presence["missing_axes"]:
        comments.append(f"Missing KDDI social-implementation axis: {axis}")
    score = presence["score"]
    if "PoC" not in text and "実証" not in text and "試作" not in text:
        comments.append("Add a concrete PoC / field validation output.")
        score = max(0.0, score - 1.0)
    if "公開" not in text and "OSS" not in text and "レポジトリ" not in text:
        comments.append("Add a shareable output: repository, demo, report, or workshop material.")
        score = max(0.0, score - 1.0)
    return {
        "score": round(score, 1),
        "missing_required_count": len(comments),
        "comments": comments,
        "axis_results": presence["axis_results"],
        "source": "KDDI Digital Innovation social-implementation proposal axes",
    }


def grant_writing_kddi_power_electronics_focus_check(text: str) -> dict:
    """Check the current KDDI power-electronics-board CAE-AI framing.

    This is intentionally domain-specific.  Use it when the KDDI Digital
    Innovation proposal is about LTspice/Radia/NGSolve/MCP for power
    electronics boards.  It checks that "companies below 1000 employees" is
    kept as the implementation target, not as the main theme.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    axis_results = {}
    missing = []
    for axis, keywords in _POWER_ELECTRONICS_FOCUS_AXES.items():
        matches = _contains_any(low, keywords)
        axis_results[axis] = {
            "ok": bool(matches),
            "matches": matches[:8],
            "keywords": keywords,
        }
        if not matches:
            missing.append(axis)

    score = round(
        10.0 * (len(_POWER_ELECTRONICS_FOCUS_AXES) - len(missing))
        / len(_POWER_ELECTRONICS_FOCUS_AXES),
        1,
    )
    comments = [f"Power-electronics focus axis missing or thin: {axis}" for axis in missing]

    sentences = [s.strip() for s in re.split(r"[。．!?！？]", text) if s.strip()]
    theme_sentences = [
        s for s in sentences
        if any(token in s for token in ("主題", "目的", "テーマ案", "テーマ", "本研究"))
    ][:5]
    if theme_sentences:
        theme_blob = "。".join(theme_sentences)
        if (
            "パワーエレクトロニクス基板" not in theme_blob
            and "パワエレ基板" not in theme_blob
            and "CAE-AI" not in theme_blob
        ):
            comments.append(
                "The opening theme sentences do not clearly state the power-electronics-board CAE-AI subject."
            )
            score = max(0.0, score - 1.0)

    risky_phrases = []
    for phrase in ("一般的なCAE導入", "すべての製造業CAE", "商用CAEの代替"):
        if phrase in text:
            risky_phrases.append(phrase)
    if "商用CAE" in text and "置換する" in text and "直ちに置換するものではない" not in text:
        risky_phrases.append("商用CAEを置換する")
    if risky_phrases:
        comments.append(
            "Avoid broad or adversarial positioning; keep commercial CAE as powerful but hard to access: "
            + ", ".join(risky_phrases)
        )
        score = max(0.0, score - 1.0)

    recommendations = [
        "Make the subject: power-electronics-board circuit/EM/thermal CAE-AI environment.",
        "Keep companies below 1000 employees as the first user and implementation field.",
        "Frame commercial CAE as powerful but expensive/difficult; the proposal creates an AI/MCP entry point.",
        "Tie the PoC to a board-level outcome: parasitics, EMC risk, heat, measurement, and handoff.",
    ]
    return {
        "score": round(score, 1),
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "comments": comments,
        "theme_sentence_examples": theme_sentences,
        "recommendations": recommendations,
        "target": "specific power-electronics-board CAE-AI theme; 1000-person firms as implementation target",
        "source": "KDDI power-electronics-board CAE-AI framing check",
    }


def grant_writing_kaken_oss_platform_check(text: str) -> dict:
    """Check KAKENHI framing for an AI-era OSS research platform proposal."""
    text = _read_text_if_path(text)
    low = text.lower()
    presence = grant_writing_section_presence(text, program="kaken_oss")
    comments = [
        f"KAKEN OSS platform axis missing or thin: {axis}"
        for axis in presence["missing_axes"]
    ]
    score = presence["score"]

    public_outputs = [
        "公開リポジトリ",
        "問題仕様",
        "参照実装",
        "ベンチマーク",
        "software.html",
        "elemag/index.php",
        "geometry.php",
        "mcp-server",
    ]
    output_hits = _contains_any(low, public_outputs)
    if len(output_hits) < 4:
        comments.append(
            "Name at least four concrete public outputs: repository, problem specs, "
            "reference implementations, benchmarks, teaching pages, or mcp-server."
        )
        score = max(0.0, score - 1.0)

    if "技術報告" not in text and "報告書" not in text:
        comments.append(
            "Treat existing Japanese technical reports as a research knowledge source, "
            "then connect them to executable assets and review history."
        )
        score = max(0.0, score - 0.5)

    silo_hits = _contains_any(
        low,
        ["属研究室化", "研究室単位", "研究室内に閉じ", "重複実装", "重複開発"],
    )
    ai_urgency_hits = _contains_any(
        low,
        ["AIが開発を加速", "AIによる開発加速", "AI時代", "高速に再生産", "増幅"],
    )
    if not silo_hits or not ai_urgency_hits:
        comments.append(
            "Make the why-now explicit: AI can accelerate duplicated, lab-siloed "
            "implementations unless shared specifications, tests, and review convert "
            "them into collaborative assets."
        )
        score = max(0.0, score - 0.5)

    reuse_hits = _contains_any(
        low,
        ["既存OSS", "再利用", "上流貢献", "上流で改良", "upstream", "meshio", "機能重複"],
    )
    if not reuse_hits:
        comments.append(
            "Add an upstream-first gate: survey existing OSS, reuse or contribute "
            "upstream, and justify new code with a tested domain-specific gap."
        )
        score = max(0.0, score - 1.0)

    quality_groups = {
        "test_corpus": ["試験データ", "テストデータ", "検証データ"],
        "expected_tolerance": ["期待値", "許容誤差", "参照値"],
        "automated_validation": ["CI", "自動テスト", "再実行"],
        "scope_and_provenance": ["既知制約", "適用範囲", "由来", "出典"],
        "independent_reproduction": ["別機関", "第三者", "独立検証", "再現確認"],
    }
    quality_results = {
        name: _contains_any(low, keywords)
        for name, keywords in quality_groups.items()
    }
    if sum(bool(hits) for hits in quality_results.values()) < 4:
        comments.append(
            "Do not treat GitHub as a repository dump. Require test provenance, "
            "expected values/tolerances, CI, documented limits, and independent re-execution."
        )
        score = max(0.0, score - 1.0)

    responsibility_phrases = ["利用者の自己責任", "使う人の自己責任"]
    responsibility_negations = [
        "自己責任を品質保証の代替にしない",
        "自己責任を理由にreleaseしない",
        "利用者への責任転嫁ではなく",
    ]
    if (
        _contains_any(low, responsibility_phrases)
        and not _contains_any(low, responsibility_negations)
    ):
        comments.append(
            "A no-warranty license does not replace scientific validation; do not "
            "shift supported-scope verification responsibility to users."
        )
        score = max(0.0, score - 1.0)

    domestic_evidence = ["伊田", "HACApK", "CEFC", "ADVENTURE", "国内共同"]
    overseas_evidence = [
        "Hollaus",
        "TU Wien",
        "TU Graz",
        "IGTE",
        "openCFS",
        "NGSolve",
    ]
    domestic_hits = _contains_any(low, domestic_evidence)
    overseas_hits = _contains_any(low, overseas_evidence)
    if not domestic_hits:
        comments.append(
            "Add named domestic preliminary evidence such as the Ida/HACApK/CEFC "
            "collaboration or an ADVENTURE connection."
        )
        score = max(0.0, score - 0.5)
    if not overseas_hits:
        comments.append(
            "Add named overseas preliminary evidence such as Hollaus/TU Wien/IGTE, "
            "TU Graz/openCFS, or the NGSolve community."
        )
        score = max(0.0, score - 0.5)

    catch_up_negations = [
        "追いつくことではない",
        "追いつく計画ではない",
        "追いつく計画ではなく",
        "追いつくことを目的にしない",
        "欧州を一方向の到達目標",
    ]
    if "追いつく" in text and not any(phrase in text for phrase in catch_up_negations):
        comments.append(
            "Avoid a catch-up frame; state reciprocal collaboration and "
            "complementarity instead of Europe as a one-way benchmark."
        )
        score = max(0.0, score - 0.5)

    environment_groups = {
        "enterprise_os": ["Windows", "Linux", "企業"],
        "shared_compute": ["mdx", "クラウド", "仮想環境"],
        "hpc": ["スパコン", "スーパーコンピュータ", "HPC"],
        "accelerator_generation": ["GPU", "A100", "H200", "ICCG"],
        "architecture_change": ["ARM", "arm64", "アーキテクチャ", "可搬", "移植"],
    }
    environment_results = {
        name: _contains_any(low, keywords)
        for name, keywords in environment_groups.items()
    }
    if sum(bool(hits) for hits in environment_results.values()) < 3:
        comments.append(
            "Frame compute as a portable execution matrix spanning enterprise "
            "Windows/Linux, mdx, HPC, GPU generations, and future CPU architectures."
        )
        score = max(0.0, score - 1.0)

    hardware_purchase_terms = ["計算機を購入", "マシンを購入", "GPUを購入"]
    portability_terms = ["可搬", "移植", "再実行", "複数環境", "アーキテクチャ"]
    if _contains_any(low, hardware_purchase_terms) and not _contains_any(low, portability_terms):
        comments.append(
            "Do not make hardware acquisition the research objective; explain how "
            "specifications, containers, CI, and AI-assisted setup survive platform change."
        )
        score = max(0.0, score - 1.0)

    platform_mentions = len(re.findall(r"JP[-‐‑–—]?MARs", text, flags=re.IGNORECASE))
    radia_mentions = len(
        re.findall(r"(?<![A-Za-z])Radia(?:/radia-mcp|-mcp)?", text, flags=re.IGNORECASE)
    )
    platform_focus = {
        "ok": platform_mentions >= max(1, radia_mentions),
        "jpmars_mentions": platform_mentions,
        "radia_mentions": radia_mentions,
    }
    if not platform_focus["ok"]:
        comments.append(
            "Keep JP-MARs as the governing research platform; position Radia and "
            "other named software as preliminary evidence or demonstrators."
        )
        score = max(0.0, score - 1.0)

    radia_integration = {
        "mentioned": bool(radia_mentions),
        "positioned_as_integration_evidence": bool(
            radia_mentions
            and _contains_any(low, ["統合", "融合", "接続", "先行実証", "インタフェース"])
        ),
    }
    if radia_integration["mentioned"] and not radia_integration["positioned_as_integration_evidence"]:
        comments.append(
            "When Radia is named, use it as preliminary evidence that specialized "
            "software and methods can be integrated, not as the proposal's main product."
        )
        score = max(0.0, score - 0.5)

    recommendations = [
        "Open with how technical reports become reviewable, executable research assets in the AI era.",
        "Use GitHub/JP-MARs to preserve issue, implementation, validation, review, and contributor history.",
        "Define lab-level siloing and explain why AI otherwise accelerates duplicated and unverified development.",
        "Survey existing OSS first; reuse it or contribute upstream before authoring a new implementation.",
        "Separate license disclaimers from scientific responsibility and require test data, tolerances, limits, CI, and independent re-execution.",
        "Show both domestic and overseas preliminary collaborations with names and outputs.",
        "Treat AI as an accelerator for search, translation, environment setup, and re-execution; researchers retain physical validation.",
        "Define a portable execution matrix across Windows/Linux, mdx, HPC, GPU generations, and future architectures such as ARM.",
        "Split platform adaptation into issues and CI jobs so institutions and companies can maintain different environments through pull requests.",
        "Keep hardware acquisition subordinate to reproducible specifications, CI, containers, and long-term governance.",
    ]
    return {
        "score": round(score, 1),
        "missing_count": presence["missing_count"],
        "missing_axes": presence["missing_axes"],
        "axis_results": presence["axis_results"],
        "public_output_hits": output_hits,
        "domestic_evidence_hits": domestic_hits,
        "overseas_evidence_hits": overseas_hits,
        "environment_results": environment_results,
        "silo_hits": silo_hits,
        "ai_urgency_hits": ai_urgency_hits,
        "reuse_hits": reuse_hits,
        "quality_results": quality_results,
        "platform_focus": platform_focus,
        "radia_integration": radia_integration,
        "comments": comments,
        "recommendations": recommendations,
        "target": "KAKENHI OSS platform: overcome lab-level siloing through upstream-first reuse, scientific quality gates, reciprocal collaboration, and AI-assisted portability",
        "source": "KAKENHI AI-era OSS research-platform framing check",
    }


def grant_writing_internal_evidence_to_external_scale_check(text: str) -> dict:
    """Check whether an internal success is evidence for external transfer.

    The check is optional: it becomes applicable only when a draft claims an
    existing pilot, operation, or preliminary result. It then asks what is
    transferred, to whom, through which route, and how independent success is
    verified. This keeps the rule useful across grant programs and domains.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    axes = {
        "internal_evidence": [
            "予備成果",
            "運用実績",
            "実証済",
            "試行済",
            "既に実現",
            "すでに実現",
            "既に運用",
            "すでに運用",
            "機能している",
            "利用している",
            "活用している",
            "導入済",
        ],
        "observed_users_or_context": [
            "利用者",
            "ユーザー",
            "現場",
            "研究室",
            "学内",
            "社内",
            "組織内",
            "教員",
            "学生",
            "担当者",
            "企業",
        ],
        "transferable_unit": [
            "問題仕様",
            "参照実装",
            "入力データ",
            "試験データ",
            "検証結果",
            "実装",
            "コード",
            "手順",
            "API",
            "教材",
            "ドキュメント",
            "ベンチマーク",
            "知識基盤",
            "データセット",
        ],
        "external_route": [
            "他機関",
            "別機関",
            "研究室間",
            "複数機関",
            "国内外",
            "国際",
            "共同研究",
            "外部展開",
            "公開",
            "社会実装",
            "技術移転",
            "普及",
            "上流貢献",
        ],
        "external_actor": [
            "第三者",
            "他機関",
            "別機関",
            "外部利用者",
            "共同研究先",
            "企業",
            "海外研究者",
        ],
        "independent_validation": [
            "再実行",
            "再現",
            "比較",
            "評価",
            "検証",
            "レビュー",
            "受入",
            "フィードバック",
        ],
    }
    axis_results = {}
    for name, keywords in axes.items():
        hits = _contains_any(low, keywords)
        axis_results[name] = {
            "ok": bool(hits),
            "matches": hits[:8],
            "keywords": keywords,
        }

    applicable = axis_results["internal_evidence"]["ok"]
    if not applicable:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": axis_results,
            "comments": [],
            "target": (
                "when internal evidence is claimed, connect it to a transferable "
                "unit, external route and actor, and independent validation"
            ),
            "source": "generic internal-evidence-to-external-scale check",
        }

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "observed_users_or_context": "内部実証を誰がどの環境で利用したかを示す。",
        "transferable_unit": (
            "外部へ渡す単位を、仕様、実装、データ、手順等として明示する。"
        ),
        "external_route": (
            "内部実証を他機関・共同研究・社会実装へ接続する経路を示す。"
        ),
        "external_actor": "内部関係者ではない利用者・検証者を明示する。",
        "independent_validation": (
            "外部での再実行、比較、評価等により移転の成否を判定する。"
        ),
    }
    comments = [comments_by_axis[name] for name in missing if name in comments_by_axis]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "comments": comments,
        "target": (
            "internal success is feasibility evidence, not external validity; state "
            "what travels, to whom, by which route, and how success is verified"
        ),
        "source": "generic internal-evidence-to-external-scale check",
    }


def grant_writing_budget_alignment_check(text: str) -> dict:
    """Check that budget items are tied to verification and implementation."""
    text = _read_text_if_path(text)
    low = text.lower()
    axes = {
        "ai_agent_costs": ["claude", "codex", "fable", "生成ai", "llm", "ai"],
        "compute_resources": ["mdx", "計算資源", "gpu", "クラウド", "サーバ"],
        "poc_experiment": ["試作", "基板", "部品", "消耗品", "計測", "評価"],
        "dissemination": ["旅費", "発表", "技術プレゼン", "ワークショップ", "報告"],
        "near_ceiling_strategy": [
            "上限",
            "助成上限",
            "上限額",
            "限度額",
            "満額",
            "ほぼ上限",
            "上限いっぱい",
            "上限近く",
        ],
        "itemized_calculation": [
            "内訳",
            "単価",
            "数量",
            "月数",
            "回数",
            "年度配分",
            "見積",
            "積算",
            "算出",
            "根拠",
            "税込",
            "税抜",
        ],
    }
    results = {}
    missing = []
    for axis, keywords in axes.items():
        matches = _contains_any(low, keywords)
        results[axis] = {"ok": bool(matches), "matches": matches, "keywords": keywords}
        if not matches:
            missing.append(axis)
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    comments = [
        _BUDGET_AXIS_COMMENTS.get(axis, f"Budget rationale missing or thin: {axis}")
        for axis in missing
    ]
    if "効率化" in text and "検証" not in text and "PoC" not in text:
        comments.append("AI費用が一般的な効率化に見える。検証ループの実行経費として説明する。")
        score = max(0.0, score - 1.0)
    return {
        "score": round(score, 1),
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": results,
        "comments": comments,
        "budget_policy": _BUDGET_POLICY,
        "target": (
            "every major cost maps to AI/tool execution, compute, PoC, or dissemination; "
            "the requested amount may be close to the ceiling when itemized and justified"
        ),
    }


def grant_writing_recommendation_letter_template(
    program: str = "kddi_digital",
    applicant: str = "菅原賢悟准教授",
) -> str:
    """Return a one-page recommendation-letter draft template."""
    if program != "kddi_digital":
        return (
            f"本学所属の{applicant}による助成申請を推薦いたします。\n\n"
            "同氏は、研究課題の遂行に必要な専門性、研究実績、研究環境を備えており、"
            "本助成によって得られる成果は当該分野の発展に資するものです。"
            "本学としても、研究実施に必要な環境整備および研究活動を支援いたします。\n\n"
            "以上の理由により、本申請を推薦いたします。"
        )
    return (
        f"本学理工学部 電気電子通信工学科 {applicant}による、"
        "KDDI財団デジタルイノベーション社会実装助成への申請を推薦いたします。\n\n"
        "同氏は、三菱電機におけるPCB EMCを含む電磁気設計・CAE実務、"
        "本学における電磁界解析・熱解析研究、ならびにRadia、NGSolve、LTspice等を"
        "MCPサーバを介してAI駆動で扱う研究開発環境の構築実績を有しています。"
        "本申請は、地域製造業における回路・電磁界・熱協調設計の属人性を低減し、"
        "PoC基板評価と公開可能な設計支援基盤を通じて社会実装へ接続するものです。\n\n"
        "本学は、同氏の研究室が有する実験・評価環境、学生を含む研究実施体制、"
        "および外部連携活動を支援し、本助成期間中の計画遂行を後押しします。"
        "以上の理由により、本申請を推薦いたします。"
    )


def grant_writing_health_report(
    text_or_path: str,
    program: str = "generic",
    skip: str = "",
) -> dict:
    """Integrated grant-writing health report.

    Args:
        text_or_path: Proposal text or an existing .md/.tex/.txt path.
        program: ``generic``, ``kddi_digital``, or ``kaken_oss``.
        skip: comma-separated tool ids to skip, e.g. ``sentence,budget``.
    """
    text = _read_text_if_path(text_or_path)
    skip_set = {s.strip().lower() for s in skip.split(",") if s.strip()}

    detailed_results: dict[str, dict] = {}
    detailed_scores: dict[str, float] = {}
    priority_issues: list[dict] = []

    if "sections" not in skip_set:
        sections = grant_writing_section_presence(text, program=program)
        detailed_results["sections"] = sections
        detailed_scores["sections"] = sections["score"]
        for axis in sections["missing_axes"]:
            priority_issues.append({
                "tool": "sections",
                "name": "section_presence",
                "severity": "HIGH",
                "score": sections["score"],
                "comments": [f"Add reviewer-visible content for: {axis}"],
            })

    if program == "kddi_digital" and "kddi" not in skip_set:
        kddi = grant_writing_kddi_digital_check(text)
        detailed_results["kddi_digital"] = kddi
        detailed_scores["kddi_digital"] = kddi["score"]
        if kddi["comments"]:
            priority_issues.append({
                "tool": "kddi",
                "name": "kddi_digital_check",
                "severity": _severity_from_score(kddi["score"]),
                "score": kddi["score"],
                "comments": kddi["comments"][:5],
            })

    if (
        program == "kddi_digital"
        and "focus" not in skip_set
        and _contains_any(text.lower(), _POWER_ELECTRONICS_FOCUS_TRIGGERS)
    ):
        focus = grant_writing_kddi_power_electronics_focus_check(text)
        detailed_results["power_electronics_focus"] = focus
        detailed_scores["power_electronics_focus"] = focus["score"]
        if focus["comments"]:
            priority_issues.append({
                "tool": "focus",
                "name": "kddi_power_electronics_focus_check",
                "severity": _severity_from_score(focus["score"]),
                "score": focus["score"],
                "comments": focus["comments"][:5],
            })

    if program in {"kaken_oss", "kaken_oss_platform"} and "kaken" not in skip_set:
        kaken = grant_writing_kaken_oss_platform_check(text)
        detailed_results["kaken_oss_platform"] = kaken
        detailed_scores["kaken_oss_platform"] = kaken["score"]
        if kaken["comments"]:
            priority_issues.append({
                "tool": "kaken",
                "name": "kaken_oss_platform_check",
                "severity": _severity_from_score(kaken["score"]),
                "score": kaken["score"],
                "comments": kaken["comments"][:5],
            })

    if "scale" not in skip_set:
        scale = grant_writing_internal_evidence_to_external_scale_check(text)
        detailed_results["internal_to_external_scale"] = scale
        if scale["applicable"]:
            detailed_scores["internal_to_external_scale"] = scale["score"]
            if scale["comments"]:
                priority_issues.append({
                    "tool": "scale",
                    "name": "internal_evidence_to_external_scale_check",
                    "severity": _severity_from_score(scale["score"]),
                    "score": scale["score"],
                    "comments": scale["comments"][:5],
                })

    if "sentence" not in skip_set:
        sent = grant_writing_analyze_sentences(text)
        detailed_results["sentence"] = sent
        if "error" not in sent:
            over = sent["over_threshold_count"]
            score = max(0.0, round(10.0 - min(over, 8) * 1.0, 1))
            detailed_scores["sentence"] = score
            if over:
                priority_issues.append({
                    "tool": "sentence",
                    "name": "analyze_sentences",
                    "severity": "MEDIUM" if over < 4 else "HIGH",
                    "score": score,
                    "comments": [f"{over} sentence(s) exceed the threshold."],
                })

    if "weak" not in skip_set:
        weak = grant_writing_count_weak_expressions(text)
        detailed_results["weak"] = weak
        score = max(0.0, round(10.0 - min(weak["total_weak_expressions"], 10), 1))
        detailed_scores["weak"] = score
        if weak["total_weak_expressions"]:
            priority_issues.append({
                "tool": "weak",
                "name": "count_weak_expressions",
                "severity": "MEDIUM" if weak["total_weak_expressions"] < 5 else "HIGH",
                "score": score,
                "comments": [f"Weak / non-committal expressions: {weak['by_pattern']}"],
            })

    if "bedrock" not in skip_set:
        bedrock = grant_writing_lint_bedrock(_prose_for_lint(text))
        detailed_results["bedrock"] = bedrock
        issue_count = bedrock.get("issue_count", 0)
        score = max(0.0, round(10.0 - issue_count * 1.5, 1))
        detailed_scores["bedrock"] = score
        if issue_count:
            priority_issues.append({
                "tool": "bedrock",
                "name": "lint_bedrock",
                "severity": _severity_from_score(score),
                "score": score,
                "comments": [i.get("rule", "issue") for i in bedrock.get("issues", [])[:5]],
            })

    if "budget" not in skip_set:
        budget = grant_writing_budget_alignment_check(text)
        detailed_results["budget"] = budget
        detailed_scores["budget"] = budget["score"]
        if budget["comments"]:
            priority_issues.append({
                "tool": "budget",
                "name": "budget_alignment_check",
                "severity": _severity_from_score(budget["score"]),
                "score": budget["score"],
                "comments": budget["comments"][:5],
            })

    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2, "UNKNOWN": 3, "LOW": 4}
    priority_issues.sort(
        key=lambda x: (sev_rank.get(x["severity"], 99), -1 * (x.get("score") or 0))
    )

    scores = list(detailed_scores.values())
    overall = round(sum(scores) / len(scores), 1) if scores else 0.0
    severity = _severity_from_score(overall)
    if overall >= 8:
        summary = "Submission logic is mostly visible; polish evidence and wording."
    elif overall >= 6:
        summary = "Core story is present, but reviewer-facing gaps remain."
    else:
        summary = "Proposal needs clearer axes, deliverables, or budget-to-verification logic."

    return {
        "overall_score": overall,
        "score_max": 10,
        "overall_severity": severity,
        "summary_comment": summary,
        "program": program,
        "detailed_scores": detailed_scores,
        "detailed_results": detailed_results,
        "priority_issues": priority_issues,
        "tools_run": sorted(detailed_results),
        "tools_skipped": sorted(skip_set),
        "total_findings": len(priority_issues),
        "hint": "Use priority_issues first; do not optimize the score mechanically.",
        "source": "radia_mcp.grant_writing public document server",
    }
