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
            payload = p.read_bytes()
            for encoding in ("utf-8-sig", "cp932"):
                try:
                    return payload.decode(encoding)
                except UnicodeDecodeError:
                    continue
            return payload.decode("utf-8", errors="replace")
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
    text = re.sub(
        r"\\begin\{(?P<figenv>figure\*?|center)\}.*?"
        r"\\includegraphics.*?\\end\{(?P=figenv)\}",
        " ",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(r"(?<!\\)\\\[.*?(?<!\\)\\\]", " 数式 ", text, flags=re.DOTALL)
    text = re.sub(
        r"\\begin\{(?:equation\*?|align\*?|gather\*?)\}.*?"
        r"\\end\{(?:equation\*?|align\*?|gather\*?)\}",
        " 数式 ",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(r"\$\$.*?\$\$", " 数式 ", text, flags=re.DOTALL)
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
        "AIにより開発が加速",
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
    "外部料金は公式料金表または見積書に基づき、提供者、URL、料金年度・改定日、参照日、"
    "税込/税抜、最低購入単位、有効期限、通貨・為替、端数処理を記録する。申請年度の料金が"
    "未公表なら現行料金による暫定積算と再確認時期を示し、公式単価、価格変動への予備幅、"
    "使用上限を区別する。費目別・年度別の総額を再計算し、直接経費上限、研究期間、間接経費"
    "の扱いと照合する。旅費は開催地の公表状況を区別し、人数、泊数、運賃、登録費、宿泊費、"
    "日当・現地交通へ分解する。"
    "最大費目を明示して中心的な研究行為へ対応付け、その他費目はサービス別に分解し、"
    "サーバ、CI、保存等の二重計上を避ける。AIと計算資源を併用する場合は、AI推論、"
    "常時稼働の低費用検証、短期集中計算を分け、各層の単位、期間、成果物を示す。"
    "GPU等の機種名を記す場合は、アクセラレータとホストCPUを公式仕様で区別し、"
    "実在する構成だけを予算化する。高速化の成功でなく、異機種間で保存すべき物理量・"
    "設計判断と不一致時の扱いを示す。"
    "科研費の基盤系種目では、採択時に申請額の約7割程度へ減額されて内定する場合が"
    "多い(充足率)。減額後も検証ループが成立する経費の優先順位を用意する。"
    "挑戦的研究は原則満額支給である。"
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
    "pricing_provenance": (
        "外部料金は公式料金表または見積書へ遡れるよう、URL、料金年度・参照日、税込/税抜、"
        "最低購入単位、有効期限、通貨・為替、端数処理を記録する。"
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


_KAKEN_REVIEW_CRITERIA_AXES = {
    "academic_importance": [
        "学術的重要性",
        "学術的意義",
        "学術的問い",
        "独自性",
        "独創性",
        "波及",
    ],
    "method_validity": [
        "研究方法",
        "研究計画",
        "妥当",
        "検証",
        "実証",
        "評価方法",
        "達成指標",
    ],
    "feasibility_environment": [
        "遂行能力",
        "研究環境",
        "研究実績",
        "準備状況",
        "実施体制",
        "研究設備",
        "予備",
    ],
}

_KAKEN_BRIEFING_NOTES = [
    "審査基準は3つ: (1)研究課題の学術的重要性、(2)研究方法の妥当性、"
    "(3)研究遂行能力及び研究環境の適切性。",
    "審査委員は約1ヶ月で多い場合100件程度の計画調書を審査する。"
    "専門外の読者でも読みやすい調書が圧倒的に採択されやすい。",
    "カラーの図・写真は審査時に白黒印刷される種目がある。"
    "色の違いだけで系列を区別しない。",
    "審査ではresearchmapが研究者番号で参照される。"
    "応募前に更新と研究者番号の登録を確認する。",
    "「人権の保護及び法令等の遵守への対応」欄は例年審査委員からの指摘が"
    "非常に多い。該当なしの場合も判断根拠を一文添える。",
    "基盤系種目は申請額の約7割程度への減額内定が多い(充足率)。"
    "挑戦的研究は原則満額支給だが採択率が低く、基盤研究との重複応募を検討する。",
]


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
        [
            "AIが開発を加速",
            "AIによる開発加速",
            "AIにより開発が加速",
            "AI時代",
            "高速に再生産",
            "増幅",
        ],
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


def grant_writing_domain_outcome_chain_check(text: str) -> dict:
    """Check that a platform or tool proposal ends in domain knowledge.

    Infrastructure, OSS, APIs, repositories, and AI interfaces can make a
    study feasible, but they are rarely the academic outcome by themselves.
    This optional check becomes applicable when such enabling technology is
    prominent and asks for the target object, measurable quantity, conditional
    knowledge product, changed decision, falsification gate, and an explicit
    statement that the technology is subordinate to the research question.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    tool_terms = [
        "プラットフォーム",
        "研究基盤",
        "ソフトウェア",
        "OSS",
        "github",
        "gitlab",
        "MCP",
        "API",
        "リポジトリ",
        "モジュール",
        "インターフェース",
        "インタフェース",
    ]
    tool_hits = _contains_any(low, tool_terms)
    if not tool_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "tool_hits": [],
            "comments": [],
            "target": "enabling technology must terminate in a field-specific, falsifiable knowledge product",
            "source": "generic tool-to-domain-outcome chain check",
        }

    axes = {
        "domain_object": [
            "機器",
            "装置",
            "材料",
            "回路",
            "磁気浮上",
            "電動機",
            "患者",
            "診断",
            "製造工程",
            "観測対象",
            "制御対象",
        ],
        "measurable_domain_quantity": [
            "設計量",
            "目的量",
            "性能",
            "損失",
            "力",
            "温度",
            "平衡",
            "剛性",
            "感度",
            "精度",
            "収率",
            "寿命",
            "制約",
        ],
        "conditional_knowledge_product": [
            "選択則",
            "設計則",
            "成立条件",
            "成立域",
            "適用範囲",
            "適用境界",
            "不能条件",
            "選択域",
            "保存条件",
            "変化する条件",
            "閾値",
        ],
        "decision_consequence": [
            "設計判断",
            "設計順位",
            "順位を確定",
            "採否",
            "選択する",
            "切り替",
            "進むべき",
            "停止条件",
            "追加解析",
            "高忠実度",
        ],
        "falsifiable_gate": [
            "合格条件",
            "不合格",
            "反証",
            "許容差",
            "包含",
            "棄却",
            "失敗条件",
            "適用不能",
            "満たさない場合",
            "主張しない",
        ],
        "tool_subordination": [
            "問いではなく",
            "手段である",
            "手段とし",
            "実装手段",
            "検証手段",
            "道具である",
            "を用いる",
            "が支援する",
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

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "domain_object": "基盤を用いて研究する具体的な対象物・現象を明示する。",
        "measurable_domain_quantity": "対象分野で測る設計量・性能量・制約を明示する。",
        "conditional_knowledge_product": (
            "構築物でなく、成立条件、境界、選択則等の新しい知見を成果にする。"
        ),
        "decision_consequence": (
            "得られる知見が、どの設計・診断・採否を変えるかまで書く。"
        ),
        "falsifiable_gate": "知見を棄却または限定する事前の判定条件を置く。",
        "tool_subordination": (
            "OSS、AI、API、MCP、GitHub等は問いでなく検証手段だと明記する。"
        ),
    }
    comments = [comments_by_axis[name] for name in missing]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "tool_hits": tool_hits[:12],
        "comments": comments,
        "target": (
            "show a complete chain from enabling technology to a measured "
            "domain quantity, conditional knowledge, changed decision, and "
            "falsification gate"
        ),
        "source": "generic tool-to-domain-outcome chain check",
    }


def grant_writing_derived_metric_validation_check(text: str) -> dict:
    """Check calibration and falsification of a proposal-specific metric.

    A newly named interval, score, index, or composite criterion can be useful,
    but it must not be tuned and judged on the same examples. This optional
    check asks for an operational definition, observable components, a bounded
    calibration set, pre-test freezing, held-out validation, an acceptance
    threshold, and a consequence for failure.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    metric_terms = [
        "判定区間",
        "許容区間",
        "信頼区間",
        "評価指標",
        "判定指標",
        "複合指標",
        "合成指標",
        "リスクスコア",
        "評価スコア",
        "選択指数",
        "安全係数",
    ]
    metric_hits = _contains_any(low, metric_terms)
    if not metric_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "metric_hits": [],
            "comments": [],
            "target": "define, calibrate, freeze, independently validate, and falsify a proposal-specific metric",
            "source": "generic derived-metric validation check",
        }

    axes = {
        "operational_definition": [
            "定義する",
            "定義した",
            "算出式",
            "計算式",
            "合成則",
            "\\Delta",
            "I_q",
            "重み付き",
        ],
        "observable_components": [
            "対計算差",
            "成分",
            "残差",
            "測定値",
            "求解差",
            "離散化差",
            "連成差",
            "入力変数",
            "説明変数",
        ],
        "bounded_calibration_set": [
            "校正集合",
            "校正ケース",
            "校正データ",
            "参照データ",
            "基準問題",
            "訓練集合",
            "学習データ",
        ],
        "pretest_freeze": [
            "事前登録",
            "事前に固定",
            "事前固定",
            "凍結",
            "検証前に固定",
            "候補比較前に固定",
            "最終候補を得る前",
        ],
        "heldout_validation": [
            "保留データ",
            "独立検証",
            "検証集合",
            "テスト集合",
            "未使用データ",
            "外部検証",
            "holdout",
            "held-out",
        ],
        "acceptance_threshold": [
            "合格条件",
            "全点包含",
            "包含率",
            "許容差以内",
            "以上を合格",
            "以下を合格",
            "閾値",
            "判定基準",
        ],
        "failure_consequence": [
            "再校正せず",
            "主張しない",
            "不合格",
            "棄却",
            "適用境界",
            "反例",
            "失敗条件",
            "満たさない場合",
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

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "operational_definition": "新しい区間・指数・スコアの算出式または手順を定義する。",
        "observable_components": "指標を構成する観測可能な成分を列挙する。",
        "bounded_calibration_set": "校正に使うケース、範囲、件数を固定する。",
        "pretest_freeze": "独立検証を見る前に式・係数・閾値を凍結する。",
        "heldout_validation": "校正に使わない保留データまたは外部データで検証する。",
        "acceptance_threshold": "包含率、許容差、件数等の合否閾値を事前に置く。",
        "failure_consequence": "不合格時に再調整で救済せず、棄却・限定・代替経路を定める。",
    }
    comments = [comments_by_axis[name] for name in missing]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "metric_hits": metric_hits[:12],
        "comments": comments,
        "target": (
            "a proposal-specific metric must have an operational definition, "
            "separate calibration and held-out validation, a frozen acceptance "
            "rule, and a failure consequence"
        ),
        "source": "generic derived-metric validation check",
    }


def grant_writing_cross_organization_pilot_check(text: str) -> dict:
    """Check whether preliminary evidence crosses an organization boundary.

    User counts, repository publication, links, and internal operation are weak
    evidence for an inter-organizational research workflow. Around a claimed
    pilot or preparation result, this check looks for a named cross-boundary
    actor, a transferred artifact, a bounded task, an observed outcome, an
    independent action, and an explicit statement of what remains unproven.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    low = text.lower()
    prior_terms = [
        "予備実証",
        "予備試験",
        "予備成果",
        "準備状況",
        "既往実績",
        "実証済",
        "再実行した",
        "統合済",
    ]
    prior_hits = _contains_any(low, prior_terms)
    if not prior_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "prior_hits": [],
            "comments": [],
            "target": "a claimed collaboration pilot must show an artifact crossing an organization boundary and producing an observed result",
            "source": "generic cross-organization preliminary-evidence check",
        }

    sentences = [
        sentence.strip()
        for sentence in re.split(r"(?<=[。．!?！？])|\n+", text)
        if sentence.strip()
    ]
    marked = [
        i
        for i, sentence in enumerate(sentences)
        if _contains_any(sentence.lower(), prior_terms)
    ]
    selected: set[int] = set()
    for index in marked:
        selected.update(range(max(0, index - 1), min(len(sentences), index + 4)))
    evidence_text = " ".join(sentences[i] for i in sorted(selected)) or text
    evidence_low = evidence_text.lower()

    axes = {
        "cross_organization_actor": [
            "大学間",
            "他機関",
            "別機関",
            "複数機関",
            "三機関",
            "四機関",
            "共同研究先",
            "開発元",
            "提供者",
        ],
        "transferred_artifact": [
            "コード",
            "実装",
            "モデル",
            "入力データ",
            "行列",
            "試験",
            "ベンチマーク",
            "変更",
            "pull request",
            "プルリクエスト",
        ],
        "bounded_task": [
            "同一問題",
            "再実行",
            "比較した",
            "統合した",
            "差し替",
            "取り込",
            "再現した",
            "検証した",
        ],
        "observed_outcome": [
            "収束",
            "反復",
            "残差",
            "誤差",
            "一致",
            "不一致",
            "採択",
            "棄却",
            "満たさなかった",
            "取り込まれ",
        ],
        "independent_action": [
            "別機関レビュー",
            "独立レビュー",
            "第三者",
            "非実装担当",
            "別機関が",
            "他機関が",
            "提供の",
            "提供した",
        ],
        "remaining_gap": [
            "未達",
            "未検証",
            "未実施",
            "主張しない",
            "に留まる",
            "残る",
            "出発点",
            "一方",
            "本研究で",
        ],
    }
    axis_results = {}
    for name, keywords in axes.items():
        hits = _contains_any(evidence_low, keywords)
        axis_results[name] = {
            "ok": bool(hits),
            "matches": hits[:8],
            "keywords": keywords,
        }

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "cross_organization_actor": "誰から誰へ資産が渡ったか、組織境界を明示する。",
        "transferred_artifact": "実際に渡したコード、モデル、データ、試験等を特定する。",
        "bounded_task": "同一問題での再実行、比較、変更等の限定した作業を示す。",
        "observed_outcome": "残差、誤差、採否、失敗等の観測結果を示す。",
        "independent_action": "資産作成者とは別の機関が行った実行・変更・レビューを示す。",
        "remaining_gap": "予備実証でまだ証明していない範囲を明記する。",
    }
    comments = [comments_by_axis[name] for name in missing]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    reviewed = bool(re.search(
        r"(?:別機関|他機関|第三者|非実装担当)[^。．]{0,50}"
        r"(?:レビューした|査読した|採択した|棄却した|変更した)",
        evidence_text,
    ))
    if axis_results["observed_outcome"]["ok"] and reviewed:
        evidence_level = "L3: independently reviewed or adopted scientific result"
    elif axis_results["observed_outcome"]["ok"]:
        evidence_level = "L2: cross-organization re-execution with an observed outcome"
    elif axis_results["transferred_artifact"]["ok"]:
        evidence_level = "L1: artifact transfer or build workflow only"
    else:
        evidence_level = "L0: publication, link, or internal-use claim only"
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "prior_hits": prior_hits[:12],
        "evidence_excerpt": evidence_text[:1000],
        "evidence_level": evidence_level,
        "independent_review_or_adoption": reviewed,
        "comments": comments,
        "target": (
            "name the cross-organization actor and artifact, show a bounded "
            "independent task and observed result, and state what remains unproven"
        ),
        "source": "generic cross-organization preliminary-evidence check",
    }


def grant_writing_named_software_abstraction_check(
    text: str,
    software_names: str = "",
) -> dict:
    """Check whether named software is used at the right proposal level.

    Titles, summaries, academic questions, aims, novelty, and impact should
    normally name the technical category (for example, OSS, a high-order
    finite-element OSS, or an OSS coupling platform), not make one software
    package the research concept.  Named software remains appropriate when it
    makes methods, preliminary evidence, collaboration, rights, or costs
    reproducible.

    Args:
        text: Proposal text or an existing .md/.tex/.txt path.
        software_names: Optional comma-separated names added to the built-in
            software list.
    """
    text = _read_text_if_path(text)
    default_names = [
        "NGSolve",
        "ONELAB",
        "openCFS",
        "FreeFEM++",
        "Gmsh",
        "GetDP",
        "preCICE",
        "OpenMDAO",
        "COMSOL",
        "JMAG",
        "ANSYS",
        "Radia",
        "MATLAB",
        "Simulink",
    ]
    extra_names = [
        name.strip() for name in software_names.split(",") if name.strip()
    ]
    names = list(dict.fromkeys(default_names + extra_names))

    core_markers = [
        "概要",
        "研究背景",
        "学術的問い",
        "中心の問い",
        "研究目的",
        "本研究の目的",
        "研究課題",
        "主題",
        "独創性",
        "新規性",
        "波及効果",
        "国際性",
        "研究動向",
        "位置付け",
        "意義",
        "タイトル",
    ]
    implementation_markers = [
        "研究方法",
        "達成指標",
        "実証",
        "検証",
        "年度計画",
        "研究計画",
        "準備状況",
        "遂行能力",
        "研究実績",
        "既往実績",
        "予備成果",
        "予備試験",
        "共同実績",
        "共同研究",
        "発表",
        "掲載",
        "リンク",
        "権利",
        "ライセンス",
        "予算",
        "経費",
        "費用",
        "リスク",
        "代替策",
        "ベンチマーク",
        "比較する",
        "用いる",
        "使用する",
        "実装する",
        "接続した",
        "統合した",
        "経験した",
        "共著",
    ]
    heading_pattern = re.compile(
        r"\\(?:sub)*section\*?\{([^{}]+)\}|^\s*#{1,6}\s+(.+?)\s*$"
    )

    def marker_hits(fragment: str, markers: list[str]) -> list[str]:
        return [marker for marker in markers if marker in fragment]

    current_section = ""
    risks = []
    allowed_mentions = []
    unclassified_mentions = []
    for line_number, line in enumerate(text.splitlines(), 1):
        heading_match = heading_pattern.search(line)
        if heading_match:
            current_section = heading_match.group(1) or heading_match.group(2) or ""

        for name in names:
            name_pattern = re.compile(
                rf"(?<![A-Za-z0-9]){re.escape(name)}(?![A-Za-z0-9])",
                re.IGNORECASE,
            )
            for match in name_pattern.finditer(line):
                sentence_start = line.rfind("。", 0, match.start()) + 1
                sentence_stop = line.find("。", match.end())
                if sentence_stop < 0:
                    sentence_stop = len(line)
                else:
                    sentence_stop += 1
                excerpt = line[sentence_start:sentence_stop].strip()
                local_core = marker_hits(excerpt, core_markers)
                section_core = marker_hits(current_section, core_markers)
                local_implementation = marker_hits(excerpt, implementation_markers)
                section_implementation = marker_hits(
                    current_section, implementation_markers
                )
                mention = {
                    "software": name,
                    "line": line_number,
                    "section": current_section,
                    "excerpt": excerpt[:320],
                }
                if local_implementation or section_implementation:
                    mention["allowed_by"] = (
                        local_implementation or section_implementation
                    )
                    allowed_mentions.append(mention)
                elif local_core or section_core:
                    mention.update({
                        "severity": "MEDIUM",
                        "framing_hits": local_core or section_core,
                        "comment": (
                            f"{name} is a named implementation, not the research "
                            "concept; state the OSS or technical category here."
                        ),
                    })
                    risks.append(mention)
                else:
                    unclassified_mentions.append(mention)

    applicable = bool(risks or allowed_mentions or unclassified_mentions)
    score = None if not applicable else max(0.0, 10.0 - 2.5 * len(risks))
    risky_names = sorted({risk["software"] for risk in risks})
    comments = [
        "問い・背景・目的・独創性・波及効果では、個別ソフト名 "
        + "、".join(risky_names)
        + " をOSSまたは技術カテゴリへ抽象化する。"
    ] if risky_names else []
    recommendations = [
        "中心概念は「OSS」「高次有限要素OSS」「既存OSS連成基盤」等の役割・技術カテゴリで記述する。",
        "固有ソフト名は、研究方法、予備成果、共同実績、権利・ライセンス、費用根拠で再現性を担保するときに限る。",
        "固有名詞を残す場合は、先に一般カテゴリと研究上の役割を示し、その一実装として名称を置く。",
    ]
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "allowed_mentions": allowed_mentions,
        "unclassified_mentions": unclassified_mentions,
        "comments": comments,
        "recommendations": recommendations if risks else [],
        "target": (
            "category-level framing in the title, summary, question, aims, "
            "novelty, and impact; named software only where it supplies "
            "reproducible implementation or feasibility evidence"
        ),
        "source": "generic named-software abstraction-level check",
    }


def grant_writing_reviewer_vocabulary_check(text: str) -> dict:
    """Check whether proposal vocabulary is accessible to the likely reviewer.

    A domain reviewer can be expected to know the field, but not every software
    acronym, foreign-university abbreviation, laboratory shorthand, or named
    benchmark. This check also keeps benchmarks in a verification role.
    """
    text = _read_text_if_path(text)
    risks: list[dict] = []
    term_results: dict[str, dict] = {}

    def context_for(match: re.Match, radius: int = 180) -> str:
        start = max(0, match.start() - radius)
        stop = min(len(text), match.end() + radius)
        return re.sub(r"\s+", " ", text[start:stop]).strip()

    def add_risk(
        risk_type: str,
        term: str,
        match: re.Match,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
    ) -> None:
        risks.append({
            "type": risk_type,
            "term": term,
            "line": text.count("\n", 0, match.start()) + 1,
            "severity": severity,
            "excerpt": context_for(match)[:360],
            "comment": comment,
            "recommendation": recommendation,
        })

    acronym_specs = {
        "OSS": {
            "pattern": r"(?<![A-Za-z0-9_])OSS(?![A-Za-z0-9_])",
            "explanations": [
                "オープンソースソフトウェア",
                "オープンソースソフトウエア",
            ],
            "recommendation": "初出を「オープンソースソフトウェア（OSS）」とする。",
        },
        "LLM": {
            "pattern": r"(?<![A-Za-z0-9_])LLM(?![A-Za-z0-9_])",
            "explanations": ["大規模言語モデル"],
            "recommendation": "略語を避けて「大規模言語モデル」と記す。",
        },
        "MCP": {
            "pattern": r"(?<![A-Za-z0-9_])MCP(?![A-Za-z0-9_])",
            "explanations": [
                "知識・実行インターフェース",
                "知識・実行インタフェース",
                "知識・実行窓口",
                "知識と実行のインターフェース",
                "知識と実行のインタフェース",
            ],
            "recommendation": (
                "初出で「AIが利用する知識・実行インターフェース」等の"
                "研究上の役割を日本語で説明する。英語名の展開だけにしない。"
            ),
        },
    }
    for term, spec in acronym_specs.items():
        matches = list(re.finditer(spec["pattern"], text))
        if not matches:
            term_results[term] = {"present": False, "count": 0, "explained": True}
            continue
        first = matches[0]
        local = context_for(first)
        explained = any(word in local for word in spec["explanations"])
        term_results[term] = {
            "present": True,
            "count": len(matches),
            "explained": explained,
            "first_line": text.count("\n", 0, first.start()) + 1,
        }
        if not explained:
            add_risk(
                "unexplained_acronym",
                term,
                first,
                f"{term}の初出に、審査者が読める日本語の意味・役割がない。",
                spec["recommendation"],
            )

    institution_aliases = {
        "TU Graz": "グラーツ工科大学",
        "TU Wien": "ウィーン工科大学",
    }
    for alias, japanese_name in institution_aliases.items():
        matches = list(re.finditer(re.escape(alias), text, flags=re.IGNORECASE))
        term_results[alias] = {
            "present": bool(matches),
            "count": len(matches),
            "preferred": japanese_name,
        }
        if matches:
            add_risk(
                "foreign_institution_alias",
                alias,
                matches[0],
                f"{alias}は分野外の審査者には大学名として即読できない。",
                f"本文では「{japanese_name}」と記し、必要な箇所だけ原語を併記する。",
            )

    shorthand_specs = {
        "MMM": {
            "pattern": r"(?<![A-Za-z0-9_])MMM(?![A-Za-z0-9_])",
            "preferred": "磁気モーメント法",
            "comment": "研究室内略称のMMMは電磁気分野でも共有を仮定できない。",
        },
        "H(curl)": {
            "pattern": (
                r"H\s*\(\s*(?:\\(?:mathrm|operatorname)\s*\{\s*)?"
                r"curl\s*\}?\s*\)"
            ),
            "preferred": "辺要素",
            "comment": "関数空間記号より、審査者が設計手法として読める用語を先に置く。",
        },
    }
    for term, spec in shorthand_specs.items():
        matches = list(re.finditer(spec["pattern"], text, flags=re.IGNORECASE))
        term_results[term] = {
            "present": bool(matches),
            "count": len(matches),
            "preferred": spec["preferred"],
        }
        if matches:
            add_risk(
                "domain_shorthand",
                term,
                matches[0],
                spec["comment"],
                (
                    f"概要・目的・意義では「{spec['preferred']}」を用い、"
                    "数理記号は必要なら方法節で補足する。"
                ),
            )

    benchmark_pattern = r"TEAM\s*(?:Problem\s*)?28"
    benchmark_matches = list(re.finditer(benchmark_pattern, text, flags=re.IGNORECASE))
    benchmark_result = {
        "present": bool(benchmark_matches),
        "count": len(benchmark_matches),
        "framed_as_verification": True,
        "limitation_stated": True,
        "engineering_value_distinguished": True,
        "negative_disclaimer_present": False,
        "used_as_significance": False,
    }
    if benchmark_matches:
        benchmark_context = " ".join(
            context_for(match, radius=260) for match in benchmark_matches
        )
        verification_terms = [
            "単純化",
            "公開基準問題",
            "初期検証",
            "整合性",
            "校正",
            "基礎検証",
        ]
        limitation_terms = [
            "工学的有用性の根拠には用いない",
            "工学的有用性の根拠にしない",
            "工学的有用性を示すものではない",
            "有用性の根拠には用いない",
            "有用性の根拠にしない",
            "有用性を示すものではない",
        ]
        real_design_terms = [
            "主実証",
            "実設計",
            "実機寸法",
            "制約付き機器設計",
            "設計判断",
            "設計量",
            "加速器電磁石",
        ]
        significance_terms = [
            "中心の問い",
            "学術的問い",
            "本研究の目的",
            "独創性",
            "新規性",
            "波及効果",
            "主成果",
        ]
        framed = any(term in benchmark_context for term in verification_terms)
        negative_disclaimer = any(term in text for term in limitation_terms)
        positive_design_bridge = framed and any(
            term in benchmark_context for term in real_design_terms
        )
        limited = positive_design_bridge
        significance = any(
            term in context_for(match, radius=120)
            for match in benchmark_matches
            for term in significance_terms
        )
        benchmark_result.update({
            "framed_as_verification": framed,
            "limitation_stated": limited,
            "engineering_value_distinguished": limited,
            "negative_disclaimer_present": negative_disclaimer,
            "used_as_significance": significance,
        })
        first = benchmark_matches[0]
        if not framed:
            add_risk(
                "benchmark_unframed",
                "TEAM Problem 28",
                first,
                "固有ベンチマークの研究計画上の役割が明示されていない。",
                "単純化した公開基準問題を初期検証または校正に用いる、と位置付ける。",
            )
        if not limited:
            add_risk(
                "benchmark_without_limit",
                "TEAM Problem 28",
                first,
                "基準問題への一致と工学的有用性が区別されていない。",
                (
                    "基準問題の役割を初期検証とし、工学的価値を示す制約付き"
                    "実設計課題と設計判断を続けて記す。予備成果を否定する"
                    "免責文は置かない。"
                ),
                severity="HIGH",
            )
        if negative_disclaimer:
            disclaimer_match = next(
                match
                for term in limitation_terms
                if (match := re.search(re.escape(term), text)) is not None
            )
            add_risk(
                "benchmark_self_negating_disclaimer",
                "TEAM Problem 28",
                disclaimer_match,
                "予備検証の成果を、直後の免責文が自ら打ち消している。",
                (
                    "基準問題で何を確認できたかを肯定形で述べ、続けて主実証の"
                    "制約付き実設計と設計判断を示す。"
                ),
            )
        if significance:
            add_risk(
                "benchmark_as_significance",
                "TEAM Problem 28",
                first,
                "固有ベンチマークが問い・目的・独創性・波及効果の代わりになっている。",
                "固有名は方法・予備結果へ下げ、工学的意義を実設計量と意思決定で示す。",
                severity="HIGH",
            )
        if len(benchmark_matches) > 2:
            add_risk(
                "benchmark_overexposure",
                "TEAM Problem 28",
                benchmark_matches[2],
                "固有ベンチマーク名の反復が、研究の主役であるように見せている。",
                "固有名は方法節で一度定義し、以後は「公開基準問題」等と記す。",
            )

    term_results["TEAM Problem 28"] = benchmark_result
    applicable = any(result.get("present") for result in term_results.values())
    deductions = sum(2.0 if risk["severity"] == "HIGH" else 1.0 for risk in risks)
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(
        dict.fromkeys(risk["recommendation"] for risk in risks)
    )
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "term_results": term_results,
        "benchmark_result": benchmark_result,
        "target": (
            "assume field expertise, not familiarity with every software/AI acronym, "
            "foreign-institution alias, laboratory shorthand, or named benchmark; "
            "use benchmarks for verification and real design tasks for engineering value"
        ),
        "source": "generic reviewer-vocabulary and benchmark-role check",
    }


def grant_writing_persuasion_quality_check(text: str) -> dict:
    """Check reviewer-facing hierarchy, equations, and defensive prose.

    The check targets a common failure mode in technically careful proposals:
    evidence is immediately negated, equations arrive before their symbols and
    decision role, or exceptions crowd out the main claim. It is a heuristic
    quality gate, not a mathematical parser.
    """
    text = _read_text_if_path(text)
    risks: list[dict] = []

    def add_risk(
        risk_type: str,
        start: int,
        excerpt: str,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
        **details,
    ) -> None:
        item = {
            "type": risk_type,
            "line": text.count("\n", 0, start) + 1,
            "severity": severity,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:360],
            "comment": comment,
            "recommendation": recommendation,
        }
        item.update(details)
        risks.append(item)

    self_negating_patterns = [
        re.compile(
            r"(?:この|その)?(?:一致|結果|成果|実績|予備実証)"
            r"[^。！？\n]{0,60}(?:根拠|有用性|意義|価値)"
            r"[^。！？\n]{0,35}(?:しない|用いない|ではない|示さない)"
        ),
        re.compile(
            r"(?:工学的|学術的)?(?:有用性|意義|価値)[^。！？\n]{0,40}"
            r"(?:示すものではない|根拠(?:に|と)(?:は)?"
            r"(?:しない|ならない|用いない))"
        ),
        re.compile(
            r"(?:結果|成果|実績)[^。！？\n]{0,50}"
            r"(?:だけを示す|にすぎない)"
        ),
    ]
    seen_self_negating: set[tuple[int, int]] = set()
    for pattern in self_negating_patterns:
        for match in pattern.finditer(text):
            span = match.span()
            if span in seen_self_negating:
                continue
            seen_self_negating.add(span)
            add_risk(
                "self_negating_evidence",
                match.start(),
                match.group(0),
                "提示した成果を同じ段落の免責文が打ち消している。",
                (
                    "成果が確認した範囲を肯定形で述べ、次の文で主実証が検証する"
                    "設計量・判断・適用範囲へ接続する。"
                ),
                severity="HIGH",
            )

    equation_pattern = re.compile(
        r"(?<!\\)\\\[(?P<bracket>.*?)(?<!\\)\\\]"
        r"|\\begin\{(?P<env>equation\*?|align\*?|gather\*?)\}"
        r"(?P<environment>.*?)\\end\{(?P=env)\}"
        r"|\$\$(?P<dollar>.*?)\$\$",
        flags=re.DOTALL,
    )
    symbol_pattern = re.compile(
        r"(?:\\[A-Za-z]+|[A-Za-z])(?:_\{[^{}\n]+\}|_[A-Za-z0-9]+)"
    )
    on_ramp_terms = [
        "定義",
        "表す",
        "求める",
        "算出",
        "評価",
        "判断",
        "区間",
        "指標",
        "関係",
        "変化",
        "比較",
        "幅",
    ]
    interpretation_terms = [
        "と表す",
        "意味",
        "示す",
        "判断",
        "順位",
        "区間",
        "係数",
        "用いる",
        "分離",
        "重な",
    ]

    def compact_latex(value: str) -> str:
        value = re.sub(r"\\(?:left|right|,|;|!|quad|qquad)", "", value)
        return re.sub(r"[\s{}$]", "", value)

    equation_count = 0
    for match in equation_pattern.finditer(text):
        equation_count += 1
        body = next(
            group
            for group in (
                match.group("bracket"),
                match.group("environment"),
                match.group("dollar"),
            )
            if group is not None
        )
        before = text[max(0, match.start() - 500):match.start()]
        after = text[match.end():min(len(text), match.end() + 500)]
        before_prose = _prose_for_lint(before)
        after_prose = _prose_for_lint(after)
        equation_excerpt = match.group(0)

        if not any(term in before_prose for term in on_ramp_terms):
            add_risk(
                "equation_without_on_ramp",
                match.start(),
                equation_excerpt,
                "数式の前に、何を判断するための式かが示されていない。",
                (
                    "式の前に対象量、比較目的、式が答える判断を一文で述べてから"
                    "数式を置く。"
                ),
                severity="HIGH",
            )

        symbols = sorted(set(symbol_pattern.findall(body)))
        surrounding = compact_latex(before + after)
        missing_symbols = [
            symbol for symbol in symbols if compact_latex(symbol) not in surrounding
        ]
        if missing_symbols:
            add_risk(
                "equation_symbols_not_introduced",
                match.start(),
                equation_excerpt,
                "数式中の記号が周辺の文章で導入されていない。",
                (
                    "各記号を物理的意味と計算操作で先に定義する。添字だけで"
                    "解析条件を推測させない。"
                ),
                severity="HIGH",
                missing_symbols=missing_symbols,
            )

        if not any(term in after_prose for term in interpretation_terms):
            add_risk(
                "equation_without_interpretation",
                match.start(),
                equation_excerpt,
                "数式の後に、値をどう研究判断へ使うかが説明されていない。",
                (
                    "式の直後に、大小・分離・閾値がどの設計判断を変えるかを"
                    "平文で述べる。"
                ),
            )

    inline_math_pattern = re.compile(
        r"(?<!\$)\$(?!\$)(?P<body>[^$\n]{1,120})(?<!\$)\$(?!\$)"
    )
    inline_condition_pattern = re.compile(
        r"^\s*(?P<symbol>(?:\\[A-Za-z]+|[A-Za-z])"
        r"(?:_\{[^{}\n]+\}|_[A-Za-z0-9]+)?)\s*"
        r"(?:=|\\leq?|\\geq?|<|>)\s*"
        r"(?:[-+]?(?:\d+(?:\.\d*)?|\.\d+)|\\infty)\s*$"
    )
    unexplained_inline_condition_count = 0
    for match in inline_math_pattern.finditer(text):
        condition = inline_condition_pattern.match(match.group("body"))
        if condition is None:
            continue
        symbol = condition.group("symbol")
        base_match = re.match(r"(?:\\[A-Za-z]+|[A-Za-z])", symbol)
        if base_match is None:
            continue
        base_symbol = base_match.group(0)
        before = text[max(0, match.start() - 400):match.start()]
        prior_symbol_pattern = re.compile(
            r"(?<!\$)\$(?!\$)[^$\n]{0,80}"
            + re.escape(base_symbol)
            + r"[^$\n]{0,80}(?<!\$)\$(?!\$)"
        )
        if prior_symbol_pattern.search(before):
            continue
        unexplained_inline_condition_count += 1
        add_risk(
            "inline_condition_before_physical_meaning",
            match.start(),
            match.group(0),
            "記号だけの条件式が、その記号の工学的意味より先に現れている。",
            (
                "現象と設計上の効果を平文で述べ、量の名称と記号を独立に"
                "定義してから条件式を置く。概要・図・成果名には工学的効果を使う。"
            ),
            severity="HIGH",
            symbol=symbol,
        )

    paragraph_negative_terms = [
        "必達成果には含めない",
        "目的としない",
        "前提にしない",
        "依存しない",
        "対象外",
        "今後の課題",
        "失敗すれば",
        "含めない",
        "ではない",
        "しない",
        "不能",
        "未達",
    ]
    defensive_exemptions = [
        "リスク",
        "研究遂行上の課題",
        "代替策",
        "失敗時",
        "不成立時",
        "倫理",
        "法令",
        "安全",
        "個人情報",
        "知的財産",
        "権利",
    ]
    negative_pattern = re.compile(
        "|".join(re.escape(term) for term in paragraph_negative_terms)
    )
    cursor = 0
    defensive_paragraph_count = 0
    for raw_paragraph in re.split(r"\n\s*\n+", text):
        start = text.find(raw_paragraph, cursor)
        cursor = max(cursor, start + len(raw_paragraph))
        paragraph = _prose_for_lint(raw_paragraph)
        if len(paragraph) < 35 or any(
            term in paragraph for term in defensive_exemptions
        ):
            continue
        negative_hits = negative_pattern.findall(paragraph)
        if len(negative_hits) >= 3:
            defensive_paragraph_count += 1
            add_risk(
                "defensive_paragraph",
                max(0, start),
                paragraph,
                "否定・除外・未達の説明が一段落に集中し、主張が見えにくい。",
                (
                    "本文では成立させる主張と判定法を先に書く。例外と代替策は"
                    "研究遂行上のリスク段落へまとめる。"
                ),
                negative_hits=negative_hits,
            )

    optional_pattern = re.compile(
        r"条件付き追加検証|必達成果には含めない|余力があれば|条件が整えば|"
        r"可能であれば[^。！？\n]{0,40}(?:追加|検証|実施)"
    )
    optional_branch_count = 0
    for match in optional_pattern.finditer(text):
        optional_branch_count += 1
        add_risk(
            "optional_branch_in_core_plan",
            match.start(),
            match.group(0),
            "任意の追加課題が中心計画に入り、必達成果の輪郭を弱めている。",
            (
                "中核計画は必達の問い・実証・成果に絞り、条件付き拡張は"
                "リスク対応または将来展開へ移す。"
            ),
        )

    memo_specs = {
        "lettered_experiment": re.compile(r"実証[AB](?![A-Za-z0-9])"),
        "letter_pair": re.compile(
            r"(?<![A-Za-z0-9])A(?:・|/)B"
            r"(?=(?:で|の|を|は|へ|と|仕様|論文|成果|実証))"
        ),
        "coded_stage": re.compile(r"(?<![A-Za-z0-9])L[1-4](?![A-Za-z0-9])"),
        "parenthetical_milestone": re.compile(
            r"[（(]\s*(?:年度末|到達点|ゴール)\s*[:：]"
        ),
        "draft_marker": re.compile(
            r"(?i:\b(?:TODO|TBD|FIXME)\b)|要確認|暫定案|仮置き"
        ),
    }
    memo_shorthand_count = 0
    for memo_type, pattern in memo_specs.items():
        matches = list(pattern.finditer(text))
        if not matches:
            continue
        memo_shorthand_count += len(matches)
        first = matches[0]
        add_risk(
            "internal_memo_shorthand",
            first.start(),
            first.group(0),
            "内部管理用の略号・進捗メモが申請本文に残っている。",
            (
                "実証は科学的内容で呼び、年度到達点は括弧メモでなく"
                "完結した文で記す。下書き状態の印は本文から除く。"
            ),
            memo_type=memo_type,
            occurrence_count=len(matches),
        )

    prose = _prose_for_lint(text)
    acronym_pattern = re.compile(
        r"(?<![A-Za-z0-9])(?:[A-Z]{2,}[A-Za-z0-9]*"
        r"(?:[-/][A-Za-z0-9]+)*|H\([a-z]+\))(?![A-Za-z0-9])"
    )
    acronym_pile_count = 0
    for sentence in re.split(r"(?<=[。．!?！？])", prose):
        acronyms = sorted(set(acronym_pattern.findall(sentence)))
        if len(acronyms) < 6:
            continue
        acronym_pile_count += 1
        start = text.find(sentence[:24])
        add_risk(
            "acronym_pile",
            max(0, start),
            sentence,
            "一文に固有略語が集中し、研究上の役割より名称の列挙が先に見える。",
            (
                "まず手法カテゴリ、比較軸、変わる設計判断を述べ、固有名は"
                "方法・実行可能性の箇所に分けて置く。"
            ),
            acronyms=acronyms,
        )

    deductions = sum(2.0 if risk["severity"] == "HIGH" else 1.0 for risk in risks)
    applicable = bool(text.strip())
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(
        dict.fromkeys(risk["recommendation"] for risk in risks)
    )
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "counts": {
            "display_equations": equation_count,
            "self_negating_evidence": len(seen_self_negating),
            "defensive_paragraphs": defensive_paragraph_count,
            "optional_branches": optional_branch_count,
            "memo_shorthand": memo_shorthand_count,
            "acronym_piles": acronym_pile_count,
            "unexplained_inline_conditions": unexplained_inline_condition_count,
        },
        "target": (
            "lead with a positive claim and reviewer decision; introduce every "
            "equation and inline condition by its physical meaning and symbols, "
            "interpret it immediately, collect exceptions in one risk section, "
            "and keep optional branches out of the core plan"
        ),
        "source": "generic reviewer-facing persuasion and equation-on-ramp check",
    }


def grant_writing_kaken_review_format_check(text: str) -> dict:
    """Check KAKENHI reviewer-format realities on a proposal draft.

    Encodes the in-house KAKENHI call briefing (R9/FY2027 call): reviewers
    judge on three criteria and read up to ~100 proposals in about a month;
    figures may be printed in monochrome for some categories; publication
    records are read through researchmap; the human-rights/legal-compliance
    box draws the most reviewer remarks; and the funding-overlap box has a
    fixed format. Fragment-level triggers gate each sub-check, so short
    excerpts stay clean; full-draft heuristics apply above ~1500 chars.
    """
    raw = _read_text_if_path(text)
    prose = _prose_for_lint(raw)
    risks: list[dict] = []

    def add_risk(
        risk_type: str,
        start: int,
        excerpt: str,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
        **details,
    ) -> None:
        item = {
            "type": risk_type,
            "line": raw.count("\n", 0, max(0, start)) + 1,
            "severity": severity,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:360],
            "comment": comment,
            "recommendation": recommendation,
        }
        item.update(details)
        risks.append(item)

    color_figure_pattern = re.compile(
        r"(?:赤|青|緑|黄|橙|紫|桃)(?:色|い)?(?:の)?"
        r"(?:実線|破線|点線|一点鎖線|線|棒|丸印|丸|矢印|領域|塗り|"
        r"プロット|曲線|マーカー|文字|枠)"
        r"|色分け|色で(?:区別|示|表)|カラーで(?:区別|示|表)"
    )
    mono_terms = [
        "白黒",
        "モノクロ",
        "グレースケール",
        "線種",
        "濃淡",
        "ハッチング",
        "マーカー形状",
    ]
    color_matches = list(color_figure_pattern.finditer(raw))
    if color_matches and not any(term in raw for term in mono_terms):
        first = color_matches[0]
        add_risk(
            "color_dependent_figure",
            first.start(),
            first.group(0),
            "色の違いだけで図の系列・領域を区別している。審査時に白黒印刷される種目がある。",
            "線種・マーカー・直接ラベル・濃淡で区別し、白黒でも判別できる図にする。",
            severity="HIGH",
            occurrence_count=len(color_matches),
        )

    subject_terms = [
        "アンケート",
        "質問紙",
        "被験者",
        "調査対象者",
        "インタビュー",
        "動物実験",
        "実験動物",
        "ヒト由来",
        "臨床",
        "個人情報",
    ]
    # NOTE: 「遵守」 alone is NOT a safeguard: the box heading itself
    # (「人権の保護及び法令等の遵守への対応」) contains it, so it would
    # suppress the check on every draft that quotes the heading.
    safeguard_terms = [
        "倫理審査",
        "倫理委員会",
        "動物実験委員会",
        "同意",
        "インフォームド",
        "承認",
        "匿名化",
        "個人情報保護",
        "適切に管理",
    ]
    subject_hits = [term for term in subject_terms if term in prose]
    if subject_hits and not any(term in prose for term in safeguard_terms):
        add_risk(
            "human_subjects_without_safeguard",
            max(0, raw.find(subject_hits[0])),
            subject_hits[0],
            "人・動物・個人情報を扱う記述があるのに、講じる対策・措置の記載が見当たらない。",
            "「人権の保護及び法令等の遵守への対応」欄に、倫理審査、同意取得、"
            "匿名化等の具体的な対策を記載する。例年、審査委員からの指摘が最も多い欄である。",
            severity="HIGH",
            subject_hits=subject_hits,
        )

    na_pattern = re.compile(r"該当\s*(?:は|事項は)?\s*(?:なし|ない|無し|ありません)")
    rationale_terms = [
        "ため",
        "ので",
        "対象としない",
        "対象とせず",
        "行わない",
        "用いない",
        "使用しない",
        "含まない",
        "扱わない",
        "のみであり",
        "のみで",
        "理由",
    ]
    bare_na_matches = []
    for match in na_pattern.finditer(raw):
        sentence_start = max(raw.rfind("。", 0, match.start()) + 1, 0)
        sentence_stop = raw.find("。", match.end())
        if sentence_stop < 0:
            sentence_stop = len(raw)
        sentence = raw[sentence_start:sentence_stop]
        if not any(term in sentence for term in rationale_terms):
            bare_na_matches.append((match, sentence))
    if bare_na_matches:
        first, sentence = bare_na_matches[0]
        add_risk(
            "not_applicable_without_rationale",
            first.start(),
            sentence,
            "「該当なし」とだけ書かれ、そう判断した根拠がない。",
            "人を対象としない数値解析のみである等、該当なしと判断した根拠を"
            "一文添える。この欄は例年審査委員からの指摘が非常に多い。",
            occurrence_count=len(bare_na_matches),
        )

    pub_triggers = ["研究遂行能力", "研究業績", "主要論文", "代表論文", "発表論文", "業績"]
    pub_hit = next((term for term in pub_triggers if term in prose), None)
    identifier_pattern = re.compile(
        r"(?:19|20)\d{2}|vol\.?\s*\d|no\.?\s*\d|pp\.?\s*\d|doi|DOI|"
        r"\d+\s*巻|\d+\s*号|第\d+巻",
        flags=re.IGNORECASE,
    )
    if pub_hit is not None and not identifier_pattern.search(raw):
        add_risk(
            "publication_not_identifiable",
            max(0, raw.find(pub_hit)),
            pub_hit,
            "業績への言及があるが、業績を特定できる情報(誌名・年・巻号等)がない。",
            "審査はresearchmapを研究者番号で参照する。調書に業績を書く場合は、"
            "特定するための十分な情報(著者、誌名、年等)を添える。",
        )

    overlap_triggers = [
        "応募中の研究費",
        "受入予定の研究費",
        "応募・受入",
        "応募中及び受入",
        "応募中および受入",
    ]
    overlap_hit = next((term for term in overlap_triggers if term in raw), None)
    if overlap_hit is not None:
        role_pattern = re.compile(
            r"(?:大学|研究所|機構|高等専門学校|高専)[^。\n]{0,15}"
            r"(?:教授|准教授|講師|助教|研究員)"
        )
        missing_parts = []
        if "相違" not in raw:
            missing_parts.append("本応募課題との相違点")
        if not re.search(r"応募(?:する)?理由", raw):
            missing_parts.append("応募する理由")
        if not role_pattern.search(raw):
            missing_parts.append("所属組織・役職(例: ○○大学教授)")
        if missing_parts:
            add_risk(
                "funding_overlap_format",
                max(0, raw.find(overlap_hit)),
                overlap_hit,
                "応募・受入状況欄に必要な記載要素が欠けている: "
                + "、".join(missing_parts),
                "国外資金・民間財団助成・受託/共同研究費も全て記載し、2件目以降は"
                "研究内容の相違点と応募する理由、所属組織・役職を添える。"
                "代表課題は分担者を含む金額、分担課題は自身の経費のみを書く。",
                missing_parts=missing_parts,
            )

    full_draft = len(prose) >= 1500
    criteria_axis_results: dict[str, dict] = {}
    if full_draft:
        low = prose.lower()
        for axis, keywords in _KAKEN_REVIEW_CRITERIA_AXES.items():
            hits = _contains_any(low, keywords)
            criteria_axis_results[axis] = {
                "ok": bool(hits),
                "matches": hits[:8],
                "keywords": keywords,
            }
        missing_criteria = [
            axis
            for axis, result in criteria_axis_results.items()
            if not result["ok"]
        ]
        if missing_criteria:
            add_risk(
                "review_criteria_axis_missing",
                0,
                "、".join(missing_criteria),
                "3つの審査基準(学術的重要性・方法の妥当性・遂行能力/環境)のうち、"
                "読み取れない軸がある: " + "、".join(missing_criteria),
                "各セクションがどの審査基準で読まれるかを意識し、3基準すべてに"
                "対応する記述を置く。",
                missing_axes=missing_criteria,
            )
        emphasis_pattern = re.compile(
            r"下線|太字|ゴシック|アンダーライン|\\underline|\\textbf|"
            r"[図表]\s*\d|Fig\.?\s*\d"
        )
        if not emphasis_pattern.search(raw):
            add_risk(
                "no_emphasis_or_figures",
                0,
                prose[:80],
                "長い本文に、強調(下線・太字・ゴシック)や図表参照が見当たらない。",
                "審査委員は約1ヶ月で最大100件程度を読む。要点への下線・太字と"
                "図表で、一読で主張が追える構成にする。",
            )

    deductions = sum(2.0 if risk["severity"] == "HIGH" else 1.0 for risk in risks)
    applicable = bool(raw.strip())
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(
        dict.fromkeys(risk["recommendation"] for risk in risks)
    )
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "full_draft_heuristics_applied": full_draft,
        "criteria_axis_results": criteria_axis_results,
        "briefing_notes": list(_KAKEN_BRIEFING_NOTES),
        "target": (
            "a proposal a reviewer can judge on the three criteria at "
            "~100-proposals-per-month reading speed: monochrome-safe figures, "
            "identifiable publications, an explicit human-rights/legal box, "
            "and a complete funding-overlap box"
        ),
        "source": "KAKENHI in-house call briefing (R9/FY2027) review-format check",
    }


def grant_writing_literature_gap_evidence_check(text: str) -> dict:
    """Check whether literature-survey evidence supports the claimed gap.

    A missing keyword or method name in a bounded corpus can help select
    examples, but it does not by itself establish field-wide non-adoption, a
    causal implementation barrier, or an academic gap.  This check looks for
    those inference jumps in a local sentence window.  It is not a systematic
    review or citation-quality validator.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    sentences = [
        sentence.strip()
        for sentence in re.split(r"(?<=[。．!?！？])|\n+", text)
        if sentence.strip()
    ]

    keyword_groups = {
        "search_method": [
            "全文検索",
            "キーワード検索",
            "検索語",
            "同一語彙",
            "検出頁",
            "目視確認",
            "出現頻度",
            "ヒット件数",
            "full-text search",
            "keyword search",
        ],
        "corpus_scope": [
            "調査対象",
            "対象文献",
            "対象資料",
            "技術報告",
            "選定した文献",
            "限定コーパス",
            "bounded corpus",
            "selected literature",
        ],
        "non_detection": [
            "確認できなかった",
            "確認されなかった",
            "見当たらなかった",
            "検出されなかった",
            "現れなかった",
            "ヒットしなかった",
            "記載がない",
            "言及がない",
            "採用例がない",
            "not found",
            "no occurrence",
            "was not identified",
        ],
        "field_scope": [
            "国内では",
            "日本では",
            "我が国では",
            "国内全",
            "当該分野",
            "分野全体",
            "コミュニティ全体",
            "海外では",
            "欧州では",
            "in japan",
            "field-wide",
            "community-wide",
        ],
        "adoption_claim": [
            "普及していない",
            "使われていない",
            "利用されていない",
            "採用されていない",
            "浸透していない",
            "未利用",
            "デファクト",
            "候補から外れ",
            "取り込まれない",
            "遅れている",
            "ガラパゴス",
            "not widely used",
            "not adopted",
        ],
        "causal_barrier": [
            "知識不足",
            "実装障壁",
            "言語障壁",
            "情報不足",
            "不足により",
            "不足から",
            "起因",
            "阻害",
            "妨げ",
            "implementation barrier",
            "language barrier",
        ],
        "academic_gap": [
            "研究障壁",
            "学術的空白",
            "研究上の空白",
            "未解決問題",
            "未解決の課題",
            "academic gap",
            "research barrier",
        ],
        "inference_link": [
            "示す",
            "裏付ける",
            "根拠",
            "意味する",
            "したがって",
            "ゆえに",
            "demonstrates",
            "therefore",
            "evidence of",
        ],
        "limited_role": [
            "背景補助",
            "補助証拠",
            "実証対象の選定",
            "実証対象選定",
            "事例選定",
            "予備調査",
            "本調査範囲",
            "網羅調査ではなく",
            "網羅的ではなく",
            "主根拠としない",
            "scope examples",
            "case selection",
        ],
    }
    count_pattern = re.compile(
        r"(?:全|計)?\s*\d+\s*(?:冊|報|件|編|誌|論文|文献|資料|報告)"
    )

    def keyword_hits(fragment: str, group: str) -> list[str]:
        low = fragment.lower()
        return [word for word in keyword_groups[group] if word.lower() in low]

    absence_indices = [
        index
        for index, sentence in enumerate(sentences)
        if keyword_hits(sentence, "non_detection")
    ]
    candidate_windows = []
    for index in absence_indices:
        start = max(0, index - 4)
        stop = min(len(sentences), index + 5)
        window = "".join(sentences[start:stop])
        search_hits = keyword_hits(window, "search_method")
        corpus_hits = keyword_hits(window, "corpus_scope")
        count_hits = count_pattern.findall(window)
        if not search_hits and not corpus_hits and not count_hits:
            continue
        candidate_windows.append({
            "sentence": index + 1,
            "excerpt": window[:360],
            "search_hits": search_hits,
            "corpus_hits": corpus_hits + count_hits,
            "non_detection_hits": keyword_hits(window, "non_detection"),
            "field_scope_hits": keyword_hits(window, "field_scope"),
            "adoption_hits": keyword_hits(window, "adoption_claim"),
            "causal_hits": keyword_hits(window, "causal_barrier"),
            "academic_gap_hits": keyword_hits(window, "academic_gap"),
            "inference_hits": keyword_hits(window, "inference_link"),
            "limited_role_hits": keyword_hits(window, "limited_role"),
        })

    risk_specs = {
        "field_generalization": (
            "限定した文献群での非検出から、国内・分野全体の未普及や採用状況へ一般化している。",
            lambda item: bool(item["field_scope_hits"] and item["adoption_hits"]),
        ),
        "unsupported_causal_inference": (
            "語の非検出から、知識・言語・実装上の障壁という原因へ飛躍している。",
            lambda item: bool(item["causal_hits"] and item["inference_hits"]),
        ),
        "absence_as_academic_gap": (
            "限定コーパスでの不在を、未解決の研究障壁・学術的空白の主根拠にしている。",
            lambda item: bool(
                item["academic_gap_hits"] and item["inference_hits"]
            ),
        ),
    }
    risks = []
    for risk_type, (comment, predicate) in risk_specs.items():
        evidence = [item for item in candidate_windows if predicate(item)]
        if evidence:
            risks.append({
                "type": risk_type,
                "severity": "HIGH",
                "comment": comment,
                "evidence": evidence[:3],
            })

    applicable = bool(candidate_windows)
    score = None if not applicable else max(0.0, 10.0 - 3.0 * len(risks))
    rewrite_strategy = [
        "限定コーパスで語を確認できないことを、分野全体の未普及・優劣・研究障壁の主根拠にしない。",
        "文献調査は背景補助、候補選定、予備調査に位置づけ、対象範囲と限界を明示する。",
        "学術的空白は、既往研究が解いていない条件、理論、比較可能性、または検証可能な仮説として独立に述べる。",
        "普及実態を主張する場合は、検索式、選定基準、対象範囲、代替語、再現可能な集計を別途示す。",
    ]
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "candidate_windows": candidate_windows[:5],
        "comments": [risk["comment"] for risk in risks],
        "rewrite_strategy": rewrite_strategy if risks else [],
        "target": (
            "use bounded literature searches for background or case selection; "
            "establish the academic gap through an unresolved scientific "
            "condition, theory, comparison, or testable hypothesis"
        ),
        "source": "generic literature-gap evidence-scope check",
    }


def grant_writing_collaborative_integration_risk_check(text: str) -> dict:
    """Check recurring risks in collaborative software-integration proposals.

    The check is domain-neutral and applies only when a draft proposes coupling,
    integration, interoperability, or cross-organization reuse. It distinguishes
    the academic question from the implementation mechanism and checks costs on
    both sides of the interface, ecosystem positioning, scope control, negative
    results, team readiness, evaluation ethics, and asset provenance.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    applicable_hits = _contains_any(
        low,
        [
            "結合",
            "統合",
            "連携",
            "連成",
            "相互運用",
            "interoperability",
            "integration",
            "coupling",
        ],
    )

    axis_groups = {
        "academic_question_vs_mechanism": [
            ["学術的問い", "中心の問い", "研究上の問い", "仮説", "research question"],
            ["実装手段", "検証手段", "研究手段", "method", "mechanism"],
        ],
        "provider_and_reuse_cost": [
            ["所有者側", "提供者側", "開発元", "provider", "maintainer"],
            ["初期整備", "初期費用", "保守", "maintenance", "setup cost"],
            ["総負担", "再利用", "損益分岐", "reuse", "total cost"],
        ],
        "existing_ecosystem_boundary": [
            ["既存基盤", "既存規格", "既存oss", "existing framework", "existing standard"],
            ["置き換えず", "置換しない", "再利用", "補完", "boundary", "reuse"],
        ],
        "core_vs_optional_scope": [
            ["中核", "必達", "成立条件", "core", "required"],
            ["独立課題", "別課題", "発展候補", "条件付き", "optional", "exploratory"],
        ],
        "negative_result_value": [
            ["結合不能", "適用境界", "不能理由", "反例", "不成立", "negative result"],
            ["成果", "判定", "同定", "明らか", "result", "outcome"],
        ],
        "team_readiness": [
            ["共著", "既往成果", "既発表", "共同遂行", "prior work", "track record"],
            ["担当", "役割", "責任", "role", "responsibility"],
            ["着手", "準備", "基礎", "既に", "ready", "readiness"],
        ],
        "evaluation_unit_and_ethics": [
            ["評価単位", "分析単位", "課題単位", "unit of analysis"],
            ["個人を評価", "個人の能力", "個人情報", "individual productivity"],
            ["倫理", "同意", "該当性", "ethics", "consent"],
        ],
        "asset_provenance_and_fallback": [
            ["権利", "保守主体", "所有者", "provenance", "maintainer"],
            ["参照実装", "公開ベンチマーク", "代替", "fallback", "reference implementation"],
        ],
    }

    axis_results = {}
    for axis, groups in axis_groups.items():
        group_hits = [_contains_any(low, group) for group in groups]
        axis_results[axis] = {
            "ok": all(bool(hits) for hits in group_hits),
            "group_matches": [hits[:8] for hits in group_hits],
            "groups": groups,
        }

    people_process_hits = _contains_any(
        low,
        [
            "学生",
            "若手",
            "工程時間",
            "手作業時間",
            "生産性",
            "被験者",
            "アンケート",
            "参加者",
        ],
    )
    if not people_process_hits:
        axis_results["evaluation_unit_and_ethics"].update(
            {"ok": True, "not_applicable": True}
        )

    if not applicable_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": axis_results,
            "comments": [],
            "target": (
                "for collaborative integration proposals, separate the research "
                "question from tooling and test lifecycle cost, scope, negative "
                "results, ethics, and provenance"
            ),
            "source": "generic collaborative-integration risk check",
        }

    missing = [axis for axis, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "academic_question_vs_mechanism": (
            "製品・プロトコル名を問いにせず、学術的問いと実装手段を分ける。"
        ),
        "provider_and_reuse_cost": (
            "利用側だけでなく、提供者側の初期整備・保守を含む総負担を再利用回数に対して評価する。"
        ),
        "existing_ecosystem_boundary": (
            "既存規格・連携基盤が担う範囲を認め、置換せず再利用する範囲と研究上の空白を示す。"
        ),
        "core_vs_optional_scope": (
            "中核実証と条件付きの発展候補を分け、発展候補を必達成果の成立条件にしない。"
        ),
        "negative_result_value": (
            "順位不変、結合不能、反例等も、条件・原因・適用境界を同定できれば成果となる成功条件を置く。"
        ),
        "team_readiness": (
            "各担当者について、既往成果、利用可能資産、役割、着手可能性を一続きで示す。"
        ),
        "evaluation_unit_and_ethics": (
            "工程を比較する場合は個人でなく課題・成果物を評価単位とし、倫理該当性と同意手続を確認する。"
        ),
        "asset_provenance_and_fallback": (
            "固有資産の権利・保守主体を確認し、利用不能時の公開ベンチマークまたは参照実装を用意する。"
        ),
    }
    comments = [comments_by_axis[axis] for axis in missing]
    score = round(10.0 * (len(axis_groups) - len(missing)) / len(axis_groups), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "people_process_hits": people_process_hits,
        "comments": comments,
        "target": (
            "a collaborative integration proposal with an academic question distinct "
            "from tooling, full lifecycle cost, ecosystem boundaries, controlled scope, "
            "valuable negative results, team readiness, ethical evaluation, and fallback assets"
        ),
        "source": "generic collaborative-integration risk check",
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
        "pricing_provenance": [
            "公式料金",
            "料金表",
            "見積書",
            "参照日",
            "改定日",
            "料金年度",
            "最低購入",
            "有効期限",
            "為替",
            "端数処理",
        ],
    }
    results = {}
    missing = []
    for axis, keywords in axes.items():
        if axis == "pricing_provenance":
            required_groups = {
                "source": ["公式料金", "料金表", "見積書", "https://", "url"],
                "vintage": ["参照日", "改定日", "料金年度", "年度料金", "2025年度", "2026年度"],
                "accounting": [
                    "税込",
                    "税抜",
                    "最低購入",
                    "有効期限",
                    "為替",
                    "端数処理",
                ],
            }
            group_matches = {
                group: _contains_any(low, group_keywords)
                for group, group_keywords in required_groups.items()
            }
            matches = [match for values in group_matches.values() for match in values]
            ok = all(group_matches.values())
            results[axis] = {
                "ok": ok,
                "matches": matches,
                "keywords": keywords,
                "required_groups": required_groups,
                "group_matches": group_matches,
            }
        else:
            matches = _contains_any(low, keywords)
            ok = bool(matches)
            results[axis] = {"ok": ok, "matches": matches, "keywords": keywords}
        if not ok:
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
            "the requested amount may be close to the ceiling when itemized, justified, "
            "and traceable to dated official prices or quotations"
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
        skip: comma-separated tool ids to skip, e.g. ``sentence,literature``.
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

    if (
        program in {"kaken_oss", "kaken_oss_platform"}
        and "abstraction" not in skip_set
    ):
        abstraction = grant_writing_named_software_abstraction_check(text)
        detailed_results["named_software_abstraction"] = abstraction
        if abstraction["applicable"]:
            detailed_scores["named_software_abstraction"] = abstraction["score"]
            if abstraction["risks"]:
                priority_issues.append({
                    "tool": "abstraction",
                    "name": "named_software_abstraction_check",
                    "severity": "MEDIUM",
                    "score": abstraction["score"],
                    "comments": abstraction["comments"][:5],
                })

    if (
        program in {"kaken_oss", "kaken_oss_platform"}
        and "vocabulary" not in skip_set
    ):
        vocabulary = grant_writing_reviewer_vocabulary_check(text)
        detailed_results["reviewer_vocabulary"] = vocabulary
        if vocabulary["applicable"]:
            detailed_scores["reviewer_vocabulary"] = vocabulary["score"]
            if vocabulary["risks"]:
                priority_issues.append({
                    "tool": "vocabulary",
                    "name": "reviewer_vocabulary_check",
                    "severity": "MEDIUM",
                    "score": vocabulary["score"],
                    "comments": vocabulary["comments"][:5],
                })

    if "persuasion" not in skip_set:
        persuasion = grant_writing_persuasion_quality_check(text)
        detailed_results["persuasion_quality"] = persuasion
        if persuasion["applicable"]:
            detailed_scores["persuasion_quality"] = persuasion["score"]
            if persuasion["risks"]:
                priority_issues.append({
                    "tool": "persuasion",
                    "name": "persuasion_quality_check",
                    "severity": _severity_from_score(persuasion["score"]),
                    "score": persuasion["score"],
                    "comments": persuasion["comments"][:5],
                })

    if "format" not in skip_set:
        review_format = grant_writing_kaken_review_format_check(text)
        detailed_results["kaken_review_format"] = review_format
        if review_format["applicable"]:
            detailed_scores["kaken_review_format"] = review_format["score"]
            if review_format["risks"]:
                priority_issues.append({
                    "tool": "format",
                    "name": "kaken_review_format_check",
                    "severity": _severity_from_score(review_format["score"]),
                    "score": review_format["score"],
                    "comments": review_format["comments"][:5],
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

    if "domain" not in skip_set:
        domain = grant_writing_domain_outcome_chain_check(text)
        detailed_results["domain_outcome_chain"] = domain
        if domain["applicable"]:
            detailed_scores["domain_outcome_chain"] = domain["score"]
            if domain["comments"]:
                priority_issues.append({
                    "tool": "domain",
                    "name": "domain_outcome_chain_check",
                    "severity": _severity_from_score(domain["score"]),
                    "score": domain["score"],
                    "comments": domain["comments"][:5],
                })

    if "metric" not in skip_set:
        metric = grant_writing_derived_metric_validation_check(text)
        detailed_results["derived_metric_validation"] = metric
        if metric["applicable"]:
            detailed_scores["derived_metric_validation"] = metric["score"]
            if metric["comments"]:
                priority_issues.append({
                    "tool": "metric",
                    "name": "derived_metric_validation_check",
                    "severity": _severity_from_score(metric["score"]),
                    "score": metric["score"],
                    "comments": metric["comments"][:5],
                })

    if "pilot" not in skip_set:
        pilot = grant_writing_cross_organization_pilot_check(text)
        detailed_results["cross_organization_pilot"] = pilot
        if pilot["applicable"]:
            detailed_scores["cross_organization_pilot"] = pilot["score"]
            if pilot["comments"]:
                priority_issues.append({
                    "tool": "pilot",
                    "name": "cross_organization_pilot_check",
                    "severity": _severity_from_score(pilot["score"]),
                    "score": pilot["score"],
                    "comments": pilot["comments"][:5],
                })

    if "literature" not in skip_set:
        literature = grant_writing_literature_gap_evidence_check(text)
        detailed_results["literature_gap_evidence"] = literature
        if literature["applicable"]:
            detailed_scores["literature_gap_evidence"] = literature["score"]
            if literature["risks"]:
                priority_issues.append({
                    "tool": "literature",
                    "name": "literature_gap_evidence_check",
                    "severity": "HIGH",
                    "score": literature["score"],
                    "comments": literature["comments"][:5],
                })

    if "integration" not in skip_set:
        integration = grant_writing_collaborative_integration_risk_check(text)
        detailed_results["collaborative_integration_risk"] = integration
        if integration["applicable"]:
            detailed_scores["collaborative_integration_risk"] = integration["score"]
            if integration["comments"]:
                priority_issues.append({
                    "tool": "integration",
                    "name": "collaborative_integration_risk_check",
                    "severity": _severity_from_score(integration["score"]),
                    "score": integration["score"],
                    "comments": integration["comments"][:5],
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
