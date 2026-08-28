"""Tool functions for the paper-writing domain of mcp-server-document.

All functions are plain callables; they get wrapped by @mcp.tool()
in the top-level server.py. No FastMCP import here (that would
accidentally register a second MCP server).
"""

from __future__ import annotations

import pathlib
import re

# Plan B Tier 1 (v0.12.0) — composite score + human-advisor comments
from .plans.T1 import paper_writing_abstract_strength  # noqa: F401
from .plans.T2 import paper_writing_contribution_clarity_score  # noqa: F401

# Plan B Tier 2 (v0.13.0)
from .plans.T3 import paper_writing_claim_quantification  # noqa: F401
from .plans.T4 import paper_writing_limitation_statement_presence  # noqa: F401
from .plans.T5 import paper_writing_related_work_density  # noqa: F401
from .plans.T6 import paper_writing_figure_referencing_coverage  # noqa: F401
from .plans.T7 import paper_writing_given_new_ordering  # noqa: F401

# Plan B Meta (v0.14.0)
from .plans.T8 import paper_writing_health_report  # noqa: F401
from .plans.T9 import paper_writing_reviewer_2_trigger_summary  # noqa: F401

# Plan B Tier 3 (v0.19.0) — paper-specific axes: triangle / reproducibility /
# stats / citation health / discussion / journal fit
from .plans.T10 import paper_writing_title_abstract_conclusion_triangle  # noqa: F401
from .plans.T11 import paper_writing_reproducibility_open_science_check  # noqa: F401
from .plans.T12 import paper_writing_statistical_reporting_compliance  # noqa: F401
from .plans.T13 import paper_writing_citation_health_4_axes  # noqa: F401
from .plans.T14 import paper_writing_discussion_structure_4_elements  # noqa: F401
from .plans.T15 import paper_writing_journal_fit_assessment  # noqa: F401

# Plan B Tier 4 / Intelligence Layer (v0.24.0) — paper 版 synthesis/orchestration
from .plans.T16 import paper_writing_root_cause_diagnosis  # noqa: F401
from .plans.T17 import paper_writing_next_5_actions  # noqa: F401
from .plans.T18 import paper_writing_run_full_workflow  # noqa: F401
from .plans.T19 import paper_writing_rewrite_suggest  # noqa: F401
from .plans.T20 import paper_writing_adaptive_health_report  # noqa: F401

# Cross-module JA-lint (v0.13.0) — grant_writing の和文 lint を paper 名前空間に re-export
from .cross_lint import (  # noqa: F401
    paper_writing_check_notation_variants,
    paper_writing_find_undefined_acronyms,
    paper_writing_acronym_usage_audit,
    paper_writing_check_kanji_ratio,
    paper_writing_lint_bedrock,
    paper_writing_check_misuse_japanese,
    paper_writing_suggest_redundancy_fixes,
    paper_writing_check_subject_predicate_distance,
)
from ._terminology_normalizer import (  # noqa: F401
    paper_writing_normalize_terminology,
    paper_writing_normalize_terminology_file,
)

# Paper-PDF download tools (v0.25.0, 2026-05-21) — IEEE Xplore + Emerald
# Cookie-seeded session pattern.  Requires caller's IP to have
# institutional access to the target publisher.
from .paper_download import (  # noqa: F401
    paper_writing_resolve_doi,
    paper_writing_ieee_doi_to_arnumber,
    paper_writing_ieee_download_pdf,
    paper_writing_emerald_download_pdf,
    paper_writing_sciencedirect_download_pdf,
    paper_writing_doi_to_bibtex,
    paper_writing_fetch_and_cite,
)


_HERE = pathlib.Path(__file__).resolve().parent
KNOWLEDGE = _HERE


def _load_skill() -> str:
    return (KNOWLEDGE / "skill.md").read_text(encoding="utf-8")


# ------------------------------------------------------------------
# Knowledge loader
# ------------------------------------------------------------------
def paper_writing_usage() -> str:
    """CAE-AI Lab の対象投稿先 (IEEJ / IEEE / 加速器系) の作文技術ガイド全体。

    IMRAD 構造 (Abstract/Intro/Method/Result/Discussion/Conclusion) の
    書き方、contribution 明示、reviewer 2 対策、citation pattern、
    figure placement、abstract 字数、英文 red-flag 検出、cover letter
    テンプレート等の総合ガイド。

    呼び出しタイミング:
    - Journal paper ドラフトのレビュー依頼
    - 「contribution が不明瞭」feedback 対応
    - Reviewer コメントへの response letter 作成
    - Abstract / Introduction の書き直し
    - 英文チェック (冠詞抜け、時制不統一、自動詞/他動詞)
    - 投稿前の最終チェックリスト

    活動ツール (実測、v0.14.0 時点):

    ## 基本 lint (existing)
    - paper_writing_count_underlines / validate_pdf_pages / check_overfull_hbox
    - paper_writing_analyze_sentences / count_weak_expressions
    - paper_writing_validate_abstract_length
    - paper_writing_check_citation_usage / check_abstract_background_ratio
    - paper_writing_check_paragraph_length / check_paragraph_opener
    - paper_writing_check_conclusion_first_use
    - paper_writing_check_passive_voice_ratio / check_tense_consistency
    - paper_writing_check_english_redflags / check_strong_adjective_budget
    - paper_writing_check_figure_caption_showing / check_figure_forward_reference
    - paper_writing_check_equation_numbering / check_self_citation_ratio
    - paper_writing_check_imrad_balance / check_word_repetition
    - paper_writing_check_sentence_ending_variety / lint_reference_format
    - paper_writing_classify_reviewer_comment / generate_cover_letter
    - paper_writing_generate_response_letter

    ## Plan B Tier 1 (v0.12.0) — score + human-advisor comments
    - paper_writing_abstract_strength (T1)
    - paper_writing_contribution_clarity_score (T2)

    ## Plan B Tier 2 (v0.13.0) — score + human-advisor comments
    - paper_writing_claim_quantification (T3)
    - paper_writing_limitation_statement_presence (T4)
    - paper_writing_related_work_density (T5)
    - paper_writing_figure_referencing_coverage (T6)
    - paper_writing_given_new_ordering (T7)

    ## Meta (v0.14.0) — 1 コールで総点検
    - paper_writing_health_report (T8) — 全 Plan B T1-T7 を束ね優先度付け
    - paper_writing_reviewer_2_trigger_summary (T9) — reviewer-2 攻撃面

    ## Advanced PDF + prose detectors (latest)
    - paper_writing_check_pdf_advanced_anomalies — reviewer-2 高頻度
      8 カテゴリ (citation order/space, acronym first use, caption position,
      equation gaps, heading capitalization, orphan/widow, image overlap).
      `check_pdf_obvious_errors` の更に深い後段。score 0-10.
    - paper_writing_check_prose_density — 圧縮 anti-pattern 検出
      (nominalisation, em-dash+semicolon chaining, jargon クラスタ,
      40 語超 long sentence). flag 率 >30% で recommendation が
      `rewrite` → `reduce_content` に切替 (圧縮の床を踏んだ signal).
    - paper_writing_suggest_concept_drops — どの概念を落とすか
      (PARALLEL_CITATION / EM_DASH_INTERPOLATION / PARENTHETICAL_ASIDE /
      TRAILING_SEMICOLON_CLAUSE / CITATION_ONLY_PARENTHETICAL).
      check_prose_density の companion。compression floor を超えた時、
      pattern match で droppable な候補を ranked list で提案。
      投稿先 (digest / journal) ごとに残す/削るは判断必要。
    - paper_writing_check_typography_hacks — 「fontを小さくするはだめ」
      の自動 enforcement。body \\small / linespread<0.95 / page-area
      override を CRITICAL/MAJOR で検出。\\vspace は INFO のみ
      (見た目のための layout tool として OK)。
      Page-limit hack hierarchy の Level 0 detector。
    - paper_writing_check_digest_human_review_triggers — 1-page digest
      human review trigger detector: abstract numeric overload, known
      Warburg/CLN framing, uncited HOIBC/Warburg first use, overloaded
      captions, opaque rank-(m,n) notation, unnumbered/uncited core
      equations, ambiguous block-matrix/N_b notation, and geometry scope
      creep.

    ## Cross-module 和文 lint (v0.13.0) — grant_writing 由来の re-export
    - paper_writing_check_notation_variants / find_undefined_acronyms
    - paper_writing_acronym_usage_audit / check_kanji_ratio
    - paper_writing_lint_bedrock / check_misuse_japanese
    - paper_writing_suggest_redundancy_fixes
    - paper_writing_check_subject_predicate_distance
    - paper_writing_normalize_terminology / normalize_terminology_file
      — 既知の表記ゆれを TeX 安全に補正。既定では `cube=>立方体` を
      日本語文脈だけに適用し、英語 caption と file path は保持。
    """
    return _load_skill()


_TARGET_VENUE_PROFILES = {
    "accelerator": {
        "label": "accelerator community",
        "aliases": (
            "ipac", "jacow", "prab", "physical review accelerators and beams",
            "pasj", "particle accelerator society of japan", "linac", "ibic",
            "napac", "particle accelerator conference", "加速器学会",
        ),
        "primary_language": "English for JACoW/PRAB; follow the official PASJ template",
        "keyword_policy": "Do not impose IEEEkeywords; follow the selected accelerator venue template.",
    },
    "ieej": {
        "label": "Institute of Electrical Engineers of Japan (IEEJ)",
        "aliases": (
            "ieej", "電気学会", "電気学会論文誌", "電気学会研究会",
            "電気学会全国大会", "電気学会部門大会",
        ),
        "primary_language": "Japanese or English body according to the venue, with bilingual metadata when required",
        "keyword_policy": "Accept IEEJ jkeyword/ekeyword blocks; do not require IEEEkeywords.",
    },
    "ieee": {
        "label": "Institute of Electrical and Electronics Engineers (IEEE)",
        "aliases": ("ieee",),
        "primary_language": "English",
        "keyword_policy": "Require the IEEE Index Terms/IEEEkeywords form used by the official template.",
    },
}


def _target_venue_alias_matches(normalized: str, alias: str) -> bool:
    if any(ord(ch) > 127 for ch in alias):
        return alias in normalized
    return re.search(
        rf"(?<![a-z0-9]){re.escape(alias)}(?![a-z0-9])",
        normalized,
    ) is not None


def paper_writing_target_venue_policy(target_venue: str = "") -> dict:
    """Enforce the CAE-AI Lab's deliberately narrow publication targets.

    Supported targets are IEEJ, IEEE, and accelerator-community venues.
    The function classifies a concrete venue name, returns its language and
    keyword lane, and rejects unrelated venue recommendations. Japanese and
    English readability are always evaluated independently and never averaged.
    """
    requested = target_venue.strip()
    supported = ["ieej", "ieee", "accelerator"]
    if not requested:
        return {
            "status": "selection_required",
            "policy_id": "cae_ai_lab_target_venues_v1",
            "target_venue": "",
            "target_category": None,
            "supported_categories": supported,
            "recommendation": (
                "Choose a concrete target in IEEJ, IEEE, or the accelerator "
                "community before applying venue-specific submission rules."
            ),
            "bilingual_policy": (
                "Japanese and English are separate review lanes; scores are "
                "not averaged and the worse applicable lane controls the gate."
            ),
        }

    normalized = requested.casefold().strip()
    matched_category = None
    # Accelerator aliases take precedence for names such as IEEE PAC.
    for category in ("accelerator", "ieej", "ieee"):
        aliases = _TARGET_VENUE_PROFILES[category]["aliases"]
        if any(_target_venue_alias_matches(normalized, alias) for alias in aliases):
            matched_category = category
            break

    if matched_category is None:
        return {
            "status": "reject",
            "policy_id": "cae_ai_lab_target_venues_v1",
            "target_venue": requested,
            "target_category": None,
            "supported_categories": supported,
            "recommendation": (
                "Do not optimize this manuscript for an unrelated venue. "
                "Retarget it to IEEJ, IEEE, or an accelerator-community venue."
            ),
            "bilingual_policy": (
                "Japanese and English are separate review lanes; scores are "
                "not averaged and the worse applicable lane controls the gate."
            ),
        }

    profile = _TARGET_VENUE_PROFILES[matched_category]
    return {
        "status": "pass",
        "policy_id": "cae_ai_lab_target_venues_v1",
        "target_venue": requested,
        "target_category": matched_category,
        "target_label": profile["label"],
        "supported_categories": supported,
        "primary_language": profile["primary_language"],
        "keyword_policy": profile["keyword_policy"],
        "template_policy": (
            "Use the current official template and verify its page, abstract, "
            "keyword, reference, and figure rules for the submission year."
        ),
        "bilingual_policy": (
            "Japanese and English are separate review lanes; scores are not "
            "averaged and the worse applicable lane controls the gate."
        ),
        "recommendation": (
            f"Apply the {matched_category} venue profile and ignore unrelated "
            "journal-fit suggestions."
        ),
    }


# ------------------------------------------------------------------
# Active diagnostic tools (shared across writing MCPs)
# ------------------------------------------------------------------
def paper_writing_count_underlines(tex_path: str) -> dict:
    """TeX/LaTeX 内の下線コマンドを実測。論文では原則ゼロを目指す。

    目標: 0 (journal paper では strong/em で十分、下線は稀)
    警告: >=5 per page
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"file not found: {tex_path}"}
    text = p.read_text(encoding="utf-8", errors="replace")
    macros = ["uline", "underline", "uwave", "uuline",
              "sout", "dashuline", "dotuline"]
    by_macro: dict[str, int] = {}
    total = 0
    for m in macros:
        n = len(re.findall(r"\\" + m + r"\{", text))
        if n:
            by_macro[m] = n
            total += n
    return {
        "file": str(p),
        "total_underlines": total,
        "by_macro": by_macro,
        "target_range": "0 (journal paper では 強調は \\emph / \\textbf)",
        "warning_threshold": ">=5 per page",
    }


def paper_writing_validate_pdf_pages(pdf_path: str, page_limit: int) -> dict:
    """PDF のページ数が投稿制限内か検証。pymupdf が必要。

    CAE-AI Lab の対象は IEEJ / IEEE / 加速器系。ページ制限は
    journal/conference、document type、投稿年度で変わるため、呼び出し側が
    公式テンプレートで確認した ``page_limit`` を明示する。
    """
    try:
        import pymupdf  # type: ignore
    except ImportError:
        return {"error": "pymupdf not installed. `pip install pymupdf`"}
    p = pathlib.Path(pdf_path)
    if not p.exists():
        return {"error": f"file not found: {pdf_path}"}
    doc = pymupdf.open(pdf_path)
    n = doc.page_count
    doc.close()
    return {
        "file": str(p),
        "page_count": n,
        "page_limit": page_limit,
        "within_limit": n <= page_limit,
        "pages_over": max(0, n - page_limit),
    }


def paper_writing_check_overfull_hbox(log_path: str) -> dict:
    """LaTeX ログ中の Overfull \\hbox 警告をカウント。journal では許容ゼロ。"""
    p = pathlib.Path(log_path)
    if not p.exists():
        return {"error": f"file not found: {log_path}"}
    text = p.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(
        r"Overfull \\hbox \(([^)]+)\) in paragraph at lines (\d+)(?:--(\d+))?",
        text,
    )
    details = [
        {"severity": m[0], "lines": m[1] + (f"-{m[2]}" if m[2] else "")}
        for m in matches[:20]
    ]
    return {
        "file": str(p),
        "overfull_count": len(matches),
        "overfull_details": details,
        "target": "0 (journal は許容なし)",
    }


def paper_writing_check_pdf_edge_overflow(
    pdf_path: str,
    margin_quantile: float = 0.99,
    tolerance_pt: float = 1.0,
    max_report: int = 30,
) -> dict:
    """論文 PDF を実際に読んで、本文/数式/図表が紙面の端からはみ出して
    いないかをチェックする (LaTeX ログを介さない直接検証)。

    手順:
      1. pymupdf で各ページの全テキスト span の bbox を取得
      2. 全 span の x_right の `margin_quantile` パーセンタイルを
         "通常の右端" とみなす (テキスト本体の右マージン位置)
      3. その通常右端より `tolerance_pt` を超えてはみ出す span を
         flag (URL の hyphenation などで 1pt 程度のはみ出しは噪音)
      4. 左端 / 上端 / 下端も同様に検査
      5. 図表 / 画像も page.get_drawings() / get_images() で抽出して
         同じ判定を適用

    LaTeX が「Overfull \\hbox warning」を出さないケース (\\sloppy +
    underfull, \\noindent + 数式が長すぎ、figure* / table* で 2 段
    幅にせず長過ぎる) も検出できる。`paper_writing_check_overfull_hbox`
    (log 経由) と組み合わせて使う。

    Args:
        pdf_path: 検査対象の PDF
        margin_quantile: 通常の右端を決める分位点 (0.95-0.99 推奨)
        tolerance_pt: この pt 未満の超過は ignore (デフォルト 1.0 pt
            ≈ 0.35 mm)
        max_report: 報告する overflow の最大件数

    Returns:
        dict: file / page_count / overflow_count / overflows (list) /
              page_bbox_pt / text_area_estimate / target / hint
    """
    try:
        import pymupdf  # type: ignore
    except ImportError:
        return {
            "error": "pymupdf not installed. `pip install pymupdf`",
            "hint": "fallback: paper_writing_check_overfull_hbox(log_path)"
        }

    p = pathlib.Path(pdf_path)
    if not p.exists():
        return {"error": f"file not found: {pdf_path}"}

    doc = pymupdf.open(pdf_path)
    n_pages = doc.page_count

    # ---- Step 1: ページ毎に text spans の bbox を集める ----
    all_right_x: list[float] = []
    all_left_x: list[float] = []
    all_top_y: list[float] = []
    all_bot_y: list[float] = []

    per_page_spans: list[list[dict]] = []  # [page][{bbox, text}]
    page_bboxes: list[tuple[float, float, float, float]] = []

    for pno in range(n_pages):
        page = doc[pno]
        rect = page.rect  # MediaBox / CropBox (pt)
        page_bboxes.append((rect.x0, rect.y0, rect.x1, rect.y1))

        spans_here: list[dict] = []
        try:
            d = page.get_text("dict")
        except Exception:
            d = {"blocks": []}
        for block in d.get("blocks", []):
            for line in block.get("lines", []):
                for span in line.get("spans", []):
                    bbox = span.get("bbox")
                    if not bbox or len(bbox) != 4:
                        continue
                    txt = span.get("text", "").strip()
                    if not txt:
                        continue
                    spans_here.append({
                        "bbox": tuple(bbox),
                        "text": txt,
                    })
                    all_left_x.append(bbox[0])
                    all_top_y.append(bbox[1])
                    all_right_x.append(bbox[2])
                    all_bot_y.append(bbox[3])
        per_page_spans.append(spans_here)

    doc.close()

    if not all_right_x:
        return {
            "file": str(p),
            "page_count": n_pages,
            "overflow_count": 0,
            "overflows": [],
            "target": "0",
            "hint": (
                "PDF に extractable text が見つからない (scan PDF か image-only?)."
                " pymupdf が text を読めないだけかもしれない — \\includepdf / \\IfFileExists"
                " で text を OCR layer なしで埋め込んでいないか確認。"
            ),
        }

    # ---- Step 2: text area の右端 / 左端を分位点で推定 ----
    def _quantile(xs: list[float], q: float) -> float:
        s = sorted(xs)
        if not s:
            return 0.0
        k = max(0, min(len(s) - 1, int(round(q * (len(s) - 1)))))
        return s[k]

    # 「通常の」textarea: 右端 = 99パーセンタイル, 左端 = 1パーセンタイル
    text_right = _quantile(all_right_x, margin_quantile)
    text_left = _quantile(all_left_x, 1.0 - margin_quantile)
    text_top = _quantile(all_top_y, 1.0 - margin_quantile)
    text_bot = _quantile(all_bot_y, margin_quantile)

    # ---- Step 3: overflow span 抽出 ----
    overflows: list[dict] = []
    for pno, spans in enumerate(per_page_spans):
        page_rect = page_bboxes[pno]
        for sp in spans:
            x0, y0, x1, y1 = sp["bbox"]
            txt = sp["text"]

            # right overflow: span 右端が text_right より大きく
            # tolerance_pt を超えていて、かつ page 境界も超えていれば critical
            right_over = x1 - text_right
            left_over = text_left - x0
            top_over = text_top - y0
            bot_over = y1 - text_bot

            # 報告対象: tolerance_pt より大きい超過
            for side, delta in (
                ("right", right_over),
                ("left", left_over),
                ("top", top_over),
                ("bottom", bot_over),
            ):
                if delta > tolerance_pt:
                    # ページ境界突破は CRITICAL
                    px0, py0, px1, py1 = page_rect
                    if side == "right":
                        page_violation = x1 > px1
                    elif side == "left":
                        page_violation = x0 < px0
                    elif side == "top":
                        page_violation = y0 < py0
                    else:
                        page_violation = y1 > py1

                    overflows.append({
                        "page": pno + 1,  # 1-indexed for humans
                        "side": side,
                        "overflow_pt": round(delta, 2),
                        "bbox_pt": [round(v, 2) for v in sp["bbox"]],
                        "text_snippet": (txt[:60] + "...") if len(txt) > 60 else txt,
                        "page_boundary_violation": page_violation,
                    })

    # 報告順: page_boundary_violation を先、その後 overflow_pt 降順
    overflows.sort(
        key=lambda r: (
            not r["page_boundary_violation"],
            -r["overflow_pt"],
        )
    )

    # ---- Step 4: 出力 ----
    page_violations = sum(1 for o in overflows if o["page_boundary_violation"])
    return {
        "file": str(p),
        "page_count": n_pages,
        "overflow_count": len(overflows),
        "page_boundary_violation_count": page_violations,
        "overflows": overflows[:max_report],
        "text_area_estimate_pt": {
            "left": round(text_left, 2),
            "right": round(text_right, 2),
            "top": round(text_top, 2),
            "bottom": round(text_bot, 2),
            "method": f"{int(margin_quantile*100)}/{int((1-margin_quantile)*100)} percentile of span bboxes",
        },
        "page_bbox_first": [round(v, 2) for v in page_bboxes[0]],
        "tolerance_pt": tolerance_pt,
        "target": (
            "page_boundary_violation_count = 0 (journal 投稿で許容ゼロ)。"
            "overflow_count > 0 でも全部 page 内ならレイアウト微調整で済む。"
        ),
        "hint": (
            "page_boundary_violation_count > 0 の場合: 該当行に "
            "\\sloppy / \\hyphenation / \\linebreak、URL は \\url + \\url-style、"
            "数式は \\resizebox / multline / split で折り返し、figure は figure* で"
            "2-column span などで対処。 paper_writing_check_overfull_hbox を併用して"
            ".log 側で line number を特定すると修正が早い。"
        ),
        "source": "pymupdf page.get_text('dict') span bboxes",
    }


def paper_writing_analyze_sentences(text: str, max_len: int = 80) -> dict:
    """文長分析 (和文)。journal では長文を避けて読みやすさ重視。

    目標: 平均 <= 60 字 / 文 (journal は grant 申請書より厳しめ)
    警告: 個別で >= 100 字 / 文

    引数:
        text: 解析対象の日本語テキスト
        max_len: 1 文の上限 (既定 80)
    """
    sentences = [s.strip() for s in re.split(r"[。．]", text) if s.strip()]
    if not sentences:
        return {"error": "no sentences found (no 。 or ．)"}
    lengths = [len(s) for s in sentences]
    avg = sum(lengths) / len(lengths)
    long_ones = [
        {"index": i, "length": ln, "head": s[:40] + ("..." if len(s) > 40 else "")}
        for i, (s, ln) in enumerate(zip(sentences, lengths))
        if ln > max_len
    ]
    return {
        "total_sentences": len(sentences),
        "avg_length": round(avg, 1),
        "max_length": max(lengths),
        "threshold": max_len,
        "over_threshold_count": len(long_ones),
        "over_threshold_examples": long_ones[:5],
        "target": f"avg <= 60 (journal), or avg <= {max_len}",
        "warning": "any sentence >= 100",
    }


def paper_writing_count_weak_expressions(text: str) -> dict:
    """弱気修飾語の出現数。journal 論文では conclusion / contribution で
    使うと査読で突っ込まれる (Reviewer 2 trigger)。

    目標: ゼロ (特に Abstract / Contribution 文)
    警告: >= 3 / ページ
    """
    # JA patterns: local _shared.hedges (was ..\\_shared.hedges in
    # mcp-server-document; now local to paper_writing per 2026-05-26
    # migration to radia-mcp).
    from ._shared.hedges import HEDGE_PATTERNS as _JA_HEDGES
    patterns = dict(_JA_HEDGES)
    patterns.update({
        # English weak hedges (ACS / IEEE にも適用)
        "it seems": r"\bit seems\b",
        "might be": r"\bmight be\b",
        "could be": r"\bcould be\b",
        "tends to": r"\btends to\b",
    })
    by_pat: dict[str, int] = {}
    total = 0
    for name, pat in patterns.items():
        n = len(re.findall(pat, text, re.IGNORECASE))
        if n:
            by_pat[name] = n
            total += n
    return {
        "total_weak_expressions": total,
        "by_pattern": by_pat,
        "target": "0 in Abstract/Contribution, <=2 per page elsewhere",
        "warning": ">=3 per page (Reviewer 2 trigger)",
    }


# ------------------------------------------------------------------
# Paper-specific diagnostic tools
# ------------------------------------------------------------------
def paper_writing_validate_abstract_length(text: str,
                                            max_words: int = 250,
                                            max_chars_ja: int = 400) -> dict:
    """Abstract 字数 / 語数が制限内か検証。言語を自動判定。

    CAE-AI Lab の対象は IEEJ / IEEE / 加速器系。既定値は診断用の一般値で、
    実投稿では選択したvenue・document type・投稿年度の公式上限を
    ``max_words`` / ``max_chars_ja`` に渡す。
    """
    # Crude language detection: if >30% of chars are CJK, treat as Japanese
    cjk = sum(1 for c in text if "぀" <= c <= "ヿ"
              or "一" <= c <= "鿿")
    ja_ratio = cjk / max(1, len(text))
    is_japanese = ja_ratio > 0.3

    if is_japanese:
        n = len(text.replace(" ", "").replace("\n", ""))
        limit = max_chars_ja
        unit = "chars (ja)"
    else:
        n = len(text.split())
        limit = max_words
        unit = "words (en)"

    return {
        "detected_language": "ja" if is_japanese else "en",
        "count": n,
        "unit": unit,
        "limit": limit,
        "within_limit": n <= limit,
        "excess": max(0, n - limit),
        "recommendation": (
            "Abstract は問題提起→手法→結果→意義 の 4 要素で構成。"
            "5-7 文が目安 (IEEE 200 words / IEEJ 400 字)。"
        ),
    }


def paper_writing_check_citation_usage(tex_path: str,
                                        bib_path: str = "") -> dict:
    """TeX 本文中の \\cite{} キーと bib file の entries を突合。
    未使用文献 / 未定義 cite を検出。

    引数:
        tex_path: 論文本体の .tex (メインファイル)
        bib_path: .bib ファイル (省略時は tex_path と同ディレクトリの *.bib)
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"tex file not found: {tex_path}"}
    text = p.read_text(encoding="utf-8", errors="replace")

    # Collect all \cite{...} keys (multiple keys in one \cite separated by ,)
    cite_keys: set[str] = set()
    for m in re.finditer(r"\\(?:cite|citep|citet|autocite|parencite)\*?"
                         r"(?:\[[^\]]*\])?\{([^}]+)\}", text):
        for k in m.group(1).split(","):
            k = k.strip()
            if k:
                cite_keys.add(k)

    # Find bib file
    if bib_path:
        bibs = [pathlib.Path(bib_path)]
    else:
        bibs = list(p.parent.glob("*.bib"))
    if not bibs:
        return {
            "error": "no .bib file found (pass bib_path= or place *.bib next to tex)"
        }

    bib_entries: set[str] = set()
    for bp in bibs:
        btext = bp.read_text(encoding="utf-8", errors="replace")
        for m in re.finditer(r"@\w+\s*\{\s*([^,\s]+)", btext):
            bib_entries.add(m.group(1))

    missing_in_bib = cite_keys - bib_entries  # cited but not in bib
    unused_in_tex = bib_entries - cite_keys   # in bib but not cited

    return {
        "tex_file": str(p),
        "bib_files": [str(b) for b in bibs],
        "total_citations": len(cite_keys),
        "total_bib_entries": len(bib_entries),
        "missing_in_bib_count": len(missing_in_bib),
        "missing_in_bib": sorted(missing_in_bib)[:20],
        "unused_in_tex_count": len(unused_in_tex),
        "unused_in_tex": sorted(unused_in_tex)[:20],
        "recommendation": (
            "missing_in_bib があれば .bib に追加。"
            "unused_in_tex が 5 件以上ある場合は refs を整理 "
            "(実際に論じていない文献の引用は査読者に嫌われる)。"
        ),
    }


def paper_writing_check_english_redflags(text: str) -> dict:
    """英文論文の典型的 red flag を検出 (冠詞、時制、自動詞/他動詞 の混同)。

    NOTE: 網羅的ではなく高頻度パターンのみ。通った後も native proof は
          推奨。
    """
    issues: list[dict] = []

    # Definite/indefinite article mis-use signals (limited heuristics)
    if re.search(r"\bIn the (paper|study|work)\b", text):
        issues.append({
            "pattern": "In the paper/study/work (without first reference)",
            "hint": "初出は 'In this paper/study/work'。'the' は既知物に限る。",
        })
    # Tense consistency (past result + present claim is fine, but all-past in
    # discussion is weak)
    # "will be" in future-work without section header
    if re.search(r"\bwe will\b", text, re.IGNORECASE):
        count = len(re.findall(r"\bwe will\b", text, re.IGNORECASE))
        issues.append({
            "pattern": f"'we will' appears {count} times",
            "hint": "Future work 以外では使わない。Method は過去形または現在形 で記述。",
        })
    # "Please refer to"
    if re.search(r"\bPlease refer to\b", text, re.IGNORECASE):
        issues.append({
            "pattern": "Please refer to",
            "hint": "論文調ではない。'See Section X' or 'Fig. Y shows ...' に置換。",
        })
    # "In order to" (can almost always be just "to")
    n_io = len(re.findall(r"\bIn order to\b", text, re.IGNORECASE))
    if n_io >= 3:
        issues.append({
            "pattern": f"'In order to' appears {n_io} times",
            "hint": "almost always replaceable with just 'to'. 冗長表現。",
        })
    # Capital after semicolon (subtle bug)
    if re.search(r";\s+[A-Z][a-z]", text):
        issues.append({
            "pattern": "Capital letter after semicolon",
            "hint": "セミコロン後は小文字が標準。'; However, ...' は '. However, ...' へ。",
        })
    # "which is" / "that is" redundancy
    n_wi = len(re.findall(r"\bwhich is\b", text, re.IGNORECASE))
    if n_wi >= 5:
        issues.append({
            "pattern": f"'which is' appears {n_wi} times",
            "hint": "'the coil, which is made of copper' -> 'the copper coil' が簡潔。",
        })
    # Passive voice over-use
    passive_count = len(re.findall(
        r"\b(is|are|was|were|been|being)\s+\w+ed\b", text))
    total_verbs_est = len(re.findall(r"\b\w+ed\b", text)) + 1
    passive_ratio = passive_count / total_verbs_est
    if passive_ratio > 0.4:
        issues.append({
            "pattern": f"High passive voice ratio ({passive_ratio:.0%})",
            "hint": "Method は passive 可だが Discussion/Contribution は active 推奨。",
        })

    return {
        "total_issues": len(issues),
        "issues": issues,
        "note": (
            "Heuristic のみ。最終的には native proof (IEEE Author Tools / "
            "DeepL Write / Grammarly Pro) を併用。"
        ),
    }


# ------------------------------------------------------------------
# Deep-learning diagnostic tools (from アクセプト術 2026-04-21 upgrade)
# ------------------------------------------------------------------

_BANNED_INTRO_OPENERS = [
    r"^\s*In this (paper|section|study|work|chapter)",
    r"^\s*This (paper|section|study|work) (presents|describes|deals)",
    r"^\s*Recent advances in",
    r"^\s*The last few years (have seen|has seen)",
    r"^\s*It is (well[- ]known|often said)",
    r"^\s*There (is|are) (many|several|a number)",
    r"^\s*First(ly)?\s*,",
    r"^\s*Several studies",
    r"^\s*Many (authors|researchers)",
]


def paper_writing_check_paragraph_opener(text: str) -> dict:
    """Introduction/Abstract の段落冒頭が禁断フレーズで始まるか検出。

    Wallwork §13.15, §14.11 のブラックリスト:
    - "In this paper, ..."
    - "Recent advances in ..."
    - "The last few years ..."
    - "It is well-known that ..."
    - "Several studies ..."

    目標: 0 (特に Abstract + Introduction の最初の段落)
    警告: >=2 / section
    """
    # Split paragraphs by blank lines
    paragraphs = [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip()]
    flagged: list[dict] = []
    for i, para in enumerate(paragraphs, start=1):
        first_line = para.split("\n", 1)[0]
        for pat in _BANNED_INTRO_OPENERS:
            if re.match(pat, first_line, re.IGNORECASE):
                flagged.append({
                    "paragraph_index": i,
                    "pattern": pat,
                    "head": first_line[:80],
                })
                break
    return {
        "total_paragraphs": len(paragraphs),
        "weak_opener_count": len(flagged),
        "flagged": flagged[:10],
        "target": "0 in Abstract/Intro opening paragraphs",
        "warning": ">=2 in a single section",
        "source": "Wallwork §13.15, §14.11",
    }


def paper_writing_check_paragraph_length(text: str,
                                          min_chars_ja: int = 150,
                                          max_chars_ja: int = 400,
                                          min_words_en: int = 75,
                                          max_words_en: int = 175,
                                          ) -> dict:
    """段落の字数/語数が適正範囲内か検出。

    目標:
        EN: 75-175 words/段落 (Wallwork §14.9)
        JP: 150-400 字/段落
    警告: >200 words / >500 字
    致命: >300 words / >800 字、または Introduction が 1 段落のみ

    引数:
        text: 複数段落の本文 (空行で区切る)
    """
    paragraphs = [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip()]

    def _is_japanese(s):
        cjk = sum(1 for c in s if "一" <= c <= "鿿" or "぀" <= c <= "ヿ")
        return cjk / max(1, len(s)) > 0.3

    lang_ja = _is_japanese(text)
    vmin = min_chars_ja if lang_ja else min_words_en
    vmax = max_chars_ja if lang_ja else max_words_en
    unit = "chars (ja)" if lang_ja else "words (en)"

    stats = []
    for i, para in enumerate(paragraphs, start=1):
        if lang_ja:
            n = len(para.replace(" ", "").replace("\n", ""))
        else:
            n = len(para.split())
        stats.append({"paragraph": i, "count": n,
                      "within": vmin <= n <= vmax})
    over_count = sum(1 for s in stats if s["count"] > vmax)
    way_over = sum(1 for s in stats if s["count"] > vmax * 1.5)
    under_count = sum(1 for s in stats if s["count"] < vmin)

    too_few_paragraphs = len(paragraphs) == 1 and sum(
        s["count"] for s in stats
    ) > vmin * 2

    return {
        "total_paragraphs": len(paragraphs),
        "language": "ja" if lang_ja else "en",
        "unit": unit,
        "target_range": f"{vmin}-{vmax} {unit}",
        "over_warn": over_count,
        "over_fatal": way_over,
        "under": under_count,
        "single_giant_paragraph": too_few_paragraphs,
        "per_paragraph": stats[:20],
        "source": "Wallwork §14.9",
    }


def _paper_writing_split_conclusion(src: str) -> tuple[str, str, str]:
    """Return text before conclusion, conclusion body, and detected heading."""
    tex_heading = re.compile(
        r"\\section\*?\s*\{([^{}]+)\}",
        re.IGNORECASE,
    )
    heading_name = re.compile(
        r"^(?:\d+(?:[.\u30fb・]\d+)*[.\u30fb　 ]*)?"
        r"(?:conclusions?|summary|closing remarks|結論|まとめ|結語)\s*$",
        re.IGNORECASE,
    )
    matches = list(tex_heading.finditer(src))
    for index, match in enumerate(matches):
        title = _paper_writing_plain_text(match.group(1)).strip()
        if not heading_name.match(title):
            continue
        end = matches[index + 1].start() if index + 1 < len(matches) else len(src)
        for boundary in (
                r"\begin{thebibliography}",
                r"\bibliographystyle",
                r"\bibliography",
                r"\end{document}"):
            boundary_at = src.find(boundary, match.end(), end)
            if boundary_at >= 0:
                end = min(end, boundary_at)
        return src[:match.start()], src[match.end():end], title

    line_heading = re.compile(
        r"(?im)^(?:#{1,4}\s*)?"
        r"((?:\d+(?:\.\d+)*[.\s]*)?"
        r"(?:conclusions?|summary|closing remarks|結論|まとめ|結語))"
        r"\s*$"
    )
    match = line_heading.search(src)
    if match:
        next_heading = re.search(r"(?m)^#{1,4}\s+.+$", src[match.end():])
        end = (match.end() + next_heading.start()
               if next_heading else len(src))
        return src[:match.start()], src[match.end():end], match.group(1).strip()
    return src, "", ""


def _paper_writing_conclusion_technical_tokens(src: str) -> set[str]:
    """Extract acronym-like and mixed-case technical names."""
    plain = _paper_writing_plain_text(src)
    patterns = (
        r"\b[A-Z]{2,}[A-Za-z0-9-]*\b",
        r"\b[A-Z][a-z]+(?:[A-Z][A-Za-z0-9]*)+\b",
        r"\b[A-Za-z]+\d+[A-Za-z0-9-]*\b",
    )
    return {
        token
        for pattern in patterns
        for token in re.findall(pattern, plain)
    }


def _paper_writing_conclusion_numeric_claims(src: str) -> set[str]:
    """Extract decimals/exponents and integers carrying a unit or ratio."""
    plain = _paper_writing_plain_text(src)
    unit = (
        r"(?:%|％|倍|×|[munpkMGTµμ]?(?:s|Hz|A|V|W|K|T|H|F|m|g|Pa|"
        r"Ω|ohm)(?:/[A-Za-z]+)?)"
    )
    pattern = re.compile(
        rf"(?<![A-Za-z0-9_.])(?:[-+]?\d+\.\d+(?:[eE][-+]?\d+)?"
        rf"|[-+]?\d+(?:[eE][-+]?\d+)|[-+]?\d+\s*{unit})"
    )
    return {re.sub(r"\s+", "", m.group(0)) for m in pattern.finditer(plain)}


def _paper_writing_conclusion_symbols(src: str) -> set[str]:
    """Extract subscripted Latin symbols from inline/display mathematics."""
    math_parts = re.findall(r"\$([^$]+)\$|\\\((.+?)\\\)|\\\[(.+?)\\\]", src,
                            flags=re.DOTALL)
    joined = " ".join(part for group in math_parts for part in group if part)
    symbols = re.findall(r"\b[A-Za-z]+_\{?[A-Za-z0-9]+\}?", joined)
    return {re.sub(r"[{}\s]", "", symbol) for symbol in symbols}


def _paper_writing_citation_keys(src: str) -> set[str]:
    keys: set[str] = set()
    for match in re.finditer(
            r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\])?\{([^{}]+)\}", src):
        keys.update(key.strip() for key in match.group(1).split(",") if key.strip())
    return keys


def paper_writing_check_conclusion_first_use(
        text_or_path: str,
        whitelist: str = "",
        max_report: int = 30) -> dict:
    """結論で初出となる技術語・変数・数値・引用を検出。

    結論には、本文を統合して得た新しい洞察・含意・展望を記載
    できる。一方、手法名、略語、数式記号、定量値、参考文献が結論で
    初めて現れる場合は、本文で説明していない証拠・方法・主張を
    結論に持ち込んだ可能性がある。本ツールはその候補を抽出する。
    """
    src, source_path = _paper_writing_text_or_path(text_or_path)
    src = _paper_writing_strip_comments(src)
    before, conclusion, heading = _paper_writing_split_conclusion(src)
    if not conclusion.strip():
        return {
            "score": None,
            "conclusion_found": False,
            "source_path": source_path,
            "comments": [
                "Conclusion/結論/まとめ section を検出できない。"
                "LaTeX \\section{} または Markdown heading を明示。"
            ],
        }

    default_whitelist = {
        "AC", "DC", "SI", "IEEE", "IEEJ", "Fig", "Figs", "Eq", "Eqs",
        "Table", "Tables", "Section", "Sections", "RMSE", "CPU", "GPU",
    }
    user_whitelist = {
        item.strip().lower() for item in re.split(r"[,\s]+", whitelist)
        if item.strip()
    }
    allowed = {item.lower() for item in default_whitelist} | user_whitelist

    before_tokens = {
        token.lower() for token in _paper_writing_conclusion_technical_tokens(before)
    }
    conclusion_tokens = _paper_writing_conclusion_technical_tokens(conclusion)
    new_terms = sorted(
        token for token in conclusion_tokens
        if token.lower() not in before_tokens
        and token.lower() not in allowed
        and token.lower().split("-", 1)[0] not in allowed
    )

    before_symbols = {
        symbol.lower() for symbol in _paper_writing_conclusion_symbols(before)
    }
    new_symbols = sorted(
        symbol for symbol in _paper_writing_conclusion_symbols(conclusion)
        if symbol.lower() not in before_symbols
    )

    before_numbers = _paper_writing_conclusion_numeric_claims(before)
    new_numbers = sorted(
        number for number in _paper_writing_conclusion_numeric_claims(conclusion)
        if number not in before_numbers
    )

    before_citations = _paper_writing_citation_keys(before)
    new_citations = sorted(_paper_writing_citation_keys(conclusion) - before_citations)

    issue_count = len(new_terms) + len(new_symbols) + len(new_numbers) + len(new_citations)
    score = max(
        0.0,
        10.0
        - 1.5 * (len(new_terms) + len(new_symbols))
        - 1.0 * (len(new_numbers) + len(new_citations)),
    )
    comments = []
    if new_terms:
        comments.append(
            f"結論で初出の技術語候補: {new_terms[:max_report]}"
        )
    if new_symbols:
        comments.append(
            f"結論で初出の数式記号候補: {new_symbols[:max_report]}"
        )
    if new_numbers:
        comments.append(
            f"結論で初出の定量値候補: {new_numbers[:max_report]}"
        )
    if new_citations:
        comments.append(
            f"結論で初出の citation key: {new_citations[:max_report]}"
        )
    if not comments:
        comments.append(
            "結論の技術語・数式記号・定量値・引用は本文で導入済み。"
        )

    return {
        "score": round(score, 1),
        "score_max": 10,
        "passed": issue_count == 0,
        "conclusion_found": True,
        "conclusion_heading": heading,
        "source_path": source_path,
        "issue_count": issue_count,
        "new_technical_terms": new_terms[:max_report],
        "new_math_symbols": new_symbols[:max_report],
        "new_numeric_claims": new_numbers[:max_report],
        "new_citation_keys": new_citations[:max_report],
        "comments": comments,
        "policy": (
            "結論で新しくすべきなのは、本文の結果から導く統合的な洞察・"
            "含意・展望。新しい手法、データ、変数、定量結果、引用は"
            "本文で先に導入。"
        ),
        "hint": (
            "文字列による補助診断。一般的な略語や意図的な将来展望は "
            "whitelist または目視で除外。"
        ),
        "source": (
            "Wallwork『日本人研究者のための論文の書き方・"
            "アクセプト術』第19章, 木下『理科系の作文技術』"
            "目標規定文"
        ),
    }


def paper_writing_check_abstract_background_ratio(abstract: str) -> dict:
    """Abstract 内で background (導入文) が占める割合を推定。

    Wallwork §13.16: >25% background は日本型 abstract で editor が reject。

    推定法: 第 1 文から順に、動詞に 'we'/'our' が無く数値結果が無い文を
    "background" として分類。数値 or 能動的動詞 (we propose/show/
    demonstrate/find/report) が出現した文から "content" に切替。

    Args:
        abstract: Abstract 本文のみ。
    """
    if not abstract.strip():
        return {"error": "abstract is empty"}
    # Split sentences: handle both .  and 。
    sentences = re.split(r"(?<=[.。])\s+", abstract.strip())
    sentences = [s for s in sentences if s.strip()]
    if not sentences:
        return {"error": "no sentences found"}

    background_chars = 0
    total_chars = 0
    break_index = None
    for i, s in enumerate(sentences):
        total_chars += len(s)
        has_digit = bool(re.search(r"\d", s))
        has_agent = bool(re.search(
            r"\b(we|our|I)\s+(propose|show|demonstrate|find|"
            r"report|present|develop|derive|measure|compare|"
            r"introduce)\b", s, re.IGNORECASE))
        if break_index is None:
            if has_digit or has_agent:
                break_index = i
            else:
                background_chars += len(s)

    ratio = background_chars / max(1, total_chars)
    return {
        "total_sentences": len(sentences),
        "background_break_index": break_index,
        "background_chars": background_chars,
        "total_chars": total_chars,
        "background_ratio": round(ratio, 3),
        "target": "<= 0.25 (Wallwork §13.16)",
        "warning": "> 0.25 (Japanese-style abstract)",
        "fatal": "> 0.40 (near-certain reject)",
    }


def paper_writing_check_abstract_no_math_no_citation(abstract: str) -> dict:
    """Abstract 内に数式 (math), citation, domain acronym が混入していないか検出。

    Editorial convention (IEEE, Elsevier, Springer, APA, ACS, RSC, AIP, AAAS,
    Nature, Science, ほぼ全ての主要 journal):

      - **数式は abstract に入れない**.  検索エンジンと indexing service が
        TeX 数式を扱えず、abstract scrape 時に化けるか欠落する.  数値結果
        (e.g. ``5.3% improvement``) は OK、TeX 数式 (``$E = mc^2$`` や
        ``\\begin{equation}``) は NG.
      - **citation は abstract に入れない**.  abstract は self-contained
        であるべきで、引用は本文で文脈付きで導入する.  例外: critical な
        先行研究を必ず指摘する必要のある一部 conference (たとえば一部の
        ML conf) — そういう場では本 lint は warning として受け取る.
      - **domain-specific acronym は abstract に入れない**.  FEM, BEM,
        MCP, LLM などの略語は abstract では避け、本文初出で
        ``Finite Element Method (FEM)`` のように展開する.  Universal acronym
        (PDF, USB, CPU など) は除外する.

    Detects:
      - Inline math: ``$...$``, ``\\(...\\)`` (但し escape `\\$` は除く)
      - Display math: ``\\[...\\]``, ``\\begin{equation}``, ``\\begin{align}``,
        ``\\begin{equation*}``, ``\\begin{align*}``, ``\\begin{eqnarray}``,
        ``\\begin{gather}``, ``\\begin{multline}``
      - Citations: ``\\cite{...}``, ``\\citet{...}``, ``\\citep{...}``,
        ``\\citeyear{...}``, ``\\citealt{...}``, ``\\citeauthor{...}``,
        ``\\nocite{...}``, ``\\fullcite{...}``, ``\\textcite{...}``,
        ``\\parencite{...}``, and bare numeric markers such as ``[1]`` or
        ``[1, 2]``
      - Domain-specific acronyms: all-caps 2-5-character tokens not in the
        universal whitelist, after math/citation/comment stripping

    Args:
        abstract: Abstract 本文のみ (``\\begin{abstract}`` ... ``\\end{abstract}``
            の中身).  ``paper_writing_extract_abstract`` で抽出した
            ``["abstract"]`` をそのまま渡す.

    Returns:
        dict with:

          - ``n_math_inline``: int
          - ``n_math_display``: int
          - ``n_citations``: int
          - ``n_acronyms``: int
          - ``math_snippets``: list of up to 5 offending math snippets (truncated)
          - ``citation_snippets``: list of up to 5 offending citation snippets
          - ``acronym_snippets``: list of up to 5 acronym contexts
          - ``status``: ``"clean"`` | ``"warning"`` | ``"fail"``
          - ``advice``: str, what to do about each finding
    """
    if not abstract or not abstract.strip():
        return {"error": "abstract is empty"}

    # ---- Math detection -------------------------------------------------
    # Inline $...$ (single dollar): NOT preceded by backslash (escaped \$).
    # The lazy match catches single inline math; we count occurrences.
    inline_math_re = re.compile(r"(?<!\\)\$[^\$\n]+?(?<!\\)\$")
    # Inline \(...\) — pair-based
    inline_paren_re = re.compile(r"\\\(.+?\\\)", flags=re.DOTALL)
    # Display math: \[...\] and the equation/align/eqnarray/gather/multline
    # environments (and their starred variants).
    display_bracket_re = re.compile(r"\\\[.+?\\\]", flags=re.DOTALL)
    display_env_re = re.compile(
        r"\\begin\{(equation\*?|align\*?|eqnarray\*?|gather\*?|multline\*?|"
        r"flalign\*?|alignat\*?)\}.+?\\end\{\1\}",
        flags=re.DOTALL,
    )

    inline_matches = inline_math_re.findall(abstract) + inline_paren_re.findall(abstract)
    # finditer + .group(0) so display_env returns the FULL match (the
    # findall above on env_re only returns the env-name group, not the
    # whole equation block).
    display_matches = (
        display_bracket_re.findall(abstract)
        + [m.group(0) for m in display_env_re.finditer(abstract)]
    )

    n_math_inline = len(inline_matches)
    n_math_display = len(display_matches)
    math_snippets = [s[:80] + ("..." if len(s) > 80 else "")
                     for s in (inline_matches + display_matches)[:5]]

    # ---- Citation detection ---------------------------------------------
    citation_re = re.compile(
        r"\\(cite|citet|citep|citeyear|citealt|citealp|citeauthor|"
        r"nocite|fullcite|textcite|parencite|footcite|autocite)"
        r"(?:\[[^\]]*\])?\{[^\}]+\}"
    )
    citation_command_spans = [m.group(0) for m in citation_re.finditer(abstract)]
    numeric_citation_re = re.compile(
        r"(?<![A-Za-z0-9])\[\s*\d+(?:\s*(?:,|-|--)\s*\d+)*\s*\](?![A-Za-z0-9])"
    )
    numeric_citation_spans = [
        m.group(0) for m in numeric_citation_re.finditer(abstract)
    ]
    citation_spans = citation_command_spans + numeric_citation_spans
    n_citation_commands = len(citation_command_spans)
    n_citation_markers = len(numeric_citation_spans)
    n_citations = len(citation_spans)
    citation_snippets = [s[:80] + ("..." if len(s) > 80 else "")
                         for s in citation_spans[:5]]

    # ---- Domain-acronym detection --------------------------------------
    # Abstracts should be readable without local acronym decoding.  Reuse
    # the same conservative whitelist and stripper as the first-use acronym
    # checker so bib keys, equations, URLs, and comments do not leak in.
    from ._undefined_acronyms import (  # local import avoids a hard cycle
        UNIVERSAL_ACRONYM_WHITELIST,
        _stripped_source,
    )

    stripped_for_acronyms = _stripped_source(abstract)
    acronym_re = re.compile(r"\b([A-Z][A-Z0-9]{1,4})\b")
    acronym_records: dict[str, dict] = {}
    for m in acronym_re.finditer(stripped_for_acronyms):
        acronym = m.group(1)
        if acronym in UNIVERSAL_ACRONYM_WHITELIST:
            continue
        if acronym in acronym_records:
            acronym_records[acronym]["n_occurrences"] += 1
            continue
        a = max(0, m.start() - 45)
        b = min(len(abstract), m.end() + 45)
        acronym_records[acronym] = {
            "acronym": acronym,
            "first_occurrence_offset": m.start(),
            "context": abstract[a:b].replace("\n", " "),
            "n_occurrences": 1,
        }
    acronyms = sorted(acronym_records)
    n_acronyms = len(acronyms)
    acronym_snippets = [
        acronym_records[a]["context"][:100]
        + ("..." if len(acronym_records[a]["context"]) > 100 else "")
        for a in acronyms[:5]
    ]

    # ---- Status + advice ------------------------------------------------
    n_total = n_math_inline + n_math_display + n_citations + n_acronyms
    if n_total == 0:
        status = "clean"
        advice = (
            "Abstract is free of math, citations, and domain acronyms — "
            "matches convention."
        )
    elif n_math_display > 0 or n_citations > 0:
        status = "fail"
        parts = []
        if n_math_display > 0:
            parts.append(
                f"{n_math_display} display-math block(s) found "
                "(\\[...\\] or equation/align environment).  Move to the body; "
                "describe the equation in prose in the abstract."
            )
        if n_citations > 0:
            parts.append(
                f"{n_citations} citation(s) found.  Remove from abstract; "
                "introduce the prior work in the Introduction / Related Work, "
                "then the abstract can just say e.g. 'building on recent work'."
            )
        if n_math_inline > 0:
            parts.append(
                f"{n_math_inline} inline-math span(s) ($...$) found.  "
                "Replace with plain text (e.g. '$E_x$' → 'the x-component "
                "of E', '$\\sigma$' → 'sigma' or 'electrical conductivity')."
            )
        if n_acronyms > 0:
            parts.append(
                f"{n_acronyms} domain acronym(s) found "
                f"({', '.join(acronyms[:8])}).  Avoid acronyms in the "
                "abstract; introduce them on first use in the body, e.g. "
                "'Finite Element Method (FEM)'."
            )
        advice = "  ".join(parts)
    else:
        # only inline math and/or acronyms, no display/citation
        status = "warning"
        parts = []
        if n_math_inline > 0:
            parts.append(
                f"{n_math_inline} inline-math span(s) ($...$) found.  Editors at "
                "most journals will ask you to replace these with plain text.  "
                "Replace e.g. '$E_x$' with 'the x-component of E'."
            )
        if n_acronyms > 0:
            parts.append(
                f"{n_acronyms} domain acronym(s) found "
                f"({', '.join(acronyms[:8])}).  Prefer plain words in the "
                "abstract; define the acronym at first body use, e.g. "
                "'Finite Element Method (FEM)'."
            )
        advice = "  ".join(parts)

    return {
        "n_math_inline": n_math_inline,
        "n_math_display": n_math_display,
        "n_citations": n_citations,
        "n_citation_commands": n_citation_commands,
        "n_citation_markers": n_citation_markers,
        "n_acronyms": n_acronyms,
        "acronyms": acronyms,
        "math_snippets": math_snippets,
        "citation_snippets": citation_snippets,
        "acronym_snippets": acronym_snippets,
        "status": status,
        "advice": advice,
        "rule": ("abstract must be self-contained: no TeX math, no \\cite{} "
                 "/ [N] citation markers, no unexplained domain acronyms "
                 "(IEEE / Elsevier / Springer / Nature / Science convention)"),
    }


def paper_writing_check_figure_caption_showing(caption: str) -> dict:
    """Figure caption が showing (describe) 形か telling (claim) 形か判定。

    Wallwork §17.9: "Figure 4 shows X" (NG) より "X was observed (Fig. 4)"
    (OK)。caption の主語が figure ではなく主張動詞。
    ただし、caption にだけ数値・解釈・結論を置かない。caption は図の同定と
    最小限の読解補助に使い、主張は本文にも必ず書く。caption が長くなる場合は
    本文を優先し、caption を短くする。

    Args:
        caption: 1 つの figure caption テキスト (Fig. 番号から末尾まで)
    """
    cap = caption.strip()
    if not cap:
        return {"error": "caption is empty"}

    # Detect "Fig. N shows/illustrates/presents ..." pattern
    tells = re.search(
        r"^(?:Fig(?:ure)?\.?\s*\d+)?\.?\s*"
        r"(shows|illustrates|presents|depicts|gives|describes)",
        cap, re.IGNORECASE
    )

    return {
        "caption_head": cap[:80],
        "is_showing_form": bool(tells),
        "recommendation": (
            "Convert to 'claim + (Fig. N)' form. "
            "NG: 'Figure 4 shows the relationship between A and B.' "
            "OK: 'The abundances of A and B were inversely related (Fig. 4).' "
            "Then restate the factual claim, numerical result, or interpretation "
            "in the body text."
        ) if tells else (
            "OK: caption opens with a claim, not a describe verb.  Still ensure "
            "that any factual claim, numerical result, or interpretation in the "
            "caption is also stated in the body text."
        ),
        "caption_body_policy": (
            "Captions identify what is plotted and help the figure stand alone, "
            "but they must not be the only place where results or conclusions "
            "appear.  If the caption becomes long, move explanation and "
            "interpretation to the body and keep the caption concise."
        ),
        "source": "Wallwork §17.9 + Sugahara caption-body policy 2026-07-07",
    }


def _paper_writing_text_or_path(text_or_path: str) -> tuple[str, str | None]:
    """Return source text and an optional path label."""
    candidate = (text_or_path or "").strip().strip("\"'")
    if candidate and "\n" not in candidate and len(candidate) < 1024:
        try:
            p = pathlib.Path(candidate)
            if p.exists() and p.is_file():
                try:
                    return p.read_text(encoding="utf-8"), str(p)
                except UnicodeDecodeError:
                    return p.read_text(encoding="cp932"), str(p)
        except (OSError, ValueError):
            pass
    return text_or_path or "", None


def _paper_writing_strip_comments(src: str) -> str:
    """Strip TeX comments conservatively, preserving escaped percent signs."""
    out = []
    for line in src.splitlines():
        cut = len(line)
        for m in re.finditer("%", line):
            if m.start() == 0 or line[m.start() - 1] != "\\":
                cut = m.start()
                break
        out.append(line[:cut])
    return "\n".join(out)


def _paper_writing_extract_environment(src: str, env: str) -> str:
    m = re.search(
        rf"\\begin\{{{re.escape(env)}\}}(.*?)\\end\{{{re.escape(env)}\}}",
        src,
        flags=re.DOTALL,
    )
    return m.group(1).strip() if m else ""


def _paper_writing_extract_braced_arg(src: str, open_brace: int) -> tuple[int, str]:
    if open_brace < 0 or open_brace >= len(src) or src[open_brace] != "{":
        return -1, ""
    depth = 1
    i = open_brace + 1
    while i < len(src):
        ch = src[i]
        if ch == "\\" and i + 1 < len(src) and src[i + 1] in "{}":
            i += 2
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i, src[open_brace + 1:i]
        i += 1
    return -1, ""


def _paper_writing_iter_command_args(src: str, command: str) -> list[dict]:
    args = []
    pat = re.compile(rf"\\{re.escape(command)}(?:\[[^\]]*\])?\s*\{{")
    for m in pat.finditer(src):
        open_brace = m.end() - 1
        end, arg = _paper_writing_extract_braced_arg(src, open_brace)
        if end >= 0:
            args.append({
                "command": command,
                "start": m.start(),
                "end": end + 1,
                "text": arg,
            })
    return args


def _paper_writing_plain_text(src: str) -> str:
    txt = _paper_writing_strip_comments(src)
    txt = re.sub(r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\])?\{[^{}]*\}", " ", txt)
    txt = re.sub(r"\\ref\{[^{}]*\}|\\eqref\{[^{}]*\}|\\label\{[^{}]*\}", " ", txt)
    txt = re.sub(r"\$[^$]*\$", " MATH ", txt)
    txt = re.sub(r"\\\(.+?\\\)|\\\[.+?\\\]", " MATH ", txt, flags=re.DOTALL)
    txt = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?", " ", txt)
    txt = re.sub(r"[{}]", " ", txt)
    txt = txt.replace(r"\%", "%")
    return re.sub(r"\s+", " ", txt).strip()


_READABILITY_SKIP_ENVIRONMENTS = (
    "equation", "equation*", "align", "align*", "alignat", "alignat*",
    "gather", "gather*", "multline", "multline*", "figure", "figure*",
    "table", "table*", "tabular", "tabular*", "thebibliography",
    "jkeyword", "ekeyword", "verbatim", "lstlisting", "tikzpicture",
)


def _paper_writing_readability_prose(src: str) -> str:
    """Return TeX prose while preserving source line positions."""

    def _keep_newlines(match: re.Match) -> str:
        return "\n" * match.group(0).count("\n")

    txt = _paper_writing_strip_comments(src)
    document = re.search(r"\\begin\{document\}", txt)
    if document:
        prefix = txt[:document.end()]
        txt = "\n" * prefix.count("\n") + txt[document.end():]
    txt = re.sub(r"\\end\{document\}[\s\S]*$", _keep_newlines, txt)

    for env in _READABILITY_SKIP_ENVIRONMENTS:
        txt = re.sub(
            rf"\\begin\{{{re.escape(env)}\}}[\s\S]*?"
            rf"\\end\{{{re.escape(env)}\}}",
            _keep_newlines,
            txt,
        )

    txt = re.sub(r"\\\[[\s\S]*?\\\]", _keep_newlines, txt)
    txt = re.sub(r"\$\$[\s\S]*?\$\$", _keep_newlines, txt)
    txt = re.sub(r"\$[^$\n]*\$", " MATH ", txt)
    txt = re.sub(r"\\\([^\n]*?\\\)", " MATH ", txt)
    txt = re.sub(
        r"\\(?:cite|citep|citet|autocite|parencite)[a-zA-Z]*\*?"
        r"(?:\[[^\]]*\])?\{[^{}]*\}",
        " REF ",
        txt,
    )
    txt = re.sub(
        r"\\(?:ref|eqref|pageref|label)\{[^{}]*\}",
        " REF ",
        txt,
    )
    txt = re.sub(
        r"\\(?:SI|SIrange|num)\*?(?:\[[^\]]*\])?"
        r"\{[^{}]*\}(?:\{[^{}]*\})?",
        " VALUE ",
        txt,
    )
    txt = re.sub(
        r"\\(?:section|subsection|subsubsection)\*?\{([^{}]*)\}",
        r"\n\n\1\n\n",
        txt,
    )
    txt = re.sub(
        r"\\(?:maketitle|newpage|clearpage|vspace|hspace)\*?"
        r"(?:\[[^\]]*\])?(?:\{[^{}]*\})?",
        " ",
        txt,
    )
    # Preserve text arguments of ordinary formatting commands. Repeating the
    # substitution handles simple nesting such as \textbf{\emph{...}}.
    for _ in range(3):
        txt = re.sub(
            r"\\[A-Za-z@]+\*?(?:\[[^\]]*\])?\{([^{}]*)\}",
            r"\1",
            txt,
        )
    txt = re.sub(r"\\[A-Za-z@]+\*?(?:\[[^\]]*\])?", " ", txt)
    txt = re.sub(r"[{}]", " ", txt)
    txt = txt.replace(r"\%", "%").replace("~", " ")
    lines = [
        re.sub(r"[ \t]+", " ", line).strip()
        for line in txt.splitlines()
    ]
    # Keep leading and repeated newlines: source_line is part of the returned
    # review contract and must continue to point at the original TeX file.
    return "\n".join(lines).rstrip()


def _paper_writing_readability_language(text: str) -> str:
    cjk = len(re.findall(r"[\u3040-\u30ff\u3400-\u9fff]", text))
    latin = len(re.findall(r"[A-Za-z]", text))
    english_words = len(re.findall(r"\b[A-Za-z]{2,}\b", text))
    if cjk >= 8 and latin >= 30:
        cjk_share = cjk / max(cjk + latin, 1)
        if 0.12 <= cjk_share <= 0.55 and english_words >= 8:
            return "mixed"
    if cjk >= 8 and cjk >= latin * 0.12:
        return "ja"
    if english_words >= 4:
        return "en"
    return "other"


def _paper_writing_english_sentences(paragraph: str) -> list[str]:
    protected = re.sub(
        r"\b(et al|e\.g|i\.e|cf|vs|Fig|Eq|Ref|Sec|Vol|No|pp|Dr|Prof)\.",
        lambda match: match.group(0).replace(".", "<DOT>"),
        paragraph,
        flags=re.IGNORECASE,
    )
    parts = re.split(
        r"(?<=[.!?])\s+(?=[A-Z\[])|(?<=[。！？])",
        protected,
    )
    return [
        part.replace("<DOT>", ".").strip()
        for part in parts
        if part.strip()
    ]


def _paper_writing_readability_result(
    language: str,
    sentences: list[dict],
    paragraphs: list[dict],
    max_report: int,
) -> dict:
    risks: list[dict] = []

    def _risk(
        kind: str,
        severity: str,
        item: dict,
        reasons: list[str],
        recommendation: str,
    ) -> None:
        risks.append({
            "type": kind,
            "severity": severity,
            "source_line": item["source_line"],
            "paragraph_index": item["paragraph_index"],
            "excerpt": re.sub(r"\s+", " ", item["text"]).strip()[:360],
            "reasons": reasons,
            "recommendation": recommendation,
        })

    if language == "ja":
        connectors = re.compile(
            r"(?:しかし|一方|さらに|また|したがって|そこで|このため|"
            r"すなわち|ところが|なお|ゆえに|ただし|まず|次に)"
        )
        predicates = re.compile(
            r"(?:である|となる|を示す|がある|できる|必要(?:である)?|"
            r"保証(?:される|する)|求める|用いる|行う|課す|得られる|"
            r"依存する|満たす|という|になる)"
        )
        for sentence in sentences:
            compact = re.sub(r"\s+", "", sentence["text"])
            n_chars = len(compact)
            n_commas = len(re.findall(r"[、，,]", sentence["text"]))
            n_connectors = len(connectors.findall(sentence["text"]))
            n_predicates = len(predicates.findall(sentence["text"]))
            terms = set(re.findall(
                r"(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9_.+-]{1,}",
                sentence["text"],
            )) - {"MATH", "REF", "VALUE"}
            reasons = []
            if n_chars > 100:
                reasons.append(f"{n_chars}字（和文目安100字以下）")
            if n_commas >= 4:
                reasons.append(f"読点・コンマ{n_commas}個")
            if n_connectors >= 3:
                reasons.append(f"接続表現{n_connectors}個")
            if n_predicates >= 4:
                reasons.append(f"述語候補{n_predicates}個")
            if len(terms) >= 5:
                reasons.append(f"英字の記号・手法名{len(terms)}個")
            high = n_chars > 180 or (
                n_chars > 140 and (n_commas >= 5 or n_predicates >= 4)
            )
            cognitively_dense = n_chars >= 80 and (
                n_commas >= 4 or n_connectors >= 2 or n_predicates >= 3
            )
            if cognitively_dense and not reasons:
                reasons.append("80字以上で節・論理接続が集中")
            if high or len(reasons) >= 2 or cognitively_dense:
                _risk(
                    "japanese_multi_claim_sentence",
                    "HIGH" if high else "MEDIUM",
                    sentence,
                    reasons,
                    (
                        "一文を『主張』『理由』『定義または条件』へ分け、"
                        "最初の文だけで要点が分かる順序にする。"
                    ),
                )

        role_patterns = {
            "background_or_problem": r"(?:従来|現状|課題|問題|ため)",
            "requirement": r"(?:要求|必要|条件|満たす)",
            "definition": r"(?:とは|すなわち|ここで|をいう)",
            "method": r"(?:本報告|本研究|提案|同定|用いる|構成)",
            "consequence": r"(?:したがって|このため|そこで|ゆえに|可能となる)",
        }
        for paragraph in paragraphs:
            n_chars = len(re.sub(r"\s+", "", paragraph["text"]))
            roles = [
                name for name, pattern in role_patterns.items()
                if re.search(pattern, paragraph["text"])
            ]
            reasons = []
            if n_chars > 500:
                reasons.append(f"段落{n_chars}字（和文目安500字以下）")
            if len(roles) >= 4 and n_chars >= 350:
                reasons.append("一段落に" + "・".join(roles) + "が混在")
            high = n_chars > 800 or (n_chars > 600 and len(roles) >= 4)
            if high or reasons:
                _risk(
                    "japanese_argument_overloaded_paragraph",
                    "HIGH" if high else "MEDIUM",
                    paragraph,
                    reasons,
                    (
                        "段落を一つの問いと一つの結論に限定する。背景、定義、"
                        "選択理由、提案内容は別段落にする。"
                    ),
                )
    else:
        nominal_suffixes = (
            "tion", "sion", "ment", "ance", "ence", "ity", "ness",
            "ship", "ization", "isation",
        )
        clause_words = re.compile(
            r"\b(?:which|that|because|although|whereas|while|thereby|"
            r"however|therefore|moreover|unless|so that)\b",
            re.IGNORECASE,
        )
        for sentence in sentences:
            words = re.findall(r"\b[A-Za-z][A-Za-z'-]*\b", sentence["text"])
            n_words = len(words)
            n_clauses = len(clause_words.findall(sentence["text"]))
            n_breaks = len(re.findall(r"[,;:]|--|[—–]", sentence["text"]))
            n_nominal = sum(
                1 for word in words
                if len(word) > 6 and word.lower().endswith(nominal_suffixes)
            )
            n_jargon = len(re.findall(
                r"\b(?:[A-Z]{2,}[A-Z0-9-]*|[A-Za-z]+-[A-Za-z-]+)\b",
                sentence["text"],
            ))
            reasons = []
            if n_words > 30:
                reasons.append(f"{n_words} words (English target <=30)")
            if n_clauses >= 3:
                reasons.append(f"{n_clauses} subordinate/logical connectors")
            if n_breaks >= 5:
                reasons.append(f"{n_breaks} comma/semicolon/dash clause breaks")
            if n_words and n_nominal / n_words >= 0.12 and n_nominal >= 4:
                reasons.append(f"{n_nominal} nominalisations")
            if n_jargon >= 5:
                reasons.append(f"{n_jargon} acronym/compound terms")
            high = n_words > 50 or (n_words > 42 and n_breaks >= 4)
            if high or len(reasons) >= 2:
                _risk(
                    "english_clause_or_noun_stack",
                    "HIGH" if high else "MEDIUM",
                    sentence,
                    reasons,
                    (
                        "Lead with one claim in an active clause. Move the reason,"
                        " qualification, or list of constraints to a new sentence."
                    ),
                )

        for paragraph in paragraphs:
            n_words = len(re.findall(
                r"\b[A-Za-z][A-Za-z'-]*\b", paragraph["text"]
            ))
            if n_words > 180:
                _risk(
                    "english_overlong_paragraph",
                    "HIGH" if n_words > 250 else "MEDIUM",
                    paragraph,
                    [f"paragraph has {n_words} words (English target <=180)"],
                    (
                        "Split the paragraph at a rhetorical boundary: problem,"
                        " method, evidence, or implication."
                    ),
                )

    high_count = sum(risk["severity"] == "HIGH" for risk in risks)
    medium_count = sum(risk["severity"] == "MEDIUM" for risk in risks)
    score = max(0, 100 - 15 * high_count - 6 * medium_count)
    if high_count or medium_count >= 4:
        status = "fail"
    elif medium_count:
        status = "warning"
    else:
        status = "pass"
    thresholds = (
        {
            "sentence_target_chars": 100,
            "sentence_high_chars": 180,
            "paragraph_target_chars": 500,
            "paragraph_high_chars": 800,
        }
        if language == "ja"
        else {
            "sentence_target_words": 30,
            "sentence_high_words": 50,
            "paragraph_target_words": 180,
            "paragraph_high_words": 250,
        }
    )
    return {
        "applicable": bool(sentences),
        "status": status,
        "heuristic_score": score,
        "score_max": 100,
        "sentence_count": len(sentences),
        "paragraph_count": len(paragraphs),
        "high_risk_count": high_count,
        "medium_risk_count": medium_count,
        "risk_count": len(risks),
        "risks": risks[:max_report],
        "thresholds": thresholds,
    }


def paper_writing_bilingual_readability_check(
    text_or_path: str,
    max_report_per_language: int = 20,
) -> dict:
    """Assess Japanese and English readability with separate criteria.

    The input may be plain prose or a TeX path. Equations, figures, tables,
    keywords, citations, and TeX commands are excluded before analysis. The
    Japanese result uses character, punctuation, predicate, and rhetorical-role
    density. The English result uses word, clause, nominalisation, and paragraph
    density. Scores are diagnostic only and are never averaged across languages.
    """
    raw, source_path = _paper_writing_text_or_path(text_or_path)
    prose = _paper_writing_readability_prose(raw)
    paragraph_matches = list(re.finditer(
        r"(?:^|\n\s*\n)(.*?)(?=\n\s*\n|\Z)",
        prose,
        flags=re.DOTALL,
    ))
    blocks: list[dict] = []
    for index, match in enumerate(paragraph_matches, start=1):
        raw_block = match.group(1)
        leading = re.match(r"\s*", raw_block).group(0)
        text = re.sub(r"\s+", " ", raw_block).strip()
        if len(text) < 20:
            continue
        language = _paper_writing_readability_language(text)
        if language == "other":
            continue
        blocks.append({
            "paragraph_index": index,
            "source_line": prose.count(
                "\n", 0, match.start(1) + len(leading)
            ) + 1,
            "language": language,
            "text": text,
        })

    by_language: dict[str, dict] = {}
    for language in ("ja", "en"):
        language_paragraphs = []
        language_sentences = []
        for block in blocks:
            block_language = block["language"]
            if block_language not in (language, "mixed"):
                continue
            language_paragraphs.append(block)
            if language == "ja":
                parts = [
                    part.strip()
                    for part in re.split(r"(?<=[。！？])", block["text"])
                    if part.strip()
                ]
            else:
                parts = _paper_writing_english_sentences(block["text"])
            for part in parts:
                part_language = _paper_writing_readability_language(part)
                if part_language not in (language, "mixed"):
                    continue
                language_sentences.append({
                    "paragraph_index": block["paragraph_index"],
                    "source_line": block["source_line"],
                    "text": part,
                })
        by_language[language] = _paper_writing_readability_result(
            language,
            language_sentences,
            language_paragraphs,
            max_report_per_language,
        )

    statuses = [
        result["status"]
        for result in by_language.values()
        if result["applicable"]
    ]
    if "fail" in statuses:
        status = "fail"
    elif "warning" in statuses:
        status = "warning"
    elif statuses:
        status = "pass"
    else:
        status = "not_applicable"

    priorities = []
    for language, result in by_language.items():
        for risk in result["risks"]:
            priorities.append({"language": language, **risk})
    priorities.sort(key=lambda item: (
        item["severity"] != "HIGH",
        item["source_line"],
    ))
    return {
        "applicable": bool(statuses),
        "status": status,
        "source_path": source_path,
        "japanese": by_language["ja"],
        "english": by_language["en"],
        "mixed_paragraph_count": sum(
            block["language"] == "mixed" for block in blocks
        ),
        "review_priorities": priorities[:max_report_per_language],
        "aggregation_policy": (
            "Japanese and English scores use different units and are not "
            "averaged. The overall status is the worse applicable language."
        ),
        "target_reader": (
            "a technically literate reviewer in the broad field who does not "
            "already know the authors' notation or argument"
        ),
        "source": "language-specific adjacent-reviewer readability heuristic",
    }


def _paper_writing_word_count(src: str) -> int:
    plain = _paper_writing_plain_text(src)
    words = re.findall(
        r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*|[\u3040-\u30ff\u3400-\u9fff]+",
        plain,
    )
    return len(words)


def _paper_writing_mask_abstract(src: str) -> str:
    return re.sub(
        r"\\begin\{abstract\}.*?\\end\{abstract\}",
        " ",
        src,
        flags=re.DOTALL,
    )


def _paper_writing_near_cite(ctx: str) -> bool:
    return bool(re.search(r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\])?\{", ctx))


def paper_writing_check_digest_human_review_triggers(
    text_or_path: str,
    max_abstract_percent_numbers: int = 1,
    max_caption_words: int = 24,
    max_geometry_families: int = 1,
    max_report: int = 20,
) -> dict:
    """Detect one-page digest issues learned from Sugahara human review.

    This checker is intentionally a warning-level review assistant, not a
    rigid grammar rule.  It catches patterns that repeatedly trigger
    reviewer/coauthor questions in short EM digests:

    - the abstract enumerates per-case exact percent errors instead of an
      aggregate claim such as "sub-percent accuracy";
    - an EM digest abstract says only "Cauer ladder models", which can
      read as a generic equivalent circuit rather than a Cauer Ladder
      Network field-reduction representation;
    - a known Warburg/Randles/CLN combination is framed as the novelty;
    - an empirical CLN--surface transition-frequency issue is described too
      narrowly as something "outside the Galerkin reduction";
    - HOIBC or Warburg is introduced in the body without a nearby citation;
    - a figure caption carries interpretation that should be in the body;
    - opaque shorthand such as rank-(1,1) appears without words;
    - a minimal block/specialization phrase can be misread as "the whole
      benchmark is represented by two basis functions";
    - an important block/Schur equation is left unnumbered, unlabeled, or
      not cited from the body;
    - block matrices such as K_bb/K_sb are shown without (s) dependence
      or without words defining bulk/surface/coupling blocks;
    - N_b is used without defining what is counted and whether the dc term
      is included;
    - the digest explains N_b by saying the uniform dc term is not counted,
      instead of using a simpler "two-rung CLN + SIBC" style count;
    - N_b and N are both used without distinguishing mixed-model bulk
      corrections from bulk-only CLN order;
    - L/R terminations are mentioned without saying they are inductive and
      resistive closures of the last CLN rung;
    - a CLN/HOIBC demonstrator uses N_b=1 without explaining why it is not
      just a dc term plus an impedance boundary condition;
    - a p_H=0 / SIBC0 example is advertised as HOIBC or high-order
      without separating the current benchmark from the higher-order
      extension;
    - Schur complement and DtN/Steklov language are equated in one sentence
      without explaining the elimination/map interpretation;
    - "SIBC scaling" is used without saying whether it means f^{-1/2};
    - CLN basis count or HOIBC order is missing from a graph/model
      description;
    - a main result figure is placed in the source before the
      numerical/verification section that explains it;
    - a numerical example is hidden as a bold "Verification." paragraph
      instead of a numbered section;
    - a figure marker such as f_N appears without saying what it marks;
    - a one-page digest expands from one benchmark to sphere/cube/polyhedra.
    - a numerical example is called proof/evidence beyond its validation scope;
    - an implementation sentence uses "the derived conditions" without
      repeating the method name whose conditions were derived.
    """
    raw_src, source_path = _paper_writing_text_or_path(text_or_path)
    looks_like_tex = bool(source_path or "\\" in raw_src)
    src = _paper_writing_strip_comments(raw_src) if looks_like_tex else raw_src
    findings: list[dict] = []

    def add(rule: str, message: str, advice: str, excerpt: str = "") -> None:
        findings.append({
            "rule": rule,
            "severity": "warning",
            "message": message,
            "advice": advice,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:240],
        })

    abstract = _paper_writing_extract_environment(src, "abstract")
    if not abstract:
        for rec in _paper_writing_iter_command_args(src, "abstract"):
            abstract = rec["text"].strip()
            break

    if abstract:
        percent_re = re.compile(
            r"(?<![A-Za-z])(?:\d+(?:\.\d+)?|"
            r"\d+\s*\\?times\s*10\^\{?[-+]?\d+\}?)\s*(?:\\%|%|percent|per cent)",
            flags=re.IGNORECASE,
        )
        percent_claims = [m.group(0) for m in percent_re.finditer(abstract)]
        if len(percent_claims) > max_abstract_percent_numbers:
            add(
                "abstract_numeric_overload",
                (f"abstract lists {len(percent_claims)} exact percent claims "
                 f"({', '.join(percent_claims[:4])})"),
                ("Use an aggregate expression in the abstract, such as "
                 "'sub-percent accuracy', and move per-case errors to the "
                 "body or figure discussion."),
                abstract,
            )

        generic_cauer_re = re.compile(
            r"\bCauer\s+ladder\s+models?\b|\bCauer\s+models?\b",
            flags=re.IGNORECASE,
        )
        m_generic_cauer = generic_cauer_re.search(abstract)
        if m_generic_cauer and not re.search(
            r"Cauer\s+Ladder\s+Network|electromagnetic\s+fields?|"
            r"field[- ]reduction|field\s+reduced|電磁場|場の縮約",
            abstract,
            flags=re.IGNORECASE,
        ):
            add(
                "abstract_generic_cauer_ladder",
                "abstract uses a generic Cauer-ladder phrase",
                ("If the contribution is an electromagnetic field-reduction "
                 "model, write 'Cauer Ladder Network representations of "
                 "electromagnetic fields' or equivalent wording.  'Cauer "
                 "ladder models' can read as a generic circuit model."),
                abstract,
            )

    body_src = _paper_writing_mask_abstract(src)
    first_section = re.search(r"\\section\*?\{", body_src)
    content_src = body_src[first_section.start():] if first_section else body_src
    full_plain = _paper_writing_plain_text(src)
    section_records = _paper_writing_iter_command_args(content_src, "section")
    numerical_section_re = re.compile(
        r"\b(numerical\s+example|numerical\s+examples|verification|"
        r"validation|results?)\b|数値例|検証|結果",
        flags=re.IGNORECASE,
    )
    numerical_sections = [
        rec for rec in section_records
        if numerical_section_re.search(_paper_writing_plain_text(rec["text"]))
    ]

    overclaim_re = re.compile(
        r"(?:基礎的|決定的|直接的|明白な|強力な)(?:な)?証拠|"
        r"\b(?:conclusive|definitive|strong)\s+evidence\b",
        flags=re.IGNORECASE,
    )
    for m in overclaim_re.finditer(full_plain):
        add(
            "overstated_numerical_evidence",
            "a result is described as evidence stronger than the stated validation scope",
            ("For a finite numerical sweep or single benchmark, prefer "
             "'supports', 'suggests', 'provides a finding', or '知見を与える'. "
             "Reserve 'proof/証拠' for a mathematical derivation or a validation "
             "program whose scope justifies that strength."),
            full_plain[max(0, m.start() - 120):m.end() + 120],
        )

    derived_subject_re = re.compile(
        r"(?:導出した|得られた)[^。．\n]{0,60}(?:条件|式|係数)(?:は|を)"
        r"[^。．\n]{0,120}[。．]|"
        r"\b(?:the\s+)?derived\s+(?:conditions?|equations?|coefficients?)\s+"
        r"(?:is|are|was|were|has|have)\b[^.]{0,120}\.",
        flags=re.IGNORECASE,
    )
    named_method_re = re.compile(
        r"(?<![A-Za-z0-9])[A-Z][A-Z0-9-]{1,}(?![A-Za-z0-9-])|"
        r"[一-龯ぁ-んァ-ン]{2,}(?:法|モデル|定式化)の"
    )
    for m in derived_subject_re.finditer(full_plain):
        sentence = m.group(0)
        if named_method_re.search(sentence):
            continue
        add(
            "ambiguous_derived_method_subject",
            "an implementation statement does not name the method that owns the derived conditions",
            ("Repeat the method name in the same sentence, especially at a paragraph "
             "boundary: for example, 'the derived MMPM moment conditions were "
             "implemented in Fortran' / '導出したMMPMのモーメント条件は...'."),
            sentence,
        )

    over_narrow_galerkin_re = re.compile(
        r"(transition\s+frequency|crossover|connection\s+frequency).{0,80}"
        r"outside.{0,40}Galerkin|"
        r"Galerkin.{0,40}outside.{0,80}"
        r"(transition\s+frequency|crossover|connection\s+frequency)|"
        r"(遷移周波数|接続周波数).{0,80}Galerkin.{0,20}外|"
        r"Galerkin.{0,20}外.{0,80}(遷移周波数|接続周波数)|"
        r"Galerkin\s*射影\s*の\s*外",
        flags=re.IGNORECASE | re.DOTALL,
    )
    m_narrow = over_narrow_galerkin_re.search(src)
    if m_narrow:
        add(
            "over_narrow_galerkin_transition_framing",
            "transition-frequency issue is framed too narrowly as outside Galerkin reduction",
            ("If the real issue is coupling a finite CLN to a Warburg/surface "
             "tail, state that the model is closed by an empirical transition "
             "frequency; do not make the claim sound specific to Galerkin "
             "projection unless that is essential."),
            src[max(0, m_narrow.start() - 120):m_narrow.end() + 120],
        )

    wall_band_re = re.compile(r"\bwall[- ]band\b|壁帯", flags=re.IGNORECASE)
    m_wall = wall_band_re.search(full_plain)
    if m_wall:
        add(
            "unclear_wall_band_term",
            "'wall band' / '壁帯' may be lab-specific or unclear",
            ("Use a more explicit phrase such as 'skin-effect transition "
             "region' or '表皮効果の遷移領域', especially in Japanese digest "
             "text."),
            full_plain[max(0, m_wall.start() - 100):m_wall.end() + 100],
        )

    # HOIBC/Warburg should be cited at first body use.  Abstract citation is
    # forbidden by a separate checker, so this intentionally scans the body.
    citation_terms = [
        ("HOIBC", r"\bHOIBC\b|higher-order\s+impedance\s+boundary\s+condition"),
        ("Warburg", r"\bWarburg\b"),
    ]
    for term, pat in citation_terms:
        m = re.search(pat, content_src, flags=re.IGNORECASE)
        if m:
            ctx = content_src[max(0, m.start() - 180):m.end() + 180]
            if not _paper_writing_near_cite(ctx):
                add(
                    "technical_term_without_citation",
                    f"{term} appears in the body without a nearby citation",
                    ("Define specialized boundary/circuit terminology at first "
                     "body use and cite the source there.  Do not put the "
                     "citation in the abstract."),
                    ctx,
                )

    # Known constructs should be acknowledged as known, then the paper should
    # shift novelty to the new coupling/derivation.
    for m in re.finditer(r"\bWarburg\b", content_src, flags=re.IGNORECASE):
        ctx = content_src[max(0, m.start() - 240):m.end() + 240]
        has_novelty = re.search(
            r"\b(propose|proposed|new|novel|first|develop|developed|"
            r"introduce|introduced)\b|新規|新しい|提案|初めて",
            ctx,
            flags=re.IGNORECASE,
        )
        has_known = re.search(
            r"\b(known|classical|classic|established|conventional|prior)\b|"
            r"既知|古典|従来",
            ctx,
            flags=re.IGNORECASE,
        )
        if has_novelty and not has_known:
            add(
                "known_construct_framing",
                "Warburg-related construction may be framed as the novelty",
                ("If Warburg/CLN combination is already known, state that it is "
                 "known and locate the novelty in the parameter-free Galerkin, "
                 "Schur-complement, or HOIBC coupling."),
                ctx,
            )
            break

    rank_re = re.compile(
        r"rank\s*[-–]?\s*(?:\$[^$]*\(\s*\d+\s*,\s*\d+\s*\)[^$]*\$|"
        r"\(?\s*\d+\s*,\s*\d+\s*\)?)",
        flags=re.IGNORECASE,
    )
    for m in rank_re.finditer(content_src):
        add(
            "opaque_rank_tuple",
            "rank-(m,n)-style notation appears without physical words",
            ("Spell out the meaning for digest readers, e.g. 'one bulk mode "
             "and one surface mode', then add the compact rank notation only "
             "if needed."),
            content_src[max(0, m.start() - 100):m.end() + 100],
        )
        break

    misleading_basis_re = re.compile(
        r"one[- ]bulk\s*/\s*one[- ]surface|"
        r"one\s+bulk\s+mode\s+and\s+one\s+surface\s+mode|"
        r"1\s*bulk\s+mode\s*\+\s*1\s*surface\s+mode|"
        r"1\s*体積モード\s*\+\s*1\s*表面モード|"
        r"基底(?:関数)?\s*(?:が|は)?\s*2\s*個",
        flags=re.IGNORECASE,
    )
    for m in misleading_basis_re.finditer(content_src):
        add(
            "misleading_minimal_basis_claim",
            "minimal bulk/surface phrase may imply the whole benchmark uses only two basis functions",
            ("Distinguish the finite bulk CLN basis from the added surface "
             "envelope/block.  Avoid wording that sounds as if a circular "
             "conductor is represented by only two spatial basis functions."),
            content_src[max(0, m.start() - 120):m.end() + 120],
        )
        break

    body_plain = _paper_writing_plain_text(content_src)
    for sent in re.split(r"(?<=[.。])\s+", body_plain):
        if (re.search(r"Schur|Schur\s+補", sent, flags=re.IGNORECASE)
                and re.search(r"Dirichlet-to-Neumann|DtN|Steklov",
                              sent, flags=re.IGNORECASE)
                and not re.search(
                    r"eliminat|Dirichlet\s+data|Neumann\s+flux|"
                    r"surface\s+Dirichlet|surface\s+potential|"
                    r"消去|Dirichlet\s*データ|Neumann\s*フラックス|"
                    r"表面.*フラックス|写す|写し",
                    sent,
                    flags=re.IGNORECASE,
                )):
            add(
                "unexplained_schur_dtn_equivalence",
                "Schur complement is equated with a DtN/Steklov map without explanation",
                ("Separate the algebraic operation from the operator "
                 "interpretation.  First say that eliminating the bulk "
                 "unknowns gives a Schur complement, then explain that the "
                 "result maps surface Dirichlet data to Neumann flux."),
                sent,
            )
            break

    sibc_scaling_re = re.compile(
        r"SIBC\s+scaling|SIBC\s*スケーリング",
        flags=re.IGNORECASE,
    )
    for m in sibc_scaling_re.finditer(content_src):
        ctx = content_src[max(0, m.start() - 120):m.end() + 120]
        if not re.search(
            r"f\s*\^\s*\{?\s*-?\s*1\s*/\s*2\s*\}?|"
            r"s\s*\^\s*\{?\s*-?\s*1\s*/\s*2\s*\}?|"
            r"square[- ]root|fractional|分数|平方根",
            ctx,
            flags=re.IGNORECASE,
        ):
            add(
                "vague_sibc_scaling",
                "SIBC scaling is mentioned without the actual exponent or tail",
                ("If 'SIBC scaling' means the surface-impedance admittance "
                 "tail, write it explicitly as f^{-1/2}, s^{-1/2}, "
                 "square-root, or fractional scaling."),
                ctx,
            )
            break

    zero_order_surface_re = re.compile(
        r"p_\{(?:\\mathrm\s*\{?H\}?|H)\}\s*=\s*0|"
        r"p_H\s*=\s*0|SIBC0|HOIBC0|zeroth[- ]order|0\s*th[- ]order|"
        r"先頭.{0,20}SIBC|0\s*次",
        flags=re.IGNORECASE,
    )
    higher_order_surface_re = re.compile(
        r"\bHOIBC\b|higher[- ]order\s+impedance|高次(?:インピーダンス|.*境界|.*表面)",
        flags=re.IGNORECASE,
    )
    if (zero_order_surface_re.search(content_src)
            and higher_order_surface_re.search(content_src)):
        separated_as_extension = re.search(
            r"higher[- ]order.{0,120}(can\s+be\s+included|can\s+be\s+needed|"
            r"needed\s+for|extension|curved|three[- ]dimensional|3D)|"
            r"benchmark.{0,120}(uses\s+only\s+the\s+leading|leading\s+SIBC|"
            r"p_\{(?:\\mathrm\s*\{?H\}?|H)\}\s*=\s*0|p_H\s*=\s*0)|"
            r"本(?:ダイジェスト|稿|例).{0,120}(先頭|p_\{(?:\\mathrm\s*\{?H\}?|H)\}\s*=\s*0|p_H\s*=\s*0)|"
            r"高次.{0,120}(曲面|3\s*次元|必要になる場合|含められる)",
            content_src,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not separated_as_extension:
            add(
                "ph_zero_higher_order_mismatch",
                "p_H=0 / SIBC0 example is framed as HOIBC or high-order",
                ("Do not sell a p_H=0 benchmark as a higher-order impedance "
                 "boundary result.  Call the current example the leading "
                 "SIBC term, and mention HOIBC only as an extension or as "
                 "needed for curved/three-dimensional applications."),
                content_src,
            )

    important_block_re = re.compile(
        r"K_\{(?:bb|bs|sb|ss)\}|Schur\s+complement|Schur\s+補",
        flags=re.IGNORECASE,
    )
    unnumbered_math_re = re.compile(
        r"\\\[(?P<bracket>.*?)\\\]|"
        r"\\begin\{(?P<env>equation\*|align\*|gather\*|multline\*)\}"
        r"(?P<starbody>.*?)\\end\{(?P=env)\}|"
        r"(?<!\\)\$\$(?P<dollar>.*?)(?<!\\)\$\$",
        flags=re.DOTALL,
    )
    for m in unnumbered_math_re.finditer(content_src):
        body = m.group("bracket") or m.group("starbody") or m.group("dollar") or ""
        if important_block_re.search(body):
            add(
                "unnumbered_core_equation",
                "important block/Schur equation is displayed without a number",
                ("Use a numbered equation environment with a label for the "
                 "main mixed-system or Schur-complement formula, then cite it "
                 "from the surrounding prose.  Unnumbered display math is fine "
                 "for side remarks, not for the core method definition."),
                body,
            )
            break

    numbered_eq_re = re.compile(
        r"\\begin\{(?P<env>equation|align|gather|multline)\}"
        r"(?P<body>.*?)\\end\{(?P=env)\}",
        flags=re.DOTALL,
    )
    for m in numbered_eq_re.finditer(content_src):
        body = m.group("body")
        if not important_block_re.search(body):
            continue
        labels = re.findall(r"\\label\{([^}]+)\}", body)
        if not labels:
            add(
                "unlabeled_core_equation",
                "important block/Schur equation has a number but no label",
                ("Add a \\label{eq:...} to the numbered core equation so the "
                 "body can refer to it explicitly."),
                body,
            )
            break
        cited = any(
            re.search(
                r"\\(?:eqref|ref|autoref|cref|Cref)\{"
                + re.escape(label)
                + r"\}",
                content_src,
            )
            for label in labels
        )
        if not cited:
            add(
                "uncited_core_equation",
                "important block/Schur equation is labeled but not cited in the body",
                ("After numbering the equation, cite it in the sentence that "
                 "defines the Schur complement or explains its role."),
                body,
            )
        break

    block_symbols = re.findall(r"K_\{(?:bb|bs|sb|ss)\}", content_src)
    if block_symbols:
        block_without_s = re.search(
            r"K_\{(?:bb|bs|sb|ss)\}(?!\s*\(\s*s\s*\))",
            content_src,
        )
        if block_without_s:
            add(
                "block_matrix_missing_s_dependence",
                "block matrix symbols appear without explicit (s) dependence",
                ("For frequency-domain reduced models, write K_bb(s), "
                 "K_bs(s), K_sb(s), and K_ss(s) unless the blocks are "
                 "truly frequency independent."),
                content_src[max(0, block_without_s.start() - 120):block_without_s.end() + 120],
            )
        has_block_words = re.search(
            r"(bulk|体積).{0,120}(surface|表面).{0,120}(coupling|結合)|"
            r"(surface|表面).{0,120}(bulk|体積).{0,120}(coupling|結合)",
            body_plain,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not has_block_words:
            add(
                "undefined_block_matrices",
                "block matrices such as K_bb/K_sb are not defined in words",
                ("Immediately before or after the equation, state which blocks "
                 "are bulk, surface, and coupling Galerkin blocks."),
                content_src,
            )

    nb_re = re.compile(r"N_\{?b\}?|N_b")
    if nb_re.search(content_src):
        m_dc_count = re.search(
            r"(uniform\s+)?d[.]?c[.]?\s+term.{0,100}(not\s+included|"
            r"not\s+counted|excluded)|"
            r"(not\s+included|not\s+counted|excluded).{0,100}"
            r"(uniform\s+)?d[.]?c[.]?\s+term|"
            r"(一様な)?直流項.{0,80}(数えず|含めない|含まない)",
            content_src,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if m_dc_count:
            add(
                "confusing_dc_count_exclusion",
                "model-order count is explained by excluding the dc term",
                ("In a short digest, avoid making the reader parse whether "
                 "the dc field is counted.  Prefer a direct physical count "
                 "such as 'two-rung CLN + SIBC0', or define N_b in one "
                 "simple sentence without negative exclusions."),
                content_src[max(0, m_dc_count.start() - 120):m_dc_count.end() + 120],
            )

        has_nb_definition = re.search(
            r"(N_\{?b\}?|N_b).{0,120}"
            r"(count|denote|number|個数|表す).{0,120}"
            r"(bulk|CLN|Krylov|basis|correction|体積|基底|補正)",
            content_src,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not has_nb_definition:
            add(
                "ambiguous_nb_definition",
                "N_b is used without defining what it counts",
                ("Define N_b in words, and say whether the uniform dc term is "
                 "included.  For example: 'N_b counts zero-boundary bulk "
                 "Krylov correction functions; the uniform dc term is not "
                 "included.'"),
                content_src,
            )

    if re.search(
        r"\b[LR][- ]term\b|\b[LR][- ]terminated\b|[LR]\s*終端",
        content_src,
        flags=re.IGNORECASE,
    ):
        has_l_def = re.search(
            r"L[- ](?:term|terminated).{0,120}induct|induct.{0,120}"
            r"L[- ](?:term|terminated)|L\s*終端.{0,80}誘導|誘導.{0,80}L\s*終端",
            content_src,
            flags=re.IGNORECASE | re.DOTALL,
        )
        has_r_def = re.search(
            r"R[- ](?:term|terminated).{0,120}resist|resist.{0,120}"
            r"R[- ](?:term|terminated)|R\s*終端.{0,80}抵抗|抵抗.{0,80}R\s*終端",
            content_src,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not (has_l_def and has_r_def):
            add(
                "undefined_lr_termination",
                "L/R CLN terminations are mentioned without defining them",
                ("Define L- and R-terminated CLNs in words: the L termination "
                 "closes the last rung inductively, and the R termination "
                 "closes it resistively.  Also state whether N is a rung, "
                 "mode, or basis count."),
                content_src,
            )

    if nb_re.search(content_src):
        has_plain_n_order = bool(re.search(
            r"(?<![A-Za-z_])\$?N\s*=\s*\d+",
            content_src,
        ))
        if has_plain_n_order:
            order_sentences = re.split(r"(?<=[.。])\s+", content_src)
            has_order_distinction = any(
                re.search(r"N_\{?b\}?|N_b", sent)
                and re.search(r"(?<![A-Za-z_])\$?N\s*=\s*\d+", sent)
                and re.search(
                    r"whereas|while|distinguish|distinct|not included|"
                    r"bulk-only|mixed model|order|ladder|counts|"
                    r"区別|混合モデル|体積のみ|次数|梯子|数える",
                    sent,
                    flags=re.IGNORECASE,
                )
                for sent in order_sentences
            )
            if not has_order_distinction:
                add(
                    "ambiguous_n_vs_nb_order",
                    "N_b and N appear without distinguishing their roles",
                    ("When both N_b and N are used, define the distinction: "
                     "N_b counts bulk correction functions in the mixed model, "
                     "whereas N is the order of the bulk-only CLN ladder."),
                    content_src,
                )
        if (re.search(r"\bCLN\b|Cauer", content_src, flags=re.IGNORECASE)
                and re.search(r"\bHOIBC\b|higher-order\s+impedance|高次.*境界",
                              content_src, flags=re.IGNORECASE)
                and re.search(r"N_\{?b\}?\s*=\s*1|N_b\s*=\s*1",
                              content_src)):
            add(
                "minimal_bulk_order_misread",
                "CLN/HOIBC demonstrator uses N_b=1, which can read as dc plus IBC only",
                ("Use at least two bulk correction functions for a mixed "
                 "CLN/HOIBC demonstration, or explicitly explain why N_b=1 "
                 "is sufficient and how it differs from the dc term plus an "
                 "impedance boundary condition."),
                content_src,
            )

    if (re.search(r"\bCLN\b|Cauer", content_src, flags=re.IGNORECASE)
            and re.search(r"\bHOIBC\b|higher-order\s+impedance|高次.*境界",
                          content_src, flags=re.IGNORECASE)):
        has_cln_order = bool(re.search(
            r"\bN\s*=\s*\d+|N_b\s*=\s*\d+|N_\{?b\}?\s*=\s*\d+|"
            r"CLN\s*(?:basis|bases|N)\s*(?:=|of)?\s*\d+|"
            r"\d+\s+bulk\s+CLN|"
            r"CLN\s*Krylov\s+basis\s+function.*?\d+|"
            r"CLN.*?基底関数.*?\d+|基底関数.*?\d+",
            content_src,
            flags=re.IGNORECASE,
        ))
        has_hoibc_order = bool(re.search(
            r"HOIBC\s*(?:order\s*)?\d|HOIBC0|HOIBC1|HOIBC2|"
            r"zeroth[- ]order|first[- ]order|second[- ]order|"
            r"p_\{(?:\\mathrm\s*\{?H\}?|H)\}\s*=\s*\d+|"
            r"p_\{?\\mathrm\{?H\}?\}?\s*=\s*\d+|p_H\s*=\s*\d+|"
            r"\d+\s*次.*HOIBC|HOIBC.*?\d+\s*次",
            content_src,
            flags=re.IGNORECASE,
        ))
        if not (has_cln_order and has_hoibc_order):
            add(
                "missing_model_order_disclosure",
                "mixed CLN/HOIBC description lacks CLN basis count and/or HOIBC order",
                ("State the model orders in the body, caption, or legend: "
                 "for example, 'bulk CLN N_b=2 with zeroth-order HOIBC "
                 "(p_H=0); bulk-only CLN references use N=10'."),
                content_src,
            )

    hidden_verification_re = re.compile(
        r"\\textbf\{\s*(Verification|Validation|Numerical\s+Example)\.?\s*\}|"
        r"\\paragraph\{\s*(Verification|Validation|Numerical\s+Example)\.?\s*\}|"
        r"\\textbf\{\s*(検証|数値例)\s*[:：.]?\s*\}",
        flags=re.IGNORECASE,
    )
    m_hidden = hidden_verification_re.search(content_src)
    if m_hidden and not numerical_sections:
        add(
            "hidden_numerical_example_section",
            "numerical example or verification is hidden in a bold paragraph",
            ("Use a numbered section such as 'Numerical Example' or "
             "'Verification' so the digest has an explicit result section.  "
             "Do not bury the main example under a bold run-in label."),
            content_src[max(0, m_hidden.start() - 160):m_hidden.end() + 160],
        )

    body_plain = _paper_writing_plain_text(content_src)
    if re.search(r"\b(circular|cylindrical|circle|cylinder)\b|円形|円柱", body_plain,
                 flags=re.IGNORECASE):
        has_radius = bool(re.search(r"\bradius\b|半径", body_plain, flags=re.IGNORECASE))
        has_reference = bool(re.search(r"\bBessel\b|ベッセル", body_plain,
                                       flags=re.IGNORECASE))
        has_conductivity = bool(re.search(
            r"\bconductivity\b|\\sigma\s*=|sigma\s*=|\bS/m\b|導電率",
            body_plain,
            flags=re.IGNORECASE,
        ))
        if not (has_radius and has_reference and has_conductivity):
            add(
                "underdescribed_circular_benchmark",
                "circular-conductor benchmark lacks radius, conductivity, and/or reference solution",
                ("For a one-page digest, put the problem setting in the body: "
                 "geometry, radius, conductivity/permeability when material "
                 "properties affect the result, excitation, and the reference "
                 "solution used for verification."),
                body_plain,
            )

    geometry_patterns = {
        "circular/cylindrical": r"\b(circular|cylindrical|circle|cylinder)\b|円形|円柱",
        "sphere": r"\b(sphere|spherical)\b|球",
        "cube": r"\b(cube|cuboid|rectangular solid)\b|立方体|直方体",
        "polyhedra": r"\b(polyhedra|polyhedral)\b|多面体",
    }
    geometry_hits = [
        name for name, pat in geometry_patterns.items()
        if re.search(pat, body_plain, flags=re.IGNORECASE)
    ]
    if len(geometry_hits) > max_geometry_families:
        add(
            "one_page_scope_creep_shapes",
            f"multiple geometry families appear in a short digest: {geometry_hits}",
            ("Keep the one-page digest on the central benchmark.  Move sphere, "
             "cube, and polyhedral extensions to the full paper unless they are "
             "the main result."),
            body_plain,
        )

    caption_records = _paper_writing_iter_command_args(content_src, "caption")
    first_figure = re.search(r"\\begin\{figure\*?\}", content_src)
    if first_figure:
        first_numerical_start = (
            min(rec["start"] for rec in numerical_sections)
            if numerical_sections else None
        )
        first_caption_text = caption_records[0]["text"] if caption_records else ""
        result_caption = re.search(
            r"\b(admittance|error|comparison|result|reference|accuracy)\b|"
            r"アドミタンス|誤差|比較|結果|参照|精度",
            first_caption_text,
            flags=re.IGNORECASE,
        )
        if result_caption and (
            first_numerical_start is None or first_figure.start() < first_numerical_start
        ):
            add(
                "result_figure_before_numerical_section",
                "main result figure is placed before the numerical/verification section",
                ("Use the standard figure environment where possible, but put "
                 "the source block in or after the section that explains the "
                 "numerical example.  Manual non-float captions, [H], "
                 "\\clearpage, or \\newpage placement are last resorts; first "
                 "adjust source position, text length, figure width, and "
                 "caption length."),
                first_caption_text,
            )

    fn_marker_re = re.compile(
        r"\$?\s*f\s*_\s*\{?N\}?\s*\$?|\\\(.*?f\s*_\s*\{?N\}?.*?\\\)",
        flags=re.IGNORECASE | re.DOTALL,
    )
    if fn_marker_re.search(content_src):
        has_fn_explanation = re.search(
            r"(highest|last|cutoff|pole\s+frequency|CLN\s+pole|"
            r"highest\s+pole|last\s+pole).{0,80}f\s*_\s*\{?N\}?|"
            r"f\s*_\s*\{?N\}?.{0,80}(highest|last|cutoff|pole\s+frequency|"
            r"CLN\s+pole|highest\s+pole|last\s+pole)|"
            r"最高極|最終極|遮断|極周波数|を示す|を表す",
            content_src,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not has_fn_explanation:
            add(
                "ambiguous_fn_marker",
                "f_N marker appears without explaining what it denotes",
                ("If a vertical line or label such as $f_N$ is shown, state "
                 "in the caption or nearby text what it marks, e.g. 'the "
                 "highest pole frequency of the bulk CLN'."),
                content_src,
            )

    caption_telling_re = re.compile(
        r"\b(whereas|while|because|therefore|hence|compared with|compares|"
        r"brackets?|miss(?:es)?|follows|agrees?|demonstrates|shows that|"
        r"error|accuracy)\b",
        flags=re.IGNORECASE,
    )
    for i, cap in enumerate(caption_records, start=1):
        n_words = _paper_writing_word_count(cap["text"])
        if n_words > max_caption_words or caption_telling_re.search(cap["text"]):
            add(
                "caption_overloaded",
                f"caption {i} carries {n_words} words or interpretive phrasing",
                ("For a one-page digest, use the caption to identify the plotted "
                 "quantity and minimal conditions.  Any factual claim, numerical "
                 "result, comparison logic, or interpretation written in the "
                 "caption must also be stated in the body.  If the caption is "
                 "long, prioritize the body and shorten the caption."),
                cap["text"],
            )

    status = "warning" if findings else "clean"
    return {
        "status": status,
        "n_findings": len(findings),
        "findings": findings[:max_report],
        "truncated": len(findings) > max_report,
        "source_path": source_path,
        "checks": {
            "max_abstract_percent_numbers": max_abstract_percent_numbers,
            "max_caption_words": max_caption_words,
            "max_geometry_families": max_geometry_families,
        },
        "summary": (
            "No digest-specific human-review triggers found."
            if not findings else
            f"{len(findings)} digest-specific human-review trigger(s) found."
        ),
        "source": "Sugahara human review, IGTE Cauer digest, 2026-06-24",
    }


def paper_writing_check_strong_adjective_budget(text: str) -> dict:
    """強調副詞/形容詞の過剰使用を検出。

    Wallwork §8.12: 2 個超で "hype" 判定、reader が信頼感を失う。

    目標: 全体で 2 個以内 (Abstract + Intro で 1-2 個 max)
    警告: 5 個超
    """
    patterns = {
        "remarkably": r"\bremarkably\b",
        "notably": r"\bnotably\b",
        "importantly": r"\bimportantly\b",
        "interestingly": r"\binterestingly\b",
        "novel": r"\bnovel\b",
        "innovative": r"\binnovative\b",
        "cutting-edge": r"\bcutting[- ]edge\b",
        "state-of-the-art": r"\bstate[- ]of[- ]the[- ]art\b",
        "crucial": r"\bcrucial\b",
        "groundbreaking": r"\bgroundbreaking\b",
        "significantly": r"\bsignificantly\b",
        "substantially": r"\bsubstantially\b",
    }
    by_pat = {}
    total = 0
    for name, pat in patterns.items():
        n = len(re.findall(pat, text, re.IGNORECASE))
        if n:
            by_pat[name] = n
            total += n
    return {
        "total_strong_adjectives": total,
        "by_adjective": by_pat,
        "target": "<=2 overall (Wallwork §8.12)",
        "warning": ">=5 overall",
        "note": (
            "Strong modifiers above threshold signal hype. Replace with "
            "quantitative claims: 'substantially faster' -> '3.2x faster "
            "(Table II)'."
        ),
    }


# ------------------------------------------------------------------
# v0.4.0: 即戦力 5 tools (IMRAD / tense / self-cite / forward-ref / eq-numbering)
# ------------------------------------------------------------------

_IMRAD_TARGETS = {
    "Abstract":    (0.04, 0.08),   # 4-8% (target ~5%)
    "Introduction":(0.12, 0.20),   # 12-20% (target ~15%)
    "Method":      (0.20, 0.30),   # 20-30% (target ~25%)
    "Result":      (0.25, 0.35),   # 25-35% (target ~30%)
    "Discussion":  (0.15, 0.25),   # 15-25% (target ~20%)
    "Conclusion":  (0.03, 0.08),   # 3-8% (target ~5%)
}


def paper_writing_check_imrad_balance(text: str) -> dict:
    """IMRAD 各セクションの字数バランスを検証する。

    Wallwork §14 の推奨: Intro 15% / Method 25% / Result 30% / Discussion 20%
    / Conclusion 5% / Abstract 5%。Japanese section headers (序論/方法/結果/
    考察/結論/Abstract) も English (Introduction/Methods/Results/Discussion/
    Conclusion/Abstract) も検出。

    Args:
        text: tex/md 本文全体。section header は `\\section{...}`, `## ...`,
              or bare 「序論」「1. 序論」等を検出。

    Returns:
        dict: section_lengths (chars/words), ratios, deviations from target.
    """
    # Section boundary patterns (bilingual)
    section_markers = {
        "Abstract":     r"(?:\\section\*?\{.*?(?:Abstract|要旨|概要)\}|^#+\s*(?:Abstract|要旨|概要))",
        "Introduction": r"(?:\\section\*?\{.*?(?:Introduction|序論|はじめに|緒言)\}|^#+\s*(?:Introduction|序論|はじめに))",
        "Method":       r"(?:\\section\*?\{.*?(?:Method|Methods|手法|方法|実験方法)\}|^#+\s*(?:Method|手法|方法))",
        "Result":       r"(?:\\section\*?\{.*?(?:Result|Results|結果|実験結果)\}|^#+\s*(?:Result|結果))",
        "Discussion":   r"(?:\\section\*?\{.*?(?:Discussion|考察|議論)\}|^#+\s*(?:Discussion|考察))",
        "Conclusion":   r"(?:\\section\*?\{.*?(?:Conclusion|結論|結言|まとめ)\}|^#+\s*(?:Conclusion|結論|まとめ))",
    }
    # Find positions
    positions = {}
    for name, pat in section_markers.items():
        m = re.search(pat, text, re.MULTILINE | re.IGNORECASE)
        if m:
            positions[name] = m.start()
    if len(positions) < 3:
        return {
            "error": "Could not identify IMRAD sections (need >=3)",
            "sections_found": list(positions.keys()),
            "hint": "Use standard section headers (\\section{...} or ## ...)",
        }

    # Sort by position
    ordered = sorted(positions.items(), key=lambda kv: kv[1])
    section_lengths = {}
    for i, (name, start) in enumerate(ordered):
        end = ordered[i+1][1] if i + 1 < len(ordered) else len(text)
        section_lengths[name] = end - start

    total = sum(section_lengths.values())
    ratios = {name: (ln / total) for name, ln in section_lengths.items()}

    deviations = []
    for name, ratio in ratios.items():
        lo, hi = _IMRAD_TARGETS.get(name, (0.0, 1.0))
        if ratio < lo:
            deviations.append({"section": name, "ratio": round(ratio, 3),
                               "target_range": [lo, hi],
                               "verdict": f"{name} too short ({ratio:.1%} < {lo:.0%})"})
        elif ratio > hi:
            deviations.append({"section": name, "ratio": round(ratio, 3),
                               "target_range": [lo, hi],
                               "verdict": f"{name} too long ({ratio:.1%} > {hi:.0%})"})
    return {
        "sections_found": [name for name, _ in ordered],
        "section_lengths_chars": section_lengths,
        "section_ratios": {k: round(v, 3) for k, v in ratios.items()},
        "deviations": deviations,
        "target": "Intro 15%, Method 25%, Result 30%, Discussion 20%, Concl 5%",
        "source": "Wallwork §14 + IMRAD conventions",
    }


def paper_writing_check_tense_consistency(text: str,
                                            section: str = "Discussion") -> dict:
    """Discussion の 3-part tense (hypothesis=現在, result=過去, background=現在完了) 混在検出。

    佐藤 Q32 + Wallwork §14.8 の標準:
      - Hypothesis: present tense ("X *is* a stimulant")
      - Your result: past tense ("the distance *was* associated")
      - Established background: present perfect ("X *has been shown* to ...")

    Args:
        text: 検査対象の英語 discussion 本文
        section: "Discussion" (default), "Method", "Result" 等。Methodの場合は
                 過去形チェックのみ。

    Returns:
        dict: mismatched sentences count + examples.
    """
    # English sentences
    sents = [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]
    if not sents:
        return {"error": "no sentences found"}

    issues = []
    # Hypothesis markers (should be present)
    hypothesis_pat = re.compile(r"\bhypothesi[sz]\w*\b|\b(?:our|the) hypothesis\b", re.I)
    # Result reference markers (Fig. N, Table N, we found/observed) → past
    result_pat = re.compile(r"\b(?:Fig\.?|Figure|Table)\s*\d+|\bwe (?:found|observed|measured|showed|demonstrated)\b", re.I)
    # Background markers
    bg_pat = re.compile(r"\b(?:has been|have been|it is (?:well[- ]known|known))\b", re.I)

    past_tense_pat = re.compile(r"\b(?:was|were|had|had been)\b", re.I)
    present_tense_pat = re.compile(r"\b(?:is|are|has|have)\b(?! been)", re.I)

    for i, s in enumerate(sents):
        has_hyp = bool(hypothesis_pat.search(s))
        has_res = bool(result_pat.search(s))
        has_past = bool(past_tense_pat.search(s))
        has_pres = bool(present_tense_pat.search(s))

        if has_hyp and has_past and not has_pres:
            issues.append({
                "index": i, "type": "hypothesis_in_past",
                "sentence": s[:120],
                "hint": "Hypothesis should be present tense.",
            })
        elif has_res and has_pres and not has_past:
            issues.append({
                "index": i, "type": "result_in_present",
                "sentence": s[:120],
                "hint": "Your specific result should be past tense.",
            })
    return {
        "section": section,
        "total_sentences": len(sents),
        "issues": issues[:15],
        "issue_count": len(issues),
        "source": "佐藤 Q32 + Wallwork §14.8",
    }


def paper_writing_check_self_citation_ratio(tex_path: str,
                                              bib_path: str = "",
                                              author_last_names: str = "") -> dict:
    """自己引用率を算出。Wallwork: <20% 推奨。

    Args:
        tex_path: .tex ファイル
        bib_path: .bib ファイル (省略時は同ディレクトリの *.bib)
        author_last_names: 自著と判定する著者の姓 (カンマ区切り、
                            例: "Sugahara,Hane,Hollaus")

    Returns:
        dict: total_citations, self_citations, ratio, flagged.
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"tex file not found: {tex_path}"}
    tex = p.read_text(encoding="utf-8", errors="replace")
    if bib_path:
        bibs = [pathlib.Path(bib_path)]
    else:
        bibs = list(p.parent.glob("*.bib"))
    if not bibs:
        return {"error": "no .bib file found"}

    authors = {a.strip() for a in author_last_names.split(",") if a.strip()}
    if not authors:
        return {"error": "provide author_last_names (comma-separated)"}

    # Parse bib entries with a brace-balanced scanner.
    # Naive regex like r"@\w+\s*\{([^@]*)" breaks on stray '@' inside fields
    # (e.g. url = {...@doi.org/...}, email = {foo@bar}), and on nested braces
    # inside author fields (e.g. {Sugahara, K. and Hane, {M.}}).
    def _extract_balanced(text: str, start: int) -> tuple[int, str]:
        """Given text[start] == '{', return (end_index_after_closing_brace, inner).
        Tracks brace depth, ignoring braces in comments is not needed for .bib.
        Returns (-1, "") if no matching close brace is found."""
        assert text[start] == "{"
        depth = 0
        i = start
        n = len(text)
        while i < n:
            c = text[i]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    return i + 1, text[start + 1 : i]
            i += 1
        return -1, ""

    bib_entries = {}
    for bp in bibs:
        btxt = bp.read_text(encoding="utf-8", errors="replace")
        # Find each @type{ header.
        for hdr in re.finditer(r"@\w+\s*\{", btxt):
            brace_start = hdr.end() - 1  # position of '{'
            entry_end, entry_body = _extract_balanced(btxt, brace_start)
            if entry_end < 0:
                continue
            # The entry body starts with:  <whitespace><cite_key>,<fields...>
            key_match = re.match(r"\s*([^,\s]+)\s*,", entry_body)
            if not key_match:
                continue
            key = key_match.group(1)
            fields_region = entry_body[key_match.end():]
            # Find author = { ... } with brace balancing (value may contain
            # nested braces like {Sugahara, K. and Hane, {M.}}).
            author_value = ""
            for am in re.finditer(r"author\s*=\s*", fields_region, re.IGNORECASE):
                vstart = am.end()
                if vstart >= len(fields_region):
                    break
                vc = fields_region[vstart]
                if vc == "{":
                    _, inner = _extract_balanced(fields_region, vstart)
                    author_value = inner
                    break
                if vc == '"':
                    # Quoted value: up to the next unescaped ".
                    end_q = fields_region.find('"', vstart + 1)
                    if end_q > 0:
                        author_value = fields_region[vstart + 1 : end_q]
                    break
                # Bare value (unusual): up to comma or end of entry.
                end_bare = re.search(r"[,\n]", fields_region[vstart:])
                author_value = (
                    fields_region[vstart : vstart + end_bare.start()]
                    if end_bare
                    else fields_region[vstart:]
                )
                break
            bib_entries[key] = author_value

    # Find all \cite{...} keys
    cite_keys = set()
    for m in re.finditer(r"\\(?:cite|citep|citet|parencite|autocite)\*?(?:\[[^\]]*\])?\{([^}]+)\}", tex):
        for k in m.group(1).split(","):
            if k.strip():
                cite_keys.add(k.strip())

    total = len(cite_keys)
    self_cites = []
    for key in cite_keys:
        author_field = bib_entries.get(key, "")
        if any(a in author_field for a in authors):
            self_cites.append(key)

    ratio = len(self_cites) / max(1, total)
    return {
        "total_citations": total,
        "self_citation_count": len(self_cites),
        "self_citation_ratio": round(ratio, 3),
        "self_citation_keys": self_cites[:20],
        "authors_checked": sorted(authors),
        "target": "<=0.20 (Wallwork §)",
        "warning": ">=0.30",
        "verdict": "OK" if ratio <= 0.20 else "WARN" if ratio < 0.30 else "HIGH",
    }


def paper_writing_check_figure_forward_reference(tex_path: str) -> dict:
    """図/表の \\label と \\ref の整合チェック (孤立ラベル / 未解決参照)。

    NOTE: ソース上での順序 (\\ref が \\label より前にあるか) は LaTeX の
    forward-reference 違反ではない。\\begin{figure}[htbp] は float なので、
    組版後の物理ページ位置はソース順序と無関係に決まる。むしろ
    「本文で先に導入 → \\begin{figure} は近くの任意位置」が推奨パターン。
    ソース順序の判定が必要なら PDF レベル (\\pageref や pdftotext) で行う。

    この関数は次の2点のみを検査する:
    - orphan_labels: \\label{fig:...}/\\label{tab:...} あるが \\ref されない
    - dangling_refs: \\ref{fig:...} あるが該当 \\label が存在しない

    Args:
        tex_path: .tex ファイル

    Returns:
        dict: labels_defined, refs_used, orphan_labels, dangling_refs.
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"tex file not found: {tex_path}"}
    tex = p.read_text(encoding="utf-8", errors="replace")

    # Collect figure/table labels and references.
    labels_defined = set()
    for m in re.finditer(r"\\label\{((?:fig|tab|table):[^}]+)\}", tex):
        labels_defined.add(m.group(1))
    refs_used = set()
    for m in re.finditer(
        r"\\(?:ref|autoref|Cref|cref|vref)\*?\{((?:fig|tab|table):[^}]+)\}",
        tex,
    ):
        refs_used.add(m.group(1))

    orphan_labels = sorted(labels_defined - refs_used)
    dangling_refs = sorted(refs_used - labels_defined)
    return {
        "labels_defined": sorted(labels_defined),
        "refs_used": sorted(refs_used),
        "orphan_labels": orphan_labels,
        "dangling_refs": dangling_refs,
        "orphan_label_count": len(orphan_labels),
        "dangling_ref_count": len(dangling_refs),
        "verdict": "OK" if not orphan_labels and not dangling_refs else "CHECK",
        "hint": (
            "ソース順序は LaTeX の float 仕様上違反ではない (PDF レベルで確認を)。"
            "ここでは orphan ラベル (定義したが参照なし) と "
            "dangling 参照 (未定義ラベルへの \\ref) のみ検出する。"
        ),
    }


def paper_writing_check_figure_uses_pdf(tex_path: str) -> dict:
    """`\\includegraphics` paths must be vector (.pdf / .eps), not raster.

    Enforces the "Figure file format: import via PDF, NOT PNG" rule
    documented in ``_em_paper_style.py`` (key ``em_paper_style`` /
    ``figure_format``).  IEEE / IEEJ Transactions, IGTE / Compumag /
    CEFC digests all expect vector line-art; embedding ``.png`` /
    ``.jpg`` causes blurry rescaling and triggers publisher pre-flight
    warnings.

    Flagged extensions: ``.png``, ``.jpg``, ``.jpeg``, ``.bmp``,
    ``.tif``, ``.tiff``.  Allowed: ``.pdf``, ``.eps``, no extension
    (lets pdflatex pick the vector copy).

    The intended pattern is ``lab_savefig(fig, "name", medium="paper",
    embed_width_cm=<actual width>)`` (emits both ``name.pdf`` and
    ``name.png`` while checking the 10 pt final-size floor) plus
    ``\\includegraphics{name.pdf}`` -- the PNG is the
    quick-preview / GitHub-README copy, the PDF is what the typesetter
    consumes.

    Args:
        tex_path: .tex file to scan.

    Returns:
        dict with ``includegraphics`` list (one entry per \\includegraphics
        invocation: ``path``, ``extension``, ``line``, ``ok``),
        ``violations`` (the ones using a raster extension), and a
        ``verdict`` of "OK" / "FAIL".
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"tex file not found: {tex_path}"}
    tex = p.read_text(encoding="utf-8", errors="replace")

    raster_exts = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff"}
    vector_exts = {".pdf", ".eps", ""}  # "" = no extension is OK
    rx = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")

    # Build (line, path) list.
    lines = tex.splitlines()
    entries = []
    line_no = 0
    pos = 0
    for line_no, line in enumerate(lines, 1):
        for m in rx.finditer(line):
            path = m.group(1).strip()
            ext = pathlib.PurePosixPath(path).suffix.lower()
            ok = ext in vector_exts and ext not in raster_exts
            entries.append({
                "path": path,
                "extension": ext or "<none>",
                "line": line_no,
                "ok": ok,
            })

    violations = [e for e in entries if not e["ok"]]
    verdict = "OK" if not violations else "FAIL"
    return {
        "includegraphics": entries,
        "n_includegraphics": len(entries),
        "violations": violations,
        "n_violations": len(violations),
        "verdict": verdict,
        "rule": ("\\includegraphics paths must be vector (.pdf / .eps) or "
                 "extension-less; .png / .jpg / .jpeg / .bmp / .tif / .tiff "
                 "are raster and rescale poorly."),
        "advice": ("Use \\includegraphics{name.pdf} (or {name} without "
                   "extension).  lab_savefig writes both name.pdf and "
                   "name.png; reference the .pdf for typesetting and keep "
                   "the .png for quick preview / GitHub README only."),
    }


def paper_writing_check_equation_numbering(tex_path: str) -> dict:
    """方程式番号 (1), (2), ... の連番欠落 / 重複をチェック。

    LaTeX の `\\begin{equation}` / `\\begin{align}` は自動採番だが、
    \\tag{} 手動指定や引用 (\\ref{eq:foo}) の欠落・不一致を検出。

    Args:
        tex_path: .tex ファイル
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"tex file not found: {tex_path}"}
    tex = p.read_text(encoding="utf-8", errors="replace")

    # Count display-math environments (approximate sequential number).
    # Includes starred (un-numbered) variants — we count them so counts reflect
    # actual math density; \eqref detection below still only matches numbered
    # \label targets so starred envs naturally drop out of cross-ref checks.
    equation_envs = len(
        re.findall(
            r"\\begin\{(?:equation|align|gather|multline|eqnarray|flalign)\*?\}"
            r"|\\\[",
            tex,
        )
    )
    # Find equation labels
    eq_labels = set()
    for m in re.finditer(r"\\label\{(eq:[^}]+)\}", tex):
        eq_labels.add(m.group(1))
    # Find equation references
    eq_refs = {}
    for m in re.finditer(r"\\(?:ref|eqref|autoref|Cref)\*?\{(eq:[^}]+)\}", tex):
        key = m.group(1)
        eq_refs[key] = eq_refs.get(key, 0) + 1

    # Dangling refs (referenced but not labeled)
    dangling = [k for k in eq_refs if k not in eq_labels]
    # Unused labels (labeled but not referenced)
    unused = [k for k in eq_labels if k not in eq_refs]

    return {
        "equation_environments": equation_envs,
        "labeled_equations": len(eq_labels),
        "unique_references": len(eq_refs),
        "total_ref_usage": sum(eq_refs.values()),
        "dangling_refs": dangling[:10],
        "unused_labels": unused[:10],
        "verdict": (
            "OK" if not dangling and not unused else "CHECK"
        ),
    }


# ------------------------------------------------------------------
# v0.4.0: 中期 5 tools (cover letter / response / classifier / passive / ref-fmt)
# ------------------------------------------------------------------

def paper_writing_generate_cover_letter(journal: str,
                                          title: str,
                                          abstract: str,
                                          contribution: str = "",
                                          corresponding_author: str = "",
                                          ) -> str:
    """投稿 cover letter の skeleton を生成。

    journal の特徴に触れ、why this journal / what's new / why now の
    3 点を押さえたテンプレを返す。ユーザーが細部を埋める前提。

    Args:
        journal: 投稿先 (例: "IEEE Transactions on Magnetics")
        title: 論文タイトル
        abstract: abstract 本文 (200-250 words)
        contribution: 箇条書きした貢献 (任意、改行区切り)
        corresponding_author: 対応著者名

    Returns:
        cover letter テンプレ (plain text)
    """
    corr = corresponding_author or "Kengo Sugahara"
    contrib_lines = ""
    if contribution:
        lines = [f"  - {c.strip()}" for c in contribution.split("\n") if c.strip()]
        contrib_lines = "\n".join(lines)

    return f"""[Date: YYYY-MM-DD]

To the Editor-in-Chief,
{journal}

Dear Editor,

We are pleased to submit our manuscript titled

  "{title}"

for consideration for publication in {journal}.

**Why this journal**
[1-2 sentences on fit between paper scope and journal's 対象分野.
Reference specific journal sections or past influential papers if relevant.]

**What is new**
{abstract[:400]}{"..." if len(abstract) > 400 else ""}

{f"**Contributions**{chr(10)}{contrib_lines}" if contribution else ""}

**Why now**
[1-2 sentences on timeliness. Why is this problem acute *now*?
Recent paradigm shifts, new measurement capabilities, recent
competing but incomplete works, policy/industry drivers.]

All authors have approved the manuscript and consent to its
submission. The work has not been published elsewhere and is not
under consideration by any other journal. We have no conflicts of
interest to declare.

Sincerely,
{corr}
(corresponding author)
[Affiliation / email / phone]
"""


def paper_writing_generate_response_letter(reviewer_comment: str,
                                             agreement_stance: str = "agree",
                                             response_draft: str = "",
                                             page_line_ref: str = "") -> str:
    """Reviewer コメントへの応答文テンプレ (佐藤 Q40 "agree first, pivot")。

    Args:
        reviewer_comment: 査読者コメント本文
        agreement_stance: "agree" (素直に受入), "partial" (一部受入+代替),
                           "disagree" (根拠を示して丁寧に不同意)
        response_draft: ユーザーが用意した応答の核 (任意)
        page_line_ref: 修正箇所 (例: "P11, L21-L25")

    Returns:
        response letter fragment (plain text)
    """
    stem = f"> {reviewer_comment.strip()}\n\n"
    if agreement_stance == "agree":
        return stem + f"""**Response**: We thank the reviewer for this valuable comment. We fully agree, and have revised the manuscript accordingly.

{response_draft}

The change has been incorporated at {page_line_ref or "[P__, L__]"}.
"""
    elif agreement_stance == "partial":
        return stem + f"""**Response**: We agree with the reviewer in that [reviewer の懸念を寛容に restate]. Unfortunately, [技術的制約 or scope 制限] makes the exact requested change infeasible in the current manuscript. Instead, we have introduced [代替 approach], which addresses the spirit of the reviewer's concern.

{response_draft}

The revised text is at {page_line_ref or "[P__, L__]"}, and Fig./Table N has been updated accordingly.
"""
    else:  # disagree
        return stem + f"""**Response**: We appreciate the reviewer's careful reading. We respectfully disagree with this point because [理由を具体データ + 引用で]. To clarify this in the manuscript, we have revised the text at {page_line_ref or "[P__, L__]"} to make our reasoning explicit.

{response_draft}

We hope this clarification addresses the reviewer's concern.
"""


def paper_writing_classify_reviewer_comment(comment: str) -> dict:
    """Reviewer コメントを佐藤 Q40 の 4-tier に自動分類。

    Difficulty:
      A: typo / 明確提案 (行動: そのまま直す)
      B: 書き直し or 引用追加要求 (行動: 同意 → 根拠 → 改訂文)
      C: 追加実験 / 追加解析要求 (行動: 範囲判定、editor に延長申請)
      D: 技術的に不可能 (行動: agree-first → 理由 → 代替提案)

    Args:
        comment: reviewer コメント 1 件

    Returns:
        dict: difficulty, confidence, triggers, suggested_stance.
    """
    c = comment.lower()
    # D-level: "impossible", "fundamental flaw", "reject" strong
    if any(k in c for k in ("fundamental flaw", "completely wrong",
                             "unable to be fixed", "intractable")):
        return {"difficulty": "D", "confidence": "medium",
                "suggested_stance": "partial-or-disagree",
                "hint": "技術的に不可能な要求。agree first で受け止め、代替提案を用意。"}
    # C-level: demands experiment/data
    if any(k in c for k in ("additional experiment", "new dataset",
                             "please measure", "please compare with",
                             "追加実験", "追加データ", "比較実験")):
        return {"difficulty": "C", "confidence": "high",
                "suggested_stance": "assess feasibility",
                "hint": ("追加実験 / 追加解析要求。実行可能性と期限を "
                          "editor に相談。"),
                "triggers": ["additional experiment / data request"]}
    # B-level: rewrite / citation
    if any(k in c for k in ("unclear", "please clarify", "please cite",
                             "missing reference", "more detail",
                             "明確に", "根拠を示", "引用", "不明瞭")):
        return {"difficulty": "B", "confidence": "high",
                "suggested_stance": "agree",
                "hint": "書き直し + 引用追加。同意して改訂版をインラインで提示。",
                "triggers": ["rewrite / citation request"]}
    # A-level: typo / simple
    if any(k in c for k in ("typo", "grammar", "spelling", "missing comma",
                             "fig.\\s*\\d+ caption", "誤植", "スペル")):
        return {"difficulty": "A", "confidence": "high",
                "suggested_stance": "agree",
                "hint": "Typo / 軽微な指摘。すぐ直して "
                         "'This has been corrected' で応答。",
                "triggers": ["typo / minor"]}
    # Default: B-level (safest)
    return {"difficulty": "B?", "confidence": "low",
            "suggested_stance": "agree-and-rewrite",
            "hint": ("明示的な keyword が無いため B-level 想定で初期応答。"
                      "内容を再読して C (実験要求) や D (不可能) に escalate "
                      "する可能性を確認。"),
            "triggers": []}


def paper_writing_check_passive_voice_ratio(text: str) -> dict:
    """英文の受動態比率を推定。Wallwork §8: Methods は passive 可、
    Discussion / Contribution は active 推奨 (ratio <= 0.40 目安)。

    Args:
        text: 検査対象の英文

    Returns:
        dict: passive_count, active_est, ratio, verdict.
    """
    # Crude heuristic: "be verb + past participle (ed|en)"
    passive_pat = re.compile(
        r"\b(?:am|is|are|was|were|be|been|being|get|gets|got|gotten)\s+\w+(?:ed|en)\b",
        re.IGNORECASE,
    )
    # Total verbs approximate: past tense verbs + "is/are" etc.
    verb_approx_pat = re.compile(
        r"\b\w+(?:ed|ing|s)\b|\b(?:is|are|was|were|be|have|has|had)\b",
        re.IGNORECASE,
    )
    passives = len(passive_pat.findall(text))
    verbs = len(verb_approx_pat.findall(text))
    ratio = passives / max(1, verbs)
    return {
        "passive_count": passives,
        "verb_count_approx": verbs,
        "passive_ratio": round(ratio, 3),
        "target": "<= 0.40 in Discussion/Contribution",
        "verdict": "OK" if ratio <= 0.40 else "HIGH",
        "hint": "Passive 多用は行為者不明でドラマ性を失う。Discussion は "
                "'We demonstrate', 'Our results show' 等 active に置換推奨。",
        "source": "Wallwork §8 + 木下 p.140",
    }


def paper_writing_lint_reference_format(bib_path: str) -> dict:
    """.bib の reference エントリーの完全性と形式を検証。

    チェック項目:
      - year フィールド欠落
      - DOI 欠落 (publisher が DOI 付ける era で重要)
      - pages 欠落 (journal 記事)
      - author フィールド形式 (lastname, firstname or firstname lastname)

    Args:
        bib_path: .bib ファイル

    Returns:
        dict: per-entry missing fields.
    """
    p = pathlib.Path(bib_path)
    if not p.exists():
        return {"error": f"bib not found: {bib_path}"}
    txt = p.read_text(encoding="utf-8", errors="replace")

    entries = []
    # Parse each @article{...} block
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,([^@]*)(?=@|\Z)",
                          txt + "@END", re.DOTALL):
        entry_type = m.group(1).lower()
        key = m.group(2)
        if entry_type == "end":
            continue
        body = m.group(3)
        fields = {}
        for fm in re.finditer(
            r"(\w+)\s*=\s*[{\"]([^}\"]*)[}\"]",
            body,
        ):
            fields[fm.group(1).lower()] = fm.group(2)
        missing = []
        if not fields.get("year"):
            missing.append("year")
        if not fields.get("author"):
            missing.append("author")
        if entry_type in ("article", "inproceedings"):
            if not fields.get("title"):
                missing.append("title")
            if entry_type == "article" and not fields.get("journal"):
                missing.append("journal")
            if entry_type == "inproceedings" and not fields.get("booktitle"):
                missing.append("booktitle")
            if not fields.get("pages"):
                missing.append("pages (optional but recommended)")
            if not fields.get("doi"):
                missing.append("doi (recommended)")
        entries.append({"key": key, "type": entry_type,
                         "missing_fields": missing})

    total = len(entries)
    with_missing = [e for e in entries if e["missing_fields"]]
    return {
        "total_entries": total,
        "entries_with_missing_fields": len(with_missing),
        "problems": with_missing[:20],
        "hint": "year/author/title は必須。IEEE 系は DOI 推奨、pages は journal で必要。",
    }


# ---------------------------------------------------------------------------
# v0.8.0 — 中島・塚本『知的な科学・技術文章の書き方』由来の 2 tools
# ---------------------------------------------------------------------------

_REPETITION_STOP_WORDS = {
    # 頻出で意味の薄い語は除外 (助動詞的)
    "こと", "これ", "それ", "この", "その", "ため", "ところ",
    "とき", "もの", "など", "ほど", "まま", "ただ", "ただし",
    "および", "または", "ならびに", "さらに", "しかし", "しかも",
    "研究", "本研究", "本稿", "本論文", "論文", "結果", "方法",
    "実験", "解析", "評価", "検討", "考察", "目的",
    "図", "表", "式", "Fig", "Table", "Eq",
}


def paper_writing_check_word_repetition(
    text: str,
    window_chars: int = 40,
    min_len: int = 2,
    max_findings: int = 40,
) -> dict:
    """同一単語の近接障害を検出する (中島・塚本『知的な科学・技術文章の書き方』§1.3.5)。

    **問題**: 同じ単語が狭い窓内で繰り返されると "幼稚な印象" を与える。
    本家の例: 「大きく変化し...影響は大きく...大きく進展している」→
    「急速に / 絶大 / 大規模」等、類語で書き換えるべき。

    **検出ロジック**:
    1. テキストから **漢字 2 字以上 / カタカナ 3 字以上 / 英単語** を token として抽出
       (ひらがな・助詞は対象外 — 繰り返し自然なため)
    2. 各 token について `window_chars` 以内で再登場があれば findings に追加
    3. `_REPETITION_STOP_WORDS` (論文特有の頻出語) は除外

    Args:
        text: 検査対象テキスト。
        window_chars: 近接とみなす最大距離 (既定 40 字、中島・塚本の例に合わせる)。
        min_len: 検出対象 token の最小長 (既定 2)。
        max_findings: 返却する findings の上限 (既定 40)。

    Returns:
        {total_findings, unique_repeated_words, findings: [{word, distance,
          first_pos, second_pos, context}, ...], hint, source}

    **注意**: 数値は hint として使う。技術用語 (ESIM 等) の反復は避け
    られないケースもあるので、全件修正は推奨しない。
    """
    import re as _re
    # Kanji / Katakana / English word tokens
    # Hiragana-only sequences are excluded (particles / conjugation endings)
    token_re = _re.compile(
        r"[一-龥]{" + str(min_len) + r",}"       # Kanji run
        r"|[ァ-ヴー]{3,}"                         # Katakana run (3+)
        r"|[A-Za-z][A-Za-z0-9\-]{1,}"             # English word
    )

    tokens: list[tuple[str, int]] = []
    for m in token_re.finditer(text):
        w = m.group(0)
        if w in _REPETITION_STOP_WORDS:
            continue
        tokens.append((w, m.start()))

    findings = []
    seen_pairs: set[tuple[str, int, int]] = set()
    for i, (w, pos) in enumerate(tokens):
        for j in range(i + 1, len(tokens)):
            w2, pos2 = tokens[j]
            if pos2 - pos > window_chars:
                break
            if w == w2:
                key = (w, pos, pos2)
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                cs = max(0, pos - 15)
                ce = min(len(text), pos2 + len(w2) + 15)
                findings.append({
                    "word": w,
                    "distance": pos2 - pos,
                    "first_pos": pos,
                    "second_pos": pos2,
                    "context": text[cs:ce].replace("\n", " "),
                })

    unique_words = sorted({f["word"] for f in findings})
    findings.sort(key=lambda f: f["distance"])

    return {
        "total_findings": len(findings),
        "unique_repeated_words": unique_words,
        "findings": findings[:max_findings],
        "window_chars": window_chars,
        "hint": (
            "同一単語の近接障害 (中島・塚本 §1.3.5) を検出。類語辞書で "
            "書き換え推奨。ただし **技術用語** (固有名詞・データ名) は "
            "反復が避けられないので全件修正は不要、hint として使うこと。"
        ),
        "source": "中島利勝・塚本真也『知的な科学・技術文章の書き方』"
                  "コロナ社 1996, §1.3.5 同一単語の近接障害",
    }


_SENTENCE_END_PATTERNS = [
    # (label, regex of ending — matches last few chars before the punctuation)
    ("である", r"である$"),
    ("だ",       r"だ$"),
    ("ある",     r"(?<!で)ある$"),
    ("ない",     r"ない$"),
    ("した",     r"した$"),
    ("された",   r"された$"),
    ("する",     r"する$"),
    ("される",   r"される$"),
    ("られる",   r"られる$"),
    ("考えられる", r"考えられる$"),
    ("思われる", r"思われる$"),
    ("示す",     r"示す$"),
    ("示した",   r"示した$"),
    ("得た",     r"得た$"),
    ("得る",     r"得る$"),
]


def paper_writing_check_sentence_ending_variety(
    text: str,
    consecutive_threshold: int = 3,
) -> dict:
    """文末表現の単調さ / 連続を検出 (中島・塚本 §1.3.2)。

    **問題**: 文末が「...である。...である。...である。」「...した。...した。」
    と続くと単調で幼稚。本家の例:
      「...である。...がある。...報告がある。...でもある。...代替フロンである。...」
    → 「...と考えられる / ...を示す / ...が報告されている」等で変化をつける。

    **ロジック**:
    1. 句点 (。．.) で文分割
    2. 各文の末尾 (句点直前 15 字) を `_SENTENCE_END_PATTERNS` で分類
    3. 同一 ending が `consecutive_threshold` 連続で発生したら "critical" に
    4. 全文末の histogram + entropy + 連続区間を返す

    Args:
        text: 検査対象テキスト。
        consecutive_threshold: 同一 ending がこの回数連続したら critical (既定 3)。

    Returns:
        {total_sentences, ending_histogram, shannon_entropy_bits,
         critical_runs: [{ending, start_index, length}], hint, source}
    """
    import math as _math
    import re as _re
    sentences = [s.strip() for s in _re.split(r"[。．\.!?！？]", text) if s.strip()]
    if not sentences:
        return {
            "total_sentences": 0,
            "ending_histogram": {},
            "shannon_entropy_bits": 0.0,
            "critical_runs": [],
            "hint": "no sentences found",
            "source": "中島・塚本 §1.3.2",
        }

    labels: list[str] = []
    for s in sentences:
        tail = s[-12:] if len(s) > 12 else s
        hit = "other"
        for lbl, pat in _SENTENCE_END_PATTERNS:
            if _re.search(pat, tail):
                hit = lbl
                break
        labels.append(hit)

    # Histogram
    hist: dict[str, int] = {}
    for lbl in labels:
        hist[lbl] = hist.get(lbl, 0) + 1

    # Shannon entropy (bits) — higher = more varied
    total = sum(hist.values())
    entropy = 0.0
    if total > 0:
        for v in hist.values():
            p = v / total
            if p > 0:
                entropy -= p * _math.log2(p)

    # Critical consecutive runs
    critical: list[dict] = []
    i = 0
    while i < len(labels):
        j = i
        while j + 1 < len(labels) and labels[j + 1] == labels[i]:
            j += 1
        run_len = j - i + 1
        if run_len >= consecutive_threshold and labels[i] != "other":
            critical.append({
                "ending": labels[i],
                "start_sentence_index": i,
                "length": run_len,
                "sample": sentences[i][-40:] if len(sentences[i]) > 40
                           else sentences[i],
            })
        i = j + 1

    return {
        "total_sentences": len(sentences),
        "ending_histogram": dict(sorted(hist.items(), key=lambda kv: -kv[1])),
        "shannon_entropy_bits": round(entropy, 3),
        "critical_runs": critical,
        "consecutive_threshold": consecutive_threshold,
        "hint": (
            "文末が `である` / `した` / `ある` 等で単調に続くと幼稚な印象。"
            "hint として使い、全件機械置換は NG — 論理的に同じ形が "
            "必要な文もある。3 連続以上の critical_runs に注目。"
        ),
        "source": "中島利勝・塚本真也『知的な科学・技術文章の書き方』"
                  "コロナ社 1996, §1.3.2 変化に富む場所へ落とす文末表現",
    }


def paper_writing_check_pdf_obvious_errors(
    pdf_path: str,
    blank_page_char_threshold: int = 50,
    max_report_per_class: int = 20,
) -> dict:
    """論文 PDF を実際に開いて「明らかな体裁エラー」を 6 種類スキャン。
    投稿前の最後の砦。`check_pdf_edge_overflow` と組み合わせて使う。

    検出項目:

    1. **未解決参照**: 本文中の literal `??` (`\\ref` が未定義のとき
       LaTeX が出力) / `[?]` (`\\cite` が未定義のとき biblatex/natbib
       が出力).  PDF 読者には「Section ??」「Eq. (??)」のように見える。

    2. **citation key の生残り**: `[Author2024]` 形式の bibtex キー
       がそのまま PDF に残っているケース (\\bibliographystyle 不一致、
       biblatex package 未読込み、\\cite{} 内のタイポ).

    3. **数式の raw LaTeX 残留**: PDF テキスト中に `\\frac`, `\\sum`,
       `\\alpha`, `\\beta`, `\\text`, `\\sqrt`, `\\mathbf` 等の生
       コマンドが現れる場合, 数式が math mode に入らずに描画された
       (` ` が抜けた、verbatim 区間、`$$` の付け忘れなど).

    4. **空白ページ**: ページ内のテキスト文字数が
       `blank_page_char_threshold` 未満.  論文の途中に空白ページが
       入るのは図の floating で text が次ページに飛んだ等のサイン。

    5. **font 埋め込み欠落 / Type 3 font**: pymupdf でフォント情報を
       取得し, embedded=False または ext='t3' (Type 3 raster) のもの
       を flag.  IEEE / Elsevier 等は Type 1/TrueType の埋め込みを必須
       要件としており, Type 3 は scaling で blur する。

    6. **画像/figure bbox が page MediaBox の外**: `page.get_images` の
       bbox を MediaBox と比較. figure が裁ち落としになっている、または
       `figure*` を忘れて 1 段 column 幅の図がはみ出している等の検出。

    Args:
        pdf_path: 検査対象 PDF
        blank_page_char_threshold: 空白ページと見なす文字数下限
        max_report_per_class: 各カテゴリの報告上限件数

    Returns:
        dict: 各カテゴリの違反一覧 + 総括スコア (0=どれかでも見つかった,
              10=完全クリーン) + hint
    """
    try:
        import pymupdf  # type: ignore
    except ImportError:
        return {"error": "pymupdf not installed. `pip install pymupdf`"}

    p = pathlib.Path(pdf_path)
    if not p.exists():
        return {"error": f"file not found: {pdf_path}"}

    doc = pymupdf.open(pdf_path)
    n_pages = doc.page_count

    # 検出パターン
    re_unresolved_ref = re.compile(
        r"(\bSection\s+\?\?|\bSec\.\s*\?\?|\bSection \(?\?\?\)?|"
        r"\bEq\.\s*\(?\?\?\)?|\bEquation\s+\?\?|"
        r"\bFig\.\s*\?\?|\bFigure\s+\?\?|"
        r"\bTable\s+\?\?|\bTab\.\s*\?\?|"
        r"\?\?\s*\.|\?\?\s*\)|\bpage\s+\?\?|"
        r"\[\?\])"
    )
    # 未展開の bibtex key — IEEE 数字スタイル [1] とは区別: `[英大文字始まりの単語+4桁年]`
    re_bibkey_residue = re.compile(
        r"\[([A-Z][A-Za-z]+\d{4}[A-Za-z]*(?:,\s*[A-Z][A-Za-z]+\d{4}[A-Za-z]*)*)\]"
    )
    # raw LaTeX command tokens
    LATEX_RAW_TOKENS = (
        r"\frac", r"\sum", r"\int", r"\prod", r"\sqrt",
        r"\alpha", r"\beta", r"\gamma", r"\delta", r"\epsilon",
        r"\theta", r"\lambda", r"\mu", r"\nu", r"\sigma", r"\omega",
        r"\Delta", r"\Sigma", r"\Omega", r"\Phi", r"\Psi",
        r"\text", r"\mathbf", r"\mathrm", r"\mathcal", r"\boldsymbol",
        r"\partial", r"\nabla", r"\infty", r"\cdot", r"\times",
        r"\leq", r"\geq", r"\neq", r"\approx",
        r"\begin{", r"\end{", r"\label{", r"\ref{", r"\cite{",
        r"\textbf", r"\textit", r"\emph",
    )
    re_raw_latex = re.compile(
        r"(\\(?:frac|sum|int|prod|sqrt|alpha|beta|gamma|delta|epsilon|"
        r"zeta|eta|theta|iota|kappa|lambda|mu|nu|xi|pi|rho|sigma|tau|"
        r"phi|chi|psi|omega|Delta|Sigma|Omega|Phi|Psi|Lambda|Theta|"
        r"text|mathbf|mathrm|mathcal|boldsymbol|partial|nabla|infty|"
        r"cdot|times|leq|geq|neq|approx|begin\{|end\{|label\{|ref\{|"
        r"cite\{|textbf|textit|emph)\b)"
    )

    # 1-3: テキスト走査
    unresolved_refs: list[dict] = []
    bibkey_residues: list[dict] = []
    raw_latex_hits: list[dict] = []
    blank_pages: list[int] = []
    image_bbox_violations: list[dict] = []

    # キャッシュ: ページごとに raw text を 1 回だけ取得
    page_texts: list[str] = []
    page_bboxes_doc: list[tuple] = []
    for pno in range(n_pages):
        page = doc[pno]
        try:
            txt = page.get_text("text") or ""
        except Exception:
            txt = ""
        page_texts.append(txt)
        page_bboxes_doc.append((page.rect.x0, page.rect.y0,
                                 page.rect.x1, page.rect.y1))

        non_ws = re.sub(r"\s+", "", txt)
        if len(non_ws) < blank_page_char_threshold:
            blank_pages.append(pno + 1)

        # 1. unresolved ref
        for m in re_unresolved_ref.finditer(txt):
            s = max(0, m.start() - 30)
            e = min(len(txt), m.end() + 30)
            context = txt[s:e].replace("\n", " ").strip()
            unresolved_refs.append({
                "page": pno + 1,
                "match": m.group(0),
                "context": context,
            })

        # 2. bibkey residue (IEEE 数字スタイルとの区別 = キー内に英大文字+数字)
        for m in re_bibkey_residue.finditer(txt):
            key = m.group(1)
            # noise reduction: "[Eq. 2024]" "[2024]" 等は数字単体なので
            # \d{4}$ にマッチしても英大文字始まりのトークンが必要
            s = max(0, m.start() - 30)
            e = min(len(txt), m.end() + 30)
            context = txt[s:e].replace("\n", " ").strip()
            bibkey_residues.append({
                "page": pno + 1,
                "key": key,
                "context": context,
            })

        # 3. raw LaTeX
        for m in re_raw_latex.finditer(txt):
            s = max(0, m.start() - 25)
            e = min(len(txt), m.end() + 25)
            context = txt[s:e].replace("\n", " ").strip()
            raw_latex_hits.append({
                "page": pno + 1,
                "match": m.group(0),
                "context": context,
            })

    # 3.5. ASCII pseudo-math residues
    # Hit symptoms in bibliography / table cells when math was authored
    # as plain ASCII like "(R/rho)^2 modulated nu" — IEEEtran bibstyle
    # silently drops the caret and collapses spaces around "rho", "nu",
    # producing visible "(R/rho)2modulatednu" in the PDF.
    #
    # Patterns we flag (in PDF-extracted text):
    #   A. literal "^[0-9a-zA-Z]" — caret in body text
    #      (math mode renders carets as superscript glyphs, not as
    #      literal "^" characters; literal carets escape filter)
    #   B. ")<digit><letter>" — broken superscript adjacent to text
    #      (e.g. "rho)2modulated" : the ")" closed a parenthesized
    #      expression, the digit was meant to be the superscript, and
    #      "modulated" is the next word that lost its leading space)
    #   C. "[a-z]<digit>[a-z]{3,}" — word-digit-word without spaces
    #      (e.g. "x2y_axis" — strong signal of lost math mode + space
    #      collapse).  Excludes 3-letter-or-shorter tail to skip
    #      legitimate alphanumeric IDs like "L1", "R2", "X3".
    pseudo_math_hits: list[dict] = []
    re_caret = re.compile(r"[A-Za-z0-9\)][\^][0-9A-Za-z]")
    re_paren_digit_word = re.compile(r"\)[0-9][a-z]{3,}")
    re_word_digit_word = re.compile(r"\b[a-z]{2,}[0-9][a-z]{4,}\b")
    for pno, txt in enumerate(page_texts):
        for label, regex in (
            ("caret_in_text", re_caret),
            ("paren_digit_word", re_paren_digit_word),
            ("word_digit_word", re_word_digit_word),
        ):
            for m in regex.finditer(txt):
                s = max(0, m.start() - 30)
                e = min(len(txt), m.end() + 30)
                context = txt[s:e].replace("\n", " ").strip()
                pseudo_math_hits.append({
                    "page": pno + 1,
                    "pattern": label,
                    "match": m.group(0),
                    "context": context,
                })

    # 4. 空白ページ -- collected above

    # 5. font embedding チェック
    font_issues: list[dict] = []
    seen_fonts: set[tuple[str, str]] = set()
    for pno in range(n_pages):
        page = doc[pno]
        try:
            fonts = page.get_fonts(full=True)
        except Exception:
            fonts = []
        for f in fonts:
            # tuple: (xref, ext, type, basefont, name, encoding, ref)
            if len(f) < 4:
                continue
            xref = f[0]
            ext = f[1] if len(f) > 1 else ""
            ftype = f[2] if len(f) > 2 else ""
            basefont = f[3] if len(f) > 3 else ""
            key = (basefont, ftype)
            if key in seen_fonts:
                continue
            seen_fonts.add(key)

            is_t3 = (ext or "").lower() in ("t3",) or "type3" in (ftype or "").lower()
            # embedded ?  pymupdf: ext == "n/a" や empty == not embedded
            is_embedded = ext not in ("", "n/a") and not is_t3
            # base14 fonts (Helvetica, Times-Roman, Courier 等) は埋め込み不要だが
            # journal によっては Type1 埋め込みが要求される — 個別 flag
            base14 = basefont in (
                "Helvetica", "Helvetica-Bold", "Helvetica-Oblique",
                "Helvetica-BoldOblique",
                "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic",
                "Courier", "Courier-Bold", "Courier-Oblique",
                "Courier-BoldOblique",
                "Symbol", "ZapfDingbats",
            )
            if is_t3:
                font_issues.append({
                    "page_first_seen": pno + 1,
                    "basefont": basefont,
                    "type": ftype,
                    "ext": ext,
                    "severity": "CRITICAL",
                    "reason": "Type 3 font (raster glyph) — scaling で blur, journal は reject 多い",
                })
            elif not is_embedded and not base14:
                font_issues.append({
                    "page_first_seen": pno + 1,
                    "basefont": basefont,
                    "type": ftype,
                    "ext": ext,
                    "severity": "WARNING",
                    "reason": "font not embedded — reader 環境で別 font に置換される可能性",
                })
            elif not is_embedded and base14:
                # base14 は技術的に OK だが IEEE 等は厳密埋め込みを要求
                font_issues.append({
                    "page_first_seen": pno + 1,
                    "basefont": basefont,
                    "type": ftype,
                    "ext": ext,
                    "severity": "INFO",
                    "reason": "PDF base-14 font (技術的に許容だが IEEE 等は埋め込み推奨)",
                })

    # 6. 画像/figure bbox vs MediaBox
    for pno in range(n_pages):
        page = doc[pno]
        page_rect = page_bboxes_doc[pno]
        px0, py0, px1, py1 = page_rect
        try:
            images = page.get_images(full=True)
        except Exception:
            images = []
        for img in images:
            xref = img[0]
            try:
                bboxes = page.get_image_rects(xref)
            except Exception:
                bboxes = []
            for bb in bboxes:
                bx0, by0, bx1, by1 = bb.x0, bb.y0, bb.x1, bb.y1
                violations = []
                tol = 0.5
                if bx0 < px0 - tol:
                    violations.append(("left", round(px0 - bx0, 2)))
                if bx1 > px1 + tol:
                    violations.append(("right", round(bx1 - px1, 2)))
                if by0 < py0 - tol:
                    violations.append(("top", round(py0 - by0, 2)))
                if by1 > py1 + tol:
                    violations.append(("bottom", round(by1 - py1, 2)))
                if violations:
                    image_bbox_violations.append({
                        "page": pno + 1,
                        "xref": xref,
                        "bbox_pt": [round(v, 2) for v in (bx0, by0, bx1, by1)],
                        "page_bbox_pt": [round(v, 2) for v in page_rect],
                        "sides_over": violations,
                    })

    doc.close()

    # スコア (各カテゴリで違反 0 ⇒ +10/7, 違反あり ⇒ +0)
    classes_present = sum(
        1 for c in (
            unresolved_refs, bibkey_residues, raw_latex_hits,
            pseudo_math_hits,
            blank_pages, font_issues, image_bbox_violations,
        ) if c
    )
    score = round(10 * (1.0 - classes_present / 7.0), 1)

    # severity 集計
    critical_count = (
        sum(1 for f in font_issues if f["severity"] == "CRITICAL")
        + len(unresolved_refs)
        + len(bibkey_residues)
        + len(image_bbox_violations)
    )
    warning_count = (
        len(raw_latex_hits)
        + len(pseudo_math_hits)
        + len(blank_pages)
        + sum(1 for f in font_issues if f["severity"] == "WARNING")
    )

    return {
        "file": str(p),
        "page_count": n_pages,
        "score": score,
        "score_max": 10,
        "critical_count": critical_count,
        "warning_count": warning_count,
        "categories": {
            "unresolved_references": {
                "count": len(unresolved_refs),
                "items": unresolved_refs[:max_report_per_class],
                "severity": "CRITICAL",
                "hint": (
                    "??/[?] が PDF に表示されたら投稿は不可。.tex を最低 2 回"
                    " pdflatex+bibtex でフルコンパイルし直す。citation は biblatex/"
                    "natbib のスタイル不一致もチェック。"
                ),
            },
            "bibkey_residues": {
                "count": len(bibkey_residues),
                "items": bibkey_residues[:max_report_per_class],
                "severity": "CRITICAL",
                "hint": (
                    "[Foo2024] 形式は未展開の bibtex キー — \\cite{} のタイポ、"
                    "bibliography 未指定、bibtex 未実行のサイン。"
                ),
            },
            "raw_latex_residues": {
                "count": len(raw_latex_hits),
                "items": raw_latex_hits[:max_report_per_class],
                "severity": "WARNING",
                "hint": (
                    "\\frac \\sum \\alpha 等が PDF にそのまま現れる場合: $...$ "
                    "の付け忘れ、verbatim 内で意図せず数式コマンド、または "
                    "fragile command in section title 等。"
                ),
            },
            "pseudo_math_residues": {
                "count": len(pseudo_math_hits),
                "items": pseudo_math_hits[:max_report_per_class],
                "severity": "WARNING",
                "hint": (
                    "ASCII pseudo-math が PDF にそのまま現れる: 典型は bib note の"
                    " plain text (e.g. '(R/rho)^2 modulated nu' が IEEEtran bibstyle "
                    "で caret 消失 + space 崩壊して '(R/rho)2modulatednu' に化ける)。"
                    "Fix: bib の note/title 内も math は $...$ 内に入れる、または該当 "
                    "note を削除して title だけ残す。pattern 種類:"
                    " caret_in_text = body text に literal '^'; "
                    " paren_digit_word = ')' の直後に digit + 単語 (superscript 消失 + "
                    "space 崩壊); "
                    " word_digit_word = 単語 + digit + 単語 (math mode 抜けで連結)。"
                ),
            },
            "blank_pages": {
                "count": len(blank_pages),
                "pages": blank_pages,
                "severity": "WARNING",
                "hint": (
                    "途中に空白ページ: float の \\clearpage 過剰、figure/table の "
                    "[H] 強制配置で text が次ページに飛んだ、\\pagebreak の置き間違い。"
                ),
            },
            "font_issues": {
                "count": len(font_issues),
                "items": font_issues[:max_report_per_class],
                "severity_mixed": True,
                "hint": (
                    "Type 3 font を含む submission は IEEE / Elsevier で reject 対象。"
                    "原因の多くは tikz/pgfplots や psfrag 等の図ファイル。 ps2pdf を "
                    "-dPDFSETTINGS=/prepress + -dEmbedAllFonts=true で書き出す、"
                    "あるいは matplotlib 図を PDF/EPS 出力時に rcParams['pdf.fonttype']=42 "
                    "を設定すると Type 1/TrueType になる。"
                ),
            },
            "image_bbox_violations": {
                "count": len(image_bbox_violations),
                "items": image_bbox_violations[:max_report_per_class],
                "severity": "CRITICAL",
                "hint": (
                    "図/画像 bbox が page MediaBox の外: figure* を忘れて 2-column 幅の "
                    "図を 1 段 column に置いた、\\includegraphics[width=...] が大きすぎ、"
                    "minipage 入れ子の幅指定ミス等。"
                ),
            },
        },
        "target": (
            "critical_count = 0 が投稿条件。warning_count は editor/reviewer の"
            "印象に影響するので可能な限りゼロに。"
        ),
        "hint": (
            "score = 10 で全カテゴリクリーン。一つでも違反があると 8.3 以下。"
            "paper_writing_check_pdf_edge_overflow / check_overfull_hbox と "
            "併用して PDF/log 両側から体裁を固めること。"
        ),
        "source": (
            "pymupdf page.get_text/get_fonts/get_image_rects + 経験則パターン "
            "(IEEE/Elsevier submission guidelines)"
        ),
    }


def paper_writing_check_pdf_advanced_anomalies(
    pdf_path: str,
    max_report_per_class: int = 20,
) -> dict:
    """論文 PDF を実際に読んで、reviewer-2 が刺してくる高頻度の体裁
    異常 8 カテゴリを検出する高度版チェッカー。

    `check_pdf_obvious_errors` (基本 7 カテゴリ) + `check_pdf_edge_overflow`
    の後に実行する。スコア 10/10 でもこの advanced 版で違反が出るケース
    がある。

    検出カテゴリ:

    1. **citation_order**: IEEE は cite の初出順に [1] [2] [3] ...
    2. **citation_space**: `word[N]` -> NG, `word [N]` -> OK (IEEE)
    3. **acronym_first_use**: (XX) 初出展開 / 二重定義 / 未定義使用
    4. **caption_position**: 表 caption は表の上 / 図 caption は図の下
    5. **equation_gaps**: (4)->(6) スキップ, 重複, 未参照
    6. **heading_capitalization**: section title が Title Case で揃うか
    7. **orphan_widow**: ページ上下に取り残された 1-2 行段落
    8. **image_overlap**: 図 vs 図 / 図 vs body text の bbox 重なり

    Args:
        pdf_path: 検査対象 PDF
        max_report_per_class: 各カテゴリの報告上限

    Returns:
        dict: file / page_count / score / categories (8 dict) +
        target / hint / source。
    """
    try:
        import pymupdf
    except ImportError:
        return {"error": "pymupdf not installed. `pip install pymupdf`"}

    p = pathlib.Path(pdf_path)
    if not p.exists():
        return {"error": f"file not found: {pdf_path}"}

    doc = pymupdf.open(pdf_path)
    n_pages = doc.page_count

    page_texts: list[str] = []
    page_blocks: list[list] = []
    page_rects: list[tuple] = []
    for pno in range(n_pages):
        page = doc[pno]
        page_rects.append((page.rect.x0, page.rect.y0,
                           page.rect.x1, page.rect.y1))
        page_texts.append(page.get_text("text") or "")
        try:
            d = page.get_text("dict")
        except Exception:
            d = {"blocks": []}
        page_blocks.append(d.get("blocks", []))

    # ============================================================
    # 1+2. Citation numerical order + space before [N]
    # ============================================================
    citation_order_issues: list[dict] = []
    citation_space_issues: list[dict] = []
    # Citation matching is delicate.  PDF text extraction collapses
    # superscripts so `10^5` becomes `105` and `[10^5, 10^7]` reads as
    # `[105, 107]`, indistinguishable from a multi-cite.  Set notation
    # like `[0, a]` looks the same as `[0, supplement]`.  We therefore
    # match a `[...]` block in full and accept it as a citation only if
    # the *entire* bracket content is one of:
    #   - a number: `[N]`
    #   - number + optional alphanumeric note: `[N, supplement]`
    #   - comma-separated numbers: `[N, M, ...]`
    #   - en-dash range: `[N-M]` / `[N–M]`
    # Numbers are restricted to ≤ 99 — IEEE TMag bibliographies almost
    # never exceed two digits, and three-digit numbers are usually math
    # artifacts (10^5 → 105).
    re_cite_bracket = re.compile(r"\[([^\[\]\n]{1,40})\]")
    re_cite_valid = re.compile(
        r"^\s*("
        r"\d{1,2}\s*(?:,\s*[A-Za-z][A-Za-z\-]+)?"   # [N] or [N, supplement]
        r"|\d{1,2}(?:\s*,\s*\d{1,2}){1,9}"           # [N, M, ...]
        r"|\d{1,2}\s*[-–—]+\s*\d{1,2}"               # [N-M]
        r")\s*$"
    )
    # Compact range inner: `[N-M]` / `[N–M]`.
    re_cite_range_inner = re.compile(
        r"^\s*(\d{1,2})\s*[-–—]+\s*(\d{1,2})\s*$")
    # Multi-cite inner: `[N, M, ...]`.
    re_cite_multi_inner = re.compile(
        r"^\s*\d{1,2}(?:\s*,\s*\d{1,2}){1,9}\s*$")
    re_cite_no_space = re.compile(r"\w\[(\d{1,2})\]")
    first_occurrence: dict[int, tuple[int, int]] = {}
    # Outer-bracket range too: `[5]–[7]` (two separate brackets).
    re_range_cite_outer = re.compile(
        r"\[(\d{1,2})\]\s*[-–—]+\s*\[(\d{1,2})\]"
    )

    def _register(k: int, pos: tuple[int, int]) -> None:
        """Register [k] preferring earliest position."""
        if k not in first_occurrence or first_occurrence[k] > pos:
            first_occurrence[k] = pos

    # Pre-compute REFERENCES heading position on each page (a `[N]` past
    # this offset is a bibliography entry; before this offset, a `[N]`
    # at line-start is just PDF text-extraction word-wrap).
    references_starts: dict[int, int] = {}
    in_references_from: "int | None" = None
    for pno, txt in enumerate(page_texts):
        m_ref = re.search(r"\b(?:REFERENCES|References)\b", txt)
        if m_ref:
            references_starts[pno] = m_ref.start()
            if in_references_from is None:
                in_references_from = pno
                break  # only the FIRST REFERENCES heading matters

    def _is_in_bibliography(pno: int, off: int) -> bool:
        if in_references_from is None:
            return False
        if pno > in_references_from:
            return True
        if pno == in_references_from and off >= references_starts[in_references_from]:
            return True
        return False

    for pno, txt in enumerate(page_texts):
        # Stage 1: walk every `[...]` block, accept only citation-shaped
        # contents.  This eliminates `[10^5, 10^7]`, `[0, a]`, etc.
        for m in re_cite_bracket.finditer(txt):
            inner = m.group(1)
            if not re_cite_valid.match(inner):
                continue
            # Skip if we are inside the bibliography region.  Do NOT
            # skip on "line starts with [N]" heuristic alone, because
            # PDF text extraction in two-column layouts often wraps
            # body-text citations like `factors [3]: alpha...` so that
            # `[3]` appears at the start of an extracted line even
            # though visually it is in the body.
            if _is_in_bibliography(pno, m.start()):
                continue
            new_pos = (pno + 1, m.start())
            # Compact range: register all integers in [lo, hi]
            mr = re_cite_range_inner.match(inner)
            if mr:
                n_lo, n_hi = int(mr.group(1)), int(mr.group(2))
                if n_hi > n_lo and n_hi - n_lo <= 30:
                    for k in range(n_lo, n_hi + 1):
                        _register(k, new_pos)
                continue
            # Multi-cite [N, M, ...]: register each
            if re_cite_multi_inner.match(inner):
                for num_str in inner.split(","):
                    _register(int(num_str.strip()), new_pos)
                continue
            # Plain [N] or [N, note]
            mn = re.match(r"\s*(\d{1,2})", inner)
            if mn:
                _register(int(mn.group(1)), new_pos)
        # Stage 2: outer-bracket range `[N]–[M]`
        for m in re_range_cite_outer.finditer(txt):
            n_lo = int(m.group(1))
            n_hi = int(m.group(2))
            if n_hi <= n_lo or n_hi - n_lo > 30:
                continue
            if _is_in_bibliography(pno, m.start()):
                continue
            new_pos = (pno + 1, m.start())
            for k in range(n_lo, n_hi + 1):
                _register(k, new_pos)

    # Sort by (page, offset, citation_number) so ties at the same
    # position (e.g. `[5]–[7]` registers [5], [6], [7] all at the same
    # offset) come out in ascending numerical order instead of
    # whatever order Python dict insertion happened to use.
    cites_in_order = sorted(first_occurrence.items(),
                            key=lambda kv: (kv[1], kv[0]))
    expected = 1
    for n, (page, _) in cites_in_order:
        if n != expected:
            citation_order_issues.append({
                "first_seen_page": page,
                "expected_number": expected,
                "actual_number": n,
            })
            expected = n + 1
        else:
            expected += 1
    # Greek letters and math symbols often appear immediately before a
    # citation in PDF text extraction even when the LaTeX source has a
    # non-breaking space (`~\cite`).  Skip those — they are a pymupdf
    # extraction artifact, not an authoring error.
    MATH_PRECEDERS = set(
        "αβγδεζηθικλμνξοπρστυφχψω"
        "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"
        "∂∇∞∝≈≠≤≥±∓×÷·⋅∈∉⊂⊆⊇∪∩∧∨¬∀∃"
        "ℝℂℕℤℚℍ𝕊𝕋"
        # Look-alike codepoints LaTeX emits in PDF for Greek letters:
        # U+2126 (Ohm sign) for Ω, U+212B (Angstrom) for Å, etc.
        "ΩÅ"
    )
    for pno, txt in enumerate(page_texts):
        for m in re_cite_no_space.finditer(txt):
            ctx_start = max(0, m.start() - 30)
            ctx = txt[ctx_start:m.end() + 5].replace("\n", " ").strip()
            # Character right before `[` is m.group(0)[0]; skip if it is
            # a math symbol (Greek or operator) — LaTeX `$\Omega$~\cite`
            # is correctly authored but pymupdf strips the visible gap.
            preceder = m.group(0)[0]
            if preceder in MATH_PRECEDERS:
                continue
            citation_space_issues.append({
                "page": pno + 1,
                "match": m.group(0),
                "context": ctx,
            })

    # ============================================================
    # 3. Acronym first-use
    # ============================================================
    acronym_issues: list[dict] = []
    re_acronym_def = re.compile(
        r"((?:[A-Z][A-Za-z\-]+\s+){1,6})\(([A-Z][A-Z0-9]{1,5})\)"
    )
    # Real acronym = at least 2 letters, all uppercase, optional trailing
    # digit suffix is allowed only if it's at the END (so `A1` matches but
    # `K0SIBC` does not).  We additionally require minimum 3 letters --
    # 2-letter all-caps tokens (MM, CU, HZ, NO, ON, OF) are too noisy to
    # treat as acronyms in IEEE-style scientific text.
    re_acronym_use = re.compile(r"\b([A-Z]{3,6}\d?|[A-Z]\d+)\b")
    # Math-context indicators near a token suggest it is a math symbol
    # extracted from LaTeX (e.g. K_{\rm SIBC} -> KSIBC, L^2 -> L2,
    # Y_{\rm CLN} -> YCLN).  When seen, demote the token.
    re_math_neighbor = re.compile(
        r"[=+\-\*/<>≤≥≈→←⇒⇐\\^_${}|]|\\[a-zA-Z]+"
    )
    acronym_defs: dict[str, list[tuple[int, int, str]]] = {}
    acronym_uses: dict[str, list[tuple[int, int]]] = {}
    acronym_math_uses: dict[str, int] = {}
    def _looks_like_real_def(phrase: str, acr: str) -> bool:
        """A `phrase (ACR)` is a real acronym definition only if the
        ACR letters roughly match the initials of phrase words.  This
        rejects math labels like `By (A1)` and table cells like
        `HZ (A1)` which structurally look like definitions but aren't.
        """
        # Strip trailing digit from acr ("A1" -> "A")
        acr_letters = "".join(c for c in acr if c.isalpha())
        if not acr_letters:
            return False
        words = [w for w in phrase.split() if w]
        if not words:
            return False
        # Require first letter of ACR matches first letter of first
        # phrase word — minimal sanity check.
        if acr_letters[0].upper() != words[0][0].upper():
            return False
        # Optionally: every ACR letter has a matching initial somewhere
        # in the phrase, in order.
        if len(acr_letters) >= 2 and len(words) >= len(acr_letters):
            initials = [w[0].upper() for w in words]
            j = 0
            for ch in acr_letters:
                while j < len(initials) and initials[j] != ch.upper():
                    j += 1
                if j == len(initials):
                    return False
                j += 1
        return True

    for pno, txt in enumerate(page_texts):
        for m in re_acronym_def.finditer(txt):
            phrase = m.group(1).strip()
            acr = m.group(2)
            if not _looks_like_real_def(phrase, acr):
                continue
            acronym_defs.setdefault(acr, []).append(
                (pno + 1, m.start(), phrase))
        for m in re_acronym_use.finditer(txt):
            acr = m.group(1)
            ctx_before = txt[max(0, m.start() - 50):m.start()]
            if ctx_before.rstrip().endswith("("):
                continue
            # Math-context demotion: if a math operator / LaTeX command
            # is within 4 chars, treat as math token, not acronym.
            tight_ctx = (
                txt[max(0, m.start() - 4):m.start()]
                + txt[m.end():m.end() + 4]
            )
            if re_math_neighbor.search(tight_ctx):
                acronym_math_uses[acr] = acronym_math_uses.get(acr, 0) + 1
                continue
            acronym_uses.setdefault(acr, []).append((pno + 1, m.start()))
    SKIP_ACRONYMS = {
        # Roman numerals
        "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
        "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII",
        # IEEE running header tokens (e.g. "IEEE TRANSACTIONS ON MAGNETICS,
        # VOL. XX, NO. Y")
        "ON", "OF", "THE", "AND", "OR", "TO", "FOR", "IN", "AT", "AS",
        "BY", "VS", "WITH", "FROM",
        "VOL", "NO", "PP", "MO", "YR",
        # Units / scales / single-letter symbols
        "HZ", "KHZ", "MHZ", "GHZ", "THZ",
        "AC", "DC", "RMS",
        "V", "A", "W", "J", "K", "T", "C", "B", "H", "E", "Q",
        # Common acronyms in EM/circuits that appear without local definition
        "IEEE", "ROM", "SPICE", "CAD", "USA", "UK", "JP", "EU",
        "GPU", "CPU", "RAM", "SI", "QD", "CLN", "PEC",
        "FEM", "BEM", "RF", "EM", "VIM", "XFEM",
        "ESIM", "SIBC", "PEEC", "PCB",
        "DOF", "OK", "CPE",
        # IEEE placeholder tokens in the cover page
        "XX", "YY", "MONTH", "VOL", "PP",
        # Math/Greek letter shortcuts that lowercase regex catches as acronyms
        "RHS", "LHS", "EMF", "MMF",
        # Body words that happen to be uppercased (figure callouts,
        # section labels, table headers, bibliographic running text)
        "ALL", "TABLE", "FIG", "FIGURE", "EQ", "REF", "SEC", "PAGE",
        "TYPE", "NAME", "CASE", "TEST", "PEAK", "BAND", "WALL",
        "SLAB", "BLOCK", "AXIS", "AXES", "BOTH", "EACH", "MORE",
        "MOST", "NEAR", "FAR", "LOW", "HIGH", "PUT", "GET", "TOP",
        "BOTTOM", "LEFT", "RIGHT", "FRONT", "BACK", "ONE", "TWO",
        "TENS", "AAA",  # AAA = method name, but appears noisy here
        "SCHUR",  # proper noun
        # Table-header / figure-caption common body words seen all-caps
        # in pymupdf extraction of \textsc / \uppercase / table headers
        "ERROR", "BASIS", "TAYLOR", "MODES", "MODE",
        # Body words capitalised in section headers / italic emphasis
        "ACROSS", "WHILE", "WHEN", "WHERE", "WITH", "FROM",
        "ABOVE", "BELOW", "NEXT", "ELSE", "THEN",
        "WHICH", "THAT", "THESE", "THOSE", "THERE",
        "INTO", "ONTO", "OVER", "UNDER", "AFTER",
        "BEFORE", "DURING", "WITHIN", "WITHOUT",
        "USING", "GIVEN", "HENCE", "THUS", "ALSO",
        # Proper nouns and standard organization acronyms
        "FOSTER", "STIELTJES", "PADE", "NEWTON", "GAUSS",
        "SIAM", "ACM", "PRL", "NATURE", "JMMM",
        "CUBOID", "COST", "MATVEC", "RATIO", "RATIOS",
        "DENSE", "SPARSE", "DIAG", "OFFDIAG", "FILL",
        "COIL", "COILS", "WIRE", "FILAMENT", "WORKPIECE",
        "STEP", "STEPS", "BAND", "BANDS", "ITER", "ITERS",
        "POLE", "POLES", "ZERO", "ZEROS", "ROOT", "ROOTS",
        "POINT", "POINTS", "INDEX", "INDICES",
        "TIME", "TIMES", "RANK", "SIZE", "MEMORY",
        "MATRIX", "MATRICES", "VECTOR", "ENERGY",
        "FORM", "FORMS", "RULES", "VALUE", "VALUES",
        "NORM", "NORMS", "FACTOR", "FACTORS", "TERM", "TERMS",
        # 2-3 letter all-caps physical / chemical symbols
        "MM", "CM", "NM", "UM", "PM", "FM",
        "CU", "FE", "AL", "AU", "AG", "PT", "PB",
        "N", "S", "U", "L", "R", "P", "M", "F", "G", "D", "O",
        # AAA = rational approximation (Adaptive Antoulas-Anderson)
        # but appears so frequently in body text it is noisy
    }
    def _is_math_subscript_collapse(acr: str) -> bool:
        """`KSIBC` = K_SIBC, `YCLN` = Y_CLN: skip if leading 1-2 letters
        stripped yields a known acronym."""
        for n_strip in (1, 2):
            tail = acr[n_strip:]
            if tail in SKIP_ACRONYMS or tail in acronym_defs:
                return True
        return False

    def _is_likely_math_token(acr: str) -> bool:
        """`L2`, `K0`, `H1` etc. — single letter + digit is almost
        always math-mode subscript/superscript collapse."""
        if len(acr) == 2 and acr[0].isalpha() and acr[1].isdigit():
            return True
        return False

    for acr, defs in acronym_defs.items():
        if acr in SKIP_ACRONYMS:
            continue
        if len(defs) >= 2:
            acronym_issues.append({
                "acronym": acr,
                "defined_count": len(defs),
                "pages": [d[0] for d in defs[:max_report_per_class]],
                "severity": "WARNING",
                "issue": "acronym defined multiple times",
            })
    for acr, uses in acronym_uses.items():
        if acr in SKIP_ACRONYMS:
            continue
        if acr in acronym_defs:
            continue
        if len(uses) < 2:
            continue
        if _is_math_subscript_collapse(acr):
            continue
        if _is_likely_math_token(acr):
            continue
        acronym_issues.append({
            "acronym": acr,
            "use_count": len(uses),
            "first_page": uses[0][0],
            "severity": "INFO",
            "issue": "acronym used without explicit definition on first use",
        })

    # ============================================================
    # 4+5. Caption position + equation number gaps
    # ============================================================
    caption_position_issues: list[dict] = []
    equation_gap_issues: list[dict] = []

    for pno in range(n_pages):
        page = doc[pno]
        try:
            images = page.get_images(full=True)
            image_rects = []
            for img in images:
                xref = img[0]
                for r in page.get_image_rects(xref):
                    image_rects.append(r)
        except Exception:
            image_rects = []
        for block in page_blocks[pno]:
            if block.get("type", 0) != 0:
                continue
            bb = block.get("bbox")
            if not bb:
                continue
            block_text = " ".join(
                span.get("text", "")
                for line in block.get("lines", [])
                for span in line.get("spans", [])
            ).strip()
            if re.match(r"^(Fig\.?|Figure)\s*\d+", block_text):
                caption_y = bb[1]
                for r in image_rects:
                    if abs(caption_y - r.y1) < 200 and caption_y < r.y0:
                        caption_position_issues.append({
                            "page": pno + 1,
                            "caption_type": "figure",
                            "caption_y": round(caption_y, 1),
                            "image_y0": round(r.y0, 1),
                            "image_y1": round(r.y1, 1),
                            "issue": "figure caption above image (IEEE: caption BELOW figure)",
                            "snippet": block_text[:60],
                        })
                        break

    eq_numbers_seen: list[int] = []
    re_eq_label = re.compile(r"\((\d{1,3})\)$", re.MULTILINE)
    for pno, txt in enumerate(page_texts):
        for m in re_eq_label.finditer(txt):
            n = int(m.group(1))
            if 1 <= n <= 99:
                eq_numbers_seen.append(n)
    sorted_eq = sorted(set(eq_numbers_seen))
    if sorted_eq:
        for i in range(len(sorted_eq) - 1):
            if sorted_eq[i + 1] - sorted_eq[i] > 1:
                gap = list(range(sorted_eq[i] + 1, sorted_eq[i + 1]))
                equation_gap_issues.append({
                    "type": "skipped",
                    "missing_numbers": gap,
                    "between": (sorted_eq[i], sorted_eq[i + 1]),
                })
        from collections import Counter
        cnt = Counter(eq_numbers_seen)
        for n, c in cnt.items():
            if c > 1:
                equation_gap_issues.append({
                    "type": "duplicate",
                    "number": n,
                    "occurrence_count": c,
                })

    # ============================================================
    # 6+7. Heading capitalization + orphan/widow
    # ============================================================
    heading_issues: list[dict] = []
    orphan_widow_issues: list[dict] = []
    # Section heading: "I. Title On One Line" — exclude any '\n' inside the
    # title (otherwise paragraph body following the title leaks into the
    # capture).
    re_section_heading = re.compile(
        r"^([IVX]+)\.\s+([A-Z][A-Za-z\-:,'’\s\$\{\}]+?)$",
        re.MULTILINE,
    )
    for txt in page_texts:
        for m in re_section_heading.finditer(txt):
            roman = m.group(1)
            title = m.group(2).strip()
            words = title.split()
            major_words = [w for w in words if len(w) >= 4]
            non_title_cased = [w for w in major_words
                                if not w[0].isupper() and w.lower() not in
                                ("and", "with", "from", "into")]
            if non_title_cased:
                heading_issues.append({
                    "section": roman,
                    "title": title,
                    "violations": non_title_cased[:5],
                    "issue": "section title not in IEEE Title Case",
                })

    for pno in range(n_pages):
        page = doc[pno]
        px0, py0, px1, py1 = page_rects[pno]
        page_h = py1 - py0
        for block in page_blocks[pno]:
            if block.get("type", 0) != 0:
                continue
            bb = block.get("bbox")
            if not bb:
                continue
            n_lines = len(block.get("lines", []))
            if n_lines >= 3:
                continue
            block_text = " ".join(
                span.get("text", "")
                for line in block.get("lines", [])
                for span in line.get("spans", [])
            ).strip()
            if len(block_text) < 30:
                continue
            # Skip IEEE/Elsevier running header & footer + page number lines.
            # Heuristic: running header is a single-line block of all-caps text
            # (IEEE: 'IEEE TRANSACTIONS ON MAGNETICS, VOL. ...') OR the paper
            # title repeated at the top.  Page numbers (just digits or
            # 'Page N') are also footers.
            blk_norm = re.sub(r"\s+", " ", block_text).strip()
            if blk_norm.isdigit():
                continue
            if re.match(r"^(IEEE|J\.|Journal|Proc|Sugahara|Page)\s", blk_norm,
                          re.IGNORECASE):
                continue
            # Skip blocks containing the running header phrase
            if "TRANSACTIONS" in blk_norm.upper() and "VOL." in blk_norm.upper():
                continue
            # Skip footer-y page number blocks
            if re.match(r"^\d+\s*$", blk_norm) or re.match(
                    r"^\w+\s+\d+$", blk_norm):
                continue
            # Skip subsection / paragraph headings: "A. Title" / "B. Title"
            # / "Step N:" / "1)" etc.  Match aggressively: a leading
            # capital letter + period is enough, even if pymupdf
            # inserts wide whitespace between the marker and the title,
            # and even if the title starts with a digit ("B. 3-Axis").
            if re.match(r"^[A-Z]\.\s*[A-Z0-9]", blk_norm):
                continue
            if re.match(r"^[A-Z]\)\s*[A-Z0-9]", blk_norm):
                continue
            if re.match(r"^Step\s+\d+", blk_norm):
                continue
            # Skip section number "I.", "II.", "III." etc. at line start
            if re.match(r"^[IVX]+\.\s*[A-Z0-9]", blk_norm):
                continue
            # Skip math-display-only blocks: lines that are >70% non-letter
            # characters (formulas in PDF text extraction look like
            # "r(f) = |Y_R(jw)|sqrt(w)/K_SIBC")
            letters = sum(1 for ch in blk_norm if ch.isalpha())
            if letters / max(len(blk_norm), 1) < 0.4:
                continue
            # Skip figure-legend text (matplotlib legend items typically
            # contain math symbols like ±1/2, ω, or specific keywords
            # such as "slope", "basis", "exact", "approximation").  These
            # are emitted as separate text blocks by matplotlib and look
            # like orphans/widows but are part of the figure.
            blk_lower = blk_norm.lower()
            if any(kw in blk_lower for kw in (
                    "slope", "basis ", "exact", "fitted", "asymptote",
                    " term ", "schur-f", "axisym", "polynomial",
                    "fit ", "ansatz", "envelope", "marker")):
                # legend / annotation block, not body text
                if any(c in blk_norm for c in "+-±×÷/∝→←≈≤≥<>"):
                    continue
                # legend item with math fraction
                if re.search(r"[+\-]\s*\d+\s*/\s*\d+", blk_norm):
                    continue
            at_top = (bb[1] - py0) < page_h * 0.12
            at_bottom = (py1 - bb[3]) < page_h * 0.08
            # Tighten "very few lines" to single-line-only for top/bottom of page
            if at_top and n_lines == 1:
                orphan_widow_issues.append({
                    "page": pno + 1,
                    "position": "top",
                    "n_lines": n_lines,
                    "snippet": block_text[:80],
                })
            elif at_bottom and n_lines == 1:
                orphan_widow_issues.append({
                    "page": pno + 1,
                    "position": "bottom",
                    "n_lines": n_lines,
                    "snippet": block_text[:80],
                })

    # ============================================================
    # 8. Image overlap + body-text overlap
    # ============================================================
    image_overlap_issues: list[dict] = []
    for pno in range(n_pages):
        page = doc[pno]
        try:
            images = page.get_images(full=True)
            image_rects = []
            for img in images:
                xref = img[0]
                for r in page.get_image_rects(xref):
                    image_rects.append((xref, r))
        except Exception:
            image_rects = []
        for i in range(len(image_rects)):
            xref_i, r_i = image_rects[i]
            for j in range(i + 1, len(image_rects)):
                xref_j, r_j = image_rects[j]
                xL = max(r_i.x0, r_j.x0)
                xR = min(r_i.x1, r_j.x1)
                yT = max(r_i.y0, r_j.y0)
                yB = min(r_i.y1, r_j.y1)
                if xR > xL and yB > yT:
                    overlap_area = (xR - xL) * (yB - yT)
                    image_overlap_issues.append({
                        "page": pno + 1,
                        "image_a_xref": xref_i,
                        "image_b_xref": xref_j,
                        "overlap_pt2": round(overlap_area, 1),
                        "issue": "two images overlap on the same page",
                    })
        for xref, r in image_rects:
            for block in page_blocks[pno]:
                if block.get("type", 0) != 0:
                    continue
                bb = block.get("bbox")
                if not bb:
                    continue
                bx0, by0, bx1, by1 = bb
                xL = max(r.x0, bx0)
                xR = min(r.x1, bx1)
                yT = max(r.y0, by0)
                yB = min(r.y1, by1)
                if xR > xL and yB > yT:
                    block_text = " ".join(
                        span.get("text", "")
                        for line in block.get("lines", [])
                        for span in line.get("spans", [])
                    ).strip()
                    if re.match(r"^(Fig\.?|Figure)\s*\d+", block_text):
                        continue
                    if len(block_text) < 15:
                        continue
                    overlap_area = (xR - xL) * (yB - yT)
                    if overlap_area < 100:
                        continue
                    image_overlap_issues.append({
                        "page": pno + 1,
                        "image_xref": xref,
                        "text_snippet": block_text[:60],
                        "overlap_pt2": round(overlap_area, 1),
                        "issue": "image bbox overlaps non-caption text block",
                    })
                    break

    doc.close()

    def _category_counts(items: list[dict], cat_severity: str) -> bool:
        """A category counts against the score if it has at least one
        item AND that item is not purely informational.  Per-item
        severity takes precedence over category-level severity; if
        items have no per-item severity, fall back to category level.
        """
        if not items:
            return False
        # If every item has an explicit severity, count only when at
        # least one is non-INFO.
        if all("severity" in it for it in items):
            return any(
                it["severity"] != "INFO" for it in items
            )
        # Otherwise rely on the category default severity.
        return cat_severity != "INFO"

    classes_present = (
        int(_category_counts(citation_order_issues, "MAJOR"))
        + int(_category_counts(citation_space_issues, "MINOR"))
        + int(_category_counts(acronym_issues, "MINOR"))
        + int(_category_counts(caption_position_issues, "MAJOR"))
        + int(_category_counts(equation_gap_issues, "WARNING"))
        + int(_category_counts(heading_issues, "MINOR"))
        + int(_category_counts(orphan_widow_issues, "INFO"))
        + int(_category_counts(image_overlap_issues, "CRITICAL"))
    )
    score = round(10 * (1.0 - classes_present / 8.0), 1)

    return {
        "file": str(p),
        "page_count": n_pages,
        "score": score,
        "score_max": 10,
        "categories": {
            "citation_order": {
                "count": len(citation_order_issues),
                "items": citation_order_issues[:max_report_per_class],
                "severity": "MAJOR",
                "hint": (
                    "IEEE は cite の初出順に [1] [2] [3] と番号付け。"
                    "cite 順序を並べ替えるか biblatex style=numeric-comp で自動順序化。"
                ),
            },
            "citation_space": {
                "count": len(citation_space_issues),
                "items": citation_space_issues[:max_report_per_class],
                "severity": "MINOR",
                "hint": (
                    "word[N] -> word [N] (IEEE 規則)。"
                    "cite 直前に意図しない `~` (non-breaking space) が無いか確認。"
                ),
            },
            "acronym_first_use": {
                "count": len(acronym_issues),
                "items": acronym_issues[:max_report_per_class],
                "severity": "MINOR",
                "hint": (
                    "略語は初出で `phrase (XXX)` 形で展開、以降は `XXX` のみ。"
                    "二重定義 -> 1 箇所に統一。"
                ),
            },
            "caption_position": {
                "count": len(caption_position_issues),
                "items": caption_position_issues[:max_report_per_class],
                "severity": "MAJOR",
                "hint": (
                    "IEEE/Elsevier: 表 caption は表の上、図 caption は図の下。"
                    "caption 命令を table 直後 / figure 直前に置く。"
                ),
            },
            "equation_gaps": {
                "count": len(equation_gap_issues),
                "items": equation_gap_issues[:max_report_per_class],
                "severity": "WARNING",
                "hint": (
                    "数式番号が飛んでいる: label を追加した後で eq を削除した "
                    "場合に発生。eqref 全部の ref を確認。"
                ),
            },
            "heading_capitalization": {
                "count": len(heading_issues),
                "items": heading_issues[:max_report_per_class],
                "severity": "MINOR",
                "hint": (
                    "IEEE は section title = Title Case (主要語 >=4 文字を頭点大)。"
                    "sentence case と混ざらないように。"
                ),
            },
            "orphan_widow": {
                "count": len(orphan_widow_issues),
                "items": orphan_widow_issues[:max_report_per_class],
                "severity": "INFO",
                "hint": (
                    "段落の 1 行だけがページ上下に取り残されている。"
                    "clubpenalty / widowpenalty を高くして抑制。"
                ),
            },
            "image_overlap": {
                "count": len(image_overlap_issues),
                "items": image_overlap_issues[:max_report_per_class],
                "severity": "CRITICAL",
                "hint": (
                    "図 vs 図、図 vs body text の重なり。"
                    "figure placement の修正か FloatBarrier 挿入で解消。"
                ),
            },
        },
        "target": "全 8 カテゴリ 0 違反。MAJOR/CRITICAL は submit 前に必修正。",
        "hint": (
            "check_pdf_obvious_errors (7 基本カテゴリ) + check_pdf_edge_overflow "
            "の後にこれを実行。score 10/10 でも advanced で違反が出るケース有。"
        ),
        "source": (
            "pymupdf page.get_text/get_images + 経験則パターン "
            "(IEEE/Elsevier submission guidelines + reviewer-2 archetype)"
        ),
    }


def paper_writing_check_prose_density(
    text: str,
    max_report: int = 10,
) -> dict:
    """Detect compression-induced prose-density anti-patterns.

    The failure mode this targets: an author hits a word/page limit and
    compresses by (a) replacing verbs with nominalisations
    ("we augment" -> "the augmentation"), (b) chaining clauses with
    em-dashes and semicolons instead of separate sentences,
    (c) clustering jargon and acronyms.  Surface detectors
    (font / citation order / overfull hbox) all pass but the prose is
    hard to read.

    Five per-sentence axes are scored; sentences with score >= 2 are
    flagged.  When >30% of sentences are flagged, the recommendation
    switches from "rewrite sentence X" to "**reduce content** --- the
    text is past its compression limit, the right move is to drop a
    concept, not compress further".

    Args:
        text: prose to analyze (LaTeX or plain).  LaTeX commands are
            roughly stripped before sentence splitting.
        max_report: max number of flagged sentences to return.

    Returns:
        dict with overall score (0-10), per-sentence issues, and a
        recommendation (rewrite | reduce_content).
    """
    # Roughly strip LaTeX so sentence splitting works on text.
    # Drop \cite{...}, \ref{...}, \eqref{...}, but keep their place
    # so sentence-break heuristics still see surrounding punctuation.
    txt = re.sub(r"\\(cite|ref|eqref|label)\{[^}]*\}", " [REF]", text)
    # Drop $...$ inline math (keep a placeholder so sentence count
    # is preserved but mathy tokens don't inflate jargon density).
    txt = re.sub(r"\$[^$]*\$", " [MATH]", txt)
    # Drop \begin{...}/\end{...} environment markers and most other
    # backslash commands, keeping their text arguments.
    txt = re.sub(r"\\(begin|end)\{[^}]*\}", " ", txt)
    txt = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", txt)
    # Collapse whitespace.
    txt = re.sub(r"\s+", " ", txt).strip()

    # Sentence split: split on `.`, `!`, `?` followed by space + capital
    # or end-of-string.  Avoid splitting on common abbreviations.
    abbrev_protected = re.sub(
        r"\b(et al|e\.g|i\.e|cf|vs|fig|eq|ref|sec|vol|no|pp|Mr|Mrs|Dr|Prof)\.",
        lambda m: m.group(0).replace(".", "<DOT>"),
        txt,
    )
    raw_sents = re.split(r"(?<=[.!?])\s+(?=[A-Z\[])", abbrev_protected)
    sentences = [s.replace("<DOT>", ".").strip() for s in raw_sents if s.strip()]

    # Nominalisation suffixes (English).  Excludes already-counted plain
    # nouns (-ity is enough; -ities and -ity merged below).
    NOM_SUFFIXES = (
        "tion", "sion", "ment", "ance", "ence", "ity",
        "ness", "ship", "ism", "ization", "isation",
    )
    # Words that ARE nominalisations but are extremely common technical
    # terms and don't carry the "stuffed prose" signal — skip them.
    NOM_SKIP = {
        "section", "equation", "construction", "function", "matrix",
        "approximation", "exception", "extension", "iteration",
        "limitation", "introduction", "conclusion", "definition",
        "verification", "compression", "completion", "computation",
        "convergence", "reference", "evidence", "experience",
        # Common physics
        "boundary", "frequency", "geometry", "permeability",
        "permittivity", "magnetism", "conductivity", "asymptote",
    }

    def _count_nominalisations(sent: str) -> int:
        words = re.findall(r"\b[A-Za-z][a-z]+\b", sent)
        n = 0
        for w in words:
            wl = w.lower()
            if wl in NOM_SKIP:
                continue
            if any(wl.endswith(suf) for suf in NOM_SUFFIXES):
                # Exclude trivial -tion words that are core domain nouns
                if len(wl) <= 5:  # "ion", "tion" alone
                    continue
                n += 1
        return n

    def _count_clause_breaks(sent: str) -> int:
        # Em-dash (—), en-dash double (--), semicolon, parenthetical em-dashes
        return (
            sent.count("—")
            + sent.count("–")
            + sent.count(";")
            + sent.count("--")
        )

    def _count_jargon(sent: str) -> int:
        # ALL-CAPS tokens (2+ letters), camelCase, technical compounds
        # (Hyphen-joined Capitalised), Greek/math placeholders
        n = 0
        n += len(re.findall(r"\b[A-Z]{2,}[A-Z0-9]*\b", sent))
        n += len(re.findall(r"\b[A-Z][a-z]+[A-Z][a-z]+\w*\b", sent))  # CamelCase
        n += len(re.findall(r"\b[A-Z][a-z]+-[A-Z][a-z]+\b", sent))    # Cauer-Taylor
        n += sent.count("[MATH]")  # math expressions
        return n

    def _word_count(sent: str) -> int:
        # Count word tokens, treating placeholders as 1
        clean = sent.replace("[MATH]", "X").replace("[REF]", "X")
        return len(re.findall(r"\b\w+\b", clean))

    flagged: list[dict] = []
    nom_total = 0
    break_total = 0
    long_sent_count = 0
    total_sents = 0

    for sent in sentences:
        n_words = _word_count(sent)
        if n_words < 8:
            continue
        total_sents += 1

        n_nom = _count_nominalisations(sent)
        n_brk = _count_clause_breaks(sent)
        n_jar = _count_jargon(sent)
        nom_per_10 = n_nom / max(n_words / 10.0, 0.1)

        nom_total += n_nom
        break_total += n_brk
        if n_words > 35:
            long_sent_count += 1

        # Score thresholds
        score = 0
        reasons = []
        if nom_per_10 >= 2.0:
            score += 1
            reasons.append(
                f"nominalisation density {n_nom}/{n_words}w "
                f"({nom_per_10:.1f}/10w)"
            )
        if n_brk >= 3:
            score += 1
            reasons.append(f"clause breaks: {n_brk} (em-dash/semicolon)")
        if n_jar >= 5 and n_words < 50:
            score += 1
            reasons.append(
                f"jargon clustering: {n_jar} terms in {n_words}w")
        if n_words > 40:
            score += 1
            reasons.append(f"long sentence: {n_words} words")

        if score >= 2:
            snippet = sent[:140] + ("..." if len(sent) > 140 else "")
            flagged.append({
                "score": score,
                "n_words": n_words,
                "nominalisations": n_nom,
                "clause_breaks": n_brk,
                "jargon": n_jar,
                "reasons": reasons,
                "snippet": snippet,
            })

    flagged.sort(key=lambda x: -x["score"])

    if total_sents == 0:
        overall_score = 10.0
        recommendation = "no text to analyze"
    else:
        ratio_flagged = len(flagged) / total_sents
        # 10 if no issues; drop linearly with flagged-ratio
        overall_score = round(max(0.0, 10.0 - 30.0 * ratio_flagged), 1)
        if ratio_flagged > 0.30:
            recommendation = (
                "REDUCE CONTENT --- the text is past its compression "
                "limit; >30% of sentences carry compression anti-patterns. "
                "The right move is to drop a concept (or section), not to "
                "compress further.  Re-evaluate which 3 messages must "
                "stay vs which are 'nice to have'."
            )
        elif ratio_flagged > 0.15:
            recommendation = (
                "REWRITE flagged sentences --- 15-30% of sentences have "
                "compression anti-patterns.  Split each into 2 with active "
                "verbs (avoid 'the augmentation' / 'the closure' / 'the "
                "construction' nominalisations)."
            )
        elif ratio_flagged > 0:
            recommendation = (
                "Minor cleanup --- a few sentences flagged.  Quick local "
                "rewrite suffices."
            )
        else:
            recommendation = "clean prose; no compression anti-patterns."

    return {
        "score": overall_score,
        "score_max": 10,
        "total_sentences": total_sents,
        "flagged_sentences": len(flagged),
        "ratio_flagged": round(len(flagged) / max(total_sents, 1), 2),
        "totals": {
            "nominalisations": nom_total,
            "clause_breaks": break_total,
            "long_sentences": long_sent_count,
        },
        "examples": flagged[:max_report],
        "recommendation": recommendation,
        "hint": (
            "Compression anti-patterns: (1) nominalisations replacing "
            "verbs ('the augmentation' vs 'we augment'); (2) em-dash + "
            "semicolon overuse to glue ideas without subordination; "
            "(3) jargon clustering (>5 technical terms in <50 words); "
            "(4) sentences over 40 words.  When >30% of sentences are "
            "flagged, the fix is REDUCE CONTENT, not further compression."
        ),
        "source": (
            "Wallwork English for Writing Research Papers §13 (sentence "
            "length); Williams Style §3 (verbs not nouns); Strunk-White "
            "rule #17 (omit needless words; corollary: when omitting "
            "would change meaning, omit the concept instead)."
        ),
    }


def paper_writing_suggest_concept_drops(
    text: str,
    max_report: int = 10,
) -> dict:
    """Suggest specific concepts to drop when prose is over the
    compression floor.

    Companion to `paper_writing_check_prose_density`: that tool flags
    *which sentences* are over-compressed; this tool suggests *what to
    drop* to recover natural prose.

    The IGTE 2026 digest debugging (2026-05-21) made the meta-rule
    concrete: when you have N concepts that won't fit in M lines, the
    fix is to drop one concept, not to compress all N harder.  This
    tool automates the manual move of looking through each
    over-compressed sentence and asking "which clause / parallel /
    citation could go without changing the core claim?"

    Patterns detected (in priority order, highest = drop-cheapest):

    1. PARALLEL_CITATION: "...the Warburg-terminated Randles cell
       \\cite{A} and Dirichlet-to-Neumann boundary closure \\cite{B}..."
       Two citations playing the same rhetorical role; drop the second.
       The IGTE digest example: dropping \\cite{GivoliKeller1989}
       (DtN) when \\cite{Randles1947} (Warburg) carries the same
       "non-rational element closes a rational interior" point.

    2. EM_DASH_INTERPOLATION: "main thought --- elaboration --- continued"
       An em-dash interpolation 3-15 words long inside an otherwise
       complete sentence.  Almost always a Nagamine-style review
       addition that doesn't fit.  Move to its own sentence or drop.

    3. PARENTHETICAL_ASIDE: "main sentence (which X)" or
       "...claim (note that Y)..."
       Non-citation parenthetical aside.  Drop if redundant with the
       surrounding prose.

    4. TRAILING_SEMICOLON_CLAUSE: "main statement; X also holds"
       A trailing clause glued via semicolon.  Often a review-driven
       addition.  Promote to its own sentence (if essential) or drop.

    5. CITATION_ONLY_PARENTHETICAL: "(see [N])", "(cf.\\cite{X})"
       A bare citation pointer.  Drop if the citation is repeated in
       the natural flow elsewhere.

    Args:
        text: prose to analyze (raw LaTeX OK; commands are roughly
            stripped before pattern matching, but `\\cite{...}` and
            similar are preserved so PARALLEL_CITATION can match).
        max_report: max candidates to return (default 10).

    Returns:
        dict with ranked drop candidates.  Each candidate carries:
        - pattern: which anti-pattern matched
        - drop_text: the exact text to remove
        - words_saved: estimated word count freed
        - rewrite_preview: surrounding sentence with the candidate
            removed (~80 chars before/after)
        - cost: rough drop-cost (1=cheap reference, 5=core claim)
        - rationale: why dropping is safe (or risky)
    """
    # We work on a LIGHTLY-stripped version: keep \cite{}, drop most
    # other LaTeX commands.  This lets PARALLEL_CITATION see the
    # \cite{} keys and pattern 2-5 work on cleanish prose.
    txt = re.sub(r"\\(begin|end)\{[^}]*\}", " ", text)
    # Keep \cite, \ref, \eqref, \label literal; collapse other commands.
    # Replace \citationcmd{X} -> [CITE:X] for easier pattern matching.
    def _cite_marker(m):
        keys = m.group(2)
        return f" [CITE:{keys}]"
    txt = re.sub(r"\\(cite|citet|citep)\{([^}]*)\}", _cite_marker, txt)
    txt = re.sub(r"\\(ref|eqref|label)\{[^}]*\}", " [REF]", txt)
    # Strip inline math
    txt = re.sub(r"\$[^$]*\$", " [MATH]", txt)
    # Strip other backslash commands
    txt = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^}]*\})?", " ", txt)
    txt = re.sub(r"\s+", " ", txt).strip()

    # Sentence split (same as prose_density)
    abbrev_protected = re.sub(
        r"\b(et al|e\.g|i\.e|cf|vs|fig|eq|ref|sec|vol|no|pp|Mr|Mrs|Dr|Prof)\.",
        lambda m: m.group(0).replace(".", "<DOT>"),
        txt,
    )
    raw_sents = re.split(r"(?<=[.!?])\s+(?=[A-Z\[])", abbrev_protected)
    sentences = [s.replace("<DOT>", ".").strip()
                 for s in raw_sents if s.strip()]

    candidates: list[dict] = []

    def _wc(s: str) -> int:
        return len(re.findall(r"\b\w+\b", s.replace("[MATH]", "X")
                                            .replace("[REF]", "X")
                                            .replace("[CITE:", "X")))

    # Pattern 1: PARALLEL_CITATION
    # "X \cite{A} and Y \cite{B}" — two same-role refs joined by "and"
    re_parallel_cite = re.compile(
        r"([\w\-]+(?:\s+[\w\-]+){1,8}\s*\[CITE:[^\]]+\])"   # X + cite
        r"\s+and\s+"
        r"([\w\-]+(?:\s+[\w\-]+){1,8}\s*\[CITE:[^\]]+\])",  # Y + cite
    )
    for sidx, sent in enumerate(sentences):
        for m in re_parallel_cite.finditer(sent):
            first = m.group(1)
            second = m.group(2)
            drop_text = " and " + second
            ws = _wc(drop_text)
            rewrite = sent[:m.start(2) - 5] + sent[m.end(2):]
            candidates.append({
                "pattern": "PARALLEL_CITATION",
                "sentence_idx": sidx,
                "drop_text": drop_text.strip(),
                "words_saved": ws,
                "rewrite_preview": rewrite[:200],
                "cost": 2,
                "rationale": (
                    f"Drop the second parallel ('{second[:60]}...'); "
                    f"the first ('{first[:60]}...') carries the same "
                    "rhetorical role."
                ),
            })

    # Pattern 2: EM_DASH_INTERPOLATION
    # "main --- elaboration --- continued"
    re_em_dash_interp = re.compile(
        r"\s+---\s+(.{8,120}?)\s+---\s+"
    )
    for sidx, sent in enumerate(sentences):
        for m in re_em_dash_interp.finditer(sent):
            inner = m.group(1)
            ws = _wc(inner)
            if ws < 3:
                continue
            rewrite = sent[:m.start()] + " " + sent[m.end():]
            candidates.append({
                "pattern": "EM_DASH_INTERPOLATION",
                "sentence_idx": sidx,
                "drop_text": "--- " + inner + " ---",
                "words_saved": ws,
                "rewrite_preview": rewrite[:200],
                "cost": 2,
                "rationale": (
                    f"Em-dash interpolation '{inner[:60]}...' splits "
                    "the main sentence; drop it or promote to its own "
                    "sentence."
                ),
            })

    # Pattern 3: PARENTHETICAL_ASIDE (non-citation only)
    re_paren_aside = re.compile(r"\s*\(([^)\[]{8,80})\)")
    for sidx, sent in enumerate(sentences):
        for m in re_paren_aside.finditer(sent):
            inner = m.group(1).strip()
            if "[CITE:" in inner or "[REF]" in inner:
                continue
            ws = _wc(inner)
            if ws < 3:
                continue
            rewrite = sent[:m.start()] + sent[m.end():]
            candidates.append({
                "pattern": "PARENTHETICAL_ASIDE",
                "sentence_idx": sidx,
                "drop_text": "(" + inner + ")",
                "words_saved": ws,
                "rewrite_preview": rewrite[:200],
                "cost": 3,
                "rationale": (
                    f"Parenthetical aside '({inner[:50]}...)'; if the "
                    "main sentence reads fine without it, drop."
                ),
            })

    # Pattern 4: TRAILING_SEMICOLON_CLAUSE
    # "main statement; X also Y"
    re_trailing_semi = re.compile(r";\s+([^;]{20,200})$")
    for sidx, sent in enumerate(sentences):
        m = re_trailing_semi.search(sent)
        if not m:
            continue
        clause = m.group(1).strip().rstrip(".!?")
        ws = _wc(clause)
        if ws < 6:
            continue
        rewrite = sent[:m.start()] + "." + sent[m.end():]
        candidates.append({
            "pattern": "TRAILING_SEMICOLON_CLAUSE",
            "sentence_idx": sidx,
            "drop_text": "; " + clause,
            "words_saved": ws,
            "rewrite_preview": rewrite[:200],
            "cost": 3,
            "rationale": (
                f"Trailing semicolon-clause '; {clause[:50]}...'; "
                "promote to its own sentence (if essential) or drop."
            ),
        })

    # Pattern 5: CITATION_ONLY_PARENTHETICAL "(see [CITE:X])", "(cf. [CITE:X])"
    re_cite_only_paren = re.compile(
        r"\s*\((?:see|cf|per|c\.f|cf\.)\s*\.?\s*\[CITE:[^\]]+\]\)"
    )
    for sidx, sent in enumerate(sentences):
        for m in re_cite_only_paren.finditer(sent):
            drop_text = m.group(0).strip()
            ws = _wc(drop_text)
            rewrite = sent[:m.start()] + sent[m.end():]
            candidates.append({
                "pattern": "CITATION_ONLY_PARENTHETICAL",
                "sentence_idx": sidx,
                "drop_text": drop_text,
                "words_saved": ws,
                "rewrite_preview": rewrite[:200],
                "cost": 1,
                "rationale": (
                    "Bare citation pointer; if the cited work is "
                    "discussed elsewhere this adds nothing."
                ),
            })

    # Sort by (cost asc, words_saved desc)
    candidates.sort(key=lambda c: (c["cost"], -c["words_saved"]))

    # Aggregate predicted savings
    total_words_saved = sum(c["words_saved"] for c in candidates[:max_report])

    if not candidates:
        recommendation = (
            "No droppable concepts found by pattern matching.  The "
            "compression issues, if any, are in sentence structure "
            "rather than removable content.  Consider rewriting "
            "flagged sentences with active verbs instead."
        )
    elif total_words_saved >= 15:
        recommendation = (
            f"Drop candidate #1 (pattern {candidates[0]['pattern']}, "
            f"saves {candidates[0]['words_saved']} words).  This often "
            "frees enough space to fit a 1-page limit AND restore "
            "natural prose flow.  If still over, drop #2."
        )
    else:
        recommendation = (
            f"Candidates exist but small savings ({total_words_saved}w "
            "total).  Combine with sentence-level rewrites for full "
            "compression-floor recovery."
        )

    return {
        "n_candidates": len(candidates),
        "candidates": candidates[:max_report],
        "total_words_saved_if_all_dropped": total_words_saved,
        "recommendation": recommendation,
        "hint": (
            "Run AFTER `check_prose_density` flags compression issues. "
            "PARALLEL_CITATION (cost 2) and CITATION_ONLY_PARENTHETICAL "
            "(cost 1) are the safest drops.  EM_DASH_INTERPOLATION and "
            "TRAILING_SEMICOLON_CLAUSE are cheap if the interpolation "
            "is a Nagamine-style late addition.  PARENTHETICAL_ASIDE "
            "needs human judgement -- review each."
        ),
        "source": (
            "Manual debugging trace, IGTE 2026 digest 2026-05-21: "
            "dropping the parallel DtN reference ([GivoliKeller1989]) "
            "next to the Warburg-Randles cell parallel ([Randles1947]) "
            "freed exactly the 2 lines needed to fit 1 page AND "
            "lifted prose_density score 6.5 -> 10.0.  This tool "
            "automates that pattern."
        ),
    }


def paper_writing_check_typography_hacks(tex_path: str) -> dict:
    """Detect typography hacks used to cram content past a page limit.

    Per the rule "fontを小さくするはだめ" (Sugahara 2026-05-21,
    feedback_compression_floor.md): when a draft cannot fit its page
    limit, shrinking font / line-spacing / margins is **STRICTLY
    FORBIDDEN**.  Reviewers spot non-standard typography immediately,
    and IEEE/IEEJ desk-reject rates are high for papers that do this.
    The correct response is to drop content (see
    `suggest_concept_drops`), not to shrink type.

    Patterns detected:

    - **CRITICAL** font-size overrides for body text:
      - `\\fontsize{X}{Y}\\selectfont` with X < 10
      - bare `\\small`, `\\footnotesize`, `\\scriptsize`, `\\tiny`
        OUTSIDE caption/table/figure/tabular environments
    - **MAJOR** line-spacing / page-area expansion:
      - `\\linespread{X}` with X < 0.95
      - `\\setlength{\\baselineskip}{X}` with X < 12pt
      - `\\setlength{\\textheight}{...}`, `\\setlength{\\textwidth}{...}`,
        `\\setlength{\\topmargin}{...}`, `\\setlength{\\oddsidemargin}{...}`
      - `\\geometry{...}` (overriding journal class margins)
    - **INFO** vspace usage (not a violation, layout tool — Sugahara
      2026-05-21 "\\vspace は、見た目が良くなるならありにしよう"):
      - `\\vspace{-Xmm}` etc. is reported for audit but does NOT count
        against the score.  Use of \\vspace for legitimate visual
        layout (template tightening, figure-text gap reduction) is
        explicitly allowed; the only concern is if it is used to cram
        content past the page limit, which is judged contextually.

    Context-aware: a `\\small` inside `\\begin{table}...\\end{table}` or
    `\\caption{...}` is NORMAL (IEEE style for tables/captions);
    flagged only when in body context.

    Args:
        tex_path: path to the .tex source file.

    Returns:
        dict with score (0-10), per-pattern issue list, severity counts,
        and a recommendation.
    """
    p = pathlib.Path(tex_path)
    if not p.exists():
        return {"error": f"file not found: {tex_path}"}
    text = p.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    issues: list[dict] = []
    critical_count = 0
    major_count = 0
    warning_count = 0

    # Track environment stack: list of environment names currently open
    env_stack: list[str] = []
    # Track if we are inside a caption (single-line {} pair)
    # Heuristic: a line containing "\caption{" until the matching "}"
    caption_open_brace_depth = 0

    # Environments that legitimise small-font commands:
    CAPTION_LIKE_ENVS = {
        "table", "table*", "figure", "figure*",
        "subfigure", "subtable",
        "tabular", "tabularx", "longtable",
        "thebibliography",  # IEEE bibliography uses smaller font
    }

    re_begin = re.compile(r"\\begin\{([a-zA-Z*]+)\}")
    re_end = re.compile(r"\\end\{([a-zA-Z*]+)\}")
    re_small_cmd = re.compile(
        r"\\(small|footnotesize|scriptsize|tiny)\b"
    )
    re_fontsize = re.compile(
        r"\\fontsize\{(\d+(?:\.\d+)?)\}\{(\d+(?:\.\d+)?)\}\s*\\selectfont"
    )
    re_linespread = re.compile(r"\\linespread\{(\d+(?:\.\d+)?)\}")
    re_baselineskip = re.compile(
        r"\\setlength\{\\baselineskip\}\{(\d+(?:\.\d+)?)(\w+)\}"
    )
    re_page_setlength = re.compile(
        r"\\setlength\{\\(textheight|textwidth|topmargin|"
        r"oddsidemargin|evensidemargin|headheight|footskip|"
        r"marginparwidth)\}"
    )
    re_geometry = re.compile(r"\\geometry\{")
    re_vspace_neg = re.compile(
        r"\\vspace\{-(\d+(?:\.\d+)?)(mm|em|pt|cm|ex)\}"
    )

    neg_vspace_count = 0

    for lineno, line in enumerate(lines, 1):
        # Strip comments (but preserve commented-out structure)
        line_nocomment = re.sub(r"(?<!\\)%.*$", "", line)

        # Update environment stack
        for m in re_begin.finditer(line_nocomment):
            env_stack.append(m.group(1))
        for m in re_end.finditer(line_nocomment):
            if env_stack and env_stack[-1] == m.group(1):
                env_stack.pop()

        # Are we inside a caption-like context?
        in_caption_like = any(e in CAPTION_LIKE_ENVS for e in env_stack)
        # \caption{ on this line (or carrying open brace from previous) ?
        if "\\caption{" in line_nocomment:
            caption_open_brace_depth += (
                line_nocomment.count("\\caption{")
            )
        # Crude brace tracking for caption: count open vs close
        if caption_open_brace_depth > 0:
            in_caption_like = True
            opens = line_nocomment.count("{")
            closes = line_nocomment.count("}")
            # If we see net close, decrement (rough)
            net = opens - closes
            if net <= 0 and "\\caption{" not in line_nocomment:
                caption_open_brace_depth = max(
                    0, caption_open_brace_depth - 1)

        # CRITICAL: fontsize override with size < 10
        for m in re_fontsize.finditer(line_nocomment):
            size = float(m.group(1))
            if size < 10:
                issues.append({
                    "line": lineno,
                    "pattern": "FONTSIZE_BELOW_10PT",
                    "severity": "CRITICAL",
                    "match": m.group(0),
                    "reason": (
                        f"\\fontsize sets body to {size}pt; IEEE/IEEJ "
                        "minimum is 10pt.  Desk-reject trigger."
                    ),
                })
                critical_count += 1

        # CRITICAL/MAJOR: small/footnotesize/scriptsize/tiny in body
        for m in re_small_cmd.finditer(line_nocomment):
            cmd = m.group(1)
            if in_caption_like:
                continue  # OK in caption/table/figure
            # Skip the \usepackage[footnotesize]{subfigure} case
            if "\\usepackage" in line_nocomment:
                continue
            # Skip if used as command option (e.g. \begin{thebibliography}{...})
            if "\\bibliography" in line_nocomment.lower():
                continue
            severity = (
                "CRITICAL" if cmd in ("scriptsize", "tiny")
                else "MAJOR"
            )
            issues.append({
                "line": lineno,
                "pattern": f"BODY_FONT_SHRINK_{cmd.upper()}",
                "severity": severity,
                "match": "\\" + cmd,
                "reason": (
                    f"\\{cmd} used in body context (no enclosing "
                    "table/figure/caption); shrinks body text below "
                    "10pt minimum."
                ),
            })
            if severity == "CRITICAL":
                critical_count += 1
            else:
                major_count += 1

        # MAJOR: linespread < 0.95
        for m in re_linespread.finditer(line_nocomment):
            ls = float(m.group(1))
            if ls < 0.95:
                issues.append({
                    "line": lineno,
                    "pattern": "LINESPREAD_BELOW_0.95",
                    "severity": "MAJOR",
                    "match": m.group(0),
                    "reason": (
                        f"\\linespread{{{ls}}} compresses line spacing; "
                        "IEEE/IEEJ require default 1.0 or larger."
                    ),
                })
                major_count += 1

        # MAJOR: baselineskip < 12pt
        for m in re_baselineskip.finditer(line_nocomment):
            val = float(m.group(1))
            unit = m.group(2)
            if unit == "pt" and val < 12:
                issues.append({
                    "line": lineno,
                    "pattern": "BASELINESKIP_BELOW_12PT",
                    "severity": "MAJOR",
                    "match": m.group(0),
                    "reason": (
                        f"baselineskip {val}{unit} below 12pt minimum."
                    ),
                })
                major_count += 1

        # MAJOR: page-area setlength overrides
        for m in re_page_setlength.finditer(line_nocomment):
            param = m.group(1)
            issues.append({
                "line": lineno,
                "pattern": f"PAGE_AREA_OVERRIDE_{param.upper()}",
                "severity": "MAJOR",
                "match": m.group(0),
                "reason": (
                    f"\\setlength{{\\{param}}} overrides journal class "
                    "default page area; reviewers spot this immediately."
                ),
            })
            major_count += 1

        # MAJOR: \geometry{...} (margin override)
        if re_geometry.search(line_nocomment):
            issues.append({
                "line": lineno,
                "pattern": "GEOMETRY_OVERRIDE",
                "severity": "MAJOR",
                "match": "\\geometry{...}",
                "reason": (
                    "\\geometry overrides journal class margins; "
                    "non-standard page area is a desk-reject trigger."
                ),
            })
            major_count += 1

        # INFO: negative vspace (audit only, NOT a score penalty)
        # Sugahara 2026-05-21: "\vspace は、見た目が良くなるならあり" --
        # \vspace is a legitimate layout tool, not a hack.  Report
        # large values for audit but do not penalize the score.
        for m in re_vspace_neg.finditer(line_nocomment):
            val = float(m.group(1))
            unit = m.group(2)
            neg_vspace_count += 1
            threshold = {
                "mm": 3.0, "em": 2.0, "pt": 8.0, "cm": 0.3, "ex": 2.0
            }.get(unit, 3.0)
            if val >= threshold:
                issues.append({
                    "line": lineno,
                    "pattern": f"VSPACE_NEG_LARGE_{unit.upper()}",
                    "severity": "INFO",
                    "match": m.group(0),
                    "reason": (
                        f"\\vspace{{-{val}{unit}}}: large negative "
                        "vspace.  OK if for visual improvement "
                        "(template tightening, figure gap reduction). "
                        "Audit if used to cram content past page limit."
                    ),
                })

    # Score: 10 - (3*critical + 1.5*major + 0.5*warning), clamped >= 0
    score = max(
        0.0,
        round(10.0 - 3.0 * critical_count - 1.5 * major_count
              - 0.5 * warning_count, 1)
    )

    if critical_count > 0:
        recommendation = (
            "CRITICAL typography hacks detected (font shrink below "
            "10pt or scriptsize/tiny in body).  Remove ALL of these "
            "before submission --- desk-reject risk is high.  The "
            "right move is to DROP CONTENT (see suggest_concept_drops)."
        )
    elif major_count > 0:
        recommendation = (
            "MAJOR typography hacks (line spacing / margin overrides / "
            "body \\small).  Reviewers will notice.  Remove these and "
            "reduce content instead."
        )
    elif warning_count > 0:
        recommendation = (
            "WARNING: significant negative vspaces detected.  Audit "
            "whether they are legitimate template tightening or "
            "ad-hoc compression."
        )
    else:
        recommendation = (
            "No typography hacks detected.  Page-limit compliance is "
            "achieved through legitimate prose and content sizing."
        )

    return {
        "file": str(p),
        "score": score,
        "score_max": 10,
        "critical_count": critical_count,
        "major_count": major_count,
        "warning_count": warning_count,
        "negative_vspace_count": neg_vspace_count,
        "issues": issues[:30],
        "n_issues": len(issues),
        "recommendation": recommendation,
        "hint": (
            "Hierarchy of bad page-limit responses (Sugahara 2026-05-21):\n"
            "  0. Typography hacks  --- STRICTLY FORBIDDEN (this tool)\n"
            "  1. Prose compression --- check_prose_density\n"
            "  2. Drop content      --- suggest_concept_drops (the right move)\n"
            "Reviewers spot non-standard typography immediately; "
            "the IEEE/IEEJ heuristic is 'if the authors had to break "
            "typography to fit, the content is too long for this venue'."
        ),
        "source": (
            "IEEE TMag author guidelines (10pt body minimum); "
            "IEEJ style requirement; Sugahara 2026-05-21 feedback "
            "'fontを小さくするはだめ'."
        ),
    }
