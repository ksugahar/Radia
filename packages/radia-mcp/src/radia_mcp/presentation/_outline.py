"""Check whether a talk announces each major section at its transition.

An outline is not merely one agenda slide near the beginning, and it is not a
retrospective statement that the problem has ended and the solution begins.
For a spoken research presentation, navigation has to appear where the topic
changes: "Motivation" before the motivation, "Proposed method" before the
method, and "Results" before the results. A repeated agenda with the current
item highlighted is equivalent to a sparse section-divider slide.

The section names are examples, not a mandatory taxonomy. Authors first
classify the actual deck into coherent acts, then place a divider at the start
of every substantial act. This checker uses the deck's 起承転結 arc to infer
the common motivation/method/results spine. A one-slide closing summary does
not need its own divider.
"""
from __future__ import annotations

import re

from ._kishotenketsu import presentation_kishotenketsu_check, read_deck


_AGENDA_TITLE = re.compile(
    r"(?i)\boutline\b|\bagenda\b|\bcontents\b|\broad ?map\b|\boverview\b"
    r"|\bwhere we are going\b|\bthe plan of the talk\b"
    r"|目次|構成|本日の流れ|発表の流れ|アウトライン|お話しする順|全体像")

_SECTION_PATTERNS = {
    "motivation": re.compile(
        r"(?i)\bmotivation\b|\bbackground\b|\bproblem\b|\bchallenge\b"
        r"|\bwhy\b|動機|背景|問題|課題|目的"),
    "method": re.compile(
        r"(?i)\bproposed? method\b|\bmethod(?:ology)?\b|\bapproach\b"
        r"|\bformulation\b|\bproposal\b|提案(?:法|手法)?|手法|方法|定式化"),
    "results": re.compile(
        r"(?i)\bresults?\b|\bvalidation\b|\bevaluation\b|\bexperiment(?:s|al)?\b"
        r"|\bbenchmark\b|結果|検証|評価|実験"),
}


def _section_labels(slide: dict) -> set[str]:
    text = f"{slide['title'] or ''}\n{slide['text']}"
    return {name for name, pattern in _SECTION_PATTERNS.items()
            if pattern.search(text)}


def _looks_like_divider(slide: dict) -> tuple[bool, str, set[str]]:
    """Recognize a section card or a repeated agenda/progress slide."""
    labels = _section_labels(slide)
    title = slide["title"] or ""
    body_lines = [line.strip() for line in (slide["text"] or "").splitlines()
                  if line.strip()]
    title_labels = {name for name, pattern in _SECTION_PATTERNS.items()
                    if pattern.search(title)}
    if len(title_labels) == 1 and len(body_lines) <= 3:
        return True, "section title", title_labels
    if _AGENDA_TITLE.search(title) and len(labels) >= 2:
        return True, "agenda/progress slide", labels
    return False, "", labels


def _near(slide_no: int, boundary: int | None, allowance: int) -> bool:
    return boundary is not None and 0 <= boundary - slide_no <= allowance


def presentation_check_outline_slide(pptx_path: str,
                                     backup_title: str = "Backup",
                                     max_slides_before_turn: int = 2) -> dict:
    """Check recurring section dividers at the talk's major transitions.

    ``max_slides_before_turn`` is retained for API compatibility. It is the
    maximum distance allowed between a method divider and the inferred turn;
    two accommodates a short goal/overview slide before the detailed method.
    """
    try:
        from pptx import Presentation  # noqa: F401
    except ImportError:
        return {"error": "python-pptx not installed."}

    _slides, _cut, main = read_deck(pptx_path, backup_title)
    if len(main) < 4:
        return {"error": f"only {len(main)} content slides; too few to need "
                         "section dividers."}

    arc = presentation_kishotenketsu_check(pptx_path, backup_title)
    arc_parts = arc.get("arc", {})
    boundaries = {
        "motivation": main[0]["slide"] if main else None,
        "method": (arc_parts.get("ten") or [None])[0],
        "results": (arc_parts.get("ketsu") or [None])[0],
    }

    dividers = []
    for slide in main:
        ok, why, labels = _looks_like_divider(slide)
        if ok:
            dividers.append({
                "slide": slide["slide"],
                "title": slide["title"],
                "detected_by": why,
                "section_labels": sorted(labels),
            })

    coverage: dict[str, list[int]] = {}
    for section, boundary in boundaries.items():
        matches = []
        for divider in dividers:
            if section not in divider["section_labels"]:
                continue
            allowance = 0 if section == "motivation" else max_slides_before_turn
            if _near(divider["slide"], boundary, allowance):
                matches.append(divider["slide"])
        coverage[section] = matches

    checks = {
        "Motivation の開始をその場で示す": bool(coverage["motivation"]),
        "Proposed method の開始をその場で示す": bool(coverage["method"]),
        "Results の開始をその場で示す": bool(coverage["results"]),
    }
    score = round(10.0 * sum(checks.values()) / len(checks), 1)
    comments = [f"{'OK  ' if value else 'FAIL'} {label}"
                for label, value in checks.items()]

    if not dividers:
        comments.append(
            "章扉が無い。スライド全体を内容で分類し、各主要章の先頭に、"
            "これから話す章名を大きく示す区切りの枚を入れる。")
    elif not all(checks.values()):
        missing = [name for name, slides in coverage.items() if not slides]
        comments.append(
            "一枚の総目次だけでは足りない。欠けている章の開始位置に章扉を置く: "
            + ", ".join(missing) + ".")

    suggested = [
        {"section": section, "insert_before_or_near_slide": boundary}
        for section, boundary in boundaries.items()
        if not coverage[section]
    ]

    return {
        "score": score,
        "score_max": 10,
        "content_slides": len(main),
        "section_boundaries": boundaries,
        "section_dividers": dividers,
        "outline_slides": dividers,
        "outline_at_boundary": sorted({n for nums in coverage.values() for n in nums}),
        "section_coverage": coverage,
        "suggested_outline": suggested,
        "checks": checks,
        "comments": comments,
        "hint": (
            "Outline は一枚の一覧ではなく、内容が切り替わる場所で現在地を知らせる"
            "章扉として使う。章名は例であり、実際の deck の内容から分類する。"
            "各回は章名だけの疎な枚、または現在章を強調した反復 agenda とする。"
        ),
        "source": "菅原による IGTE'26 デッキ査読の明確化 2026-09-11。",
    }
