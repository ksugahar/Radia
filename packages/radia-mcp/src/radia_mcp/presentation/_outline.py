"""Check whether a talk announces each major section at its transition.

An outline is not merely one agenda slide near the beginning, and it is not a
retrospective statement that the problem has ended and the solution begins.
For a spoken research presentation, navigation has to appear where the topic
changes: "Motivation" before the motivation, "Proposed method" before the
method, and "Results" before the results. Each divider repeats the complete
agenda and highlights only the section that starts now, so the audience sees
both the route and its current position.

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
    """Recognize a repeated agenda/progress slide."""
    labels = _section_labels(slide)
    title = slide["title"] or ""
    body_lines = [line.strip() for line in (slide["text"] or "").splitlines()
                  if line.strip()]
    complete = set(_SECTION_PATTERNS).issubset(labels)
    if _AGENDA_TITLE.search(title) and complete:
        return True, "agenda/progress slide", labels
    if complete and len(body_lines) <= 6:
        return True, "repeated full agenda", labels
    return False, "", labels


def _run_signature(run) -> tuple[str, bool, float | None]:
    """Return style signals that can distinguish the active agenda item."""
    try:
        colour = str(run.font.color.rgb or "")
    except (AttributeError, TypeError):
        colour = ""
    size = run.font.size.pt if run.font.size is not None else None
    return colour, bool(run.font.bold), size


def _highlighted_section_labels(pptx_slide) -> set[str]:
    """Find the uniquely styled section label in a complete agenda."""
    from collections import Counter

    signatures: dict[str, list[tuple[str, bool, float | None]]] = {}
    for shape in pptx_slide.shapes:
        if not getattr(shape, "has_text_frame", False):
            continue
        shape_labels = {name for name, pattern in _SECTION_PATTERNS.items()
                        if pattern.search(shape.text or "")}
        if len(shape_labels) < 2:
            continue
        for paragraph in shape.text_frame.paragraphs:
            for run in paragraph.runs:
                for name, pattern in _SECTION_PATTERNS.items():
                    if pattern.search(run.text or ""):
                        signatures.setdefault(name, []).append(_run_signature(run))

    primary = {name: Counter(values).most_common(1)[0][0]
               for name, values in signatures.items() if values}
    counts = Counter(primary.values())
    if len(primary) < 2 or not counts:
        return set()
    ordinary, ordinary_count = counts.most_common(1)[0]
    if ordinary_count < 2:
        return set()
    return {name for name, signature in primary.items() if signature != ordinary}


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
        from pptx import Presentation
    except ImportError:
        return {"error": "python-pptx not installed."}

    prs = Presentation(pptx_path)
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
            highlighted = _highlighted_section_labels(prs.slides[slide["slide"] - 1])
            dividers.append({
                "slide": slide["slide"],
                "title": slide["title"],
                "detected_by": why,
                "section_labels": sorted(labels),
                "highlighted_sections": sorted(highlighted),
            })

    coverage: dict[str, list[int]] = {}
    for section, boundary in boundaries.items():
        matches = []
        for divider in dividers:
            if not set(boundaries).issubset(divider["section_labels"]):
                continue
            if divider["highlighted_sections"] != [section]:
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
            "反復 Outline が無い。スライド全体を内容で分類し、各主要章の"
            "先頭で全章を再掲し、今から始まる章だけを強調する。")
    elif not all(checks.values()):
        missing = [name for name, slides in coverage.items() if not slides]
        comments.append(
            "各章の先頭で全項目を再掲し、現在章だけを色・太字・大きさ等で"
            "一意に強調する。満たしていない章: "
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
            "Outline は内容が切り替わる場所で繰り返す。毎回すべての章を同じ順で"
            "示し、今から始まる章だけを赤字などで強調する。章名は例であり、"
            "実際の deck の内容から分類する。"
        ),
        "source": "菅原による IGTE'26 デッキ査読の明確化 2026-09-11。",
    }
