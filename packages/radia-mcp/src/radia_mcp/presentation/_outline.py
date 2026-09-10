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
of every substantial act. This checker discovers the author's ordered section
list from the repeated agenda itself; names such as Theory, Implementation or
Discussion therefore work without a fixed IMRAD vocabulary. A one-slide
closing summary does not need its own divider unless the author lists it as a
section.
"""
from __future__ import annotations

import re
from collections import Counter

from ._kishotenketsu import read_deck

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


def _clean_item(text: str) -> str:
    """Remove list decoration while preserving the author's section name."""
    return re.sub(r"^\s*(?:[•●▪◦-]\s*|\(?\d{1,2}\)?[.):、-]?\s*)", "", text).strip()


def _normalise_item(text: str) -> str:
    return re.sub(r"\s+", " ", _clean_item(text)).casefold()


def _section_key(label: str) -> str:
    """Keep familiar keys for compatibility; accept arbitrary author labels."""
    matches = [name for name, pattern in _SECTION_PATTERNS.items()
               if pattern.search(label)]
    if len(matches) == 1:
        return matches[0]
    ascii_key = re.sub(r"[^a-z0-9]+", "_", label.casefold()).strip("_")
    return ascii_key or _normalise_item(label)


def _run_signature(run) -> tuple[str, bool, float | None]:
    """Return style signals that can distinguish the active agenda item."""
    try:
        colour = str(run.font.color.rgb or "")
    except (AttributeError, TypeError):
        colour = ""
    size = run.font.size.pt if run.font.size is not None else None
    return colour, bool(run.font.bold), size


def _agenda_block(
        pptx_slide) -> list[tuple[str, tuple[tuple[str, bool, float | None], ...]]]:
    """Return the most plausible 2--8 item agenda block on one slide."""
    blocks = []
    for shape in pptx_slide.shapes:
        if not getattr(shape, "has_text_frame", False):
            continue
        items = []
        for paragraph in shape.text_frame.paragraphs:
            text = "".join(run.text for run in paragraph.runs).strip()
            if not text:
                continue
            signatures = tuple(dict.fromkeys(
                _run_signature(run) for run in paragraph.runs if run.text.strip()))
            items.append((_clean_item(text), signatures))
        if not 2 <= len(items) <= 8:
            continue
        blocks.append(items)
    return max(blocks, key=len) if blocks else []


def _highlighted_item_indices(
        block: list[tuple[str, tuple[tuple[str, bool, float | None], ...]]]) -> list[int]:
    """Find the single agenda line whose style differs from the other lines."""
    line_signatures = [set(signatures) for _text, signatures in block]
    occurrences = Counter(signature for signatures in line_signatures
                          for signature in signatures)
    unique_style = [index for index, signatures in enumerate(line_signatures)
                    if any(occurrences[signature] == 1 for signature in signatures)]
    if len(unique_style) == 1:
        return unique_style
    bold = [index for index, signatures in enumerate(line_signatures)
            if any(signature[1] for signature in signatures)]
    if len(bold) == 1 and len(bold) < len(block):
        return bold
    primary = [next(iter(signatures), ("", False, None))
               for signatures in line_signatures]
    counts = Counter(primary)
    if not counts:
        return []
    ordinary, ordinary_count = counts.most_common(1)[0]
    if ordinary_count < 2:
        return []
    return [index for index, signature in enumerate(primary)
            if signature != ordinary]


def presentation_check_outline_slide(pptx_path: str,
                                     backup_title: str = "Backup",
                                     max_slides_before_turn: int = 2) -> dict:
    """Check recurring section dividers at the talk's major transitions.

    ``max_slides_before_turn`` is retained for API compatibility. Dynamic
    section discovery no longer assumes a fixed method/results taxonomy.
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

    raw = []
    for slide in main:
        block = _agenda_block(prs.slides[slide["slide"] - 1])
        if not block:
            continue
        agenda = tuple(_normalise_item(text) for text, _signature in block)
        raw.append((slide, block, agenda))

    repeated = Counter(agenda for _slide, _block, agenda in raw)
    eligible = [entry for entry in raw
                if _AGENDA_TITLE.search(entry[0]["title"] or "")
                or repeated[entry[2]] >= 2]
    reference = (max(eligible, key=lambda entry: (repeated[entry[2]], len(entry[2])))[2]
                 if eligible else ())
    reference_entry = next((entry for entry in eligible if entry[2] == reference), None)
    labels = [text for text, _signature in reference_entry[1]] if reference_entry else []
    keys = []
    for label in labels:
        base = _section_key(label)
        key = base
        suffix = 2
        while key in keys:
            key = f"{base}_{suffix}"
            suffix += 1
        keys.append(key)

    dividers = []
    for slide, block, agenda in eligible:
        if agenda != reference:
            continue
        active_indices = _highlighted_item_indices(block)
        active = [keys[index] for index in active_indices if index < len(keys)]
        dividers.append({
            "slide": slide["slide"],
            "title": slide["title"],
            "detected_by": ("agenda/progress slide" if _AGENDA_TITLE.search(slide["title"] or "")
                            else "repeated full agenda"),
            "section_labels": labels,
            "highlighted_sections": active,
        })

    coverage: dict[str, list[int]] = {key: [] for key in keys}
    expected = 0
    for divider in dividers:
        active = divider["highlighted_sections"]
        if expected < len(keys) and active == [keys[expected]]:
            coverage[keys[expected]].append(divider["slide"])
            expected += 1
        elif active:
            break

    checks = {f"{label} の開始をその場で示す": bool(coverage[key])
              for key, label in zip(keys, labels)}
    score = (round(10.0 * sum(checks.values()) / len(checks), 1)
             if checks else 0.0)
    comments = [f"{'OK  ' if value else 'FAIL'} {label}"
                for label, value in checks.items()]

    if not dividers:
        comments.append(
            "反復 Outline が無い。スライド全体を内容で分類し、各主要章の"
            "先頭で全章を再掲し、今から始まる章だけを強調する。")
    elif not all(checks.values()):
        missing = [label for key, label in zip(keys, labels) if not coverage[key]]
        comments.append(
            "各章の先頭で全項目を再掲し、現在章だけを色・太字・大きさ等で"
            "一意に強調する。満たしていない章: "
            + ", ".join(missing) + ".")

    suggested = [{"section": label, "insert_before_or_near_slide": None}
                 for key, label in zip(keys, labels) if not coverage[key]]
    boundaries = {key: (coverage[key][0] if coverage[key] else None)
                  for key in keys}

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
