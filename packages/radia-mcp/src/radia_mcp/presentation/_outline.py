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

_ITEM_PREFIX = re.compile(
    r"^\s*(?:(?:\d+|[ivxlcdm]+)[\s.)\-:：]+|[-–—•●○▪▫]\s*)",
    re.IGNORECASE)


def _section_labels(slide: dict) -> set[str]:
    text = f"{slide['title'] or ''}\n{slide['text']}"
    return {name for name, pattern in _SECTION_PATTERNS.items()
            if pattern.search(text)}


def _normalise_agenda_item(text: str) -> str:
    """Return a comparison key without numbering or display whitespace."""
    return re.sub(r"\s+", " ", _ITEM_PREFIX.sub("", text or "")).strip().casefold()


def _agenda_items(slide: dict) -> list[str]:
    """Extract the ordered rows of a titled agenda/progress slide."""
    if not _AGENDA_TITLE.search(slide["title"] or ""):
        return []
    lines = [line.strip() for line in (slide["text"] or "").splitlines()
             if line.strip()]
    # A progress divider needs enough context to show a route, while a long
    # prose overview is not an agenda merely because its title says Overview.
    if not 3 <= len(lines) <= 8:
        return []
    keys = [_normalise_agenda_item(line) for line in lines]
    if any(not key for key in keys) or len(set(keys)) != len(keys):
        return []
    return lines


def _looks_like_divider(slide: dict) -> tuple[bool, str, set[str], list[str]]:
    """Recognize a repeated agenda/progress slide."""
    labels = _section_labels(slide)
    title = slide["title"] or ""
    body_lines = [line.strip() for line in (slide["text"] or "").splitlines()
                  if line.strip()]
    agenda_items = _agenda_items(slide)
    if agenda_items:
        return True, "agenda/progress slide", labels, agenda_items
    complete = set(_SECTION_PATTERNS).issubset(labels)
    if _AGENDA_TITLE.search(title) and complete:
        return True, "agenda/progress slide", labels, body_lines
    if complete and len(body_lines) <= 6:
        return True, "repeated full agenda", labels, body_lines
    return False, "", labels, []


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


def _highlighted_agenda_items(pptx_slide, agenda_items: list[str]) -> list[str]:
    """Return agenda rows whose dominant run style differs from the majority."""
    from collections import Counter

    by_key = {_normalise_agenda_item(item): item for item in agenda_items}
    signatures: dict[str, list[tuple[str, bool, float | None]]] = {
        key: [] for key in by_key
    }
    for shape in pptx_slide.shapes:
        if not getattr(shape, "has_text_frame", False):
            continue
        for paragraph in shape.text_frame.paragraphs:
            key = _normalise_agenda_item(paragraph.text or "")
            if key not in signatures:
                continue
            signatures[key].extend(
                _run_signature(run) for run in paragraph.runs
                if (run.text or "").strip())

    primary = {
        key: Counter(values).most_common(1)[0][0]
        for key, values in signatures.items() if values
    }
    counts = Counter(primary.values())
    if len(primary) != len(agenda_items) or not counts:
        return []
    ordinary, ordinary_count = counts.most_common(1)[0]
    if ordinary_count < 2:
        return []
    highlighted = {key for key, signature in primary.items()
                   if signature != ordinary}
    return [item for item in agenda_items
            if _normalise_agenda_item(item) in highlighted]


def _repeated_agenda_coverage(dividers: list[dict]) -> dict:
    """Check an arbitrary taxonomy by repeated order and unique emphasis."""
    from collections import Counter

    sequences = [tuple(_normalise_agenda_item(item)
                       for item in divider["agenda_items"])
                 for divider in dividers if divider["agenda_items"]]
    if not sequences:
        return {"applicable": False}
    canonical, repeats = Counter(sequences).most_common(1)[0]
    if repeats < 2:
        return {"applicable": False}

    matching = [divider for divider, sequence in zip(
        [d for d in dividers if d["agenda_items"]], sequences)
        if sequence == canonical]
    display = matching[0]["agenda_items"]
    coverage = {key: [] for key in canonical}
    observed = []
    for divider in matching:
        highlighted = divider["highlighted_items"]
        if len(highlighted) != 1:
            continue
        key = _normalise_agenda_item(highlighted[0])
        if key not in coverage:
            continue
        coverage[key].append(divider["slide"])
        observed.append(key)

    # Extra dividers within one chapter are harmless, but the first occurrence
    # of every chapter must follow the agenda order without skipping backward.
    first_sequence = [key for key in canonical if coverage[key]]
    observed_unique = list(dict.fromkeys(observed))
    order_ok = observed_unique == first_sequence
    complete = all(coverage.values()) and observed_unique == list(canonical)
    return {
        "applicable": True,
        "items": display,
        "normalised_items": list(canonical),
        "coverage": {
            display[index]: coverage[key]
            for index, key in enumerate(canonical)
        },
        "observed_highlight_order": [
            display[list(canonical).index(key)] for key in observed_unique
        ],
        "same_order": all(sequence == canonical for sequence in sequences),
        "highlight_order_ok": order_ok,
        "complete": complete and all(sequence == canonical for sequence in sequences),
    }


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
        ok, why, labels, agenda_items = _looks_like_divider(slide)
        if ok:
            highlighted = _highlighted_section_labels(prs.slides[slide["slide"] - 1])
            highlighted_items = _highlighted_agenda_items(
                prs.slides[slide["slide"] - 1], agenda_items)
            dividers.append({
                "slide": slide["slide"],
                "title": slide["title"],
                "detected_by": why,
                "section_labels": sorted(labels),
                "highlighted_sections": sorted(highlighted),
                "agenda_items": agenda_items,
                "highlighted_items": highlighted_items,
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

    semantic_checks = {
        "Motivation の開始をその場で示す": bool(coverage["motivation"]),
        "Proposed method の開始をその場で示す": bool(coverage["method"]),
        "Results の開始をその場で示す": bool(coverage["results"]),
    }
    recurring = _repeated_agenda_coverage(dividers)
    custom_taxonomy = recurring.get("applicable", False) and not (
        len(recurring.get("normalised_items", [])) == 3
        and all(len({name for name, pattern in _SECTION_PATTERNS.items()
                     if pattern.search(item)}) == 1
                for item in recurring.get("normalised_items", [])))
    if custom_taxonomy:
        checks = {
            f"{item} の開始を反復 Outline で示す": bool(slides)
            for item, slides in recurring["coverage"].items()
        }
        checks["全章を毎回同じ順序で再掲する"] = recurring["same_order"]
        checks["強調章が章順に進む"] = recurring["highlight_order_ok"]
        score = (10.0 if recurring["complete"]
                 else round(10.0 * sum(checks.values()) / len(checks), 1))
        scoring_mode = "repeated-agenda"
    else:
        checks = semantic_checks
        score = round(10.0 * sum(checks.values()) / len(checks), 1)
        scoring_mode = "semantic-boundaries"
    comments = [f"{'OK  ' if value else 'FAIL'} {label}"
                for label, value in checks.items()]

    if not dividers:
        comments.append(
            "反復 Outline が無い。スライド全体を内容で分類し、各主要章の"
            "先頭で全章を再掲し、今から始まる章だけを強調する。")
    elif not all(checks.values()):
        if custom_taxonomy:
            missing = [name for name, slides
                       in recurring["coverage"].items() if not slides]
            problem = (("満たしていない章: " + ", ".join(missing) + ".")
                       if missing else "章の再掲順または強調順が一致しない。")
        else:
            missing = [name for name, slides in coverage.items() if not slides]
            problem = "満たしていない章: " + ", ".join(missing) + "."
        comments.append(
            "各章の先頭で全項目を再掲し、現在章だけを色・太字・大きさ等で"
            "一意に強調する。" + problem)

    if custom_taxonomy:
        suggested = [
            {"section": section, "insert_before_or_near_slide": None}
            for section, slides in recurring["coverage"].items() if not slides
        ]
    else:
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
        "recurring_agenda": recurring,
        "scoring_mode": scoring_mode,
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
