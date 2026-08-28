"""Classify main, closing, and backup slides in a presentation deck."""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence


_BACKUP_TITLE_RE = re.compile(
    r"(?:^\s*(?:補足|質問対策|予備|付録|補遺|appendix|backup)\s*[：:]?"
    r"|[［\[]\s*(?:質疑用|質問用|backup|appendix)\s*[］\]])",
    re.IGNORECASE,
)

_CLOSING_TITLE_RE = re.compile(
    r"(?:謝辞|ご清聴|q\s*&\s*a|q\s*and\s*a|質疑|questions?|"
    r"acknowledg(?:e)?ments?|thank\s*you|thanks)",
    re.IGNORECASE,
)


def is_backup_title(title: str) -> bool:
    """Return True for an explicitly labelled appendix or Q&A slide."""
    return bool(_BACKUP_TITLE_RE.search(str(title or "")))


def is_closing_title(title: str) -> bool:
    """Return True for thanks, acknowledgement, or Q&A closing slides."""
    return bool(_CLOSING_TITLE_RE.search(str(title or "")))


def is_hidden_slide(slide) -> bool:
    """Return whether the OOXML slide is hidden from the slide show."""
    try:
        return bool(slide._element.get("show") == "0")
    except Exception:
        return False


def classify_deck_sections(
    slides: Sequence,
    title_getter: Callable[[object], str],
) -> dict:
    """Split slides into audience-facing main slides and backup slides.

    Explicitly hidden slides and slides labelled ``補足`` / ``質疑用`` /
    ``Appendix`` / ``Backup`` are backup slides.  A closing slide remains in
    the main presentation because its spoken acknowledgement consumes time;
    a trailing labelled appendix after that closing slide does not.
    """
    rows = []
    main_slide_numbers = []
    backup_slide_numbers = []
    closing_slide_numbers = []
    for slide_no, slide in enumerate(slides, 1):
        title = str(title_getter(slide) or "").strip()
        hidden = is_hidden_slide(slide)
        backup = hidden or is_backup_title(title)
        closing = is_closing_title(title) and not backup
        section = "backup" if backup else "closing" if closing else "main"
        rows.append({
            "slide": slide_no,
            "title": title,
            "section": section,
            "hidden": hidden,
        })
        if backup:
            backup_slide_numbers.append(slide_no)
        else:
            main_slide_numbers.append(slide_no)
        if closing:
            closing_slide_numbers.append(slide_no)

    return {
        "main_slide_numbers": main_slide_numbers,
        "backup_slide_numbers": backup_slide_numbers,
        "closing_slide_numbers": closing_slide_numbers,
        "main_slide_count": len(main_slide_numbers),
        "backup_slide_count": len(backup_slide_numbers),
        "slides": rows,
    }

