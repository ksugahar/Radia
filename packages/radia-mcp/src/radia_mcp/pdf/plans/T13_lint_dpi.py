"""Tier 1 — pdf_lint_image_detection_dpi.

Scan Python source for the non-DPI-invariant image-content-detection
anti-pattern that bit us 2026-05-20 in the (now-retired) 2-up PDF
split path's centered-cover detection.  The same anti-pattern lives
on in any image-projection heuristic that lacks a DPI-aware
threshold (Kindle scanner, OCR pre-trim, etc.) so the lint stays
useful even after the original caller was deleted.

Anti-pattern signature
----------------------
A boolean reduction along an axis of an ink/dark-pixel mask:

  ``( <array> < <intensity> ).any(axis=0)``
  ``( <array> < <intensity> ).any(axis=1)``

At low render DPI (≤150) this works because each page edge has few
anti-aliased pixels.  At higher DPI (≥200) it silently produces
false-positive "content" detections in columns / rows that contain
only a handful of anti-alias artifacts, ruining margin / bbox
heuristics that depend on the column / row being honestly blank.

DPI-invariant alternative
-------------------------
A count threshold scaled to image dimension:

  ``col_dark = (arr < 200).sum(axis=0)``
  ``is_content_col = col_dark >= max(5, int(H * 0.01))``

The threshold ``max(5, int(H * 0.01))`` rejects sparse anti-alias edge
noise (a few pixels per column) while still triggering on any real
text content (tens to hundreds of dark pixels per column when there's
a glyph).

The lint reports occurrences with file:line + matched expression and
a one-line suggestion.  Useful when AI-authoring image-projection
heuristics in the lab.
"""
from __future__ import annotations

import os
import pathlib
import re
from typing import List, Optional

# Match  ``(<expr> < <num>).any(axis=0)``  or ``.any(axis=1)``.
# Allow whitespace.  Don't require any specific variable name.
_PATTERN = re.compile(
    r"\(\s*[^()]+?\s*<\s*\d+\s*\)\s*\.\s*any\s*\(\s*axis\s*=\s*[01]\s*\)"
)

# Less strict: ``.any(axis=0)``  alone, after any expression, useful for
# catching variant formulations where the comparison is named.
_LOOSE_PATTERN = re.compile(r"\.any\s*\(\s*axis\s*=\s*[01]\s*\)")


def _scan_file(path: pathlib.Path, strict: bool = True) -> List[dict]:
    hits: list[dict] = []
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return hits
    pattern = _PATTERN if strict else _LOOSE_PATTERN
    for m in pattern.finditer(text):
        # 1-based line number
        line_no = text.count("\n", 0, m.start()) + 1
        line_start = text.rfind("\n", 0, m.start()) + 1
        line_end = text.find("\n", m.end())
        line_text = text[line_start:line_end if line_end != -1 else None].strip()
        # Heuristic context filter: ignore lines that look unrelated (no
        # array/image variable names).  This is best-effort.
        snippet = m.group(0)
        hits.append({
            "file": str(path),
            "line": line_no,
            "match": snippet,
            "line_text": line_text[:200],
        })
    return hits


def pdf_lint_image_detection_dpi(
    target: str,
    strict: bool = True,
    extensions: tuple[str, ...] = (".py",),
) -> dict:
    """Scan a file or directory for the DPI-not-invariant anti-pattern.

    Args:
        target: file or directory path.
        strict: when True (default), require the comparison ``(_ < N)``
            be immediately followed by ``.any(axis=0|1)``.  When False,
            flag every ``.any(axis=0|1)`` regardless of context (more
            false positives, useful for codebase-wide sweep).
        extensions: file extensions to scan (default ``(".py",)``).

    Returns:
        dict with::

            {
              "target":   <path>,
              "n_files":  int,
              "n_hits":   int,
              "items": [
                {"file": str, "line": int, "match": str, "line_text": str},
                ...
              ],
              "suggestion": str,
            }
    """
    src = pathlib.Path(target)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return {"target": str(src), "error": f"not found: {src}"}

    files: list[pathlib.Path] = []
    if src.is_file():
        if src.suffix.lower() in extensions:
            files = [src]
    else:
        for ext in extensions:
            files.extend(src.rglob(f"*{ext}"))

    hits: list[dict] = []
    for f in files:
        hits.extend(_scan_file(f, strict=strict))

    suggestion = (
        "Replace `(arr < N).any(axis=0)` with a count-threshold:\n"
        "  col_dark = (arr < 200).sum(axis=0)\n"
        "  is_content_col = col_dark >= max(5, int(H * 0.01))\n"
        "See memory/feedback_image_content_detection_dpi.md for the rationale."
    )
    return {
        "target": str(src),
        "n_files": len(files),
        "n_hits": len(hits),
        "items": hits,
        "suggestion": suggestion,
    }
