"""Tier 2 — pdf_bookmarks_from_outline / pdf_add_bookmark.

Generate PDF bookmarks (outline) from a list, or extract existing
bookmarks. Useful for long submissions / theses where reviewers
navigate via the outline pane.
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_list_bookmarks(pdf_path: str) -> str:
    """List existing bookmarks (outline) in a PDF."""
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    doc = fitz.open(p)
    toc = doc.get_toc(simple=True)
    doc.close()
    if not toc:
        return f"pdf_list_bookmarks: {p}\n  no bookmarks"
    lines = [f"pdf_list_bookmarks: {p}", f"  {len(toc)} bookmark(s):"]
    for level, title, page in toc:
        lines.append(f"  {'  '*(level-1)}- p.{page:3d}: {title}")
    return "\n".join(lines)


def pdf_set_bookmarks(pdf_path: str,
                     outline: list[tuple[int, str, int]],
                     out_path: str | None = None) -> str:
    """Replace the PDF's outline with the given list.

    Args:
        pdf_path: Source PDF.
        outline: list of ``(level, title, page)`` — level starts at 1.
        out_path: Destination. Defaults to in-place rewrite.

    Returns:
        Status with bookmark count.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    out = pathlib.Path(out_path) if out_path else p
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    doc = fitz.open(p)
    doc.set_toc([list(t) for t in outline])
    doc.save(out, incremental=(out == p), encryption=fitz.PDF_ENCRYPT_KEEP)
    doc.close()
    return f"pdf_set_bookmarks: {p} -> {out}  ({len(outline)} bookmark(s))"
