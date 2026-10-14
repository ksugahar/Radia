"""Tier 2 — pdf_extract_text.

Extract text content of a PDF for further analysis (e.g. word count,
reading-grade-level, search). Includes a layout-preserving mode.
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_extract_text(pdf_path: str,
                     mode: str = "plain",
                     page_range: str | None = None) -> str:
    """Extract text from a PDF.

    Args:
        pdf_path: Source PDF.
        mode: ``plain`` (raw text, default), ``layout`` (preserves
              line / column structure), ``dict`` (structured page dict
              — returns a JSON-like string).
        page_range: Optional 1-based range like ``"3-7"``. Defaults to
                    all pages.

    Returns:
        The extracted text, or a structured representation.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    doc = fitz.open(p)
    n = len(doc)
    pages = set(range(1, n + 1))
    if page_range:
        pages = set()
        for chunk in page_range.replace(" ", "").split(","):
            if "-" in chunk:
                a, b = chunk.split("-", 1)
                pages.update(range(int(a), int(b) + 1))
            elif chunk:
                pages.add(int(chunk))
        pages &= set(range(1, n + 1))
    mode_kwarg = {"plain": "text", "layout": "blocks", "dict": "dict"}.get(mode, "text")
    out = []
    for i in sorted(pages):
        page = doc[i - 1]
        out.append(f"=== page {i} ===")
        out.append(str(page.get_text(mode_kwarg)))
    doc.close()
    return "\n".join(out)
