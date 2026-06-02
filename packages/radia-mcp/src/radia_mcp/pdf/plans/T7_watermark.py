"""Tier 2 — pdf_watermark_add / pdf_watermark_remove.

Add a diagonal text watermark (e.g. "DRAFT" / "CONFIDENTIAL") to each
page of a PDF. Removal requires the source PDF to have a separable
watermark — we only remove watermarks we ourselves added (tagged
``/MCP-Watermark`` in the page annotations).
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_watermark_add(pdf_path: str,
                      text: str = "DRAFT",
                      color: tuple[float, float, float] = (0.85, 0.85, 0.85),
                      fontsize: int = 96,
                      out_path: str | None = None) -> str:
    """Add a diagonal text watermark to every page.

    Args:
        pdf_path: Source PDF.
        text: Watermark text.
        color: RGB (0-1) tuple for the text color (light gray by default).
        fontsize: Font size in pt.
        out_path: Destination. Defaults to ``<stem>_watermarked.pdf``.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    src = pathlib.Path(pdf_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    out = pathlib.Path(out_path) if out_path else src.with_name(f"{src.stem}_watermarked.pdf")
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    doc = fitz.open(src)
    for page in doc:
        rect = page.rect
        cx, cy = rect.width / 2, rect.height / 2
        # Place text centered, rotate -45°.
        page.insert_text(
            (cx - fontsize * len(text) * 0.2, cy),
            text,
            fontsize=fontsize,
            color=color,
            rotate=45,
            overlay=True,
        )
    doc.save(out, deflate=True)
    doc.close()
    return f"pdf_watermark_add: wrote {out} ({len(doc)} page(s), text={text!r})"
