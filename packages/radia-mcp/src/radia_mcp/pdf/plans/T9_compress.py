"""Tier 2 — pdf_compress.

Reduce PDF file size by rewriting with deflate compression and
optionally downscaling embedded images. Many journal submission
systems have a hard 10 MB / 25 MB cap.
"""
from __future__ import annotations

import io
import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_compress(pdf_path: str,
                 out_path: str | None = None,
                 image_dpi: int = 150,
                 image_quality: int = 75) -> str:
    """Recompress a PDF.

    Strategy:
        - Re-save with ``deflate=True``, ``garbage=4`` (drop unused).
        - Optionally downscale rasters above ``image_dpi`` and re-encode
          JPEG at ``image_quality``.

    Args:
        pdf_path: Source PDF.
        out_path: Destination. Defaults to ``<stem>_compressed.pdf``.
        image_dpi: Target DPI for rasters (default 150 — Web-friendly).
                   Set to 0 to skip image downscaling.
        image_quality: JPEG quality 1-100 (default 75).

    Returns:
        Before/after size + savings ratio.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    src = pathlib.Path(pdf_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    out = pathlib.Path(out_path) if out_path else src.with_name(f"{src.stem}_compressed.pdf")
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    before = src.stat().st_size

    doc = fitz.open(src)
    if image_dpi and image_dpi > 0:
        for page in doc:
            for img in page.get_images(full=True):
                xref = img[0]
                try:
                    pix = fitz.Pixmap(doc, xref)
                except Exception:
                    continue
                if pix.n - pix.alpha > 3:  # CMYK or weird; skip
                    pix = None
                    continue
                # Compute current DPI from the smallest dimension on page.
                rects = page.get_image_rects(xref)
                if not rects:
                    pix = None
                    continue
                rect = rects[0]
                cur_dpi = pix.width / (rect.width / 72) if rect.width > 0 else 9999
                if cur_dpi <= image_dpi:
                    pix = None
                    continue
                scale = image_dpi / cur_dpi
                new_w = int(pix.width * scale)
                new_h = int(pix.height * scale)
                if new_w < 8 or new_h < 8:
                    pix = None
                    continue
                small = fitz.Pixmap(pix, new_w, new_h)
                # Re-embed via raw JPEG.
                try:
                    jpg = small.tobytes(output="jpeg", jpg_quality=image_quality)
                    doc._update_object(xref, b"<<>>")  # noqa: SLF001
                    doc.update_stream(xref, jpg)
                except Exception:
                    pass
                pix = None
                small = None

    doc.save(out, deflate=True, garbage=4, clean=True)
    doc.close()
    after = out.stat().st_size
    saved = (before - after) / before * 100 if before else 0
    return (f"pdf_compress: {src} ({before/1024:.1f} KB) -> "
            f"{out} ({after/1024:.1f} KB, saved {saved:.0f}%)")
