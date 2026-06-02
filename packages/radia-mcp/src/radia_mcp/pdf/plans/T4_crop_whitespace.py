"""Tier 2 — pdf_crop_whitespace.

Auto-crop white margins from each page of a PDF. Useful for figures
that were exported with extra whitespace (PowerPoint exports tend to
have ~0.5cm margins around content).

Uses pymupdf's pixmap analysis to find the content bounding box.
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_crop_whitespace(pdf_path: str,
                        out_path: str | None = None,
                        margin_pt: float = 5.0,
                        threshold: int = 250) -> str:
    """Auto-crop whitespace from each page.

    Algorithm:
      1. Render the page at 72 DPI.
      2. Find the bounding box of pixels darker than ``threshold``.
      3. Add ``margin_pt`` padding and set the page's ``CropBox``.
      4. Save the cropped PDF.

    Args:
        pdf_path: Source PDF.
        out_path: Destination. Defaults to ``<stem>_cropped.pdf``.
        margin_pt: Padding around content (default 5pt ≈ 1.8mm).
        threshold: Pixel value < this counts as "content" (default 250
                   — anti-aliased text against white).

    Returns:
        Status string with per-page crop deltas.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    src = pathlib.Path(pdf_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    out = pathlib.Path(out_path) if out_path else src.with_name(f"{src.stem}_cropped.pdf")
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    doc = fitz.open(src)
    lines = [f"pdf_crop_whitespace: {src} -> {out}"]
    for i, page in enumerate(doc, 1):
        pix = page.get_pixmap(matrix=fitz.Matrix(1, 1), alpha=False)
        w, h = pix.width, pix.height
        samples = bytes(pix.samples) if hasattr(pix.samples, "__bytes__") else pix.samples
        stride = pix.stride
        comps = pix.n
        x_min, y_min, x_max, y_max = w, h, 0, 0
        for y in range(h):
            row = samples[y * stride:y * stride + w * comps]
            for x in range(w):
                # darkest component < threshold = content
                px = row[x * comps:x * comps + comps]
                if min(px) < threshold:
                    if x < x_min:
                        x_min = x
                    if x > x_max:
                        x_max = x
                    if y < y_min:
                        y_min = y
                    if y > y_max:
                        y_max = y
        if x_min >= x_max or y_min >= y_max:
            lines.append(f"  page {i}: blank, skipping crop")
            continue
        rect = fitz.Rect(
            max(0, x_min - margin_pt),
            max(0, y_min - margin_pt),
            min(w, x_max + margin_pt),
            min(h, y_max + margin_pt),
        )
        # CropBox is in PDF user units (default 1pt = 1px at 72 DPI rendering).
        page.set_cropbox(rect)
        delta_w = w - rect.width
        delta_h = h - rect.height
        lines.append(f"  page {i}: cropped {delta_w:.0f}pt W + {delta_h:.0f}pt H")
    doc.save(out, deflate=True)
    doc.close()
    return "\n".join(lines)
