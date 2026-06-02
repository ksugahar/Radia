"""Tier 1 — pdf_inspect.

Page-by-page inspection: dimensions, font list, image / vector
breakdown, total file size. Useful at submission time to verify the
PDF has the expected character (e.g., 1-column journal vs 2-column).
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_inspect(pdf_path: str, verbose: bool = False) -> str:
    """Inspect a PDF and report dimensions / structure per page.

    Args:
        pdf_path: Path to the PDF.
        verbose: If True, list every page's dimensions; otherwise show
                 a summary (page count, modal page size, total fonts).

    Returns:
        Multi-line report.
    """
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    fitz = _load_pymupdf()
    if fitz is None:
        return ("Error: pymupdf not installed. Install with "
                "`pip install -e .[pdf]`.")
    doc = fitz.open(p)
    n = len(doc)
    sizes: dict[tuple[int, int], int] = {}
    img_total = 0
    drawings_total = 0
    fonts: set[str] = set()
    per_page: list[str] = []
    for i, page in enumerate(doc, 1):
        w_pt, h_pt = page.rect.width, page.rect.height
        w_mm = w_pt * 25.4 / 72
        h_mm = h_pt * 25.4 / 72
        key = (round(w_mm), round(h_mm))
        sizes[key] = sizes.get(key, 0) + 1
        img_n = len(page.get_images(full=True))
        drawings_n = len(page.get_drawings())
        img_total += img_n
        drawings_total += drawings_n
        if verbose:
            per_page.append(f"  page {i}: {w_mm:.0f}x{h_mm:.0f} mm, "
                            f"images={img_n}, drawings={drawings_n}")
        for f_ in page.get_fonts(full=True):
            fonts.add(f_[3] if len(f_) > 3 else str(f_))
    doc.close()
    modal = max(sizes.items(), key=lambda kv: kv[1])
    lines = [f"pdf_inspect: {p}",
             f"  pages: {n}",
             f"  modal page size: {modal[0][0]}x{modal[0][1]} mm ({modal[1]}/{n} pages)",
             f"  total images: {img_total}",
             f"  total drawing ops: {drawings_total}",
             f"  unique fonts: {len(fonts)}"]
    if verbose:
        lines.extend(per_page)
    return "\n".join(lines)
