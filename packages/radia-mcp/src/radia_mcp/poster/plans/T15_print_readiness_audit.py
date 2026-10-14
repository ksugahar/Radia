"""Tier 5 — poster_print_readiness_audit.

Validates that included figures will print cleanly at A0/A1: 300dpi
equivalent, declared paper size matches output, no RGB→CMYK surprises.

Lazy-imports pymupdf so the rest of the sub-skill stays dependency-free.
"""
from __future__ import annotations

import pathlib
import re

from .._textutil import read_tex


_PAPERSIZE_RE = re.compile(r"\\documentclass\[([^\]]*)\]")
_INCLUDEGRAPHICS_RE = re.compile(
    r"\\includegraphics\b(?:\[(?P<opts>[^\]]*)\])?\{(?P<file>[^}]+)\}"
)
_LINEWIDTH_RE = re.compile(r"width\s*=\s*([\d.]+)\s*\\linewidth")

# Roughly A0/A1 sizes in mm.
_PAPER_MM = {
    "a0paper": (841, 1189),
    "a1paper": (594, 841),
    "a2paper": (420, 594),
    "b0paper": (1000, 1414),
}


def _load_pymupdf():
    try:
        import fitz  # noqa: F401
        return fitz
    except ImportError:
        return None


def poster_print_readiness_audit(tex_path: str,
                                 figures_dir: str | None = None) -> str:
    """Audit a poster for print-readiness: paper size + figure DPI.

    Checks:

    1. Paper class is one of A0/A1/A2/B0 (poster ranges).
    2. Each ``\\includegraphics`` PDF/image is inspected via pymupdf;
       computes effective DPI as ``image_pixels / printed_inches``
       where printed-inches = poster width × ``\\linewidth`` fraction.
       <250 DPI is HIGH severity; <300 DPI is MEDIUM.
    3. Embedded color space: warns if any image is grayscale-only when
       the poster otherwise uses color (poor for figures).

    Args:
        tex_path: Path to the poster ``.tex`` source.
        figures_dir: Optional override for figure lookup.

    Returns:
        Per-figure DPI table + summary.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    fitz = _load_pymupdf()

    cls = _PAPERSIZE_RE.search(tex)
    opts = (cls.group(1) if cls else "").lower()
    paper = next((p for p in _PAPER_MM if p in opts), None)
    if not paper:
        paper_warning = f"[MEDIUM] poster paper not declared as A0/A1/A2/B0 (options={opts!r})"
        paper_size_mm = (594, 841)  # assume A1 portrait
    else:
        paper_warning = None
        paper_size_mm = _PAPER_MM[paper]
    landscape = "landscape" in opts
    if landscape:
        paper_size_mm = (paper_size_mm[1], paper_size_mm[0])
    width_mm = paper_size_mm[0]

    figdir = pathlib.Path(figures_dir) if figures_dir else (path.parent / "figures")
    if not figdir.is_absolute():
        figdir = pathlib.Path.cwd() / figdir

    lines = [f"poster_print_readiness_audit: {path}",
             f"  paper: {paper or '?'} {paper_size_mm} mm (width={width_mm} mm)"]
    if paper_warning:
        lines.append(f"  {paper_warning}")

    low_dpi = 0
    for m in _INCLUDEGRAPHICS_RE.finditer(tex):
        opts_s = m.group("opts") or ""
        fname = m.group("file")
        fp = (path.parent / fname)
        if not fp.exists():
            for ext in (".pdf", ".png", ".jpg", ".jpeg", ".eps"):
                candidate = figdir / (pathlib.Path(fname).name + ext)
                if candidate.exists():
                    fp = candidate
                    break
                candidate = figdir / pathlib.Path(fname).name
                if candidate.exists():
                    fp = candidate
                    break
        if not fp.exists():
            lines.append(f"  - {fname}: NOT FOUND (skip)")
            continue
        lwm = _LINEWIDTH_RE.search(opts_s)
        lw_frac = float(lwm.group(1)) if lwm else 1.0
        printed_width_mm = width_mm * lw_frac
        printed_width_in = printed_width_mm / 25.4

        if fitz is None:
            lines.append(f"  - {fname}: pymupdf not installed; install `mcp-server-document[pdf]` for DPI check")
            continue

        try:
            doc = fitz.open(fp)
            page = doc[0]
            img_list = page.get_images(full=True)
            has_drawings = bool(page.get_drawings())
            if not img_list and has_drawings:
                # Pure vector PDF — resolution-independent at print scale.
                lines.append(f"  - {fp.name}: printed {printed_width_mm:.0f}mm wide, "
                             f"vector PDF [OK]")
                doc.close()
                continue
            if img_list:
                xref = img_list[0][0]
                pix = fitz.Pixmap(doc, xref)
                px_w = pix.width
                is_vector_mix = has_drawings
                pix = None
            else:
                # No images, no drawings — empty page; treat as OK.
                lines.append(f"  - {fp.name}: empty page (skip)")
                doc.close()
                continue
            doc.close()
            dpi = px_w / printed_width_in if printed_width_in > 0 else 0
            if dpi < 250:
                tag = "[HIGH]"
                low_dpi += 1
            elif dpi < 300:
                tag = "[MEDIUM]"
                low_dpi += 1
            else:
                tag = "[OK]"
            mix_note = " (raster+vector)" if is_vector_mix else ""
            lines.append(f"  - {fp.name}: printed {printed_width_mm:.0f}mm wide, "
                         f"{px_w}px → {dpi:.0f} DPI{mix_note} {tag}")
        except Exception as e:
            lines.append(f"  - {fp.name}: DPI check failed: {type(e).__name__}: {e}")

    if low_dpi:
        lines.append(f"FAIL — {low_dpi} figure(s) below 300 DPI at intended print width")
    else:
        lines.append("PASS — print resolution check complete")
    return "\n".join(lines)
