"""Tier 3 — pdf_health_report.

Composite check for journal submission readiness:
    - All fonts embedded (HIGH if not)
    - Metadata set (LOW if Title/Author empty)
    - File size within typical journal cap (10 MB warn, 25 MB fail)
    - Page count reasonable (>100 = thesis, not paper)
"""
from __future__ import annotations

import pathlib

from .T1_font_embed import pdf_font_embed_check
from .T2_inspect import pdf_inspect
from .T6_metadata import pdf_get_metadata


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_health_report(pdf_path: str,
                      max_size_mb: float = 10.0) -> str:
    """Submission-readiness audit for journal / conference PDF upload.

    Args:
        pdf_path: Source PDF.
        max_size_mb: Hard size cap typical for the target venue
                     (default 10 MB — many IEEE / Elsevier journals).

    Returns:
        Per-axis verdict + overall score.
    """
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."

    issues: list[tuple[str, str]] = []

    # Size
    size_mb = p.stat().st_size / (1024 * 1024)
    if size_mb > max_size_mb:
        issues.append(("HIGH", f"file size {size_mb:.1f} MB > cap {max_size_mb} MB"))
    elif size_mb > 0.8 * max_size_mb:
        issues.append(("LOW", f"file size {size_mb:.1f} MB is close to {max_size_mb} MB cap"))

    # Font embed
    fe = pdf_font_embed_check(str(p))
    if "FAIL" in fe:
        issues.append(("HIGH", "fonts not fully embedded — many submission systems reject"))

    # Metadata
    md = pdf_get_metadata(str(p))
    md_low = md.lower()
    if "title          :" in md_low or "title           :" in md_low:
        for line in md.splitlines():
            if line.strip().startswith("title"):
                if line.split(":", 1)[-1].strip() == "":
                    issues.append(("LOW", "metadata Title is empty"))
            if line.strip().startswith("author"):
                if line.split(":", 1)[-1].strip() == "":
                    issues.append(("LOW", "metadata Author is empty"))

    # Page count sanity
    doc = fitz.open(p)
    n = len(doc)
    doc.close()
    if n > 100:
        issues.append(("LOW", f"{n} pages — looks like a thesis, not a paper"))
    if n == 0:
        issues.append(("HIGH", "PDF has 0 pages"))

    severities = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
    for sev, _ in issues:
        severities[sev] += 1
    penalty = severities["HIGH"] * 15 + severities["MEDIUM"] * 6 + severities["LOW"] * 2
    score = max(0, 100 - penalty)

    lines = [f"pdf_health_report: {p}",
             f"  size: {size_mb:.2f} MB, pages: {n}",
             f"  composite score: {score}/100"]
    if not issues:
        lines.append("PASS — ready for journal submission")
    else:
        lines.append(f"FAIL — {len(issues)} issue(s):")
        for sev, msg in issues:
            lines.append(f"  [{sev}] {msg}")
    return "\n".join(lines)
