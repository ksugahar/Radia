"""Tier 5 — poster_font_embed_check.

Inspects a compiled PDF to verify all fonts are embedded (subset is
acceptable). Print shops reject PDFs with non-embedded fonts.
"""
from __future__ import annotations

import pathlib
import shutil
import subprocess


def poster_font_embed_check(pdf_path: str) -> str:
    """Run ``pdffonts`` and report any non-embedded font.

    ``pdffonts`` ships with TeXLive's poppler bundle and Xpdf. The
    ``emb`` and ``sub`` columns must both be ``yes`` for print
    safety.

    Args:
        pdf_path: Path to the compiled poster PDF.

    Returns:
        Tabular pdffonts output + verdict.
    """
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    if not shutil.which("pdffonts"):
        return ("Error: pdffonts not on PATH. Install TeXLive (which "
                "bundles poppler / pdffonts) or xpdf-utils.")
    try:
        proc = subprocess.run(
            ["pdffonts", str(p)], capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=60,
        )
    except subprocess.TimeoutExpired:
        return "Error: pdffonts timed out after 60s."
    if proc.returncode != 0:
        return f"Error: pdffonts returned {proc.returncode}: {proc.stderr.strip()}"

    out = proc.stdout
    lines_out = out.splitlines()
    if len(lines_out) < 3:
        return f"poster_font_embed_check: {p}\n  pdffonts produced no font rows; PDF may be empty"

    # First 2 lines are header + separator. Each row: name type encoding emb sub uni object ID
    not_embedded = []
    for row in lines_out[2:]:
        cols = row.split()
        if len(cols) < 5:
            continue
        emb = cols[3]
        if emb.lower() != "yes":
            not_embedded.append(row)

    rep = [f"poster_font_embed_check: {p}", "  --- pdffonts output ---"]
    rep.extend("  " + ln for ln in lines_out[:2 + len(lines_out)])
    if not_embedded:
        rep.append(f"FAIL — {len(not_embedded)} font(s) NOT embedded:")
        rep.extend(f"  {r}" for r in not_embedded)
        rep.append("  remedy: pdflatex -dPDFSETTINGS=/prepress, or set EmbedAllFonts=true in updmap.cfg")
    else:
        rep.append("PASS — all fonts embedded")
    return "\n".join(rep)
