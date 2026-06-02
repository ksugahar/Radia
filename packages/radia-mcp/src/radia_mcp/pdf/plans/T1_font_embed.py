"""Tier 1 — pdf_font_embed_check.

Wrapper around ``pdffonts`` (poppler / xpdf-utils) — verifies all
fonts are embedded (subset is acceptable). Used at submission time for
journals and at print time for posters.

Generalization of the poster-specific ``poster_font_embed_check``.
"""
from __future__ import annotations

import pathlib
import shutil
import subprocess


def pdf_font_embed_check(pdf_path: str) -> str:
    """Verify all fonts are embedded in the PDF.

    Args:
        pdf_path: Path to the PDF to check.

    Returns:
        pdffonts output + verdict.
    """
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    if not shutil.which("pdffonts"):
        return ("Error: pdffonts not on PATH. Install TeXLive (bundles "
                "poppler / pdffonts) or xpdf-utils.")
    try:
        proc = subprocess.run(
            ["pdffonts", str(p)], capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=60,
        )
    except subprocess.TimeoutExpired:
        return "Error: pdffonts timed out after 60s."
    if proc.returncode != 0:
        return f"Error: pdffonts returned {proc.returncode}: {proc.stderr.strip()}"
    lines_out = proc.stdout.splitlines()
    if len(lines_out) < 3:
        return f"pdf_font_embed_check: {p}\n  no fonts found in PDF (empty?)"
    # pdffonts output uses fixed columns delimited by '---' on the separator
    # line. Parse the separator to get exact column spans, since the "type"
    # field can contain spaces (e.g. "Type 1C") and break naive split().
    sep_line = lines_out[1] if len(lines_out) > 1 else ""
    # Find span ranges of '---' groups.
    spans = []
    i = 0
    while i < len(sep_line):
        if sep_line[i] == "-":
            j = i
            while j < len(sep_line) and sep_line[j] == "-":
                j += 1
            spans.append((i, j))
            i = j
        else:
            i += 1
    header = lines_out[0] if lines_out else ""
    columns = [header[a:b].strip() for a, b in spans]
    try:
        emb_idx = columns.index("emb")
    except ValueError:
        # Header layout changed; fall back to position-based guess.
        emb_idx = 3

    not_embedded = []
    for row in lines_out[2:]:
        if len(spans) <= emb_idx:
            continue
        a, b = spans[emb_idx]
        emb = row[a:b].strip().lower() if len(row) >= b else ""
        if emb != "yes":
            not_embedded.append(row)
    rep = [f"pdf_font_embed_check: {p}"]
    rep.extend("  " + ln for ln in lines_out[:min(20, len(lines_out))])
    if not_embedded:
        rep.append(f"FAIL — {len(not_embedded)} font(s) NOT embedded:")
        rep.extend(f"  {r}" for r in not_embedded)
        rep.append("  remedy: set EmbedAllFonts=true in updmap.cfg / -dPDFSETTINGS=/prepress")
    else:
        rep.append("PASS — all fonts embedded")
    return "\n".join(rep)
