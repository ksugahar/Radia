"""Tier 2 — pdf_metadata get/set.

PDF document metadata (Title / Author / Subject / Keywords / Producer)
matters for journal submission systems — many require non-empty fields.
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_get_metadata(pdf_path: str) -> str:
    """Print the PDF's metadata dictionary."""
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    p = pathlib.Path(pdf_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    doc = fitz.open(p)
    meta = doc.metadata or {}
    doc.close()
    lines = [f"pdf_get_metadata: {p}"]
    for k in ("title", "author", "subject", "keywords", "creator",
              "producer", "creationDate", "modDate", "format", "encryption"):
        v = meta.get(k, "")
        lines.append(f"  {k:15s}: {v}")
    return "\n".join(lines)


def pdf_set_metadata(pdf_path: str,
                     title: str | None = None,
                     author: str | None = None,
                     subject: str | None = None,
                     keywords: str | None = None,
                     out_path: str | None = None) -> str:
    """Set PDF document metadata fields.

    Args:
        pdf_path: Source PDF.
        title / author / subject / keywords: Optional overrides.
        out_path: Destination. Defaults to in-place rewrite.

    Returns:
        Status string with the new metadata.
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
    meta = dict(doc.metadata or {})
    if title is not None:
        meta["title"] = title
    if author is not None:
        meta["author"] = author
    if subject is not None:
        meta["subject"] = subject
    if keywords is not None:
        meta["keywords"] = keywords
    doc.set_metadata(meta)
    doc.save(out, incremental=(out == p), encryption=fitz.PDF_ENCRYPT_KEEP)
    doc.close()
    lines = [f"pdf_set_metadata: wrote {out}"]
    for k in ("title", "author", "subject", "keywords"):
        v = meta.get(k, "")
        lines.append(f"  {k}: {v}")
    return "\n".join(lines)
