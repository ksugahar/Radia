"""Locate the canonical bibliography bundled with ``radia-mcp``.

The public package keeps one ``references.bib`` data file so paper-writing
checks have a deterministic default in source trees and installed wheels.
Callers may still pass an explicit bibliography when reviewing another
project.
"""
from __future__ import annotations

from pathlib import Path

from .._bibparse import read_bib_file

CANONICAL = Path(__file__).resolve().parents[1] / "data" / "references.bib"


def bibliography_canonical_path() -> str:
    """Return the bundled canonical bibliography path and entry count."""
    if not CANONICAL.is_file():
        return f"canonical bibliography not found at {CANONICAL}"
    entries = [entry for entry in read_bib_file(CANONICAL) if not entry.kind.startswith("@")]
    return (
        "bibliography_canonical_path\n"
        f"  {CANONICAL}\n"
        f"  {len(entries)} entries, {CANONICAL.stat().st_size / 1024:.0f} KB"
    )
