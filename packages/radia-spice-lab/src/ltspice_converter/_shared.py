"""Shared paths and helpers for ltspice_converter.

The training datasets (textbook DB, LTspice bundled examples, third-party
GitHub repos) are intentionally NOT shipped with this package -- only the
trained converter itself and a handful of author-authored test fixtures.
See README.md for the copyright rationale.
"""
from __future__ import annotations

import os
from pathlib import Path

_PACKAGE_DIR = Path(__file__).parent              # src/ltspice_converter/
_PARSER_DIR = _PACKAGE_DIR / "parser"


def decode_asc_bytes(raw: bytes) -> str:
    """Decode .asc / .net bytes to text, covering every encoding LTspice
    has shipped over the years.

    LTspice files appear in the wild as:
      - UTF-16 LE with BOM (LTspice 17+, the common modern case)
      - UTF-16 LE/BE *without* a BOM (some bundled symbol-test schematics
        such as ``UniversalOpAmp*.asc`` -- the case that motivated this)
      - UTF-8 (with or without BOM)
      - Windows-1252 / Latin-1 (older 8-bit; e.g. a bare 0xB5 for µ)

    The historical `raw[:2] == b"\\xff\\xfe"` check only caught BOM-marked
    UTF-16 LE; a BOM-less UTF-16 file decoded as UTF-8 yields text riddled
    with NUL characters, so `Version` / `SHEET` / `WIRE` never match and
    the schematic reads as empty.  This helper detects that case by the
    tell-tale density of NUL bytes (ASCII text in UTF-16 has a NUL in
    every other byte) and is shared by the CLI reader and the parser so
    both agree.
    """
    if raw[:2] == b"\xff\xfe":
        return raw[2:].decode("utf-16-le", errors="replace")
    if raw[:2] == b"\xfe\xff":
        return raw[2:].decode("utf-16-be", errors="replace")
    if raw[:3] == b"\xef\xbb\xbf":
        return raw[3:].decode("utf-8", errors="replace")
    # BOM-less UTF-16: ASCII-range content carries a NUL in every other
    # byte, so a quarter-or-more NUL density is a reliable tell.
    head = raw[:400]
    if head.count(0) > len(head) // 4:
        if head[1::2].count(0) >= head[0::2].count(0):
            return raw.decode("utf-16-le", errors="replace")
        return raw.decode("utf-16-be", errors="replace")
    for enc in ("utf-8", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def _find_ltspice() -> str | None:
    """Auto-detect LTspice executable (Windows). Returns None on non-Windows
    or if LTspice is not installed.
    """
    candidates = [
        Path(os.environ.get("PROGRAMFILES", "")) / "ADI" / "LTspice" / "LTspice.exe",
        Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "ADI" / "LTspice" / "LTspice.exe",
        Path(os.environ.get("PROGRAMFILES", "")) / "LTC" / "LTspiceXVII" / "XVIIx64.exe",
    ]
    for p in candidates:
        if p.exists():
            return str(p)
    return None
