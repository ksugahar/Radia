"""Tier 2 — poster_color_contrast_wcag.

Source: WCAG 2.1 AA — 3:1 for ≥18pt text or ≥14pt bold; 4.5:1 for
smaller. We scan for ``colback=X`` blocks and find the ``coltext=Y``
(or ``\\color{Y}``) declared in the same tcolorbox — that's the actual
text/bg pair, no brute-force cross product.
"""
from __future__ import annotations

import re

from .._textutil import read_tex
from .._color import parse_definecolor, contrast_ratio


_TCBSET_PAIR_RE = re.compile(
    r"colback\s*=\s*(?P<bg>[A-Za-z0-9!]+)[^}]*?coltext\s*=\s*(?P<text>[A-Za-z0-9!]+)",
    re.DOTALL,
)
_TCB_INLINE_RE = re.compile(
    r"\\begin\{tcolorbox\}\[[^\]]*colback\s*=\s*(?P<bg>[A-Za-z0-9!]+)[^\]]*\]"
    r"[\s\S]{0,400}?\\color\{(?P<text>[A-Za-z0-9!]+)\}",
)

_AA_LARGE = 3.0
_AA_SMALL = 4.5


def _normalize(name: str) -> str:
    return name.split("!")[0].strip()


def poster_color_contrast_wcag(tex_path: str) -> str:
    """Check WCAG 2.1 AA contrast for *actually-paired* text/bg combinations.

    Pair-discovery rules:
      - ``tcbset{...}`` blocks with both ``colback=`` and ``coltext=``
      - ``\\begin{tcolorbox}[colback=X]`` followed within 400 chars by
        ``\\color{Y}``

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Per-pair contrast ratio + verdict.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    palette = parse_definecolor(tex)
    palette.setdefault("black", (0, 0, 0))
    palette.setdefault("white", (1, 1, 1))

    pairs: set[tuple[str, str]] = set()
    for m in _TCBSET_PAIR_RE.finditer(tex):
        pairs.add((_normalize(m.group("text")), _normalize(m.group("bg"))))
    for m in _TCB_INLINE_RE.finditer(tex):
        pairs.add((_normalize(m.group("text")), _normalize(m.group("bg"))))

    # Always include default: black-on-white body text.
    pairs.add(("black", "white"))

    lines = [f"poster_color_contrast_wcag: {path}",
             f"  palette: {len(palette)} colors; paired text/bg combinations: {len(pairs)}"]
    failures = 0
    for tc, bg in sorted(pairs):
        if tc not in palette or bg not in palette:
            lines.append(f"  - {tc} on {bg}: color not defined (skip)")
            continue
        if tc == bg:
            continue
        cr = contrast_ratio(palette[tc], palette[bg])
        verdict = "OK " if cr >= _AA_LARGE else "FAIL"
        if verdict == "FAIL":
            failures += 1
        note = ""
        if cr < _AA_LARGE:
            note = f" (<{_AA_LARGE}, fails AA-Large)"
        elif cr < _AA_SMALL:
            note = f" (≥{_AA_LARGE} OK for ≥18pt but <{_AA_SMALL} for body)"
        lines.append(f"  [{verdict}] {tc} on {bg}: ratio={cr:.2f}{note}")
    if failures:
        lines.append(f"[HIGH] {failures} text/bg pair(s) fail WCAG AA-Large 3:1")
    else:
        lines.append("PASS — all paired text/bg combinations meet WCAG AA-Large")
    return "\n".join(lines)
