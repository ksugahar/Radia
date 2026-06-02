"""Tier 2 — poster_color_count_321.

Source: Editage『学会ポスター発表完全ガイド6』— 3-color rule with
70% base / 25% main / 5% accent split. Posters with > 5 unique
``\\definecolor`` declarations tend to look chaotic.
"""
from __future__ import annotations

from .._textutil import read_tex
from .._color import parse_definecolor


_MAX_COLORS = 5
_WARN_COLORS = 3


def poster_color_count_321(tex_path: str) -> str:
    """Count unique colors and warn if the palette violates the 3-color rule.

    The 3-color rule from Editage targets visual coherence at 1m reading
    distance: too many hues compete for attention.

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Palette list + verdict.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    palette = parse_definecolor(tex)
    n = len(palette)

    lines = [f"poster_color_count_321: {path}",
             f"  unique \\definecolor declarations: {n} (Editage target ≤{_WARN_COLORS}, hard ceiling {_MAX_COLORS})"]
    for name, (r, g, b) in sorted(palette.items()):
        lines.append(f"  - {name}: RGB({r:.2f}, {g:.2f}, {b:.2f})")

    if n > _MAX_COLORS:
        lines.append(f"[MEDIUM] {n} colors > {_MAX_COLORS}; palette risks visual chaos at 1m")
    elif n > _WARN_COLORS:
        lines.append(f"[LOW] {n} colors > {_WARN_COLORS}; consider consolidating to base/main/accent")
    else:
        lines.append("PASS — palette within 3-color rule")
    return "\n".join(lines)
