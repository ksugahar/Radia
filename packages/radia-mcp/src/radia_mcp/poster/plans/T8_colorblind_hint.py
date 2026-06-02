"""Tier 2 — poster_colorblind_hint.

Source: WCAG 2.1 — don't rely on color alone; ~5% of adults have CVD.
Editage and ACS Chemical Health & Safety both flag red/green pairs as
the canonical hazard.
"""
from __future__ import annotations

from .._textutil import read_tex
from .._color import parse_definecolor, simulate_deuteranopia, color_distance


# Threshold below which two colors collapse under deuteranopia
# simulation; chosen so red-vs-green pairs fail and blue-vs-orange pass.
_COLLAPSE_THRESHOLD = 0.15


def poster_colorblind_hint(tex_path: str) -> str:
    """Simulate deuteranopia and flag pairs that collapse to similar colors.

    For each pair of named colors in ``\\definecolor`` palette, apply a
    Brettel-Viénot-style channel mix and measure the Euclidean distance.
    Pairs with simulated distance < 0.15 are flagged as a CVD hazard.

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Report of collapsed pairs + recommendations.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    palette = parse_definecolor(tex)
    if len(palette) < 2:
        return f"poster_colorblind_hint: {path}\n  palette too small to check CVD pairs"

    items = sorted(palette.items())
    hazards = []
    for i, (n1, c1) in enumerate(items):
        for n2, c2 in items[i + 1:]:
            normal_d = color_distance(c1, c2)
            if normal_d < 0.1:
                # Pair already nearly identical without CVD; skip.
                continue
            sim_d = color_distance(simulate_deuteranopia(c1), simulate_deuteranopia(c2))
            if sim_d < _COLLAPSE_THRESHOLD:
                hazards.append((n1, n2, normal_d, sim_d))

    lines = [f"poster_colorblind_hint: {path}",
             f"  palette: {len(palette)} colors checked pairwise"]
    if not hazards:
        lines.append("PASS — no hazardous CVD pairs detected")
        return "\n".join(lines)

    lines.append(f"[MEDIUM] {len(hazards)} pair(s) collapse under deuteranopia simulation:")
    for n1, n2, nd, sd in hazards:
        lines.append(f"  - {n1} ↔ {n2}: normal-d={nd:.2f}, CVD-d={sd:.2f} (<{_COLLAPSE_THRESHOLD})")
    lines.append("  remedy: pair colors with shape/pattern/dash style, not color alone")
    return "\n".join(lines)
