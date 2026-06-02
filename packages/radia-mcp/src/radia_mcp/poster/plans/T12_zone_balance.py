"""Tier 3 — poster_zone_balance_check.

Validates the three-zone (left bar / center billboard / right bar) area
split of a Morrison-style betterposter. Target proportions: 25 / 50 / 25
(±10% tolerance).
"""
from __future__ import annotations

import re

from .._textutil import read_tex


_MINIPAGE_RE = re.compile(r"\\begin\{minipage\}\[t\]\{([\d.]+)\\linewidth\}")


def poster_zone_balance_check(tex_path: str) -> str:
    """Check that minipage column widths approximate 0.25 / 0.50 / 0.25.

    Args:
        tex_path: Path to a betterposter ``.tex`` source.

    Returns:
        Detected zone widths + verdict.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    widths = [float(m.group(1)) for m in _MINIPAGE_RE.finditer(tex)]
    lines = [f"poster_zone_balance_check: {path}",
             f"  detected minipage widths: {widths}"]
    if len(widths) != 3:
        lines.append(f"[HIGH] expected 3 minipage zones, found {len(widths)}")
        return "\n".join(lines)

    target = (0.25, 0.50, 0.25)
    diffs = [abs(w - t) for w, t in zip(widths, target)]
    bad = [(i, w, t, d) for i, (w, t, d) in enumerate(zip(widths, target, diffs)) if d > 0.05]
    if bad:
        lines.append("[MEDIUM] zone(s) more than ±0.05 off Morrison 25/50/25 target:")
        for i, w, t, d in bad:
            zone = ["left", "center", "right"][i]
            lines.append(f"  - {zone}: {w:.2f} vs target {t:.2f} (Δ={d:+.2f})")
    else:
        lines.append("PASS — three zones approximate Morrison 25/50/25 split")
    return "\n".join(lines)
