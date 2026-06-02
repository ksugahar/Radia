"""Tier 1 — poster_fontsize_by_distance.

Source: Editage『学会ポスター発表完全ガイド6』— viewing-distance formula:

    character_height_mm = viewing_distance_m / 250 * 1000   (≡ d/0.25)

so 5m → 20mm → ~57pt; 3m → 12mm → ~34pt; 2m → 8mm → ~23pt; 1m → 4mm → ~11pt.
Since posters are read at 1-3m (body) and 5m+ (title), we apply the
formula per element class.

The tool ingests the user's intended viewing distance per element and
checks that ``\\fontsize`` declarations are large enough.
"""
from __future__ import annotations

import re

from .._textutil import read_tex


_FONTSIZE_RE = re.compile(r"\\fontsize\{(\d+)(?:pt)?\}\{(\d+)(?:pt)?\}")

# 1 pt ≈ 0.3528 mm; required pt = required_mm / 0.3528
_PT_PER_MM = 1 / 0.3528


def _required_pt(distance_m: float, safety: float = 1.0) -> int:
    """Required font height in pt for the given viewing distance.

    Default safety factor 1.0 follows the Editage formula. Use 1.2 if
    you want a comfortable margin.
    """
    h_mm = distance_m / 0.25 * safety
    return int(round(h_mm * _PT_PER_MM))


# Default expected viewing distances per element class.
_DEFAULTS = {
    "title": 5.0,
    "subtitle": 3.0,
    "heading": 3.0,
    "body": 2.0,
    "caption": 1.5,
    "refs": 1.0,
}


def poster_fontsize_by_distance(tex_path: str,
                                distance_overrides: dict[str, float] | None = None) -> str:
    """Verify that ``\\fontsize{X}`` values are large enough for their role.

    Element class is inferred from the surrounding context:
    - largest font ≥50pt → title
    - 35-49pt → subtitle / heading
    - 24-34pt → body / caption
    - 18-23pt → refs

    Then each class is checked against the required pt for its viewing
    distance (Editage formula: ``h_mm = d_m / 0.25``).

    Args:
        tex_path: Path to the poster ``.tex`` source.
        distance_overrides: Optional ``{"body": 1.5, ...}`` to override
                            default viewing distances.

    Returns:
        Per-class required vs. observed pt, with HIGH severity if observed < required.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    distances = dict(_DEFAULTS)
    if distance_overrides:
        distances.update(distance_overrides)

    sizes = [int(a) for a, _ in _FONTSIZE_RE.findall(tex)]
    if not sizes:
        return f"poster_fontsize_by_distance: {path}\n  no \\fontsize declarations found"

    title = max(sizes)
    body_candidates = [s for s in sizes if 22 <= s <= 30]
    body = min(body_candidates) if body_candidates else None
    ref_candidates = [s for s in sizes if 16 <= s <= 22]
    refs = min(ref_candidates) if ref_candidates else None
    heading_candidates = [s for s in sizes if 30 < s < title]
    heading = min(heading_candidates) if heading_candidates else None

    rows = [
        ("title", title, distances["title"]),
        ("heading", heading, distances["heading"]),
        ("body", body, distances["body"]),
        ("refs", refs, distances["refs"]),
    ]

    lines = [f"poster_fontsize_by_distance: {path}",
             "  Editage formula: required_pt = (viewing_m / 0.25) * 2.835"]
    issues = 0
    for role, observed, dist in rows:
        if observed is None:
            lines.append(f"  - {role}: no candidate font detected (skipped)")
            continue
        required = _required_pt(dist)
        verdict = "OK " if observed >= required else "FAIL"
        if verdict == "FAIL":
            issues += 1
        lines.append(f"  - {role}: observed={observed}pt, viewing={dist}m → required≥{required}pt [{verdict}]")

    if issues:
        lines.append(f"[HIGH] {issues} class(es) below required pt for declared viewing distance")
    else:
        lines.append("PASS — all element classes meet viewing-distance requirement")
    return "\n".join(lines)
