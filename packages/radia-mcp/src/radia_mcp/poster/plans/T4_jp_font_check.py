"""Tier 1 — poster_jp_font_check.

Source: Editage. Gothic for Japanese body, sans-serif for English.
Mincho/serif body text is unreadable at 1m. Emphasis should use bold,
not font switching.
"""
from __future__ import annotations

import re

from .._textutil import read_tex


# Family selectors that indicate Gothic / sans-serif (acceptable).
_GOTHIC_TOKENS = (
    "\\gtfamily", "\\sffamily", "\\sfdefault",
    "meiryo", "メイリオ", "yu gothic", "游ゴシック", "ヒラギノ角ゴ",
    "noto sans cjk jp", "noto sans jp", "ipaexgothic", "ipag",
    "ms gothic", "ms ゴシック",
)
# Family selectors that indicate Mincho / serif (problematic for body).
_MINCHO_TOKENS = (
    "\\mcfamily", "\\rmfamily",  # rmfamily defaults to serif under bxjsarticle
    "ipaexmincho", "ipam", "yu mincho", "游明朝", "ヒラギノ明朝",
    "ms mincho", "ms 明朝", "noto serif cjk jp",
)
# Acceptable English serif fonts that often appear in newtxtext setup
# but should NOT be the body family if the poster is mostly Japanese.
_EN_SERIF_PACKAGES = ("newtxtext", "mathptmx", "times")


def poster_jp_font_check(tex_path: str) -> str:
    """Inspect Japanese font family declarations in a poster .tex.

    The check follows Editage's guidance:

    - Body text should resolve to a Gothic / sans-serif Japanese family.
    - Mincho / serif Japanese fonts are unreadable for body text at 1m.
    - Emphasis should rely on ``\\textbf`` / ``\\bfseries``, not font
      switching to a different family.

    Heuristics:
        - bxjsarticle defaults to Mincho unless ``\\renewcommand{\\familydefault}{\\sfdefault}``
          or ``\\renewcommand{\\kanjifamilydefault}{\\gtdefault}`` is set.
        - Explicit ``{\\mcfamily ...}`` blocks in body region are warned.

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Status with family default detection + per-emphasis findings.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    low = tex.lower()
    has_sf_default = ("\\renewcommand{\\familydefault}{\\sfdefault}" in tex or
                       "\\renewcommand\\familydefault{\\sfdefault}" in tex)
    has_gt_kanji = ("\\renewcommand{\\kanjifamilydefault}{\\gtdefault}" in tex or
                     "\\renewcommand\\kanjifamilydefault{\\gtdefault}" in tex)

    gothic_hits = [t for t in _GOTHIC_TOKENS if t.lower() in low]
    mincho_hits = [t for t in _MINCHO_TOKENS if t.lower() in low]

    issues: list[str] = []
    if not (has_sf_default or has_gt_kanji or gothic_hits):
        issues.append("[HIGH] no Gothic / sans-serif default detected; "
                      "bxjsarticle defaults to Mincho/Serif body — unreadable at 1m")
    if mincho_hits:
        # Allow rmfamily only inside narrow regions (math, references).
        body_mincho = [t for t in mincho_hits if t not in ("\\rmfamily",)]
        if body_mincho:
            issues.append(f"[MEDIUM] explicit Mincho/serif tokens used in body: {body_mincho}")

    # Emphasis-by-font-switch: count {\rmfamily ...} or {\mcfamily ...} occurrences
    # near short text (heuristic for emphasis).
    switch_pat = re.compile(r"\{\\(rmfamily|mcfamily)\s+([^{}]{1,40})\}")
    switches = switch_pat.findall(tex)
    if switches:
        issues.append(f"[LOW] {len(switches)} font-family switch(es) used for emphasis; prefer \\textbf")

    lines = [f"poster_jp_font_check: {path}",
             f"  sfdefault set: {has_sf_default}",
             f"  gtdefault for kanji: {has_gt_kanji}",
             f"  gothic tokens seen: {gothic_hits or 'none'}",
             f"  mincho tokens seen: {mincho_hits or 'none'}"]
    if issues:
        lines.append("issues:")
        lines.extend(f"  {i}" for i in issues)
    else:
        lines.append("PASS — Japanese body resolves to Gothic/sans-serif")
    return "\n".join(lines)
