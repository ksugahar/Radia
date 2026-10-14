"""Tier 6 — poster_next_5_actions.

Re-sort the lint findings by impact × effort and return the top 5
concrete action items with rewrite examples.
"""
from __future__ import annotations

import re

from .T1_word_count import poster_word_count
from .T2_line_length import poster_line_length
from .T3_caption_self_contained import poster_caption_self_contained
from .T4_jp_font_check import poster_jp_font_check
from .T5_fontsize_by_distance import poster_fontsize_by_distance
from .T6_typography_lints import poster_typography_lints
from .T7_color_contrast_wcag import poster_color_contrast_wcag
from .T8_colorblind_hint import poster_colorblind_hint
from .T9_color_count_321 import poster_color_count_321


# impact (1-5), effort (1-5: 1=trivial, 5=heavy rework), example fix
_ACTION_LIBRARY = {
    "fontsize_by_distance.HIGH":     (5, 2, "全 \\fontsize を Editage 式 (h_mm = d_m/0.25) で逆算し直し、タイトル 60pt / 本文 24pt 以上に揃える"),
    "color_contrast.HIGH":           (5, 1, "text/bg ペアで WCAG <3:1 のものを黒+白 or 明色+濃色に置換 (\\definecolor)"),
    "word_count.HIGH":               (4, 3, "本文 1000 語超過 → 各 contentbox を 200 語に圧縮、長文は caption に逃がす"),
    "caption_self_contained.MEDIUM": (4, 1, "図 caption に「軸+単位+凡例」を 1 行で書き足す"),
    "jp_font.HIGH":                  (4, 1, "\\renewcommand{\\familydefault}{\\sfdefault} と {\\kanjifamilydefault}{\\gtdefault} を preamble に追加"),
    "colorblind.MEDIUM":             (3, 2, "赤×緑ペアを 青×橙 に置換、または線種 (実線/破線) で 2 重符号化"),
    "color_count.MEDIUM":            (3, 1, "\\definecolor を 3 色に削減 (base / main / accent)"),
    "line_length.HIGH":              (3, 3, "multicols の段数を 1 列増やすか \\columnsep を広げて 1 行 22-32 字に整える"),
    "typography.LOW":                (2, 1, "下線→\\textbf、ALL CAPS→sentence case、ダブルスペース → 単一スペースに置換"),
}


def _scan(report: str, axis_name: str):
    sev = {sev: len(re.findall(rf"\[{sev}\]", report)) for sev in ("HIGH", "MEDIUM", "LOW")}
    return [(axis_name, lvl) for lvl, n in sev.items() for _ in range(n)]


def poster_next_5_actions(tex_path: str) -> str:
    """Return the top 5 actions ranked by (impact - 0.5*effort).

    Each entry includes a concrete LaTeX-level fix example so the user
    can act immediately. Same pattern as
    ``grant_writing_next_5_actions``.
    """
    reports = [
        ("fontsize_by_distance", poster_fontsize_by_distance(tex_path)),
        ("color_contrast", poster_color_contrast_wcag(tex_path)),
        ("word_count", poster_word_count(tex_path)),
        ("caption_self_contained", poster_caption_self_contained(tex_path)),
        ("jp_font", poster_jp_font_check(tex_path)),
        ("colorblind", poster_colorblind_hint(tex_path)),
        ("color_count", poster_color_count_321(tex_path)),
        ("line_length", poster_line_length(tex_path)),
        ("typography", poster_typography_lints(tex_path)),
    ]

    issues: list[tuple[str, str]] = []
    for axis, rep in reports:
        issues.extend(_scan(rep, axis))

    # Build scored items.
    scored = []
    for axis, level in issues:
        key = f"{axis}.{level}"
        if key in _ACTION_LIBRARY:
            impact, effort, example = _ACTION_LIBRARY[key]
            score = impact - 0.5 * effort
            scored.append((score, axis, level, impact, effort, example))

    # Sort by score, then take unique (axis,level) up to 5.
    scored.sort(key=lambda x: -x[0])
    seen = set()
    top: list[tuple] = []
    for entry in scored:
        sig = (entry[1], entry[2])
        if sig in seen:
            continue
        seen.add(sig)
        top.append(entry)
        if len(top) == 5:
            break

    lines = [f"poster_next_5_actions: {tex_path}"]
    if not top:
        lines.append("  no impactful actions — health is high")
        return "\n".join(lines)
    for i, (score, axis, level, impact, effort, example) in enumerate(top, 1):
        lines.append(f"  {i}. [{level}] {axis}  (impact={impact}, effort={effort}, score={score:.1f})")
        lines.append(f"     → {example}")
    return "\n".join(lines)
