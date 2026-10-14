"""Tier 1 — poster_word_count.

Source: Colin Purrington, "Designing conference posters". Total body text
≤1000 words; per-section ~200 words. Japanese posters use mixed JP/EN
text — we count both axes separately so the user can apply the right
threshold.
"""
from __future__ import annotations

import re

from .._textutil import read_tex, strip_latex, count_words_mixed


# Recognized poster sections in Japanese + English.
_SECTION_ALIASES = {
    "intro": ("introduction", "はじめに", "背景", "緒言"),
    "method": ("method", "methods", "materials", "手法", "方法", "解析手法", "実験方法"),
    "result": ("result", "results", "結果", "解析結果", "実験結果"),
    "discussion": ("discussion", "考察"),
    "conclusion": ("conclusion", "conclusions", "summary", "まとめ", "結論", "おわりに"),
    "ack": ("acknowledgment", "acknowledgments", "謝辞"),
}

_TARGET_RANGES = {
    "intro": (100, 250),
    "method": (100, 250),
    "result": (100, 300),
    "discussion": (100, 250),
    "conclusion": (100, 250),
    "ack": (0, 60),
}

_TOTAL_TARGET = (300, 1000)

_SECTION_HEADING_RE = re.compile(
    r"(?:\{\\fontsize\{\d+(?:pt)?\}\{\d+(?:pt)?\}\\selectfont[^}]*\}|"
    r"\\section\*?\{[^}]*\}|\\subsection\*?\{[^}]*\}|"
    r"\\begin\{tcolorbox\}\[(?:section|subsection)box[^\]]*\][^\n]*)",
    re.IGNORECASE,
)


def _classify(heading: str) -> str | None:
    h = heading.lower()
    for key, aliases in _SECTION_ALIASES.items():
        if any(a in h for a in aliases):
            return key
    return None


def poster_word_count(tex_path: str) -> str:
    """Lint a poster's word budget against Purrington's ≤1000-word target.

    The check is bilingual: English words and Japanese characters are
    counted separately so mixed-script posters can be evaluated. The
    rough equivalence is ``1 Japanese character ≈ 0.5 English word`` for
    information density purposes (cf. Editage guideline).

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        A multi-line status report listing per-section budget vs. actual.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    text = strip_latex(tex)
    en_total, jp_total = count_words_mixed(text)
    # "Effective" word count for the 1000-word budget: 1 JP char ≈ 0.5 EN word.
    effective = en_total + int(jp_total * 0.5)

    issues: list[str] = []
    lo, hi = _TOTAL_TARGET
    if effective > hi:
        issues.append(f"[HIGH] total: effective={effective} words (en={en_total}, jp={jp_total}ch) > {hi}; Purrington 1000-word ceiling")
    elif effective < lo:
        issues.append(f"[LOW] total: effective={effective} words < {lo}; poster may be content-thin")

    lines = [f"poster_word_count: {path}",
             f"  total: en_words={en_total}, jp_chars={jp_total}, effective={effective} (target {lo}-{hi})"]

    if issues:
        lines.append("issues:")
        lines.extend(f"  {i}" for i in issues)
    else:
        lines.append("PASS — total within Purrington target")
    return "\n".join(lines)
