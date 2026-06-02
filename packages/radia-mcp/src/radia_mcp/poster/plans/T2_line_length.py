"""Tier 1 — poster_line_length.

In a LaTeX source the typeset line length is determined by multicol /
column-width, not by line breaks in the source. So the *meaningful*
lint here is sentence length, which tracks comprehension load at 1m
reading distance.

Source: Purrington (45-65 char English lines) translated to sentence
length, plus Editage observation that Japanese sentences over ~80
characters cause comprehension drop on posters.
"""
from __future__ import annotations

import re

from .._textutil import read_tex, strip_latex, count_japanese_chars


_REFERENCE_BLOCK = re.compile(
    r"\\begin\{enumerate\}\[label=\{\[\\arabic\*\]\}.*?\\end\{enumerate\}",
    re.DOTALL,
)

_EN_HI = 35  # too-long English sentence in WORDS (≈ 200 chars)
_JP_HI = 80  # too-long Japanese sentence in CHARACTERS


def poster_line_length(tex_path: str) -> str:
    """Flag sentences that are too long for poster reading.

    Splits on Japanese 。 and English ``. !? ``. Each fragment longer
    than 35 English words or 80 Japanese characters is flagged. Short
    fragments (≤10 chars) are ignored — they're titles / captions /
    bullet labels.

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Sentence-length distribution + flagged examples.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    text = strip_latex(_REFERENCE_BLOCK.sub("", tex))
    # Split into sentences.
    sentences = [s.strip() for s in re.split(r"(?<=[。.！？!?])\s*", text)]
    sentences = [s for s in sentences if len(s) > 10]

    en_too_long = []
    jp_too_long = []
    for s in sentences:
        jp = count_japanese_chars(s)
        if jp >= len(s) / 3:
            if jp > _JP_HI:
                jp_too_long.append(s)
        else:
            words = len([w for w in s.split() if any(c.isalpha() for c in w)])
            if words > _EN_HI:
                en_too_long.append(s)

    total = len(sentences)
    bad = len(en_too_long) + len(jp_too_long)
    lines = [f"poster_line_length: {path}",
             f"  sentences analyzed: {total}",
             f"  too long: en={len(en_too_long)}, jp={len(jp_too_long)} "
             f"(thresholds: >{_EN_HI} en words / >{_JP_HI} jp chars)"]

    if total == 0:
        lines.append("  no analyzable sentences — template skeleton")
        return "\n".join(lines)

    ratio = bad / total
    if ratio >= 0.3:
        lines.append(f"[HIGH] {ratio:.0%} of sentences run long — split or compress")
    elif ratio >= 0.15:
        lines.append(f"[MEDIUM] {ratio:.0%} of sentences run long")
    else:
        lines.append("PASS — sentence-length distribution is healthy")

    for s in (en_too_long + jp_too_long)[:3]:
        lines.append(f"  - {s[:80]}{'…' if len(s) > 80 else ''}")
    return "\n".join(lines)
