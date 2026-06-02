"""Tier 3 — poster_betterposter_billboard_lint.

For a Morrison-style poster, the central billboard's *plain-language*
quality determines the whole format's payoff (eye-tracking pilot showed
~10× faster comprehension when the billboard is genuinely plain).
"""
from __future__ import annotations

import re

from .._textutil import read_tex


_JARGON_HINTS = (
    # Engineering/physics jargon that hints at over-technicality.
    "PDE", "FEM", "FEA", "BVP", "FFT", "DFT", "ODE",
    "tensor", "manifold", "eigenvalue", "covariant", "contravariant",
    "テンソル", "多様体", "固有値", "共形変換", "ヤコビアン", "ラプラシアン",
    "DCT", "STFT", "ICA", "PCA", "VAE", "GAN", "MCMC", "Transformer",
)


def poster_betterposter_billboard_lint(tex_path: str) -> str:
    """Lint the central billboard text for plain-language compliance.

    Rules (from Morrison 2019 + ScienceUX eye-tracking findings):

    - Billboard ≥120pt (else not a "billboard")
    - ≤15 words / ≤30 Japanese characters
    - ≤2 jargon tokens (graduate-level terms)
    - No bullet markers (it's a sentence, not a list)
    - No reference markers like ``[1]``

    Args:
        tex_path: Path to a betterposter ``.tex`` source.

    Returns:
        Per-rule verdict + actionable fixes.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    # Find the largest \fontsize block — that's the billboard. We want
    # the actual main finding, not a sidebar heading.
    candidates = list(re.finditer(
        r"\\fontsize\{(?P<size>\d{2,3})\}\{\d+\}\\selectfont\\bfseries\s+(?P<body>[\s\S]+?)(?:\\\\\[|\\end\{)",
        tex,
    ))
    if not candidates:
        return (f"poster_betterposter_billboard_lint: {path}\n"
                f"  no large \\fontsize\\bfseries block detected — is this a betterposter?")
    m = max(candidates, key=lambda mm: int(mm.group("size")))

    size = int(m.group("size"))
    body_raw = m.group("body")
    # Strip simple LaTeX commands.
    body = re.sub(r"\\[A-Za-z@]+\*?\{?([^}]*)\}?", lambda mm: mm.group(1), body_raw)
    body = re.sub(r"[{}]", "", body).strip()

    en_words = len([w for w in re.findall(r"[A-Za-z]+", body) if len(w) > 1])
    jp_chars = sum(1 for ch in body if "぀" <= ch <= "ヿ" or "一" <= ch <= "鿿")
    jargon_hits = [j for j in _JARGON_HINTS if j.lower() in body.lower()]
    has_bullets = bool(re.search(r"[•・\-]", body))
    has_ref = bool(re.search(r"\[\d+\]", body))

    issues = []
    if size < 120:
        issues.append(f"[HIGH] billboard font {size}pt < 120pt; not large enough to dominate the poster")
    if en_words > 15:
        issues.append(f"[HIGH] {en_words} English words > 15; trim to a single plain-language sentence")
    if jp_chars > 30:
        issues.append(f"[HIGH] {jp_chars} JP chars > 30; trim to a single plain-language sentence")
    if len(jargon_hits) > 2:
        issues.append(f"[MEDIUM] {len(jargon_hits)} jargon tokens ({jargon_hits}); use lay-language equivalents")
    if has_bullets:
        issues.append("[MEDIUM] bullet glyphs in billboard; rewrite as a sentence")
    if has_ref:
        issues.append("[LOW] reference marker in billboard; move [N] to side bar")

    lines = [f"poster_betterposter_billboard_lint: {path}",
             f"  billboard size: {size}pt; words: en={en_words}, jp={jp_chars}",
             f"  text: {body[:80]}{'…' if len(body) > 80 else ''}"]
    if issues:
        lines.append(f"FAIL — {len(issues)} issue(s):")
        lines.extend(f"  {i}" for i in issues)
    else:
        lines.append("PASS — billboard meets Morrison plain-language criteria")
    return "\n".join(lines)
