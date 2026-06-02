"""Tier 1 — poster_typography_lints.

Source: Colin Purrington's automated-checkable rules. Cheap regex
hygiene that catches the most common typographic errors on posters.
"""
from __future__ import annotations

import re

from .._textutil import read_tex


_RULES = []


def _rule(name, severity, hint):
    def deco(fn):
        _RULES.append((name, severity, hint, fn))
        return fn
    return deco


@_rule("double_space", "LOW", "use single space between sentences (Purrington)")
def _double_space(tex: str):
    # Only count whitespace inside text-like context (skip leading indent).
    return len(re.findall(r"[A-Za-z\.,;:!?]  +(?=[A-Z])", tex))


_KNOWN_ACRONYMS = {
    "IEEE", "IEEJ", "IEICE", "IEEE-TPEL", "COMSOL", "MATLAB", "TIKZ",
    "ANSYS", "CST", "JSAP", "ABSTRACT", "INTRODUCTION", "METHODS",
    "RESULTS", "DISCUSSION", "CONCLUSION", "REFERENCES",
    "FEM", "PDE", "BEM", "PML", "WCAG", "CSS", "HTML",
    "TOKYO", "OSAKA", "JAPAN", "TIPS",
}


@_rule("all_caps_run", "MEDIUM", "avoid runs of ALL CAPS in body / titles (Purrington)")
def _all_caps(tex: str):
    """Count ALL CAPS runs while ignoring widely-known acronyms."""
    hits = re.findall(r"\b[A-Z]{5,}\b(?!\))", tex)
    return sum(1 for h in hits if h not in _KNOWN_ACRONYMS)


@_rule("title_case_section", "LOW", "use sentence case, not Title Case (Purrington)")
def _title_case_section(tex: str):
    n = 0
    for m in re.finditer(r"\\section\*?\{([^}]+)\}", tex):
        title = m.group(1)
        words = [w for w in re.split(r"\s+", title) if re.match(r"^[A-Z][a-z]+$", w)]
        if len(words) >= 3:
            # 3+ capitalized lowercase words = Title Case
            n += 1
    return n


@_rule("underline_for_emphasis", "MEDIUM", "use \\textit / \\textbf instead of \\underline (Purrington)")
def _underline(tex: str):
    return len(re.findall(r"\\underline\b", tex))


@_rule("bad_quotation_marks", "LOW", "use `` '' for curly quotes, not \" or '' (Purrington)")
def _quotes(tex: str):
    # Straight ASCII quotes inside text mode
    return len(re.findall(r'(?<![\\=\w])"(?!`)', tex))


@_rule("data_is", "LOW", "'data are' (plural), not 'data is' (Purrington)")
def _data_is(tex: str):
    return len(re.findall(r"\bdata is\b", tex, re.IGNORECASE))


@_rule("trailing_colon_title", "LOW", "avoid colons in \\title (Purrington)")
def _colon_title(tex: str):
    m = re.search(r"\\title\{([^}]+)\}", tex)
    return 1 if m and ":" in m.group(1) else 0


@_rule("bullets_in_heading", "LOW", "no bullets on section headings (Purrington)")
def _bullets_in_heading(tex: str):
    # heuristic: \section{• something} or \section{- foo}
    return len(re.findall(r"\\(sub)?section\*?\{\s*[•・\-\*]", tex))


@_rule("naked_double_prime", "LOW", "use `` '' for quotes, not double-prime glyph (Purrington)")
def _double_prime(tex: str):
    return len(re.findall(r"[″〞]", tex))


def poster_typography_lints(tex_path: str) -> str:
    """Run cheap regex-based typography hygiene checks.

    Each rule is sourced from Purrington's "Designing conference posters"
    list of mistakes; severity defaults to LOW unless the rule materially
    hurts readability (e.g. underline-for-emphasis, ALL CAPS run).

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Tally per rule + total issue count.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    lines = [f"poster_typography_lints: {path}"]
    total = 0
    for name, severity, hint, fn in _RULES:
        try:
            n = fn(tex)
        except re.error:
            n = 0
        if n:
            total += 1
            lines.append(f"  [{severity}] {name}: {n} hit(s) — {hint}")
    if total == 0:
        lines.append("PASS — no typography hygiene issues")
    else:
        lines.append(f"FAIL — {total} rule(s) triggered")
    return "\n".join(lines)
