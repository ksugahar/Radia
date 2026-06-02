"""Tier 7 — poster_elevator_pitch_generate.

Generate a 3-minute spoken script from the poster .tex. Allocation:
30s hook → 90s main result → 60s implication/limitation/next.

Heuristic — no LLM call. Pulls the title, abstract paragraph (first
contentbox), and the まとめ / Conclusion contentbox, and renders them
into a speakable script with timing cues.
"""
from __future__ import annotations

import re

from .._textutil import read_tex


_TITLE_RE = re.compile(
    r"\\fontsize\{60\}\{72\}\\selectfont\\bfseries\\color\{[^}]+\}\s+([^}]+?)\s*\}"
)
_ABSTRACT_RE = re.compile(
    r"はじめに[\s\S]+?\\end\{tcolorbox\}\s*\{\\fontsize\{24pt\}\{33pt\}\\selectfont\\par\s*([\s\S]+?)\\par\s*\}",
)
_SUMMARY_RE = re.compile(
    r"(?:まとめ|conclusion|summary)[\s\S]+?\\begin\{tcolorbox\}\[contentbox[\s\S]+?\{\\fontsize\{24pt\}\{33pt\}\\selectfont\\par\s*([\s\S]+?)\\par",
    re.IGNORECASE,
)


def _clean(s: str) -> str:
    s = re.sub(r"\\[A-Za-z@]+\*?\{([^}]*)\}", lambda m: m.group(1), s)
    s = re.sub(r"\\[A-Za-z@]+\*?", "", s)
    s = re.sub(r"[{}]", "", s)
    return re.sub(r"\s+", " ", s).strip()


def poster_elevator_pitch_generate(tex_path: str) -> str:
    """Render a 3-minute speakable script from the poster source.

    Sections:
    - 0:00-0:30 Hook (title + 1-sentence motivation)
    - 0:30-2:00 Main result (most of the abstract body)
    - 2:00-3:00 Implication / next step (summary box if present)

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        A scripted pitch with timing cues, ready to rehearse.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    tm = _TITLE_RE.search(tex)
    title = _clean(tm.group(1)) if tm else "(タイトル未検出)"

    am = _ABSTRACT_RE.search(tex)
    abstract = _clean(am.group(1)) if am else ""
    sm = _SUMMARY_RE.search(tex)
    summary = _clean(sm.group(1)) if sm else ""

    # Split abstract into sentences (Japanese 。 / English . without
    # requiring whitespace separator).
    chunks = [c.strip() for c in re.split(r"(?<=[。.！？!?])", abstract) if c.strip()]
    hook = chunks[0] if chunks else "(導入文がありません)"
    body = "".join(chunks[1:]) if len(chunks) > 1 else "(本文を contentbox に書き足してください)"

    lines = [
        f"poster_elevator_pitch_generate: {path}",
        "",
        "[0:00-0:30] Hook ----------------------------------",
        f"  ご覧いただきありがとうございます。私は{title}について発表しています。",
        f"  {hook}",
        "",
        "[0:30-2:00] Main result ---------------------------",
        f"  {body if body else '(本文を contentbox に書き足してください)'}",
        "",
        "[2:00-3:00] Implication / next --------------------",
        f"  {summary if summary else '(まとめ contentbox がありません — 結論と次の一歩を口頭で補足)'}",
        "  ご質問ありますでしょうか?",
    ]
    return "\n".join(lines)
