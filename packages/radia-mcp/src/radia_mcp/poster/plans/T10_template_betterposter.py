"""Tier 3 — poster_template_betterposter.

A0 landscape, Morrison 2019 three-zone format (left TL;DR / center
billboard / right details) adapted for Japanese with bxjsarticle +
tcolorbox.

Source: Mike Morrison, "How to design a better research poster" (2019);
ScienceUX eye-tracking pilot (2020).
"""
from __future__ import annotations

import pathlib
import re

from .._textutil import read_tex  # noqa: F401 — kept for symmetry


_HERE = pathlib.Path(__file__).resolve().parent.parent
_TEMPLATE = _HERE / "templates" / "betterposter_a0_landscape.tex"


_FIELD_PATTERNS: dict[str, re.Pattern[str]] = {
    "billboard": re.compile(
        r"(\{\\fontsize\{120\}\{140\}\\selectfont\\bfseries\s*)([\s\S]+?)(\s*\})"
    ),
    "billboard_subtitle": re.compile(
        r"(\{\\fontsize\{50\}\{60\}\\selectfont\s+)([\s\S]+?)(\s*\}\s*\\end\{center\})"
    ),
    "author": re.compile(
        r"(\{\\fontsize\{22\}\{28\}\\selectfont\s+)([^}]+?)(\s*\})"
    ),
    "title": re.compile(
        r"(\\bfseries\s+タイトル / Title\}\\\\\[4pt\]\s*\{\\fontsize\{22\}\{28\}\\selectfont\s+)([^}]+?)(\s*\})"
    ),
}


def poster_template_betterposter(out_path: str | None = None,
                                 replace: dict[str, str] | None = None) -> str:
    """Return (or write) the A0-landscape #betterposter template.

    Single-finding billboard in the center, TL;DR sidebar on the left,
    detail sidebar on the right. The billboard text MUST be in plain
    language — see ``poster_betterposter_billboard_lint``.

    ``replace`` keys:
        - ``billboard``: the main finding (≤15 words, plain language)
        - ``billboard_subtitle``: 1-sentence elaboration
        - ``author``: byline
        - ``title``: small-print title in the right bar

    Args:
        out_path: Optional output ``.tex`` path. If None, the source is
                  returned as a string.
        replace: Optional field substitutions.

    Returns:
        Status line when ``out_path`` is given; rendered .tex otherwise.
    """
    tex_src = _TEMPLATE.read_text(encoding="utf-8")

    warnings: list[str] = []
    n_subs = 0

    if replace:
        unknown = set(replace) - set(_FIELD_PATTERNS)
        if unknown:
            warnings.append(f"unknown replace keys ignored: {sorted(unknown)}")
        for field, pat in _FIELD_PATTERNS.items():
            if field not in replace:
                continue
            new_val = replace[field]
            new_src, n = pat.subn(
                lambda m, v=new_val: f"{m.group(1)}{v}{m.group(3)}",
                tex_src, count=1,
            )
            if n == 0:
                warnings.append(f"field {field!r}: pattern did not match")
            else:
                tex_src = new_src
                n_subs += 1

    if out_path is None:
        head = ""
        if warnings:
            head = "% poster_template_betterposter warnings:\n" + "".join(
                f"%   - {w}\n" for w in warnings
            ) + "\n"
        return head + tex_src

    out = pathlib.Path(out_path)
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tex_src, encoding="utf-8")
    lines = [f"Wrote betterposter template to: {out}", f"  substitutions: {n_subs}"]
    if warnings:
        lines.append("  warnings:")
        lines.extend(f"    - {w}" for w in warnings)
    return "\n".join(lines)
