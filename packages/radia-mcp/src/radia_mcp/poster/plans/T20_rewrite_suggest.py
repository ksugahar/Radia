"""Tier 6 — poster_rewrite_suggest.

Template-based rewrite candidates (no LLM call needed). Mirrors the
``grant_writing_rewrite_suggest`` / ``presentation_rewrite_suggest``
pattern: each target offers 3-4 alternative phrasings keyed to typical
weaknesses.
"""
from __future__ import annotations


_TEMPLATES: dict[str, list[str]] = {
    "title": [
        "{topic}: A {method} Approach for {problem}",
        "{problem}を{method}で解く: {topic}",
        "{method}による{topic}の{benefit}",
        "{topic}の{benefit}を実現する{method}",
    ],
    "subtitle": [
        "Maxwell 方程式の{transformation}を用いた{domain}の定式化",
        "{transformation}による{domain}の統一的取扱い",
        "{domain}を{transformation}で{benefit}する",
    ],
    "billboard": [
        "{method}で{problem}が{benefit}できる。",
        "{problem}は{cause}で起こる。{method}で解決。",
        "{benefit}を{x_fold}倍。{method}による。",
    ],
    "abstract": [
        "本研究では、{topic}に対して{method}を提案する。{problem}を{benefit}し、"
        "{validation}で妥当性を示した。{audience}にとって{implication}が得られる。",
        "{problem}を解くため、{method}を開発した。{key_result}を {validation} で確認。"
        "従来の{baseline}に比べ {benefit}。",
    ],
    "takeaway": [
        "1. {method}により{benefit}\n2. {limitation}に注意\n3. {next_step}が今後の課題",
        "・{benefit}を{x_fold}倍\n・{generalization}が可能\n・{availability}",
    ],
    "qr_label": [
        "Scan for full paper",
        "Scan for preprint + data",
        "Scan for code + dataset",
        "論文 + データ ↓",
    ],
}


def poster_rewrite_suggest(target: str, slot_values: dict[str, str] | None = None) -> str:
    """Return 3-4 candidate phrasings for a target poster element.

    ``target`` is one of: title, subtitle, billboard, abstract, takeaway,
    qr_label. ``slot_values`` fills in template placeholders like
    ``{method}``, ``{problem}``, ``{benefit}``.

    Args:
        target: One of the supported target keys.
        slot_values: Optional dict to fill ``{slot}`` placeholders.

    Returns:
        Numbered candidates with placeholders kept verbatim when no slot
        value is supplied, so the user can see where to insert content.
    """
    if target not in _TEMPLATES:
        return f"Error: unknown target {target!r}. Supported: {sorted(_TEMPLATES)}"
    cands = _TEMPLATES[target]
    sv = slot_values or {}
    lines = [f"poster_rewrite_suggest target={target}"]
    for i, tmpl in enumerate(cands, 1):
        try:
            rendered = tmpl.format_map(_Defaulting(sv))
        except KeyError as e:
            rendered = f"(missing slot: {e})"
        lines.append(f"  {i}. {rendered}")
    return "\n".join(lines)


class _Defaulting(dict):
    """Dict subclass that returns ``{slot_name}`` for missing keys."""

    def __missing__(self, key):
        return "{" + key + "}"
