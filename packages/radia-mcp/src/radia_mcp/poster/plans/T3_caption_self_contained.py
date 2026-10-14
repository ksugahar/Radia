"""Tier 1 — poster_caption_self_contained.

Source: Editage『学会ポスター発表完全ガイド6』— captions should let a
viewer understand the figure without reading body text. Required tokens:
axis label / unit / legend / statistical info (when relevant).
"""
from __future__ import annotations

import re

from .._textutil import read_tex


def _extract_captions(tex: str) -> list[str]:
    """Find ``\\caption`` / ``\\captionof{figure}`` blocks with balanced braces."""
    out: list[str] = []
    for marker in ("\\captionof{figure}{", "\\caption{"):
        i = 0
        while True:
            j = tex.find(marker, i)
            if j < 0:
                break
            start = j + len(marker)
            depth = 1
            k = start
            while k < len(tex) and depth > 0:
                if tex[k] == "{":
                    depth += 1
                elif tex[k] == "}":
                    depth -= 1
                k += 1
            if depth == 0:
                out.append(tex[start:k - 1])
            i = k
    return out

# Keyword inventories (JP + EN).
_AXIS_KW = (
    "axis", "x-axis", "y-axis", "horizontal", "vertical",
    "縦軸", "横軸", "X軸", "Y軸", "z軸", "Z軸",
)
_UNIT_KW = (
    "[Hz]", "[s]", "[m]", "[mm]", "[A]", "[V]", "[T]", "[W]", "[Pa]",
    "[°]", "[%]", "[dB]", "/m", "/s",
    "ヘルツ", "メートル", "ボルト", "アンペア", "テスラ",
)
_LEGEND_KW = (
    "legend", "key", "color", "dashed", "solid", "blue", "red", "black",
    "凡例", "実線", "破線", "青", "赤", "黒",
)
_STAT_KW = (
    "n=", "n =", "p<", "p =", "p<0", "mean", "std", "sd",
    "n＝", "p＝", "平均", "標準偏差", "誤差", "誤差棒",
)


def _hit(s: str, kw_list) -> bool:
    low = s.lower()
    return any(k.lower() in low for k in kw_list)


def poster_caption_self_contained(tex_path: str) -> str:
    """Score each ``\\caption`` / ``\\captionof{figure}`` for self-sufficiency.

    Each caption is checked for four token classes — axis label, unit,
    legend hint, statistical info. The score is the count of *classes
    present* (max 4). Captions scoring < 2 are flagged. Statistical
    keywords are only expected for "result" figures and don't count
    against captions in setup / concept diagrams.

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Per-caption report + summary.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    captions = _extract_captions(tex)
    if not captions:
        return f"poster_caption_self_contained: {path}\n  no captions detected — figures lack self-contained labels"

    lines = [f"poster_caption_self_contained: {path}", f"  captions analyzed: {len(captions)}"]
    flagged = 0
    for i, cap in enumerate(captions, 1):
        clean = re.sub(r"\\[A-Za-z@]+\*?\{[^}]*\}", "", cap)
        clean = re.sub(r"\\[A-Za-z@]+\*?", "", clean)
        clean = re.sub(r"[{}]", "", clean).strip()
        # Likely-result figure if it mentions distribution/comparison/etc.
        is_result = any(k in clean for k in (
            "distribution", "comparison", "result", "結果", "比較", "分布", "誤差",
        ))
        # Concept figures (short, conceptual title) get a relaxed bar.
        is_concept = (len(clean) <= 20 and not is_result)
        score = 0
        present: list[str] = []
        for kw_list, name in (
            (_AXIS_KW, "axis"),
            (_UNIT_KW, "unit"),
            (_LEGEND_KW, "legend"),
        ):
            if _hit(clean, kw_list):
                score += 1
                present.append(name)
        if is_result:
            if _hit(clean, _STAT_KW):
                score += 1
                present.append("stat")
        else:
            # Concept figures get a free pass on stat.
            score += 1

        if is_concept:
            verdict = "OK*"
            present.append("concept-figure")
        elif score >= 2:
            verdict = "OK "
        else:
            verdict = "WEAK"
            flagged += 1
        snip = clean[:50] + ("…" if len(clean) > 50 else "")
        lines.append(f"  [{verdict}] fig {i}: {score}/4 ({','.join(present) or '-'}) — {snip}")

    if flagged:
        lines.append(f"[MEDIUM] {flagged} caption(s) below 2/4 — viewers cannot grasp figure without body text")
    else:
        lines.append("PASS — all captions self-contained")
    return "\n".join(lines)
