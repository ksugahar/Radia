"""Tier 6 — poster_root_cause_diagnosis.

Map the per-axis severity profile to one of 5 typical root-cause
patterns, so the user gets "this poster is title-heavy" rather than a
flat list of 30 lint hits.
"""
from __future__ import annotations

import re

from .T17_health_report import poster_health_report


_PATTERNS = [
    {
        "name": "title-heavy",
        "trigger": lambda axis: axis["word_count"]["HIGH"] >= 1 and axis["fontsize_by_distance"].get("HIGH", 0) == 0,
        "summary": "本文が長すぎて図やまとめが圧倒されている (Purrington 1000-word 超え)",
        "remedy": "本文を Hook / Method / Result / Conclusion ≤200語ずつに分解し、図の caption に逃がす",
    },
    {
        "name": "figure-starved",
        "trigger": lambda axis: axis["caption_self_contained"].get("MEDIUM", 0) >= 1,
        "summary": "図の caption が自己完結していない (軸+単位+凡例 不足)",
        "remedy": "各図の caption に「軸ラベル+単位+凡例+(統計情報)」を必ず入れる",
    },
    {
        "name": "color-chaotic",
        "trigger": lambda axis: (axis["color_count"].get("LOW", 0) + axis["color_count"].get("MEDIUM", 0) >= 1
                                  or axis["colorblind"].get("MEDIUM", 0) >= 1
                                  or axis["color_contrast"].get("HIGH", 0) >= 1),
        "summary": "配色がまとまっておらず、WCAG / CVD で読みにくい",
        "remedy": "70/25/5 ルールに従って 3 色に絞り、赤×緑ペアを避け、text/bg の WCAG AA 3:1 を確保",
    },
    {
        "name": "billboard-buried",
        "trigger": lambda axis: axis["fontsize_by_distance"].get("HIGH", 0) >= 1,
        "summary": "視認距離に対してフォントが小さく、1m / 3m / 5m での読み取りが破綻",
        "remedy": "Editage 公式 (h_mm = d_m / 0.25) で必要 pt を逆算し、タイトル≥50pt / 本文≥24pt を厳守",
    },
    {
        "name": "typography-noisy",
        "trigger": lambda axis: sum(axis["typography"].values()) >= 3,
        "summary": "全角/半角・引用符・大文字運用などのタイポグラフィに古い癖",
        "remedy": "下線禁止、ALL CAPS 回避、sentence case、ダブルスペース禁止 (Purrington)",
    },
]


def poster_root_cause_diagnosis(tex_path: str) -> str:
    """Diagnose which of the 5 typical poster failure patterns apply.

    Builds on ``poster_health_report``'s severity counts; matches each
    pattern by a hand-crafted predicate and prints a one-line summary +
    remedy for every triggered pattern.
    """
    report = poster_health_report(tex_path)
    # Parse the per-axis severity counts back out of the health report.
    axis: dict[str, dict[str, int]] = {}
    for m in re.finditer(r"  - (\w+): H=(\d+) M=(\d+) L=(\d+)", report):
        axis[m.group(1)] = {"HIGH": int(m.group(2)), "MEDIUM": int(m.group(3)), "LOW": int(m.group(4))}
    if not axis:
        return f"poster_root_cause_diagnosis: {tex_path}\n  health report could not be parsed (file unreadable?)"

    lines = [f"poster_root_cause_diagnosis: {tex_path}"]
    triggered = []
    for pat in _PATTERNS:
        try:
            if pat["trigger"](axis):
                triggered.append(pat)
        except KeyError:
            continue

    if not triggered:
        lines.append("  no typical root-cause pattern matched — poster looks structurally sound")
        return "\n".join(lines)

    for pat in triggered:
        lines.append(f"  [{pat['name']}]")
        lines.append(f"    症状: {pat['summary']}")
        lines.append(f"    対策: {pat['remedy']}")
    return "\n".join(lines)
