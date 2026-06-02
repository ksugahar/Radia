"""Tier 6 — poster_adaptive_health_report.

Adjust severity per conference preset. Same pattern as
``grant_writing_adaptive_health_report``.
"""
from __future__ import annotations

from .T17_health_report import poster_health_report


_PRESETS = {
    "COMSOL": {
        "name": "COMSOL Conference",
        "boost": ["figure_dominance", "color_contrast"],
        "demote": ["typography", "color_count"],
        "note": "図中心 (COMSOL の数値例が見せ場)。タイポ系は重視されない",
    },
    "IEEJ": {
        "name": "電気学会全国大会 / 部門大会",
        "boost": ["jp_font", "fontsize_by_distance", "caption_self_contained"],
        "demote": [],
        "note": "日本語ポスター。Gothic 必須、視認距離公式厳守。最大 3 分でコア質問が来る",
    },
    "IEEE_Magnetics": {
        "name": "IEEE Magnetics Society (Intermag / Magnetism & Magnetic Materials)",
        "boost": ["caption_self_contained", "color_contrast", "color_count"],
        "demote": ["jp_font"],
        "note": "英語、図優位。CVD 配慮を強く期待される (Magnetics Letter 投稿規程と同等)",
    },
    "Compumag": {
        "name": "Compumag / IGTE",
        "boost": ["caption_self_contained", "fontsize_by_distance"],
        "demote": ["color_count"],
        "note": "計算電磁気の国際会議。grayscale ok、equation slide compliance 重視",
    },
    "応物": {
        "name": "応用物理学会 (JSAP)",
        "boost": ["jp_font", "fontsize_by_distance"],
        "demote": [],
        "note": "日英混在。技術文書らしい簡潔さが好まれる",
    },
}


def poster_adaptive_health_report(tex_path: str,
                                  conference: str = "IEEJ") -> str:
    """Health report tuned for a target conference's review culture.

    Args:
        tex_path: Path to the poster ``.tex`` source.
        conference: One of the keys in ``_PRESETS``. Falls back to
                    ``IEEJ`` (electrical-engineering Japan domestic).

    Returns:
        The base health report + conference-specific commentary.
    """
    if conference not in _PRESETS:
        return f"Error: unknown conference {conference!r}. Supported: {sorted(_PRESETS)}"

    base = poster_health_report(tex_path)
    cfg = _PRESETS[conference]

    lines = [f"poster_adaptive_health_report: {tex_path} (preset={conference})",
             f"  venue: {cfg['name']}",
             f"  note: {cfg['note']}",
             "",
             "--- base health report ---",
             base,
             "",
             "--- preset adjustments ---"]
    if cfg["boost"]:
        lines.append(f"  boost (重要視される軸): {', '.join(cfg['boost'])}")
    if cfg["demote"]:
        lines.append(f"  demote (相対的に許容される軸): {', '.join(cfg['demote'])}")
    return "\n".join(lines)
