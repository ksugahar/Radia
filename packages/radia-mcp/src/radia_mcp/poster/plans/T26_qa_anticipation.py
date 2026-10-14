"""Tier 7 — poster_qa_anticipation_list.

Heuristic Q&A anticipation based on poster content patterns. Each
question is annotated with the typical motivation behind it
(reviewer-type taxonomy) so the user can rehearse adequate answers.
"""
from __future__ import annotations

import re

from .._textutil import read_tex


_TAXONOMY = [
    {
        "trigger": lambda tex: "PML" in tex or "ABC" in tex.upper(),
        "q": "PML / 吸収境界条件と比較してどう違うのか？",
        "why": "比較ベンチマーク要求 (reviewer 1 型)",
    },
    {
        "trigger": lambda tex: re.search(r"渦電流|eddy current|Eddy", tex, re.IGNORECASE),
        "q": "周波数依存性はどう扱っているのか？高周波と低周波で同じ手法か？",
        "why": "適用範囲の質問 (reviewer 1 型)",
    },
    {
        "trigger": lambda tex: re.search(r"COMSOL|FEM|FEA", tex, re.IGNORECASE),
        "q": "メッシュ密度や次数依存性を確認したか？",
        "why": "数値検証要求 (reviewer 2 型)",
    },
    {
        "trigger": lambda tex: re.search(r"hertz|antenna|アンテナ", tex, re.IGNORECASE),
        "q": "解析解との誤差を定量化しているか？dB 単位で示してほしい",
        "why": "定量評価要求",
    },
    {
        "trigger": lambda tex: re.search(r"3D|three-dimensional|立体|3次元", tex, re.IGNORECASE),
        "q": "3D 化したときの計算時間はどれくらい増えるのか？",
        "why": "実用性 / スケーラビリティ",
    },
    {
        "trigger": lambda tex: "ケルビン変換" in tex or "Kelvin" in tex,
        "q": "外部領域の物体は実際にどのような構造を想定しているのか？",
        "why": "応用イメージ要求 (新規分野からの聴衆)",
    },
    {
        "trigger": lambda tex: "PDE" in tex or "係数型" in tex,
        "q": "COMSOL の係数型 PDE での実装で使った具体的なオプションを教えてほしい",
        "why": "再現性要求 (実装者からの聴衆)",
    },
    {
        "trigger": lambda tex: re.search(r"ベンチマーク|benchmark", tex, re.IGNORECASE) is None,
        "q": "他のオープンソースのベンチマーク (TEAM 等) で検証する予定はあるか？",
        "why": "汎用性確認 (常套質問)",
    },
    {
        "trigger": lambda tex: True,  # 常に追加
        "q": "本研究のオリジナリティを 30 秒で説明してほしい",
        "why": "コア質問 (必須準備)",
    },
    {
        "trigger": lambda tex: True,
        "q": "今後の展開は何か？最も近い応用ターゲットは？",
        "why": "次の一歩 (必須準備)",
    },
]


def poster_qa_anticipation_list(tex_path: str) -> str:
    """Anticipate likely poster Q&A and tag each with reviewer-type motivation.

    Triggers are content-based (e.g. presence of PML, COMSOL, 3D
    keywords). Output is a numbered list — directly usable as a
    rehearsal prompt.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    qs = []
    for entry in _TAXONOMY:
        try:
            if entry["trigger"](tex):
                qs.append((entry["q"], entry["why"]))
        except Exception:
            continue

    lines = [f"poster_qa_anticipation_list: {path}",
             f"  {len(qs)} anticipated question(s):"]
    for i, (q, why) in enumerate(qs, 1):
        lines.append(f"  Q{i}. {q}")
        lines.append(f"      (why: {why})")
    return "\n".join(lines)
