"""Plan B T6: Figure / Table referencing coverage 診断.

Tufte / Wallwork §17 / reviewer-2 archetype「図は話題に上ったか」に基づく粗い温度計。
\\label{fig:X} / \\label{tab:X} が本文の \\ref で何回参照されているかを集計し、

    - unreferenced (count == 0) : page space の無駄 → 削除か \\ref 追加
    - singleton  (count == 1)   : redundant の疑い → 2 回以上議論するか削る
    - well_referenced (count>=2): 複数箇所で議論 → 図が仕事をしている

を測定する。sec:/eq: など別 prefix は対象外 (figure/table に絞って評価する)。

設計方針:
- score は 0-10 の粗い温度計。9 以上は誤差範囲、6 以下は改善サイン。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動化 numeric 判定には使わない。
- paper_writing_check_figure_forward_reference (orphan / dangling) と補完関係。
  (あちらは 0/1 存在チェック、本ツールは カウントで redundant 疑いを flag)
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Patterns
# ------------------------------------------------------------------
# figure / table 系の label のみ対象 (sec:/eq: は除外)
_LABEL_PATTERN = re.compile(
    r"\\label\{((?:fig|tab|table):[^}]+)\}"
)

# \ref / \autoref / \Cref / \cref / \vref / \pageref / \eqref (eqref は対象外)
_REF_PATTERN = re.compile(
    r"\\(?:ref|autoref|Cref|cref|vref|pageref)\*?"
    r"\{((?:fig|tab|table):[^}]+)\}"
)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_figure_referencing_coverage(tex: str) -> dict:
    """Every \\label{fig:X}/\\label{tab:X} が本文で \\ref される回数を集計。

    unreferenced = 無駄な page space、referenced_once = redundant 疑いを flag。

    Args:
        tex: 論文原稿 (LaTeX)。

    Returns:
        dict: score / score_max / comments / observations / hint / source
    """
    if tex is None:
        tex = ""

    # ---- Step 1: label 収集 (順序を保持、重複定義は最初の位置を記録) ----
    labels_defined: dict[str, int] = {}
    for m in _LABEL_PATTERN.finditer(tex):
        name = m.group(1)
        if name not in labels_defined:
            labels_defined[name] = m.start()

    # ---- Step 2: ref カウント ----
    ref_counts: dict[str, int] = {name: 0 for name in labels_defined}
    for m in _REF_PATTERN.finditer(tex):
        name = m.group(1)
        if name in ref_counts:
            ref_counts[name] += 1
        # dangling refs (label 未定義への ref) は T6 では扱わない
        # (forward_reference ツール側で扱う)

    # ---- Step 3: カテゴリ分類 ----
    unreferenced_labels = sorted(
        [name for name, c in ref_counts.items() if c == 0]
    )
    singleton_labels = sorted(
        [name for name, c in ref_counts.items() if c == 1]
    )
    well_referenced = sorted(
        [name for name, c in ref_counts.items() if c >= 2]
    )

    total_labels = len(labels_defined)
    if total_labels > 0:
        coverage_ratio = len(well_referenced) / total_labels
    else:
        coverage_ratio = 1.0

    observations = {
        "labels_defined": dict(labels_defined),
        "ref_counts": dict(ref_counts),
        "unreferenced_labels": unreferenced_labels,
        "singleton_labels": singleton_labels,
        "well_referenced": well_referenced,
        "total_labels": total_labels,
        "coverage_ratio": round(coverage_ratio, 3),
    }

    # ---- Step 4: score ----
    if total_labels == 0:
        score = 10.0
    else:
        raw = 10 - (len(unreferenced_labels) * 2) - len(singleton_labels)
        score = float(max(0, min(10, raw)))

    # ---- Step 5: comments (PRIMARY) ----
    comments: list[str] = []

    if total_labels == 0:
        comments.append(
            "figure/table labels が検出されない (文章のみの論文?)。該当なし。"
        )
        comments.append(
            "図表を入れる予定なら \\label{fig:...} / \\label{tab:...} を"
            "付けてから再診断を。"
        )
    else:
        if unreferenced_labels:
            n = len(unreferenced_labels)
            shown = ", ".join(unreferenced_labels[:5])
            comments.append(
                f"参照されない figure/table が {n} 件: {shown}。"
                "本文で議論されない図は page space の無駄 — "
                "\\ref を追加 or 削除。"
            )
        if singleton_labels:
            n = len(singleton_labels)
            shown = ", ".join(singleton_labels[:5])
            comments.append(
                f"1 回しか参照されない figure/table が {n} 件: {shown}。"
                "redundant の疑い。2 回以上議論するか、思い切って削る "
                "(Tufte / Wallwork §17)。"
            )
        if (
            not unreferenced_labels
            and not singleton_labels
            and well_referenced
        ):
            comments.append(
                f"全 {len(well_referenced)} figure/table が複数箇所で"
                "参照されている。good coverage。"
            )
        # 追加の補助 comment
        if unreferenced_labels or singleton_labels:
            comments.append(
                "図表は「本文で 2 回以上議論される」ことが推奨 "
                "(Tufte)。導入・結果考察・まとめ等で複数回 \\ref を挟む。"
            )

    # 最少 2 comment 確保
    if len(comments) < 2:
        comments.append(
            "score は粗い指標。skill.md §図表 / Wallwork §17 を直接参照して"
            "人間判断で補完すること。"
        )
    comments = comments[:5]

    hint = (
        "\\label 検出は fig:/tab: prefix で絞っている。"
        "sec:/eq: など別 prefix は対象外。"
    )
    source = (
        "Tufte / Wallwork §17 / reviewer-2 archetype "
        "'図は話題に上ったか'"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": source,
    }
