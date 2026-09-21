"""Shared policy for editing documents under a hard page limit.

The policy is intentionally exposed through both the paper-writing and
grant-writing MCP servers. Page pressure must change *which content is kept*,
not make the prose that survives harder to understand.
"""

from __future__ import annotations


_DOCUMENT_LABELS = {
    "paper": "論文・研究会原稿",
    "grant": "助成金申請書",
}


def page_limit_content_selection_policy(document_type: str) -> dict:
    """Return the mandatory cross-document page-limit editing protocol."""
    if document_type not in _DOCUMENT_LABELS:
        allowed = ", ".join(sorted(_DOCUMENT_LABELS))
        raise ValueError(
            f"document_type must be one of: {allowed}; got {document_type!r}"
        )

    return {
        "policy_id": "page_limit_select_content_do_not_compress_prose",
        "document_type": document_type,
        "document_label": _DOCUMENT_LABELS[document_type],
        "mandatory": True,
        "principle": (
            "ページ制約は文章圧縮で解決しない。残す文章は自然で自己完結したままにし、"
            "収まらない分は優先度の低い内容を論点単位で削除または適切な欄へ移す。"
            "意味を明確に保てない一文は、圧縮版を残さず、その一文・論点ごと削除する。"
        ),
        "forbidden_actions": [
            {
                "id": "compress_surviving_prose",
                "rule": (
                    "行数や文字数を保つために、残す文の主語、対象、条件、因果、"
                    "指示対象を省略しない。"
                ),
            },
            {
                "id": "pack_relations_into_labels",
                "rule": (
                    "複数の課題・人物・手法の関係を、名詞列、括弧、スラッシュ、"
                    "『その』だけで一文へ押し込まない。"
                ),
            },
            {
                "id": "typography_for_content_packing",
                "rule": (
                    "追加内容を残す目的で本文フォント、行間、余白、見出し間隔を"
                    "縮めない。"
                ),
            },
            {
                "id": "reward_shortness_without_semantic_review",
                "rule": (
                    "短くなったこと、ページに収まったこと、機械点が上がったことだけで"
                    "改稿を改善と判定しない。"
                ),
            },
        ],
        "required_sequence": [
            {
                "step": 1,
                "action": "choose_semantically_sound_baseline",
                "rule": (
                    "機械点が最も高い版ではなく、人が読んで意味関係を理解できる版を"
                    "基準稿にする。圧縮で意味が壊れた高得点版を出発点にしない。"
                ),
            },
            {
                "step": 2,
                "action": "improve_content_and_write_clear_version",
                "rule": (
                    "基準稿の内容と日本語を、行数を気にせず改善する。誰が、何を、"
                    "どの対象で、何のために行うかが一読で分かる自然な文として確定する。"
                ),
            },
            {
                "step": 3,
                "action": "build_then_confirm_actual_overflow",
                "rule": (
                    "内容改善後の稿をPDFへ組版し、実際の超過箇所と必要削減量を確認する。"
                    "超過がなければ削除しない。"
                ),
            },
            {
                "step": 4,
                "action": "inventory_retained_meaning",
                "rule": (
                    "残す主張ごとに、主体、対象、操作、条件、課題間の関係、得られる知見を"
                    "確認し、圧縮による意味落ちを許さない。"
                ),
            },
            {
                "step": 5,
                "action": "rank_sentence_and_claim_importance",
                "rule": (
                    "各一文・論点が、中心的な問い、審査項目、主要な方法・検証、"
                    "固有の証拠のどれを担うかで重要度を評価する。短さや削減文字数を"
                    "重要度の代用にしない。"
                ),
            },
            {
                "step": 6,
                "action": "drop_whole_low_priority_content_unit",
                "rule": (
                    "最低重要度の一文を丸ごと削除する。必要なら、重複説明、副次的な例、"
                    "補助的な証拠・結果・背景を一つの論点として丸ごと削除または"
                    "適切な欄へ移す。"
                ),
            },
            {
                "step": 7,
                "action": "repair_transition_without_compression",
                "rule": (
                    "削除後の接続だけを自然に直す。残した複数の意味を短い一文へ"
                    "再充填しない。"
                ),
            },
            {
                "step": 8,
                "action": "rebuild_and_human_read",
                "rule": (
                    "再組版してページ数を確認し、修正箇所を前後の文脈込みで音読・通読する。"
                    "意味が弱くなった場合は採用しない。"
                ),
            },
        ],
        "protected_content": [
            "中心的な問い・目的",
            "主体と対象の対応",
            "比較・検証条件",
            "課題間の区別と関係",
            "得られる新しい知見",
            "結論を支える主要な証拠",
        ],
        "importance_rubric": [
            {
                "level": 4,
                "meaning": "必須",
                "criterion": "中心的な問い・目的または公式審査項目へ直接答える。",
            },
            {
                "level": 3,
                "meaning": "中核",
                "criterion": "主要な方法、比較条件、検証、成果・知見を成立させる。",
            },
            {
                "level": 2,
                "meaning": "有力な裏付け",
                "criterion": "実現可能性や主張を支える固有で代替しにくい証拠である。",
            },
            {
                "level": 1,
                "meaning": "補助",
                "criterion": "副次例、追加背景、補助的な結果で、なくても主筋が通る。",
            },
            {
                "level": 0,
                "meaning": "削除候補",
                "criterion": "重複、脱線、役割不明、または他の一文で既に回収されている。",
            },
        ],
        "deletion_rule": (
            "最低重要度の完全な一文を最初の削除候補とする。同点なら、中心的な問いとの"
            "距離が遠く、他文と証拠が重複する一文を先に削る。削除後に指示語、接続、"
            "番号参照、論理の飛躍が生じないことを確認する。"
        ),
        "baseline_rule": (
            "改稿前後の機械点が、意味の通る基準稿の選択を上書きしてはならない。"
            "まず意味の通る旧稿を内容面から改善し、組版後に必要な削除だけを行う。"
        ),
        "preferred_drop_order": [
            "同じ主張を繰り返す説明",
            "中心結論に不要な固有名詞や副次例",
            "補助的な背景説明",
            "副次的な結果または将来課題",
            "他欄・補足資料へ移せる証拠",
        ],
        "stop_condition": (
            "必要な意味を保ったまま丸ごと落とせる論点がない場合は、"
            "ページ制約と内容範囲が衝突していると報告して著者判断へ戻す。"
            "文を圧縮して自動解決しない。"
        ),
        "acceptance_rule": (
            "採用条件は、ページ内に収まり、かつ残した各文の意味関係が改稿前以上に"
            "明確であること。ページ適合や機械点の上昇だけでは採用しない。"
            "圧縮して意味不明な一文を残すより、その一文・論点がない版を選ぶ。"
        ),
    }
