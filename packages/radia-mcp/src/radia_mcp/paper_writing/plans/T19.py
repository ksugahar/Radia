"""Paper T19 (Intelligence Layer S5): Rewrite suggest (paper 用 target)
(paper_writing_rewrite_suggest).

paper 特有の rewrite target 11 種 × 3-5 candidate を template-based で提供。
"""

from __future__ import annotations


REWRITE_TEMPLATES: dict[str, dict] = {
    "abstract_structure": {
        "label": "Abstract 構造 (Bg→Methods→Results→Contribution)",
        "templates": [
            "【Background】近年 X 領域では Y が課題である。"
            "【Methods】本研究は Z を用いて W を評価する。"
            "【Results】M 条件下で N% の improvement を確認。"
            "【Contribution】本研究の貢献は (1)...(2)...(3)... である。",
            "【1 文問題提起】X で P の未解決問題がある。"
            "【2 文手法】我々は Q で R を行う。"
            "【3 文結果】S の結果、T を達成。"
            "【4 文 contribution】本研究の core contribution は U である。",
        ],
        "instruction": "Abstract は 150-250 words (英語) / 300-400 字 (日本語) の 4 要素構成。background は 25% 以内。",
    },
    "contribution_bullets": {
        "label": "Contribution 3 点列挙 (Intro 末尾)",
        "templates": [
            "The contributions of this paper are as follows:\n"
            "(1) We propose X, which achieves Y with Z% improvement over state-of-the-art.\n"
            "(2) We establish a theoretical foundation for W through M.\n"
            "(3) We release the dataset and code at https://github.com/...",
            "本論文の貢献は以下の 3 点である:\n"
            "(1) X 手法を提案し、Y 精度で Z を達成した。\n"
            "(2) W 理論基盤を M により確立した。\n"
            "(3) データセットとコードを公開した。",
        ],
        "instruction": "Introduction の最終段落に 3 点列挙で placing。各 contribution は 1-2 文 + 数値裏付け。",
    },
    "hedging_to_assertion": {
        "label": "Hedging (弱気) → Assertion (断定) (5 例)",
        "templates": [
            "「It seems that X」 → 「X」",
            "「may be」 → 「is」 (data で裏付けあれば)",
            "「could be」 → 「we found」",
            "「might indicate」 → 「indicates」",
            "「probably」 → 削除 or specific probability",
        ],
        "instruction": "Results は assertion、Discussion は moderated hedging (Wallwork §8.6)。純 results 文で hedging は使わない。",
    },
    "passive_to_active_paper": {
        "label": "Passive → Active (paper 版 5 例)",
        "templates": [
            "「X was measured」 → 「We measured X」",
            "「The results are shown in Fig. 2」 → 「Fig. 2 shows the results」",
            "「It is concluded that...」 → 「We conclude that...」",
            "「A method was developed」 → 「We developed a method」",
            "「The data were analyzed」 → 「We analyzed the data」",
        ],
        "instruction": "paper でも active voice が基本 (IMRAD 内各 section で異なる、Methods は受動 許容)。",
    },
    "statistical_reporting_fix": {
        "label": "p 値単独 → 4 要素完全 (5 例)",
        "templates": [
            "「p < 0.05」 → 「(t(18) = 2.41, p = 0.027, d = 0.54, 95% CI [0.06, 1.02], n = 20)」",
            "「significant」 → 「significantly higher (M = 3.4, SD = 0.8) vs control (M = 2.7, SD = 0.6), t(38) = 3.12, p = 0.003, Cohen's d = 0.87, n = 40」",
            "「no difference」 → 「no significant difference (p = 0.31, 95% CI [−0.15, 0.42])」",
            "「highly significant」 → 「(F(2, 57) = 14.3, p < 0.001, η² = 0.33)」",
            "OR 形式: 「X倍」 → 「(OR = 2.45, 95% CI [1.32, 4.56], p = 0.004, n = 250)」",
        ],
        "instruction": "APA 7th + CONSORT + STROBE: 全 p 値に effect size + 95% CI + n を添付。",
    },
    "discussion_4_elements": {
        "label": "Discussion 4 要素 skeleton",
        "templates": [
            "¶1 Interpretation: Our findings suggest that X through mechanism Y.\n"
            "¶2 Implications: These results have implications for Z field.\n"
            "¶3 Limitations: The study is limited by W (sample size / specific condition).\n"
            "¶4 Future work: Future studies should extend to V and examine U.\n"
            "¶5 Generalization: Our results may generalize to T settings, with caveats about S.",
        ],
        "instruction": "Discussion を Results 再記述で終わらせない。4 要素を明示的に articulate (Day & Gastel, Wallwork Ch.13)。",
    },
    "limitation_paragraph": {
        "label": "Limitation 段落 (3 テンプレ)",
        "templates": [
            "Our study has several limitations. First, the sample size (n = X) limits statistical power for Y subgroup analyses. Second, we used Z metric which may not generalize to W contexts. Third, the experimental setup assumes V which holds for U but not T.",
            "本研究には複数の限界がある。第一に、サンプルサイズ (n = X) が小さく Y 分析の統計的検出力が限定的である。第二に、Z 条件に限定しており W 設定での汎化は未検証である。第三に、V 仮定が成り立つ範囲に限られる。",
            "A limitation of this work is the assumption that X is constant, which may not hold in Y scenarios. While our results are robust under Z, caution is warranted when extending to W.",
        ],
        "instruction": "書かないと reviewer #2 に掘られる。3-5 件が標準、多すぎると自信欠如印象。",
    },
    "reviewer_2_response_prep": {
        "label": "Reviewer #2 想定反論テンプレ (3 種)",
        "templates": [
            "Anticipated: 「Why only n = X samples?」\n"
            "Pre-emptive response in paper: 「While larger samples would strengthen generalization, n = X is sufficient to detect effects of Cohen's d ≥ 0.5 at α = 0.05 with power 0.8 (see power analysis in §Methods).」",
            "Anticipated: 「Related work incomplete」\n"
            "Pre-emptive response: Add 「Our approach differs from [Smith 2023] in X, from [Jones 2024] in Y」 in Related Work.",
            "Anticipated: 「No comparison with baseline Z」\n"
            "Pre-emptive response: If Z is impractical, add 「We did not compare with Z because of W (unavailability / orthogonal goal / etc.)」 in Limitations.",
        ],
        "instruction": "paper T9 reviewer_2_trigger_summary で flag された弱点に対し、pre-emptive な defense を paper 内に仕込む。",
    },
    "title_claim_form": {
        "label": "Title を claim form に (5 案)",
        "templates": [
            "名詞句タイトル → claim form:\n"
            "「Analysis of X」 → 「X Achieves Y with Z% Improvement」",
            "「A Study of X」 → 「X Enables Y via Z」",
            "「Method for X」 → 「X: A Novel Approach for Y with Z Features」",
            "「Towards X」 → 「X Reduces Y by Z%」",
            "Results-first: 「N-fold Speedup in X via Y」",
        ],
        "instruction": "Nature/Science は claim form、IEEE は名詞句も許容。journal の慣習を事前 check。",
    },
    "intro_banned_opener_fix": {
        "label": "Intro 冒頭禁句 → 問題提起型 (3 案)",
        "templates": [
            "❌ 「In this paper, we propose...」 (Wallwork §14.11 banned)\n"
            "✅ 「Despite significant progress in X, the problem of Y remains unsolved due to Z.」",
            "❌ 「Recently, X has attracted much attention.」\n"
            "✅ 「X technology has reached Y, but a critical gap remains in Z.」",
            "❌ 「This paper studies X.」\n"
            "✅ 「Understanding X is critical for Y, yet current approaches suffer from Z.」",
        ],
        "instruction": "Intro 第 1 文は問題提起 or gap identification。「In this paper」「Recently」「This paper」で始めない。",
    },
    "cover_letter_hook": {
        "label": "Cover letter opening 3 テンプレ",
        "templates": [
            "Dear Editor,\n\nWe are pleased to submit our manuscript titled \"X\" for consideration in [Journal]. This work advances the field of Y by Z, which directly aligns with the journal's scope on [scope keyword].",
            "Dear Editor,\n\nPlease find enclosed our manuscript \"X\" for publication in [Journal]. We believe this submission is suitable for the [special issue/section on Y] because Z.",
            "Dear Prof. [Editor Name],\n\nWe submit the enclosed manuscript \"X\" which we believe contributes to [Journal]'s readership interested in Y. The core novelty is Z, supported by W evidence.",
        ],
        "instruction": "Journal の aims/scope keyword を必ず含める (T15 で特定)。editor 名が分かれば呼称つき。",
    },
}


def paper_writing_rewrite_suggest(
    target: str = "",
    context_text: str = "",
    n_candidates: int = 3,
) -> dict:
    """paper 用 rewrite candidate 生成 (11 target × 3-5 candidate)。"""
    if not target:
        return {
            "available_targets": list(REWRITE_TEMPLATES.keys()),
            "comments": [f"target 指定必須。available: {list(REWRITE_TEMPLATES.keys())}"],
            "source": "paper T19 rewrite suggest",
        }
    if target not in REWRITE_TEMPLATES:
        return {
            "available_targets": list(REWRITE_TEMPLATES.keys()),
            "error": f"unknown target: {target}",
            "source": "paper T19 rewrite suggest",
        }

    spec = REWRITE_TEMPLATES[target]
    candidates = spec["templates"][:n_candidates]

    comments = [
        f"Target: {target} ({spec['label']}) — {len(candidates)} 案",
        f"Instruction: {spec['instruction']}",
    ]
    if context_text:
        comments.append(f"Context: {context_text[:120]}…")
    comments.append(
        "candidates は template。placeholder を context に応じて置換後、"
        "文脈整合を人間が最終確認。"
    )

    return {
        "target": target,
        "label": spec["label"],
        "candidates": candidates,
        "n_candidates": len(candidates),
        "instruction": spec["instruction"],
        "context_provided": bool(context_text),
        "available_targets": list(REWRITE_TEMPLATES.keys()),
        "comments": comments,
        "hint": (
            "paper T19 は grant T50 の paper 版 (11 target)。"
            "IMRAD 各 section で典型的な rewrite 需要をカバー。"
            "T17 next_5_actions と組合わせて、priority action に対応。"
        ),
        "source": (
            "Intelligence Layer S5 (paper 版)。Wallwork + Day & Gastel + "
            "APA 7th + NEJM statistical reporting の標準を template 化。"
            "(paper-writing 2026-04 採用、v0.24.0)。"
        ),
    }
