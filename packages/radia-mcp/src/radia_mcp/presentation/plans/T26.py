"""Presentation T26 (Intelligence Layer S5): Rewrite suggest (presentation 用)
(presentation_rewrite_suggest).
"""

from __future__ import annotations


REWRITE_TEMPLATES: dict[str, dict] = {
    "slide_title_claim_form": {
        "label": "Slide title 名詞句 → claim form (5 例)",
        "templates": [
            "「Results」 → 「PEEC Achieves 4× Faster Simulation」",
            "「Method」 → 「Surface Impedance Method: 10% Accuracy Improvement」",
            "「Discussion」 → 「Key Insight: Domain Y Dominates Above 100kHz」",
            "「Introduction」 → 「Why EM Analysis Matters for AI Era」",
            "「Conclusions」 → 「Take-home: PEEC is 4× Faster at Same Accuracy」",
        ],
        "instruction": "slide title は動詞 + 数値 の claim form で書く (skill.md L98)。",
    },
    "take_home_one_sentence": {
        "label": "Take-home 1 文 (5 案)",
        "templates": [
            "本研究は X を Y で実現し、Z に波及する。",
            "Take-home: X is Y times faster with Z% better accuracy.",
            "X を、Y で、N 倍に。本研究の核はこの 1 行。",
            "X 問題 → Y 解決 → Z 展開。",
            "If you remember one thing: X solves Y with Z.",
        ],
        "instruction": "末尾 slide または Hook 直後に 1 文で take-home を明示。",
    },
    "opening_hook_3sec": {
        "label": "3 秒 Hook (5 案)",
        "templates": [
            "【問題提起】電磁界解析は AI 時代に...の課題に直面している。",
            "【意外性】これまで不可能と考えられてきた X が、実は可能。",
            "【数字先行】N 倍速い電磁界解析を、OSS で実現。",
            "【場面】ある現場で、3 日かけて 1 通りしか試行できなかった。本研究は...",
            "【質問型】なぜ電磁界解析は国産化されないのか? 答えは...",
        ],
        "instruction": "冒頭 2 slide で reviewer の注意を掴む。TED 風ではなく理系プレゼン向け。",
    },
    "qa_backup_slide_template": {
        "label": "Q&A backup slide テンプレ (3 案)",
        "templates": [
            "Backup: 前提条件詳細\n"
            "- 使用装置: X 社 Y モデル (校正日 YYYY-MM-DD)\n"
            "- 環境: 温度 T°C、湿度 H%\n"
            "- 被験者 / サンプル: n = N, selection criteria: ...\n"
            "- 計算機 spec: CPU, RAM, GPU",
            "Backup: 統計詳細\n"
            "- 検定手法: X-test, α = 0.05, Bonferroni correction\n"
            "- Power analysis: 1-β = 0.80 at d = 0.5\n"
            "- Effect size: Cohen's d = X, 95% CI = [Y, Z]",
            "Backup: Related work 比較表\n"
            "| Method | Year | Metric 1 | Metric 2 | Note |\n"
            "| Smith  | 2023 | 85%      | 100ms    | ... |\n"
            "| Ours   | 2026 | 92%      | 30ms     | +7% @ 3.3× |",
        ],
        "instruction": "Q&A で頻出する質問 3-5 種に対応する backup slide を事前準備。",
    },
    "mini_imrad_phase_skeleton": {
        "label": "Mini-IMRAD 7 phase 骨格",
        "templates": [
            "Slide 1: Title (animated with author + date)\n"
            "Slide 2-3: Background (既存手法の限界)\n"
            "Slide 4: Problem / Research Question (1 主張 slide)\n"
            "Slide 5-8: Methods (4 slide 以内)\n"
            "Slide 9-14: Results (各 slide 1 figure + 1 claim)\n"
            "Slide 15-16: Discussion (interpretation + limitations)\n"
            "Slide 17: Conclusion + Take-home\n"
            "Slide 18-20: Backup (Q&A 用)",
        ],
        "instruction": "20 分発表の標準配分。15 分なら Methods / Results を圧縮。",
    },
    "figure_slide_compliance_template": {
        "label": "Figure slide 5 要素 skeleton",
        "templates": [
            "Slide title (claim form): 「X achieves Y at Z conditions」\n"
            "[Figure with axis labels, units, legend, error bars]\n"
            "  - x 軸: 条件 [単位]\n"
            "  - y 軸: 結果 [単位]\n"
            "  - 凡例: Method A, Method B, Baseline\n"
            "  - Error bar: ±SD (n=N)\n"
            "Caption (下): Fig. X: 説明文 (出典 [N] or 自作)",
        ],
        "instruction": "理系プレゼンの figure slide 5 要素 (軸 + 単位 + 凡例 + caption + 出典) を全部含める。",
    },
    "equation_slide_compliance_template": {
        "label": "Equation slide 4 要素 skeleton",
        "templates": [
            "Slide title: 「Surface Impedance Method: Governing Equation」\n\n"
            "    ∂E/∂t + (1/σ) ∇·J = 0   ... (1)\n\n"
            "Where (記号定義):\n"
            "  - E: 電場 [V/m]\n"
            "  - σ: 電気伝導率 [S/m]\n"
            "  - J: 電流密度 [A/m²]\n"
            "Physical meaning: (物理意味)\n"
            "  この式は電磁誘導による電場変化と電流の関係を表す。\n"
            "Reference: Stratton (1941) Eq. (3.42) [citation]",
        ],
        "instruction": "数式 slide に 記号定義 + 物理意味 + 次元 + 出典 の 4 要素必須。",
    },
    "results_slide_stat_template": {
        "label": "Results slide 4 要素 skeleton",
        "templates": [
            "Slide title (claim form): 「PEEC is 4× Faster with Equivalent Accuracy」\n\n"
            "[Bar chart: PEEC vs FEM vs Baseline]\n"
            "  - y 軸: Simulation time [s]\n"
            "  - Error bars: ±SD (n=10 runs per condition)\n"
            "  - Annotation: *** p<0.001 (paired t-test)\n\n"
            "Results: PEEC 125 ± 12 s vs FEM 512 ± 45 s\n"
            "  - 4.1× speedup (95% CI: 3.8 to 4.4), Cohen's d = 12.4\n"
            "  - Accuracy RMSE < 0.2% (within measurement uncertainty)",
        ],
        "instruction": "Results slide に error bar + n + 統計検定 + effect size の 4 要素 (paper T12 相当)。",
    },
    "discussion_slide_4elements": {
        "label": "Discussion slide 4 要素 skeleton",
        "templates": [
            "Slide: Discussion (2 枚に分けて)\n\n"
            "Slide 1 (Interpretation + Implications):\n"
            "  - Our result shows X because of mechanism Y\n"
            "  - Broader implications: Z field benefits from W\n\n"
            "Slide 2 (Limitations + Future):\n"
            "  - Limitations: P (sample size / specific condition)\n"
            "  - Future work: Q (extend to R)\n"
            "  - Generalization caveats: only S if T holds",
        ],
        "instruction": "Discussion は interpretation / limitations / future / generalization の 4 要素。Results 再記述で終わらせない。",
    },
}


def presentation_rewrite_suggest(
    target: str = "",
    context_text: str = "",
    n_candidates: int = 3,
) -> dict:
    if not target:
        return {
            "available_targets": list(REWRITE_TEMPLATES.keys()),
            "comments": [f"target 指定必須。available: {list(REWRITE_TEMPLATES.keys())}"],
            "source": "presentation T26 rewrite suggest",
        }
    if target not in REWRITE_TEMPLATES:
        return {
            "available_targets": list(REWRITE_TEMPLATES.keys()),
            "error": f"unknown target: {target}",
            "source": "presentation T26 rewrite suggest",
        }
    spec = REWRITE_TEMPLATES[target]
    candidates = spec["templates"][:n_candidates]
    comments = [
        f"Target: {target} ({spec['label']}) — {len(candidates)} 案",
        f"Instruction: {spec['instruction']}",
    ]
    if context_text:
        comments.append(f"Context: {context_text[:120]}…")
    return {
        "target": target,
        "label": spec["label"],
        "candidates": candidates,
        "n_candidates": len(candidates),
        "instruction": spec["instruction"],
        "context_provided": bool(context_text),
        "available_targets": list(REWRITE_TEMPLATES.keys()),
        "comments": comments,
        "hint": "presentation T26 は grant T50 / paper T19 の pptx 版。9 target。",
        "source": "Intelligence Layer S5 (presentation 版)。(v0.25.0)",
    }
