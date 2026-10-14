"""Plan B T12: Statistical reporting compliance check
(paper_writing_statistical_reporting_compliance).

実験系・臨床系論文では p 値単独報告が長年問題視されており、近年の
APA 7th / NEJM / BMJ guidelines では **p 値 + effect size + 信頼区間
+ サンプルサイズ** の 4 点セット報告が標準。CONSORT/STROBE/PRISMA も
同様の要請。本ツールは draft 中の統計報告箇所をスキャンし、これら 4 要素
の充足度を診断する。Reviewer #2 トリガー筆頭。

検出する 4 要素 (per p 値 occurrence):
  1. **p 値** — p = 0.03 / p < 0.05 / p > .05
  2. **effect size** — Cohen's d / OR / RR / η² / r / β / Δ など
  3. **confidence interval** — 95% CI [x, y] / 95% CI: x to y
  4. **sample size** — n = X / N = X / per group n / total n

設計方針:
- 各 p 値検出ごとに「半径 200 字以内」に他の 3 要素があるかを確認
- score は 「効果量 + CI + n が揃った p 値 / 全 p 値」の比率 × 10
- PRIMARY OUTPUT は per_p_value findings + missing per occurrence
- 純理論系では p 値が出ないので score=10 (該当なし) を返す
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# パターン定義
# ------------------------------------------------------------------
_P_VALUE_RE = re.compile(
    r"\b(?:p|P)\s*"
    r"(?:[<>=≤≥]|＝|＜|＞)\s*"
    r"(0?\.\d+|\d\.?\d*[eE][-+]?\d+|0\.\d+)(?![\d.])"
    r"(?!\s*(?:W|V|A|Hz|J|T|K|H|F|Pa|Ω|ohm)\b)"
)

_EFFECT_SIZE_PATTERNS = [
    r"\b(?:Cohen'?s?\s*d|d\s*=)",
    r"\b(?:OR|odds\s+ratio)\s*[=:]",
    r"\b(?:RR|relative\s+risk)\s*[=:]",
    r"\b(?:HR|hazard\s+ratio)\s*[=:]",
    r"\bη\s*[2²]?\s*[=:]",  # eta squared
    r"\b(?:partial\s+)?eta(?:[\s\-]squared)?\s*[=:]",
    r"\bω\s*[2²]?\s*[=:]",
    r"\b(?:Cohen'?s?\s*)?f\s*[2²]?\s*[=:]",
    r"\bφ\s*[=:]",  # phi
    r"\br\s*[=:]\s*[\-]?0?\.\d",  # correlation coef
    r"\bρ\s*[=:]",  # spearman rho
    r"\bτ\s*[=:]",  # kendall tau
    r"\bβ\s*[=:]\s*[\-]?\d",  # standardized beta
    r"\bΔ\s*[=:]",  # delta (often effect)
    r"\beffect\s+size",
    r"効果量",
    r"効果サイズ",
    # Mean difference / SMD
    r"\b(?:mean\s+diff|MD)\s*[=:]",
    r"\bSMD\s*[=:]",
]

_CI_PATTERNS = [
    r"\b9[05]\s*%\s*CI\s*[\[(]\s*[\-]?\d",
    r"\b9[05]\s*%\s*CI\s*[:：]\s*[\-]?\d",
    r"\bconfidence\s+interval",
    r"\b信頼区間",
    r"\bCI\s*=?\s*\[",
    r"\b9[05]\s*%\s*confidence",
    # bracketed numeric range often denotes CI: [1.2, 3.4]
    r"\[\s*[\-]?\d+\.?\d*\s*,\s*[\-]?\d+\.?\d*\s*\]",
]

_SAMPLE_SIZE_PATTERNS = [
    r"\bn\s*=\s*\d+",
    r"\bN\s*=\s*\d+",
    r"\bn\s*=\s*\d+\s*per\s+group",
    r"\btotal\s+(?:n|N)\s*=?\s*\d+",
    r"\bsample\s+size\s+(?:of|=|:|：)\s*\d+",
    r"被験者\s*\d+\s*名",
    r"サンプル\s*(?:数|サイズ)\s*[=:：]\s*\d+",
    r"参加者\s*\d+\s*名",
    r"\b\d+\s+participants?",
    r"\b\d+\s+subjects?",
]


def _has_pattern_in_window(
    patterns: list[str], window: str
) -> tuple[bool, list[str]]:
    hits: list[str] = []
    for pat in patterns:
        for m in re.finditer(pat, window, re.IGNORECASE):
            hits.append(m.group(0))
    return (len(hits) > 0, hits[:3])


def paper_writing_statistical_reporting_compliance(
    text: str, window_chars: int = 250
) -> dict:
    """各 p 値の周辺で effect size / CI / sample size の有無を診断。

    APA 7th / NEJM / CONSORT が要請する p 値 + 3 要素の compliance を
    各 p 値 occurrence ごとに評価。

    Args:
        text: 論文全文
        window_chars: p 値の前後 N 字を「同 statement」とみなす範囲 (default 250)

    Returns:
        dict: score / total_p_values / per_p_value / summary / comments
    """
    if text is None:
        text = ""

    p_value_matches = list(_P_VALUE_RE.finditer(text))
    n_p = len(p_value_matches)

    if n_p == 0:
        return {
            "score": 10.0,
            "score_max": 10,
            "total_p_values": 0,
            "per_p_value": [],
            "summary": {
                "all_4_complete": 0,
                "p_only": 0,
                "p_with_some": 0,
            },
            "comments": [
                "p 値の使用を検出せず。純理論系・解析系論文の可能性。"
                "本ツールは実験/臨床系のための診断なので適用対象外。"
            ],
            "hint": (
                "p 値が無ければ統計報告 compliance は不要。"
                "解析系なら理論的境界・収束性などの別軸 (T3 claim_quantification) を。"
            ),
            "source": "T12 statistical reporting (paper-writing 2026-04)",
        }

    per_p_value: list[dict] = []
    summary = {"all_4_complete": 0, "p_only": 0, "p_with_some": 0}

    for m in p_value_matches:
        center = m.start()
        win_start = max(0, center - window_chars)
        win_end = min(len(text), center + window_chars)
        window = text[win_start:win_end]

        has_effect, effect_hits = _has_pattern_in_window(
            _EFFECT_SIZE_PATTERNS, window
        )
        has_ci, ci_hits = _has_pattern_in_window(_CI_PATTERNS, window)
        has_n, n_hits = _has_pattern_in_window(_SAMPLE_SIZE_PATTERNS, window)

        components_present = sum([has_effect, has_ci, has_n])
        completeness = (1 + components_present) / 4  # p value itself counts

        if components_present == 3:
            summary["all_4_complete"] += 1
            verdict = "COMPLETE"
        elif components_present == 0:
            summary["p_only"] += 1
            verdict = "P_ONLY"
        else:
            summary["p_with_some"] += 1
            verdict = "PARTIAL"

        per_p_value.append({
            "pos": m.start(),
            "p_value_text": m.group(0),
            "context": text[
                max(0, m.start() - 60):
                min(len(text), m.end() + 60)
            ].replace("\n", " "),
            "has_effect_size": has_effect,
            "effect_size_hits": effect_hits,
            "has_ci": has_ci,
            "ci_hits": ci_hits,
            "has_sample_size": has_n,
            "sample_size_hits": n_hits,
            "completeness_ratio": round(completeness, 2),
            "verdict": verdict,
        })

    # Score: average completeness × 10
    if per_p_value:
        avg = sum(p["completeness_ratio"] for p in per_p_value) / len(
            per_p_value
        )
        score = round(avg * 10, 1)
    else:
        score = 10.0

    # ---- comments ----
    comments: list[str] = []

    comments.append(
        f"検出した p 値: {n_p} 件。"
        f"COMPLETE (p+effect+CI+n) {summary['all_4_complete']} / "
        f"PARTIAL {summary['p_with_some']} / "
        f"P_ONLY {summary['p_only']}。"
    )

    if summary["p_only"] >= 1:
        comments.append(
            f"P_ONLY (p 値単独報告) が {summary['p_only']} 件。"
            "APA 7th / NEJM / CONSORT 違反。reviewer #2 トリガー筆頭。"
            "「p < 0.05」だけでは effect の大きさが不明 — Cohen's d / "
            "OR / RR / 95% CI / n を必ず添える。"
        )

    if summary["p_with_some"] >= 1:
        comments.append(
            f"PARTIAL (3 要素のうち一部欠落) が {summary['p_with_some']} 件。"
            "per_p_value の verdict='PARTIAL' を確認し、欠落要素を補強。"
        )

    if summary["all_4_complete"] == n_p and n_p > 0:
        comments.append(
            f"全 {n_p} 件の p 値で 4 要素 (p / effect / CI / n) が揃って "
            "報告されている。CONSORT/STROBE 準拠。"
        )

    # 多重比較補正の hint (>= 5 p 値だと multiple comparison が話題に)
    if n_p >= 5:
        if not re.search(
            r"Bonferroni|Holm|Benjamini|FDR|adjust(?:ed)?\s+p|"
            r"multiple\s+(?:test|comparison)|多重比較|多重検定",
            text, re.IGNORECASE,
        ):
            comments.append(
                f"p 値 {n_p} 件と多めだが、多重比較補正 "
                "(Bonferroni / Holm / FDR 等) の言及なし。"
                "5 件以上の検定で補正なしは reviewer #2 の指摘対象。"
            )

    if score < 5:
        comments.append(
            f"compliance score {score}/10 — desk reject or major revision "
            "リスク高。statistical reporting を全面的に補強推奨。"
        )

    comments = comments[:7]

    hint = (
        "window_chars (default 250) の範囲で「同 statement」を判定するため、"
        "離れた表中に CI/n が書かれていても拾えないことがある。"
        "verdict='P_ONLY' は false positive 可能性 (表参照の場合) があるので "
        "context を確認してから判断。理論系論文では適用外。"
    )

    return {
        "score": score,
        "score_max": 10,
        "total_p_values": n_p,
        "per_p_value": per_p_value[:30],  # Cap at 30 to keep response size sane
        "summary": summary,
        "comments": comments,
        "hint": hint,
        "source": (
            "APA Publication Manual 7th ed. (2020) Chapter 6 + "
            "NEJM Statistical Reporting Guidelines + "
            "CONSORT 2010 + STROBE Statement (2007) — "
            "p value alone insufficient; effect size + CI + n required。"
            "(paper-writing 2026-04 採用)。"
        ),
    }
