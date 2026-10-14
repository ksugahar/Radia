"""Plan B T14: Discussion structure 4 elements
(paper_writing_discussion_structure_4_elements).

論文の Discussion section は **Results 再記述で終わる** ことが頻発し、
reviewer の不満理由筆頭の 1 つ。標準的な Discussion は以下 4 要素を含む
(Day & Gastel 7th ed. Chapter 13, Wallwork Ch.13):

  1. **interpretation_implications** — 結果の意味、機構的解釈、broader implications
  2. **limitations** — 制約 / scope / バイアス (既存 T4 と相補)
  3. **future_work** — 次の研究の方向性 / open questions
  4. **generalization_caveats** — 結果の一般化可能性とその限界

設計方針:
- Discussion section を抽出 (regex で section heading)
- 各要素の存在と質 (どの程度書かれているか) を評価
- score は 4 軸の coverage + 各軸の文字数重み
- PRIMARY OUTPUT は per_element findings + missing_elements
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Discussion section extraction
# ------------------------------------------------------------------
_DISCUSSION_PATTERNS = [
    re.compile(
        r"\\section\*?\{(?:Discussion|考察)\}(.*?)"
        r"(?=\\section\b|\\bibliography|\Z)",
        re.DOTALL | re.IGNORECASE,
    ),
    re.compile(
        r"^#{1,3}\s*(?:Discussion|考察)\s*\n+(.+?)(?=\n#{1,3}\s|\Z)",
        re.MULTILINE | re.DOTALL | re.IGNORECASE,
    ),
]


# ------------------------------------------------------------------
# Element keyword patterns
# ------------------------------------------------------------------
_INTERPRETATION_PATTERNS = [
    r"(?:our|the|these)\s+(?:results?|findings?|data)\s+"
    r"(?:suggest|indicate|imply|demonstrate|show)",
    r"(?:this|the)\s+(?:result|finding)\s+can\s+be\s+(?:explained|interpreted)",
    r"mechanism|underlying|explanation|insight",
    r"これらの(?:結果|知見|データ)は[^。]{5,40}(?:示唆|意味|示し)",
    r"機構(?:的)?(?:解釈|考察)",
    r"考えられる(?:理由|要因|機構)",
    r"implication[s]?\s+(?:for|of)",
    r"broader\s+implications?",
    r"impact\s+(?:on|for)",
]

_LIMITATIONS_PATTERNS = [
    r"\b(?:limitation|caveat|shortcoming|weakness|drawback)s?\b",
    r"\bhowever,?\s+(?:our|the|this)\s+(?:study|approach|method)",
    r"\bnotably,?\s+(?:our|we|this)",
    r"制限|限界|制約|不足|問題点",
    r"本研究の(?:限界|課題|制約)",
    r"考慮すべき(?:点|事項)",
    r"未解決の(?:問題|課題)",
]

_FUTURE_WORK_PATTERNS = [
    r"future\s+(?:work|research|direction|study|studies|investigation)",
    r"further\s+(?:research|investigation|study|work)",
    r"next\s+step[s]?",
    r"open\s+question[s]?",
    r"今後の(?:課題|研究|展開|方向|展望)",
    r"次の(?:段階|ステップ|展開)",
    r"将来(?:的)?(?:に|な)\s*[^。]{5,30}(?:検討|展開|拡張)",
    r"今後[^。]{0,30}(?:検討|拡張|発展)",
]

_GENERALIZATION_PATTERNS = [
    r"\b(?:may|might|could)\s+(?:not\s+)?generaliz",
    r"\bgeneraliz(?:ability|able|ation)",
    r"\bextend(?:s|ed|ing)?\s+to\s+other",
    r"\bapplicable\s+to\s+other",
    r"\bsetting[s]?\s+(?:beyond|other than)",
    r"一般化|汎化|外挿|適用範囲|適用可能性",
    r"他の(?:設定|条件|分野|領域)で",
    r"より広い(?:文脈|範囲|条件)",
    r"caveat[s]?\s+(?:for|in)\s+generaliz",
]


def _extract_discussion(text: str) -> str:
    for pat in _DISCUSSION_PATTERNS:
        m = pat.search(text)
        if m:
            return m.group(1).strip()
    return ""


def _check_element(text: str, patterns: list[str]) -> dict:
    hits: list[str] = []
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            hits.append(m.group(0))
            if len(hits) >= 6:
                break
    return {
        "hit_count": len(hits),
        "examples": hits[:3],
        "covered": len(hits) >= 1,
        "well_developed": len(hits) >= 2,
    }


def paper_writing_discussion_structure_4_elements(text: str) -> dict:
    """Discussion section の 4 要素 (interpretation / limitations /
    future_work / generalization) の存在と質を診断。

    Args:
        text: 論文全文 (LaTeX or markdown)

    Returns:
        dict: score / discussion_chars / per_element / comments
    """
    if text is None:
        text = ""

    discussion = _extract_discussion(text)

    if not discussion:
        return {
            "score": None,
            "discussion_detected": False,
            "comments": [
                "Discussion / 考察 section が検出できなかった。"
                "\\section{Discussion} or '## Discussion' / '## 考察' を "
                "明示するか、Conclusion と統合されている場合は別の "
                "tool (T10 triangle 等) で確認。"
            ],
            "hint": "section 抽出失敗。",
            "source": "T14 discussion structure (paper-writing 2026-04)",
        }

    elements = {
        "interpretation_implications": _check_element(
            discussion, _INTERPRETATION_PATTERNS
        ),
        "limitations": _check_element(discussion, _LIMITATIONS_PATTERNS),
        "future_work": _check_element(discussion, _FUTURE_WORK_PATTERNS),
        "generalization_caveats": _check_element(
            discussion, _GENERALIZATION_PATTERNS
        ),
    }

    # ---- score ----
    n_covered = sum(1 for e in elements.values() if e["covered"])
    n_well_developed = sum(
        1 for e in elements.values() if e["well_developed"]
    )
    # 4 covered, well_developed が bonus
    score = (n_covered / 4 * 7) + (n_well_developed / 4 * 3)
    score = round(min(10.0, score), 1)

    # ---- comments ----
    comments: list[str] = []

    discussion_chars = len(discussion)
    comments.append(
        f"Discussion 検出: {discussion_chars} 文字。"
        f"4 要素中 {n_covered} 軸 covered, {n_well_developed} 軸 well_developed。"
    )

    not_covered = [k for k, v in elements.items() if not v["covered"]]

    if "interpretation_implications" in not_covered:
        comments.append(
            "⚠ Interpretation/Implications が無い。"
            "Discussion を Results 再記述で終わらせがちな典型 NG。"
            "「結果は X を示唆する (mechanism は Y)」「broader implications: Z」"
            "の段落を 1 つ加える。"
        )

    if "limitations" in not_covered:
        comments.append(
            "⚠ Limitations 言及なし。"
            "reviewer は『書き手が認識している limitation』を期待する。"
            "未記載は『盲点なのか隠してるのか』疑問符が付く。"
            "(既存 T4 limitation_statement_presence と整合確認)。"
        )

    if "future_work" in not_covered:
        comments.append(
            "⚠ Future work 言及なし。"
            "「次にやるべきこと」「open questions」の 1 段落で論文の "
            "学術的意義を広げる。Conclusion 直前が定位置。"
        )

    if "generalization_caveats" in not_covered:
        comments.append(
            "⚠ Generalization caveats なし。"
            "「本結果は X 条件でのみ成立」「他の Y 設定では未検証」のような "
            "適用範囲の明示が不足。reviewer の『どこまで一般化できるか』"
            "への先回り回答。"
        )

    if discussion_chars < 1500:
        comments.append(
            f"Discussion が {discussion_chars} 字と短い。"
            "標準的な journal 論文では 1500-4000 字程度。"
            "短い場合は 4 要素の各々が薄い可能性。"
        )

    if n_covered == 4 and n_well_developed >= 3:
        comments.append(
            f"Discussion 4 要素揃い、score {score}/10 — "
            "標準的な論文 Discussion 構造を充足。"
        )

    comments = comments[:7]

    hint = (
        "Discussion を抽出して 4 要素の keyword detection。"
        "well_developed = 同要素の hit が 2 件以上 (1 件だけだと触れた程度)。"
        "既存 T4 (limitation_statement_presence) と相補的 — T4 は paper 全体、"
        "T14 は Discussion section 内に限定して厳格に評価。"
        "純技法系論文 (proof-of-concept, theorem paper) では generalization "
        "が薄くてもよい場合あり — context に応じて判断。"
    )

    return {
        "score": score,
        "score_max": 10,
        "discussion_detected": True,
        "discussion_chars": discussion_chars,
        "discussion_head": discussion[:300],
        "n_covered": n_covered,
        "n_well_developed": n_well_developed,
        "per_element": elements,
        "comments": comments,
        "hint": hint,
        "source": (
            "Day & Gastel 'How to Write and Publish a Scientific Paper' "
            "(7th ed.) Chapter 13 + Wallwork 'English for Writing Research "
            "Papers' Chapter 13 — Discussion section standard 4-element "
            "structure (interpretation / limitations / future / generalization)。"
            "(paper-writing 2026-04 採用)。"
        ),
    }
