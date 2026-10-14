"""直訳調・AI調の日本語を見つける共通検査。

2026-09-02、基盤C計画調書を人手で通読したところ、機械検査（bedrock、誤用、
表記ゆれ）をすべて通過した原稿に、英語から直訳したような言い回しと文法上の
誤りが 9 か所残っていた。代表例は次のとおり。

- 「判定則を加速器電磁石設計へ発展する」  -- 「発展する」は自動詞。他動詞は
  「発展させる」。英語 develop をそのまま写した形。
- 「接着層（glue）である」  -- 日本語の術語に英単語を括弧で添える注記。
  審査者に英語を読ませる理由がない。略語（MCP 等）の定義とは別。
- 「劇的な差」「設計判断を保存できる」「判定則を移転できる」「区間へ戻す」
  -- dramatic / preserve / transfer / feed back の直訳。

同じ検査を paper-writing、grant-writing、presentation の三つが共有する。
点数は機械的な密度であって文章の良し悪しではない。指摘は候補であり、
定義語（「反証」「凍結」等）や数学の用語（「写す」）は残してよい。
"""
from __future__ import annotations

import re

from .language import is_japanese

# 自動詞を他動詞のように使った誤り。develop / improve / spread を「〜を発展する」
# と写すと生じる。「〜を発展させる」「〜を向上させる」が正しい。目的語と動詞の
# 間に「〜へ」「〜に」の到達点が挟まる形（「判定則を…設計へ発展する」）も拾う。
# 「が」「は」「で」を挟むと別の主語の自動詞になり得るので、そこで打ち切る。
_INTRANSITIVE_AS_TRANSITIVE = re.compile(
    r"を(?:[^、。をがはでも]{0,20}[へに])?"
    r"(発展|向上|進化|成立|充実|普及|増加|減少|低下|上昇|定着|浸透|拡散|安定|深化|進展)"
    r"(?:する|した|し(?=[、て]))"
)

# 「接着層（glue）」のように、日本語の語の直後に小文字英単語を括弧で添えた注記。
# 大文字の略語定義「知識・実行インターフェース（MCP）」は対象にしない。
_ENGLISH_GLOSS = re.compile(
    r"([一-鿿ぁ-ゟ゠-ヿー]{2,})[（(]([a-z][a-z \-]{2,})[）)]"
)

# 直訳定型句。左が検出パターン、中が短い理由、右が置き換えの候補。
_CALQUES: list[tuple[str, str, str]] = [
    (r"劇的(?:な|に)", "dramatic の直訳", "大きな／著しく"),
    (r"を可能に(?:する|し)", "enable の直訳", "〜できるようにする／〜できる"),
    (r"に焦点を当て", "focus on の直訳", "〜を対象にする／〜に絞る"),
    (r"に対処(?:する|し)", "address の直訳", "〜に対応する／〜を解決する"),
    (r"重要な役割を(?:果た|担)", "play an important role の直訳", "〜に欠かせない／〜を担う"),
    (r"という事実", "the fact that の直訳", "削除して事実をそのまま書く"),
    (r"であることに注意", "note that の直訳", "〜に注意する"),
    (r"することが重要である", "it is important to の直訳", "〜する必要がある／〜する"),
    (r"の一つである", "one of the の直訳", "削除するか「〜である」"),
    (r"包括的(?:な|に)", "comprehensive の直訳", "全体を／幅広く"),
    (r"潜在的(?:な|に)", "potential の直訳", "〜し得る／将来の"),
    (r"を確実に(?:する|し)", "ensure の直訳", "〜を保証する／必ず〜する"),
    (r"することが可能(?:となる|である|になる)", "it becomes possible の直訳", "〜できる"),
    (r"を活用(?:する|し)", "leverage の直訳", "〜を使う／用いる"),
    (r"を最大限に", "maximize の直訳", "できるだけ／十分に"),
    (r"(?:シームレス|ロバスト|スケーラブル|エコシステム|レバレッジ|インサイト|ベストプラクティス)",
     "英語をカタカナにしただけ", "継ぎ目なく／頑健／拡張できる／生態系／知見／定石"),
]

# 生成AIの日本語に多い、内容を伴わない強調語。数だけ返し、点数には軽く入れる。
_AI_TICS: list[tuple[str, str]] = [
    (r"本質的に", "本質的に"),
    (r"根本的に", "根本的に"),
    (r"多角的|多面的", "多角的／多面的"),
    (r"を促進(?:する|し)", "〜を促進する"),
    (r"に寄与(?:する|し)", "〜に寄与する"),
    (r"を強化(?:する|し)", "〜を強化する"),
    (r"効果的に", "効果的に"),
]


def _split_sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[。．！？!?])|\n+", text) if s.strip()]


def _excerpt(sentence: str, match: re.Match[str], width: int = 40) -> str:
    start = max(0, match.start() - width)
    end = min(len(sentence), match.end() + width)
    head = "…" if start > 0 else ""
    tail = "…" if end < len(sentence) else ""
    return head + sentence[start:end] + tail


def check_translationese(text: str) -> dict:
    """直訳調・AI調の候補を、文法誤り・英語注記・定型句・強調語の順に返す。"""
    if not text or not is_japanese(text):
        return {
            "applicable": False,
            "score": None,
            "risk_count": 0,
            "risks": [],
            "comments": [],
            "recommendations": [],
            "target": "日本語として書かれた文（直訳と分かる語順・語彙・自他動詞の誤りがない）",
            "source": "translationese check (2026-09-02 基盤C通読)",
        }

    sentences = _split_sentences(text)
    risks: list[dict] = []

    grammar_examples = []
    for sentence in sentences:
        for match in _INTRANSITIVE_AS_TRANSITIVE.finditer(sentence):
            verb = match.group(1)
            grammar_examples.append({
                "excerpt": _excerpt(sentence, match),
                "found": match.group(0),
                "suggestion": f"を{verb}させる",
            })
    if grammar_examples:
        risks.append({
            "type": "intransitive_verb_used_transitively",
            "severity": "HIGH",
            "count": len(grammar_examples),
            "comment": "自動詞を「〜を…する」の形で使っている。他動詞は「〜させる」。",
            "recommendation": "「発展する」「向上する」は自動詞。「発展させる」「向上させる」にする。",
            "examples": grammar_examples[:8],
        })

    gloss_examples = []
    for sentence in sentences:
        for match in _ENGLISH_GLOSS.finditer(sentence):
            gloss_examples.append({
                "excerpt": _excerpt(sentence, match),
                "term": match.group(1),
                "gloss": match.group(2),
            })
    if gloss_examples:
        risks.append({
            "type": "english_gloss_after_japanese_term",
            "severity": "MEDIUM",
            "count": len(gloss_examples),
            "comment": "日本語の術語に英単語を括弧で添えている。略語定義ではないので読者に英語を読ませる理由がない。",
            "recommendation": "括弧の英語を外す。訳語が定着していないなら日本語で一度説明する。",
            "examples": gloss_examples[:8],
        })

    calque_hits: list[dict] = []
    for pattern, why, better in _CALQUES:
        regex = re.compile(pattern)
        examples = []
        for sentence in sentences:
            for match in regex.finditer(sentence):
                examples.append(_excerpt(sentence, match))
        if examples:
            calque_hits.append({
                "pattern": pattern,
                "why": why,
                "better": better,
                "count": len(examples),
                "examples": examples[:4],
            })
    if calque_hits:
        risks.append({
            "type": "calque_phrase",
            "severity": "LOW",
            "count": sum(hit["count"] for hit in calque_hits),
            "comment": "英語の定型句をそのまま写した言い回し。",
            "recommendation": "各候補の「better」欄の言い回しに置き換えるか、具体的な動作や量で書く。",
            "examples": calque_hits[:12],
        })

    tic_hits: list[dict] = []
    for pattern, label in _AI_TICS:
        count = len(re.findall(pattern, text))
        if count:
            tic_hits.append({"label": label, "count": count})
    if tic_hits:
        risks.append({
            "type": "ai_emphasis_vocabulary",
            "severity": "LOW",
            "count": sum(hit["count"] for hit in tic_hits),
            "comment": "内容を伴わない強調語。生成AIの日本語に多い。",
            "recommendation": "何がどう変わるかを量や動作で書けるなら、その語は外す。",
            "examples": tic_hits,
        })

    score = 10.0
    score -= 2.0 * len(grammar_examples)
    score -= 1.0 * len(gloss_examples)
    score -= 0.5 * len(calque_hits)
    score -= 0.25 * len(tic_hits)
    score = round(max(0.0, score), 1)

    return {
        "applicable": True,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": [risk["comment"] for risk in risks],
        "recommendations": [risk["recommendation"] for risk in risks],
        "sentence_count": len(sentences),
        "target": "日本語として書かれた文（直訳と分かる語順・語彙・自他動詞の誤りがない）",
        "source": "translationese check (2026-09-02 基盤C通読)",
    }
