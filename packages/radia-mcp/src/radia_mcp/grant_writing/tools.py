"""Tool functions for grant proposal writing.

All functions are plain callables; ``server.py`` wraps them as MCP tools.
The Japanese technical-prose lint helpers are re-exported from the
grant-writing implementation that was already preserved inside
``paper_writing._ja_lint`` during the public radia-mcp migration.
"""
from __future__ import annotations

import pathlib
import re

from radia_mcp.paper_writing._ja_lint import (
    grant_writing_acronym_usage_audit as _ja_acronym_usage_audit,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_check_kanji_ratio as _ja_check_kanji_ratio,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_check_misuse_japanese as _ja_check_misuse_japanese,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_check_notation_variants as _ja_check_notation_variants,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_check_subject_predicate_distance as _ja_subject_predicate_distance,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_find_undefined_acronyms as _ja_find_undefined_acronyms,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_lint_bedrock as _ja_lint_bedrock,
)
from radia_mcp.paper_writing._ja_lint import (
    grant_writing_suggest_redundancy_fixes as _ja_suggest_redundancy_fixes,
)

from .._shared.hedges import HEDGE_PATTERNS, scan_hedges

_HERE = pathlib.Path(__file__).resolve().parent


def _load_skill() -> str:
    return (_HERE / "skill.md").read_text(encoding="utf-8")


def _read_text_if_path(text_or_path: str) -> str:
    """Treat a short existing .md/.tex/.txt path as a file, else as text."""
    if text_or_path is None:
        return ""
    s = str(text_or_path)
    if "\n" in s or "\r" in s:
        return s
    if len(s) > 320:
        return s
    try:
        p = pathlib.Path(s)
        if p.exists() and p.suffix.lower() in {".md", ".tex", ".txt"}:
            payload = p.read_bytes()
            for encoding in ("utf-8-sig", "cp932"):
                try:
                    return payload.decode(encoding)
                except UnicodeDecodeError:
                    continue
            return payload.decode("utf-8", errors="replace")
    except OSError:
        pass
    return s


# Sentences the form wrote, not the applicant. Every marker was taken from the
# instruction paragraphs of two real 科研費 forms (令和2年度 S-14 and 平成28年度
# S-1-8); together they match all of them and none of the applicant prose in
# either document. Instructions are 10% of one form and 40% of the other.
_FORM_INSTRUCTION = re.compile(
    r"本欄に[はも]|本欄は|"
    r"記述すること|記入すること|記入の上|"
    r"(?:記述|記入|参照|選択|要約|確認|留意|注意)して\s*(?:ください|下さい)|"
    r"公募要領|記入要領|作成・記入要領|審査されます|"
    r"て\s*も\s*可(?![能])|てもよい|ても構いません|"
    r"空欄のまま|記述欄を削除|記入欄を削除|"
    # Older PDF forms wrap this ethics-box example across several physical
    # lines. The fragments do not retain the polite sentence ending, so the
    # general form-voice rule below cannot identify them on its own.
    r"例えば[、，,]?\s*個人情報を伴う|"
    r"個人情報を伴うアンケート調査[・･]インタビュー調査|"
    r"(?:研究)?助成に関するアンケート|"
    r"承認手続[きがの]?必要となる(?:調査|研究|実験)"
)

# ですます調. The form speaks this way; a proposal body does not.
_POLITE_ENDING = re.compile(r"(?:ます|です|ません|ましょう|ください)[。．]?$")

_CITATION_YEAR = re.compile(r"(?:19|20)\d{2}")


def _strip_citation_items(match: re.Match[str]) -> str:
    """Drop citation entries from a list and keep the remaining items apart.

    An item is a citation when it carries a four-digit year and at least two
    commas, which is how a publication list reads in every form this suite has
    seen. The year is not required to be parenthesised: an accepted paper is
    listed as ``IGTE Symposium 2026 (accepted)``.

    Surviving items are terminated with a Japanese full stop so the sentence
    splitter cannot fuse consecutive bullets into one pseudo-sentence. An
    English period does not end a sentence for that splitter, so it does not
    count as an existing terminator here.
    """
    kept: list[str] = []
    for item in re.split(r"\\item\b", match.group("body"))[1:]:
        stripped = item.strip()
        if not stripped:
            continue
        if _CITATION_YEAR.search(stripped) and stripped.count(",") >= 2:
            continue
        if not stripped.endswith(("。", "．", "！", "？")):
            stripped += "。"
        kept.append(stripped)
    return " " + " ".join(kept) + " " if kept else " "


def _prose_for_lint(text: str) -> str:
    """Remove common LaTeX scaffolding before prose-oriented checks.

    Section-presence and program checks still inspect the original source. This
    normalization is only for sentence length, hedge, and Japanese prose lint;
    otherwise template comments and commands are reported as applicant prose.

    A form's printed instruction text is stripped whether or not the document
    is LaTeX. A proposal extracted from Word carries no backslashes at all, and
    returning early for such text skipped the very normalization it needed
    most: instructions are 40% of one real 科研費 form.
    """
    if "\\" in text or re.search(r"(?m)^\s*%", text):
        text = _strip_latex_scaffolding(text)
    return _finish_prose(text)


def _strip_latex_scaffolding(text: str) -> str:
    text = re.sub(r"(?m)^\s*%.*$", " ", text)
    text = re.sub(
        r"\\begin\{(?P<figenv>figure\*?|center)\}.*?"
        r"\\includegraphics.*?\\end\{(?P=figenv)\}",
        " ",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(r"(?<!\\)\\\[.*?(?<!\\)\\\]", " 数式 ", text, flags=re.DOTALL)
    text = re.sub(
        r"\\begin\{(?:equation\*?|align\*?|gather\*?)\}.*?"
        r"\\end\{(?:equation\*?|align\*?|gather\*?)\}",
        " 数式 ",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(r"\$\$.*?\$\$", " 数式 ", text, flags=re.DOTALL)
    text = re.sub(r"\$[^$]*\$", " 数式 ", text)
    # Form fields often use Japanese-named commands whose values are metadata,
    # not proposal prose (for example ``\newcommand{\研究種目名}{該当なし}``).
    # Remove the complete definition before the generic command stripper can
    # leave both arguments behind as ordinary text.
    text = re.sub(
        r"\\(?:re)?newcommand\*?\{[^{}]*\}(?:\[[^\]]*\])?\{[^{}]*\}",
        " ",
        text,
    )
    # A 研究業績リスト is citations, not prose. Left alone, the \item markers are
    # stripped as commands and the whole list merges into one sentence of many
    # hundred characters, so every proposal trips the sentence-length check on
    # the publication list its form requires.
    text = re.sub(
        r"\\begin\{(?P<listenv>enumerate|itemize)\}(?P<body>.*?)\\end\{(?P=listenv)\}",
        _strip_citation_items,
        text,
        flags=re.DOTALL,
    )
    # A heading is its own line. Rendered inline it fused with the paragraph
    # below it, and the field title plus the opening sentence was reported as
    # one 93-character sentence.
    text = re.sub(
        r"\\(?:section|subsection|subsubsection)\*?\{([^{}]*)\}",
        r"\n\1\n",
        text,
    )
    text = re.sub(
        r"\\(?:textbf|textit|emph|underline)\*?\{([^{}]*)\}",
        r" \1 ",
        text,
    )
    text = re.sub(
        r"\\(?:begin|end|input|include|includegraphics|bibliography|bibliographystyle)"
        r"\*?(?:\[[^\]]*\])?\{[^{}]*\}",
        " ",
        text,
    )
    text = re.sub(r"\\[A-Za-z@]+\*?(?:\[[^\]]*\])?", " ", text)
    text = text.replace("\\\\", " ").replace("{", " ").replace("}", " ")
    # Stripping the command but keeping its braces left the arguments behind:
    # \vspace*{0zw} and \rule{\linewidth}{1pt} became the tokens 0zw and 1pt,
    # which then joined the following sentence and were linted as prose.
    text = re.sub(
        r"(?<![0-9A-Za-z])-?[0-9]*\.?[0-9]+\s*"
        r"(?:zw|zh|pt|mm|cm|in|em|ex|bp|dd|sp|truept|truemm)(?![0-9A-Za-z])",
        " ",
        text,
    )
    return text


def _finish_prose(text: str) -> str:
    """Drop the form's own instruction text, then normalise whitespace.

    Every Japanese application form carries several hundred characters of
    instructions, and linting them reports the funder's writing as the
    applicant's defects: 「冒頭にその概要を簡潔にまとめて記述し、本文には、
    (1)…」 was reported as a 逆茂木 sentence in an adopted proposal.
    """
    lines = [re.split(r"(?<=[。．!?！？])", line) for line in text.split("\n")]
    every = [s for line in lines for s in line if s.strip()]
    polite = [s for s in every if _POLITE_ENDING.search(s.strip())]
    # A Japanese proposal body is written in である調 and the form's own
    # instructions in ですます調. Measured across five real documents, polite
    # endings are 0-3% of the text and every one of them belongs to the form:
    # 「…承認手続が必要となる調査・研究・実験などが対象となります。」 was
    # linted as the applicant's 逆茂木 sentence. The ratio guard leaves a
    # proposal that is genuinely written in ですます調 alone.
    drop_polite = bool(every) and len(polite) / len(every) < 0.3

    kept: list[str] = []
    for line in lines:
        sentences = [
            s for s in line
            if s.strip()
            and not _FORM_INSTRUCTION.search(s)
            and not (drop_polite and _POLITE_ENDING.search(s.strip()))
        ]
        if sentences:
            kept.append("".join(sentences))
    text = "\n".join(kept)

    # Collapse runs, but keep newlines: they are the segment boundaries the
    # sentence and co-occurrence checks rely on, and flattening them fused
    # every heading into the paragraph below it.
    text = re.sub(r"[^\S\n]+", " ", text)
    return re.sub(r"\n\s*\n+", "\n", text).strip()


# The shared Japanese helpers also serve paper-writing and intentionally accept
# already-extracted prose. The grant-writing server accepts whole .tex files,
# so every public wrapper must apply the same form/scaffolding filter as the
# integrated health report. Otherwise a direct MCP call and the health report
# disagree about the same proposal.
def grant_writing_lint_bedrock(text: str) -> dict:
    """Lint applicant prose after removing form and LaTeX scaffolding."""
    return _ja_lint_bedrock(_prose_for_lint(_read_text_if_path(text)))


def grant_writing_suggest_redundancy_fixes(text: str) -> dict:
    """Suggest redundancy fixes only in applicant prose."""
    return _ja_suggest_redundancy_fixes(_prose_for_lint(_read_text_if_path(text)))


def grant_writing_check_misuse_japanese(text: str) -> dict:
    """Check Japanese misuse only in applicant prose."""
    return _ja_check_misuse_japanese(_prose_for_lint(_read_text_if_path(text)))


def grant_writing_check_kanji_ratio(
    text: str,
    min_ratio: float = 0.18,
    max_ratio: float = 0.40,
) -> dict:
    """Measure the kanji ratio of applicant prose, not the application form."""
    return _ja_check_kanji_ratio(
        _prose_for_lint(_read_text_if_path(text)),
        min_ratio=min_ratio,
        max_ratio=max_ratio,
    )


def grant_writing_check_subject_predicate_distance(
    text: str,
    max_chars: int = 40,
) -> dict:
    """Check subject-predicate distance only in applicant prose."""
    # Word drafts in the private corpus use both Japanese commas, ``、`` and
    # ``，``. The shared checker keys on 「は、」「が、」; without normalising
    # punctuation, an otherwise identical Word draft has no analyzable topic.
    prose = _prose_for_lint(_read_text_if_path(text)).replace("，", "、")
    return _ja_subject_predicate_distance(
        prose,
        max_chars=max_chars,
    )


def grant_writing_find_undefined_acronyms(
    text: str,
    whitelist: str = "",
    min_len: int = 2,
    max_len: int = 6,
    context_window: int = 200,
) -> dict:
    """Find undefined acronyms in applicant prose only."""
    return _ja_find_undefined_acronyms(
        _prose_for_lint(_read_text_if_path(text)),
        whitelist=whitelist,
        min_len=min_len,
        max_len=max_len,
        context_window=context_window,
    )


def grant_writing_acronym_usage_audit(
    text: str,
    whitelist: str = "",
    min_uses_for_abbrev: int = 3,
    min_len: int = 2,
    max_len: int = 6,
) -> dict:
    """Audit acronym use in applicant prose only."""
    return _ja_acronym_usage_audit(
        _prose_for_lint(_read_text_if_path(text)),
        whitelist=whitelist,
        min_uses_for_abbrev=min_uses_for_abbrev,
        min_len=min_len,
        max_len=max_len,
    )


def grant_writing_check_notation_variants(text: str) -> dict:
    """Check notation variants in applicant prose only."""
    return _ja_check_notation_variants(_prose_for_lint(_read_text_if_path(text)))


def _contains_any(text_lower: str, keywords: list[str]) -> list[str]:
    return [kw for kw in keywords if kw.lower() in text_lower]


def _severity_from_score(score: float | None) -> str:
    if score is None:
        return "UNKNOWN"
    if score < 4:
        return "CRITICAL"
    if score < 6:
        return "HIGH"
    if score < 8:
        return "MEDIUM"
    return "LOW"


def grant_writing_usage() -> str:
    """Return the grant-writing guide."""
    return _load_skill()


def grant_writing_kaken_review_axes() -> dict:
    """Return the current official review axes for KAKENHI B/C (General).

    This is a source-grounded reference map, not a prediction or score. It
    keeps the three research-plan elements, the separately rated
    internationality element, and the additional budget-validity assessment
    distinct because they affect a proposal in different ways.
    """
    review_regulations = (
        "https://www.jsps.go.jp/file/storage/kaken_0103_shinsakitei_g_4984/"
        "hyoukakitei260622.pdf"
    )
    web_input_guide = (
        "https://www.jsps.go.jp/file/storage/kaken_kiban_2026_g_4978/"
        "web_yoryo_kiban.pdf"
    )
    review_page = (
        "https://www.jsps.go.jp/j-grantsinaid/01_seido/03_shinsa/index.html"
    )
    reviewer_pamphlet = (
        "https://www.jsps.go.jp/file/storage/kaken_pamph_j2026/"
        "kakenhi2026.pdf"
    )
    return {
        "scheme": "科研費 基盤研究(B・C)(一般)",
        "verified_on": "2026-08-21",
        "source_revision": "審査規程: 2026-06-22改正",
        "purpose": (
            "公式基準を申請書のレビュー観点へ写す参照表。"
            "採否予測やキーワード採点には用いない。"
        ),
        "review_process": {
            "method": "同一審査委員による2段階書面審査",
            "reviewers_for_scientific_research_c": 3,
            "first_stage": (
                "研究計画の3要素と国際性を個別に絶対評価し、"
                "総合評点は審査区分内の相対評価で付す。"
            ),
            "second_stage": (
                "ボーダーゾーン等を、他の審査委員の1段階目意見も"
                "参照して再評価する。"
            ),
        },
        "research_plan_axes": [
            {
                "id": "academic_importance",
                "label": "研究課題の学術的重要性",
                "review_questions": [
                    "学術的に推進すべき重要な課題か。",
                    "核心となる学術的問いが明確で、独自性・創造性があるか。",
                    "着想、国内外の研究動向、研究上の位置づけが明確か。",
                    "より広い学術、科学技術、社会への波及が期待できるか。",
                ],
            },
            {
                "id": "method_validity",
                "label": "研究方法の妥当性",
                "review_questions": [
                    "目的に対する研究方法が具体的かつ適切か。",
                    "研究経費が研究計画と整合しているか。",
                    "目的達成に必要な準備が整っているか。",
                ],
            },
            {
                "id": "capability_environment",
                "label": "研究遂行能力及び研究環境の適切性",
                "review_questions": [
                    "これまでの研究活動から十分な遂行能力を確認できるか。",
                    "必要な施設、設備、資料等の研究環境が整っているか。",
                ],
            },
        ],
        "separate_scored_axis": {
            "id": "internationality",
            "label": "研究課題の国際性",
            "review_question": (
                "世界の研究を将来けん引する、協同により世界の研究へ貢献する、"
                "又は日本独自の研究として高い価値を創出することが期待できるか。"
            ),
            "allocation_note": (
                "国際性の高い課題は、若手研究者等への助成調整や、"
                "応募額を尊重した配分の対象になり得る。"
            ),
        },
        "budget_validity": {
            "role": (
                "研究方法では計画との整合性を評価し、これとは別に配分額の"
                "判断材料として経費の妥当性・必要性を確認する。"
                "学術的重要性と同格の独立した総合評点軸ではない。"
            ),
            "review_questions": [
                "経費内容が妥当で、有効利用が見込まれるか。",
                "設備備品が研究計画の遂行に真に必要か。",
                "設備、旅費、人件費・謝金のいずれかが90%を超える場合も有効利用できるか。",
            ],
            "consequence": (
                "基盤研究(B・C)では、複数の審査委員が経費に問題ありとした場合、"
                "平均より低い充足率となる。"
            ),
            "persuasiveness_test": [
                "各費目がどの研究行為、年度、成果物に対応するか。",
                "単価、数量、月数・回数と根拠資料から再計算できるか。",
                "最大費目が研究の中心作業と一致しているか。",
                "減額時も核心の検証ループを維持できる優先順位があるか。",
            ],
        },
        "budget_entry_requirements": [
            "機械器具を単なる『一式』とせず、内訳を示す。",
            "設備備品費・消耗品費は必要性と積算根拠を示す。",
            "旅費は出張の事項・目的ごとに示す。",
            "人件費・謝金は用途を分け、判明していれば身分、人数、月数も示す。",
            "その他経費も事項ごとに分け、必要性と積算根拠を示す。",
            "年度内で特定費目が90%を超える場合や大きな割合を占める場合は、研究上の必要性を明記する。",
            "研究代表者・研究分担者本人の人件費・謝金は直接経費に計上しない。",
        ],
        "interpretation": [
            "申請書を『4軸で均等採点』する制度ではない。研究計画3要素を中心に総合評価される。",
            "国際性は別に評定されるため、共同研究だけでなく世界への価値の出方を明示する。",
            "予算の説得力は金額の大小ではなく、研究計画との対応、積算可能性、必要性で作る。",
            "公式基準に論文数の足切りはない。業績は各担当を遂行できる証拠として結びつける。",
        ],
        "sources": [
            {
                "title": "科研費 審査及び評価に関する規程（2026-06-22改正）",
                "url": review_regulations,
                "supports": "審査方式、評定要素、国際性、研究経費の妥当性",
            },
            {
                "title": "令和9年度 基盤研究等 Web入力要領",
                "url": web_input_guide,
                "supports": "経費明細、必要性、積算根拠、90%超の説明、対象外経費",
            },
            {
                "title": "日本学術振興会 審査・評価について",
                "url": review_page,
                "supports": "現行の審査規程・評定基準への公式入口",
            },
            {
                "title": "科研費パンフレット2026",
                "url": reviewer_pamphlet,
                "supports": "基盤研究(C)の1課題当たり審査委員数",
            },
        ],
    }


def grant_writing_analyze_sentences(text: str, max_len: int = 90) -> dict:
    """Analyze Japanese sentence length for grant proposals.

    Grant drafts can tolerate denser prose than slides, but application
    reviewers still need a clear one-claim-per-sentence rhythm.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    # A newline ends a segment too. A 年度計画 matrix is a run of short cells
    # with no full stop, and joining them reported a Gantt table in an adopted
    # proposal as a single 455-character sentence.
    sentences = [s.strip() for s in re.split(r"[。．!?！？\n]", text) if s.strip()]
    if not sentences:
        return {"error": "no sentences found"}
    lengths = [len(s) for s in sentences]
    long_ones = [
        {"index": i, "length": n, "head": s[:50] + ("..." if len(s) > 50 else "")}
        for i, (s, n) in enumerate(zip(sentences, lengths))
        if n > max_len
    ]
    return {
        "total_sentences": len(sentences),
        "avg_length": round(sum(lengths) / len(lengths), 1),
        "max_length": max(lengths),
        "threshold": max_len,
        "over_threshold_count": len(long_ones),
        "over_threshold_examples": long_ones[:8],
        "target": "avg <= 70-80 chars; one reviewer-relevant claim per sentence",
        "warning": "sentences above 90 chars usually hide two claims",
    }


_ADJACENT_REVIEWER_INFRA_TERMS = (
    "OSS", "API", "MCP", "GitHub", "CI", "AI", "リポジトリ", "サーバ",
    "インターフェース", "インタフェース", "モジュール", "版管理",
)
_ADJACENT_REVIEWER_SCIENCE_TERMS = (
    "解析", "設計", "最適化", "離散化", "求解", "測定", "実験", "計算",
    "モデル", "物理量", "磁場", "損失", "効率", "軌道", "インピーダンス",
)
_ADJACENT_REVIEWER_DECISION_TERMS = (
    "条件", "判定", "順位", "再現", "反証", "検証", "許容差", "適用限界",
    "採否", "凍結", "高忠実度", "確定",
)
_ADJACENT_REVIEWER_AUDIT_TERMS = (
    "再現", "再構成", "反証", "判定", "採否", "凍結", "版管理", "許容差",
    "適用限界", "変更履歴",
)
_ADJACENT_REVIEWER_EVIDENCE_PATTERNS = (
    re.compile(r"(?:実装|導入|再実行|完了|確認|収束|採択|発表|共著|取り込)"),
    re.compile(r"\d[\d,.]*(?:\\?%|~?k?Hz|自由度|反復|件|ケース|回)"),
    re.compile(r"(?:共同研究|予備実証|回帰試験|解析解|有限差分一致)"),
)
_ADJACENT_REVIEWER_TAKEAWAY_PATTERN = re.compile(
    r"(?:これら(?:の[^。！？\n]{0,24})?|この(?:実績|成果|準備|予備実証)|以上)"
    r"(?:により|から)|"
    r"研究項目[0-9０-９]+[^。！？\n]{0,64}(?:開始|着手|実行)できる"
)


def grant_writing_adjacent_reviewer_readability_check(text: str) -> dict:
    """Find prose that is short but cognitively dense for an adjacent reviewer.

    Sentence length and phrase-level redundancy do not explain every reading
    failure. A sentence can be under 60 characters and still ask the reader to
    unpack technical nouns, an infrastructure layer, a scientific operation,
    and a decision rule at once. This diagnostic intentionally has no score.
    """
    prose = _prose_for_lint(_read_text_if_path(text))
    sentences = [
        segment.strip()
        for segment in re.split(r"[。．!?！？\n]", prose)
        if segment.strip()
    ]
    paragraphs = [
        line.strip()
        for line in prose.splitlines()
        if len(line.strip()) >= 30
        and not re.fullmatch(r"[（(【\[].{0,45}[）)】\]]", line.strip())
    ]
    if not sentences:
        return {"applicable": False, "risk_count": 0, "risks": []}

    jp_chars = re.findall(r"[一-龥々〆ヵヶぁ-んァ-ヶー]", prose)
    kanji_chars = re.findall(r"[一-龥々〆ヵヶ]", prose)
    kanji_ratio = len(kanji_chars) / len(jp_chars) if jp_chars else 0.0
    risks: list[dict] = []

    def add_risk(
        risk_type: str,
        excerpt: str,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
        **details,
    ) -> None:
        item = {
            "type": risk_type,
            "severity": severity,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:360],
            "comment": comment,
            "recommendation": recommendation,
        }
        item.update(details)
        risks.append(item)

    dense_sentences: list[dict] = []
    notation_piles: list[dict] = []
    representation_mismatches: list[dict] = []
    ambiguous_relation_phrases: list[dict] = []
    scope_without_deliverables: list[dict] = []
    takeaways_after_evidence: list[dict] = []
    representation_pattern = re.compile(
        r"(?P<answer>設計則|選択則|指針|適用条件|成立条件|知見|成果)"
        r"を[、，,\s]*(?P<representation>[^。！？\n]{0,36}"
        r"(?:区間|指標|コード|実装|リポジトリ))として"
        r"(?:与える|示す|提示する)"
    )
    ambiguous_patterns = (
        re.compile(r"(?:一方|双方)へ(?:統一|還流|移行|集約)"),
        re.compile(r"異なるコード系譜(?:間|を|の)"),
        re.compile(r"(?:同一|共通)の?設計量で採否(?:する|を判断する)"),
    )
    scope_markers = ("必達範囲", "必須範囲", "最低限の達成範囲")
    deliverable_pattern = re.compile(
        r"(?:実証|検証|解明|確立|同定|導出|提示|実装|構築|開発|評価|"
        r"明らかに|条件を示|法則を示|知見を得)"
    )
    for index, sentence in enumerate(sentences):
        local_jp = re.findall(r"[一-龥々〆ヵヶぁ-んァ-ヶー]", sentence)
        local_kanji = re.findall(r"[一-龥々〆ヵヶ]", sentence)
        local_ratio = len(local_kanji) / len(local_jp) if local_jp else 0.0
        comma_count = len(re.findall(r"[、，,]", sentence))
        latin_terms = list(dict.fromkeys(re.findall(
            r"(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9_.()+/-]{1,}",
            sentence,
        )))
        if len(sentence) >= 36 and local_ratio >= 0.60 and comma_count >= 3:
            dense_sentences.append({
                "index": index,
                "length": len(sentence),
                "kanji_ratio": round(local_ratio, 3),
                "comma_count": comma_count,
                "excerpt": sentence[:240],
            })
        has_prose_predicate = bool(re.search(
            r"(?:する|した|用いる|使う|示す|求める|比べる|比較|評価|解析|"
            r"設計|接続|判定|構築|開発|明らか|である|となる)",
            sentence,
        ))
        if len(sentence) >= 34 and len(latin_terms) >= 3 and has_prose_predicate:
            notation_piles.append({
                "index": index,
                "terms": latin_terms[:8],
                "excerpt": sentence[:240],
            })
        for match in representation_pattern.finditer(sentence):
            representation_mismatches.append({
                "index": index,
                "answer": match.group("answer"),
                "representation": match.group("representation"),
                "excerpt": sentence[:240],
            })
        vague_hits = [
            match.group(0)
            for pattern in ambiguous_patterns
            for match in pattern.finditer(sentence)
        ]
        if vague_hits:
            ambiguous_relation_phrases.append({
                "index": index,
                "phrases": list(dict.fromkeys(vague_hits)),
                "excerpt": sentence[:240],
            })
        if (
            any(marker in sentence for marker in scope_markers)
            and not deliverable_pattern.search(sentence)
        ):
            scope_without_deliverables.append({
                "index": index,
                "excerpt": sentence[:240],
            })

    if dense_sentences:
        first = dense_sentences[0]
        add_risk(
            "compressed_concept_density",
            first["excerpt"],
            (
                "短い文でも漢語と列挙が密集し、隣接分野の審査者は一文の中で"
                "複数の概念を展開する必要がある。"
            ),
            (
                "最初の文は対象・困りごと・得る答えの一つに絞る。手法名、"
                "条件、評価量は次の文へ一段ずつ展開する。"
            ),
            examples=dense_sentences[:6],
        )

    if notation_piles:
        first = notation_piles[0]
        add_risk(
            "notation_or_method_pile",
            first["excerpt"],
            "一文に三つ以上の英字略語・手法名があり、関係より名称が先に見える。",
            (
                "一般名と役割を先に述べ、固有名・略語は一文に一つを目安に"
                "導入する。複数手法の対応は図表又は別文へ分ける。"
            ),
            examples=notation_piles[:6],
        )

    if representation_mismatches:
        first = representation_mismatches[0]
        add_risk(
            "result_representation_type_mismatch",
            first["excerpt"],
            (
                "研究上の答えと、その数理・実装上の表現を『として』で直結し、"
                "何を明らかにする研究かが読みにくい。"
            ),
            (
                "まず条件・設計則・知見を平易に述べ、区間・指標・コード等は"
                "それをどのように表現又は検証するかとして次の文へ分ける。"
            ),
            examples=representation_mismatches[:5],
        )

    if ambiguous_relation_phrases:
        first = ambiguous_relation_phrases[0]
        add_risk(
            "ambiguous_relation_or_decision_object",
            first["excerpt"],
            (
                "統一・還流・採否の対象を内部略語や指示語へ預けており、"
                "隣接分野の審査者が対象を確定しにくい。"
            ),
            (
                "開発母体、内部形式、解析手法、設計候補、適用可否など、"
                "何を統一せず、何を判断するのかを名詞で明示する。"
            ),
            examples=ambiguous_relation_phrases[:5],
        )

    if scope_without_deliverables:
        first = scope_without_deliverables[0]
        add_risk(
            "required_scope_without_deliverable",
            first["excerpt"],
            "必達範囲が対象課題の列挙だけで、研究期間内に何を達成するかがない。",
            (
                "『二課題で結合条件を実証する』のように、対象とともに"
                "検証・解明・確立する成果を動詞で示す。"
            ),
            severity="HIGH",
            examples=scope_without_deliverables[:5],
        )

    layer_paragraphs: list[dict] = []
    for index, paragraph in enumerate(paragraphs):
        infra = [term for term in _ADJACENT_REVIEWER_INFRA_TERMS if term in paragraph]
        science = [
            term for term in _ADJACENT_REVIEWER_SCIENCE_TERMS if term in paragraph
        ]
        decision = [
            term for term in _ADJACENT_REVIEWER_DECISION_TERMS if term in paragraph
        ]
        if len(infra) >= 2 and len(science) >= 2 and len(decision) >= 2:
            layer_paragraphs.append({
                "index": index,
                "infrastructure": infra[:6],
                "science": science[:6],
                "decision": decision[:6],
                "excerpt": paragraph[:300],
            })

        takeaway_match = _ADJACENT_REVIEWER_TAKEAWAY_PATTERN.search(paragraph)
        if takeaway_match:
            evidence_prefix = paragraph[:takeaway_match.start()]
            evidence_signals = sum(
                1
                for pattern in _ADJACENT_REVIEWER_EVIDENCE_PATTERNS
                if pattern.search(evidence_prefix)
            )
            latin_terms = list(dict.fromkeys(re.findall(
                r"(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9_.()+/-]{1,}",
                evidence_prefix,
            )))
            if len(latin_terms) >= 2:
                evidence_signals += 1
            relative_position = takeaway_match.start() / max(len(paragraph), 1)
            if (
                evidence_signals >= 2
                and takeaway_match.start() >= 48
                and relative_position >= 0.28
            ):
                takeaways_after_evidence.append({
                    "index": index,
                    "takeaway": takeaway_match.group(0),
                    "takeaway_relative_position": round(relative_position, 2),
                    "evidence_signals_before_takeaway": evidence_signals,
                    "terms_before_takeaway": latin_terms[:8],
                    "excerpt": paragraph[:360],
                })
    if layer_paragraphs:
        first = layer_paragraphs[0]
        add_risk(
            "three_layer_paragraph",
            first["excerpt"],
            (
                "同じ段落で研究基盤、科学的操作、設計判断を同時に説明している。"
                "審査者は何が主張で何が手段かを保持し続けなければならない。"
            ),
            (
                "段落を『対象分野の問題』『何を比較・測定するか』『再現基盤が"
                "どう支えるか』の順に分け、基盤名は最後へ下げる。"
            ),
            severity="HIGH",
            examples=layer_paragraphs[:5],
        )

    if takeaways_after_evidence:
        first = takeaways_after_evidence[0]
        add_risk(
            "takeaway_after_evidence",
            first["excerpt"],
            (
                "実績が研究計画に何を可能にするかが段落後半まで現れず、"
                "審査者は固有名・手法名・数値を保持してから結論を逆算する必要がある。"
            ),
            (
                "段落冒頭で研究項目の開始点又は得られる判断を述べる。次に一般語で"
                "技術的な役割を説明し、手法名・数値・発表先はその根拠として後へ置く。"
            ),
            severity="HIGH",
            examples=takeaways_after_evidence[:5],
            rewrite_order=[
                "reviewer_takeaway",
                "plain_language_role",
                "specific_method_or_evidence",
                "limit_or_remaining_question",
            ],
        )

    audit_hits = [
        {"term": term, "count": prose.count(term)}
        for term in _ADJACENT_REVIEWER_AUDIT_TERMS
        if prose.count(term)
    ]
    audit_count = sum(item["count"] for item in audit_hits)
    audit_density = 1000.0 * audit_count / max(len(prose), 1)
    audit_paragraph_count = sum(
        1
        for paragraph in paragraphs
        if any(term in paragraph for term in _ADJACENT_REVIEWER_AUDIT_TERMS)
    )
    if audit_density >= 5.0 and audit_paragraph_count >= 4:
        add_risk(
            "distributed_assurance_repetition",
            "、".join(f"{item['term']} {item['count']}回" for item in audit_hits),
            (
                "再現・反証・判定・版管理などの保証語が多くの段落へ分散し、"
                "科学的な発見より監査手順が前面に出ている。"
            ),
            (
                "各研究項目では得る知見を先に書き、共通の再現・反証手順は"
                "検証設計の一段落へ集約する。"
            ),
            audit_density_per_1000=round(audit_density, 2),
            paragraph_count=audit_paragraph_count,
        )

    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(dict.fromkeys(
        risk["recommendation"] for risk in risks
    ))
    return {
        "applicable": True,
        "score": None,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "metrics": {
            "sentence_count": len(sentences),
            "average_sentence_length": round(
                sum(map(len, sentences)) / len(sentences), 1
            ),
            "kanji_ratio": round(kanji_ratio, 3),
            "compressed_dense_sentence_count": len(dense_sentences),
            "notation_or_method_pile_count": len(notation_piles),
            "result_representation_mismatch_count": len(
                representation_mismatches
            ),
            "ambiguous_relation_phrase_count": len(
                ambiguous_relation_phrases
            ),
            "required_scope_without_deliverable_count": len(
                scope_without_deliverables
            ),
            "three_layer_paragraph_count": len(layer_paragraphs),
            "takeaway_after_evidence_count": len(takeaways_after_evidence),
            "assurance_term_density_per_1000_chars": round(audit_density, 2),
        },
        "revision_protocol": {
            "sequence": [
                "reviewer_takeaway",
                "plain_language_role",
                "specific_method_or_evidence",
                "limit_or_remaining_question",
            ],
            "instruction": (
                "Preserve necessary technical terms, but make the reviewer-facing "
                "claim readable before asking the reader to unpack method names, "
                "numbers, publications, or infrastructure."
            ),
        },
        "diagnosis": (
            "Sentence length alone is insufficient. Short sentences can remain hard "
            "when abstract nouns, method names, infrastructure, and decision rules "
            "are compressed into the same reading unit."
        ),
        "target_reader": (
            "a reviewer who knows the broad field but not the applicant's software, "
            "laboratory shorthand, or exact numerical formulation"
        ),
        "source": (
            "generic adjacent-domain reviewer readability diagnostic; non-scoring"
        ),
    }


_GRANT_DOCUMENT_TYPES = {
    "grant",
    "grant_application",
    "grant_proposal",
    "funding_application",
    "japanese_grant_proposal",
    "助成金申請",
    "科研費申請",
}
_MANUSCRIPT_DOCUMENT_TYPES = {
    "conference_manuscript",
    "paper",
    "research_manuscript",
    "research_meeting",
    "research_meeting_manuscript",
    "technical_paper",
    "研究会原稿",
    "論文",
}


def grant_writing_japanese_genre_contract(document_type: str) -> dict:
    """Route Japanese prose to grant or research-manuscript review criteria.

    Grant applications and research-meeting manuscripts share low-level
    Japanese mechanics, but they do not share a scoring objective. This
    contract prevents readable completed-work prose from being mistaken for a
    persuasive and feasible funding proposal.
    """
    requested = str(document_type).strip()
    normalized = requested.casefold().replace("-", "_").replace(" ", "_")
    shared_foundation = [
        "clear_sentence_boundaries",
        "short_modifier_scope",
        "subject_predicate_proximity",
        "notation_and_term_consistency",
    ]
    grant_criteria = [
        "reviewer_visible_problem_and_why_now",
        "academic_question_and_proposed_move",
        "feasibility_team_and_preliminary_evidence",
        "verifiable_deliverables_schedule_and_budget",
        "committed_future_wording",
    ]
    manuscript_criteria = [
        "definitions_assumptions_and_symbol_introduction",
        "equation_and_method_reproducibility",
        "result_figure_table_and_citation_traceability",
        "evidence_bounded_claims_and_limitations",
        "completed_work_reported_in_appropriate_tense",
    ]
    base = {
        "requested_document_type": requested,
        "shared_foundation": shared_foundation,
        "grant_proposal_criteria": grant_criteria,
        "research_manuscript_criteria": manuscript_criteria,
        "policy": (
            "Share only foundational Japanese lint. Never average or reuse the "
            "genre-specific score across grant proposals and research manuscripts."
        ),
    }
    if normalized in _GRANT_DOCUMENT_TYPES:
        return {
            **base,
            "applicable": True,
            "status": "supported",
            "canonical_document_type": "grant_proposal",
            "review_owner": "grant-writing",
            "review_goal": (
                "help a time-limited funding reviewer decide why the work matters "
                "now and whether the proposed team can deliver it"
            ),
        }
    if normalized in _MANUSCRIPT_DOCUMENT_TYPES:
        return {
            **base,
            "applicable": False,
            "status": "wrong_genre",
            "canonical_document_type": "research_meeting_manuscript",
            "expected_document_type": "grant_proposal",
            "review_owner": "paper-writing",
            "review_goal": (
                "make definitions, methods, equations, evidence, figures, and "
                "citations traceable as a completed scientific argument"
            ),
            "route_to": {
                "server": "mcp-server-paper-writing",
                "tools": [
                    "paper_writing_bilingual_readability_check",
                    "paper_writing_em_submission_gate",
                ],
            },
            "reason": (
                "Research-meeting manuscripts and grant proposals require "
                "different Japanese review criteria."
            ),
        }
    return {
        **base,
        "applicable": False,
        "status": "unsupported_document_type",
        "canonical_document_type": None,
        "expected_document_type": "grant_proposal",
        "review_owner": None,
        "reason": (
            "Declare a grant-proposal or research-manuscript document type "
            "before applying a genre-specific Japanese score."
        ),
    }


def grant_writing_japanese_readability_score(
    text: str,
    document_type: str,
) -> dict:
    """Score Japanese grant prose with Japanese-specific writing criteria.

    This is a 100-point readability diagnostic, not a funding prediction or a
    scientific-merit score. It combines six observable axes: one-claim
    sentence rhythm, Japanese logical order, subject-predicate proximity,
    lexical load, notation consistency, and committed proposal wording.
    English prose is deliberately not scored or averaged into this result.
    ``document_type`` is required so research-meeting manuscripts cannot be
    silently evaluated against grant-proposal criteria.
    """
    genre = grant_writing_japanese_genre_contract(document_type)
    if not genre["applicable"]:
        return {
            "applicable": False,
            "status": genre["status"],
            "score": None,
            "score_max": 100,
            "document_type": genre.get("canonical_document_type"),
            "reason": genre["reason"],
            "genre_contract": genre,
            "scoring_policy": (
                "No score was produced because Japanese readability scores "
                "are genre-specific."
            ),
        }

    prose = _prose_for_lint(_read_text_if_path(text))
    japanese_chars = re.findall(r"[\u3040-\u30ff\u3400-\u9fff]", prose)
    if len(japanese_chars) < 20:
        return {
            "applicable": False,
            "status": "not_applicable",
            "score": None,
            "score_max": 100,
            "document_type": "grant_proposal",
            "reason": "at least 20 Japanese characters are required",
            "genre_contract": genre,
            "scoring_policy": (
                "Japanese grant prose only. English is neither scored nor "
                "averaged into this diagnostic."
            ),
        }

    sentence = grant_writing_analyze_sentences(prose)
    kanji = grant_writing_check_kanji_ratio(
        prose,
        min_ratio=0.30,
        max_ratio=0.60,
    )
    subject = grant_writing_check_subject_predicate_distance(prose)
    bedrock = grant_writing_lint_bedrock(prose)
    adjacent = grant_writing_adjacent_reviewer_readability_check(prose)
    misuse = grant_writing_check_misuse_japanese(prose)
    notation = grant_writing_check_notation_variants(prose)
    weak = grant_writing_count_weak_expressions(prose)

    overlong = sentence.get("over_threshold_count", 0)
    max_length = sentence.get("max_length", 0)
    average_length = sentence.get("avg_length", 0.0)
    sentence_score = 25 - min(18, 6 * overlong)
    if max_length > 140:
        sentence_score -= 10
    elif max_length > 110:
        sentence_score -= 5
    if average_length > 80:
        sentence_score -= 3
    sentence_score = max(0, sentence_score)

    bedrock_count = bedrock.get("issue_count", 0)
    logical_order_score = max(0, 20 - min(20, 6 * bedrock_count))

    subject_violations = subject.get("violation_count", 0)
    subject_score = max(0, 15 - min(15, 5 * subject_violations))

    adjacent_high = sum(
        risk.get("severity") == "HIGH" for risk in adjacent.get("risks", [])
    )
    adjacent_medium = sum(
        risk.get("severity") == "MEDIUM" for risk in adjacent.get("risks", [])
    )
    lexical_penalty = 8 * adjacent_high + 5 * adjacent_medium
    kanji_ratio = kanji.get("kanji_ratio", 0.0)
    if kanji_ratio < 0.20 or kanji_ratio > 0.70:
        lexical_penalty += 8
    elif kanji_ratio < 0.30 or kanji_ratio > 0.60:
        lexical_penalty += 4
    lexical_score = max(0, 20 - min(20, lexical_penalty))

    misuse_count = misuse.get("total_matches", 0)
    notation_count = notation.get("total_findings", 0)
    consistency_score = max(
        0,
        10 - min(10, 3 * misuse_count + 2 * notation_count),
    )

    weak_count = weak.get("total_weak_expressions", 0)
    commitment_score = max(0, 10 - min(10, 2 * weak_count))

    axes = {
        "one_claim_sentence_rhythm": {
            "score": sentence_score,
            "score_max": 25,
            "evidence": sentence,
        },
        "japanese_logical_order": {
            "score": logical_order_score,
            "score_max": 20,
            "evidence": bedrock,
        },
        "subject_predicate_proximity": {
            "score": subject_score,
            "score_max": 15,
            "evidence": subject,
        },
        "lexical_and_concept_load": {
            "score": lexical_score,
            "score_max": 20,
            "evidence": {
                "kanji": kanji,
                "adjacent_reviewer": adjacent,
            },
        },
        "notation_and_usage_consistency": {
            "score": consistency_score,
            "score_max": 10,
            "evidence": {
                "misuse": misuse,
                "notation": notation,
            },
        },
        "committed_proposal_wording": {
            "score": commitment_score,
            "score_max": 10,
            "evidence": weak,
        },
    }
    score = sum(axis["score"] for axis in axes.values())
    status = "pass" if score >= 85 else "warning" if score >= 70 else "fail"
    priorities = [
        {
            "axis": name,
            "lost_points": axis["score_max"] - axis["score"],
            "score": axis["score"],
            "score_max": axis["score_max"],
        }
        for name, axis in axes.items()
        if axis["score"] < axis["score_max"]
    ]
    priorities.sort(key=lambda item: (-item["lost_points"], item["axis"]))
    return {
        "applicable": True,
        "status": status,
        "score": score,
        "score_max": 100,
        "document_type": "grant_proposal",
        "japanese_character_count": len(japanese_chars),
        "genre_contract": genre,
        "scoring_axes": axes,
        "revision_priorities": priorities,
        "thresholds": {
            "pass": "85-100",
            "warning": "70-84",
            "fail": "0-69",
        },
        "scoring_policy": (
            "Japanese grant prose only. English is neither scored nor "
            "averaged. Kanji ratio is a lightly weighted technical-prose "
            "signal; sentence structure and reviewer load carry more weight."
        ),
        "interpretation": (
            "This score estimates Japanese reading load and writing mechanics. "
            "It does not assess novelty, feasibility, scientific merit, or the "
            "probability of funding."
        ),
        "sources": [
            "Kinoshita: Japanese technical-writing principles",
            "Honda: modifier order, punctuation, and kanji balance",
            "Kitahara: Japanese usage diagnostics",
            "CAE-AI Lab adjacent-domain reviewer readability evidence",
        ],
    }


_REVIEWER_MOMENTUM_STAKE_PATTERN = re.compile(
    r"(?:患者|治療|電力消費|省エネルギー|環境|安全|市場|コスト|期間|"
    r"熟練|設計(?:判断|期間|候補|者)|候補選択|比較対象|解析経路|"
    r"研究選択|性能|品質|効率|損失|"
    r"普及|運用|負担|採用|社会|産業|製造)[^。！？\n]{0,44}"
    r"(?:求め|望まれ|必要|左右|削減|短縮|向上|改善|支える|妨げ|課題|"
    r"価値|影響|依存|困難|偏る)|"
    r"(?:求め|望まれ|必要|左右|削減|短縮|向上|改善|妨げ|課題|依存|"
    r"困難|偏る)[^。！？\n]{0,44}(?:患者|治療|電力消費|省エネルギー|環境|"
    r"安全|市場|コスト|期間|熟練|設計|性能|品質|効率|損失|普及|運用)"
)
_REVIEWER_MOMENTUM_BOTTLENECK_PATTERN = re.compile(
    r"(?:しかし|一方|ところが|にもかかわらず|課題|障壁|難所|困難|"
    r"未解決|未確立|未実現|実現されていない|できない|難しい|依存|"
    r"試行錯誤|要する|負担|限界|妨げ|左右|偏る|再実装|実装し直)"
)
_REVIEWER_MOMENTUM_SPECIFICITY_PATTERN = re.compile(
    r"(?:\d[\d,.]*(?:\\?%|年|月|日|時間|件|回|倍|円|人|自由度)?|"
    r"半年|手作業|試行錯誤|既存コード|接続|再実装|初期化運転|人件費|"
    r"誤差|未実現|未確立|未検証|設計時間|作業負荷|適用例|候補選択)"
)
_REVIEWER_MOMENTUM_OPPORTUNITY_PATTERN = re.compile(
    r"(?:近年|急速|進展|蓄積|新た|自由度|利用可能|可能にな|実用化|"
    r"使える|到来|発展|拡大|解放)"
)
_REVIEWER_MOMENTUM_MOVE_PATTERN = re.compile(
    r"(?:そこで)?本(?:研究|提案|事業|研究開発)(?:では|は|の目的)|"
    r"(?:中心|学術的)の問い|目的は|"
    r"そこで[^。！？\n]{0,100}(?:定義|構築|提案|確立|検証|明らか|"
    r"開発|導出|比較|定量化)"
)
_REVIEWER_MOMENTUM_PAYOFF_PATTERN = re.compile(
    r"(?:明らかに|確立|実現|算出|短縮|削減|向上|改善|可能に|開拓|"
    r"普及|解放|判断|選択|確定|再現|反証|変える|取り戻す|断つ|"
    r"定量化|体系化|示す|条件を求め|適用条件)"
)
_REVIEWER_MOMENTUM_METHOD_PATTERN = re.compile(
    r"(?:FEM|RNA|CLN|MCP|OSS|API|AIエージェント|有限要素|"
    r"数値解析|解析手法|モデル|アルゴリズム|ソルバー|メッシュ|"
    r"等価回路|インターフェース|モジュール)"
)
_REVIEWER_MOMENTUM_HYPE_PATTERN = re.compile(
    r"(?:世界初|世界で初めて|全く新しい|画期的|革新的|圧倒的|劇的|"
    r"パラダイムシフト|夢の|究極の)"
)
_REVIEWER_MOMENTUM_EVIDENCE_PATTERN = re.compile(
    r"(?:\d[\d,.]*(?:\\?%|年|月|日|時間|件|回|倍|円|人)?|文献|"
    r"報告|実測|実証|実装|確認|採択|査読|共同研究|市場|既に|済み)"
)


def grant_writing_reviewer_momentum_check(text: str) -> dict:
    """Check whether the opening makes a reviewer want to keep reading.

    Readability and interest are different. A readable opening can still feel
    like an inventory, while an exciting opening can be empty hype. This
    non-scoring diagnostic looks for the reusable arc found in strong grant
    prose: a reviewer-visible stake, a concrete bottleneck, the research move,
    and the observable change that move would unlock. It also keeps method
    names subordinate to that arc.
    """
    prose = _prose_for_lint(_read_text_if_path(text))
    all_sentences = [
        segment.strip()
        for segment in re.split(r"[。．!?！？\n]", prose)
        if segment.strip()
    ]
    if len(all_sentences) < 3:
        return {
            "applicable": False,
            "score": None,
            "risk_count": 0,
            "risks": [],
            "reason": "at least three opening sentences are needed",
        }

    lead_sentences: list[str] = []
    lead_chars = 0
    for sentence in all_sentences[:24]:
        if lead_sentences and lead_chars + len(sentence) > 2400:
            break
        lead_sentences.append(sentence)
        lead_chars += len(sentence)

    def first_index(pattern: re.Pattern[str]) -> int | None:
        return next(
            (index for index, sentence in enumerate(lead_sentences)
             if pattern.search(sentence)),
            None,
        )

    stake_index = first_index(_REVIEWER_MOMENTUM_STAKE_PATTERN)
    bottleneck_index = first_index(_REVIEWER_MOMENTUM_BOTTLENECK_PATTERN)
    opportunity_index = first_index(_REVIEWER_MOMENTUM_OPPORTUNITY_PATTERN)
    move_index = first_index(_REVIEWER_MOMENTUM_MOVE_PATTERN)
    method_index = first_index(_REVIEWER_MOMENTUM_METHOD_PATTERN)

    bottleneck_specific = False
    if bottleneck_index is not None:
        nearby = "。".join(
            lead_sentences[max(0, bottleneck_index - 1):bottleneck_index + 2]
        )
        bottleneck_specific = bool(
            _REVIEWER_MOMENTUM_SPECIFICITY_PATTERN.search(nearby)
        )

    payoff_index = None
    if move_index is not None:
        payoff_index = next(
            (
                index
                for index in range(move_index, min(len(lead_sentences), move_index + 9))
                if _REVIEWER_MOMENTUM_PAYOFF_PATTERN.search(
                    lead_sentences[index]
                )
            ),
            None,
        )

    method_terms_before_move: list[str] = []
    method_scope_end = move_index if move_index is not None else min(
        len(lead_sentences), 8
    )
    generic_method_terms = {
        "有限要素", "数値解析", "解析手法", "モデル", "アルゴリズム",
        "ソルバー", "メッシュ", "等価回路", "インターフェース",
        "モジュール",
    }
    for sentence in lead_sentences[:method_scope_end]:
        method_terms_before_move.extend(
            term
            for term in _REVIEWER_MOMENTUM_METHOD_PATTERN.findall(sentence)
            if term not in generic_method_terms
        )
        method_terms_before_move.extend(
            re.findall(
                r"(?<![A-Za-z0-9_])[A-Z][A-Za-z0-9_.()+/-]{1,}", sentence
            )
        )
    method_terms_before_move = list(dict.fromkeys(method_terms_before_move))

    risks: list[dict] = []

    def add_risk(
        risk_type: str,
        excerpt: str,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
        **details,
    ) -> None:
        item = {
            "type": risk_type,
            "severity": severity,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:360],
            "comment": comment,
            "recommendation": recommendation,
        }
        item.update(details)
        risks.append(item)

    if stake_index is None or stake_index > 3:
        excerpt = "。".join(lead_sentences[:4])
        add_risk(
            "reviewer_stakes_delayed",
            excerpt,
            "導入部で、誰のどの判断・性能・負担が変わる研究かが見えない。",
            (
                "最初の一〜二文で、対象分野の読者が既に価値を理解できる"
                "設計判断、時間、精度、損失、人への効果のいずれかを示す。"
            ),
        )

    if bottleneck_index is None or not bottleneck_specific:
        excerpt = (
            "。".join(lead_sentences[:6])
            if bottleneck_index is None
            else lead_sentences[bottleneck_index]
        )
        add_risk(
            "bottleneck_not_concrete",
            excerpt,
            "課題が抽象語に留まり、なぜ今の方法では前へ進めないかを想像しにくい。",
            (
                "未解決の操作、手作業、時間、誤差、適用不能、試行錯誤など、"
                "研究で取り除く具体的な詰まりを一つ示す。"
            ),
        )

    if (
        method_index is not None
        and method_index <= 1
        and (stake_index is None or stake_index > method_index)
        and (bottleneck_index is None or bottleneck_index > method_index)
    ):
        add_risk(
            "method_before_problem",
            lead_sentences[method_index],
            "対象分野の困りごとより先に手法・モデルが現れ、研究の必要性が後付けに見える。",
            (
                "方法名の前に、既存の選択や運用が何に阻まれているかを述べる。"
                "手法はその障壁を外す一手として導入する。"
            ),
            severity="HIGH",
        )

    if len(method_terms_before_move) >= 5:
        add_risk(
            "lead_method_inventory",
            "、".join(method_terms_before_move),
            "導入部の主役が増え、研究課題より手法名の在庫表として読まれる。",
            (
                "導入部の中核概念は二〜三個に絞る。残りの固有名・略語は、"
                "研究項目、図、実行可能性の根拠へ下げる。"
            ),
            method_terms=method_terms_before_move,
        )

    if move_index is None:
        add_risk(
            "research_move_missing_from_lead",
            "。".join(lead_sentences[-4:]),
            "問題を示した後、申請者が何を変えるかが導入部で宣言されていない。",
            (
                "『そこで本研究では』に続け、研究期間内の一手を一文で述べる。"
                "技術要素の列挙ではなく、障壁をどう扱う研究かを主語述語で示す。"
            ),
            severity="HIGH",
        )
    elif payoff_index is None:
        add_risk(
            "payoff_missing_after_research_move",
            "。".join(lead_sentences[move_index:move_index + 4]),
            "研究の一手はあるが、それにより何が判断・実現できるかが続かない。",
            (
                "方法の直後に、短縮・精度・設計判断・新たに可能になる検証など、"
                "審査者が成果を思い描ける変化を置く。"
            ),
        )

    unsupported_hype: list[dict] = []
    for index, sentence in enumerate(lead_sentences):
        hype = _REVIEWER_MOMENTUM_HYPE_PATTERN.search(sentence)
        if hype is None:
            continue
        nearby = "。".join(
            lead_sentences[max(0, index - 2):min(len(lead_sentences), index + 3)]
        )
        if _REVIEWER_MOMENTUM_EVIDENCE_PATTERN.search(nearby):
            continue
        unsupported_hype.append({
            "index": index,
            "phrase": hype.group(0),
            "excerpt": sentence,
        })
    if unsupported_hype:
        first = unsupported_hype[0]
        add_risk(
            "unsupported_excitement_language",
            first["excerpt"],
            "強い形容だけで期待を作っており、研究上の緊張や根拠がない。",
            (
                "『画期的』『世界初』を足す代わりに、新たな可能性、まだ使えない"
                "理由、本研究が外す障壁を具体化する。"
            ),
            examples=unsupported_hype[:5],
        )

    arc_complete = (
        stake_index is not None
        and stake_index <= 3
        and bottleneck_index is not None
        and bottleneck_specific
        and move_index is not None
        and payoff_index is not None
        and stake_index <= bottleneck_index <= move_index <= payoff_index
    )
    productive_tension = (
        opportunity_index is not None
        and bottleneck_index is not None
        and opportunity_index <= bottleneck_index
    )
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(dict.fromkeys(
        risk["recommendation"] for risk in risks
    ))
    return {
        "applicable": True,
        "score": None,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "metrics": {
            "lead_sentence_count": len(lead_sentences),
            "lead_character_count": lead_chars,
            "stake_sentence": stake_index,
            "opportunity_sentence": opportunity_index,
            "bottleneck_sentence": bottleneck_index,
            "bottleneck_is_specific": bottleneck_specific,
            "research_move_sentence": move_index,
            "payoff_sentence": payoff_index,
            "method_sentence": method_index,
            "method_terms_before_move": method_terms_before_move,
            "arc_complete": arc_complete,
            "productive_tension": productive_tension,
        },
        "revision_protocol": {
            "sequence": [
                "reviewer_visible_stakes",
                "concrete_bottleneck_or_unused_opportunity",
                "research_move",
                "observable_payoff",
                "bounded_evidence",
            ],
            "instruction": (
                "Create interest from a concrete unresolved tension, not from "
                "adjectives. Keep the opening's core concepts few, and introduce "
                "method names only after the reviewer understands the problem."
            ),
        },
        "diagnosis": (
            "Strong adopted prose is easy to enter and gives the reviewer a reason "
            "to continue: a familiar value is blocked by a specific obstacle, and "
            "the proposed research unlocks an observable change."
        ),
        "source": (
            "generic reviewer-momentum diagnostic informed by adopted grant prose "
            "and adjacent-reader feedback; non-scoring"
        ),
    }


def grant_writing_count_weak_expressions(text: str) -> dict:
    """Count hedges and grant-specific non-commitment phrases."""
    text = _prose_for_lint(_read_text_if_path(text))
    total, by_pattern = scan_hedges(text)
    patterns = dict(HEDGE_PATTERNS)
    patterns.update({
        "検討する": r"検討(?:する|します|を行う|を進める)",
        "目指す": r"目指(?:す|します)",
        "努める": r"努め(?:る|ます)",
        "今後": r"今後",
    })
    # Recompute with the grant-specific additions included.
    total = 0
    by_pattern = {}
    for name, pat in patterns.items():
        n = len(re.findall(pat, text))
        if n:
            by_pattern[name] = n
            total += n

    # 「など」 hedges only when it trails a single item. After an enumeration
    # -- 「高次要素、補助空間前処理など」 -- it is ordinary Japanese for "and
    # the like", and counting it told an applicant to damage correct prose.
    # A fixed-width lookbehind cannot express "no comma earlier in the
    # clause", so it is counted here instead of in the pattern table.
    # The form's own section headings are not the applicant's prose.
    form_headings = ("研究目的、研究方法など", "研究目的、研究方法等")
    nado = 0
    for match in re.finditer("など", text):
        window = text[max(0, match.start() - 12):match.end() + 2]
        if any(h in window for h in form_headings):
            continue
        clause_start = max(
            text.rfind(mark, 0, match.start())
            for mark in ("。", "．", "\n", "、", "，")
        )
        clause = text[clause_start + 1:match.start()]
        if "、" in clause or "，" in clause:
            continue  # enumeration
        # An enumeration whose last separator is the nearest mark also ends
        # up here, so require the clause to look like a lone noun.
        if len(clause) > 12:
            continue
        preceding = text[max(0, clause_start - 30):clause_start]
        if "、" in preceding or "，" in preceding:
            continue
        nado += 1
    if nado:
        by_pattern["など"] = nado
        total += nado
    return {
        "total_weak_expressions": total,
        "by_pattern": by_pattern,
        "target": "0 in aims / outcomes; replace with deliverables and dates",
        "warning": ">=5 suggests the plan reads as intent rather than execution",
    }


_GENERIC_AXES = {
    "social_or_scientific_problem": ["課題", "問題", "ニーズ", "必要", "社会", "産業"],
    "objective": ["目的", "本研究", "本提案", "実現", "構築", "開発"],
    "novelty": ["新規", "独自", "先行", "特徴", "差別化", "初"],
    "method": ["手法", "方法", "解析", "実験", "検証", "評価"],
    "feasibility": ["実績", "経験", "環境", "体制", "遂行", "可能"],
    "schedule": ["年度", "計画", "スケジュール", "1年目", "2年目", "3年目"],
    "outcomes": ["成果", "公開", "発表", "論文", "レポート", "デモ"],
    "budget": ["予算", "経費", "費用", "助成", "使途", "計算資源"],
}

_KDDI_DIGITAL_AXES = {
    "social_issue": ["社会的課題", "地域", "中小企業", "製造業", "人材", "現場"],
    "social_implementation": ["社会実装", "実装", "PoC", "実証", "試作", "技術プレゼン"],
    "digital_use": ["デジタル", "AI", "生成AI", "MCP", "LTspice", "SPICE", "CAE"],
    "field_validation": ["基板", "計測", "評価", "実験", "熱", "EMC", "厚銅"],
    "open_outputs": ["公開", "OSS", "レポジトリ", "デモ", "教材", "報告書"],
    "schedule": ["年度", "1年目", "2年目", "3年目", "スケジュール", "マイルストーン"],
    "budget_alignment": ["予算", "Claude", "Codex", "Fable", "MDX", "計算資源"],
    "applicant_fit": ["三菱電機", "EMC", "IH", "Radia", "NGSolve", "LTspice", "実績"],
}

_KAKEN_OSS_PLATFORM_AXES = {
    "academic_question": [
        "学術的な問い",
        "学術的問い",
        "研究課題",
        "問い",
        "共同研究サイクル",
    ],
    "technical_reports_as_source": [
        "技術報告",
        "報告書",
        "日本語",
        "知識源",
        "既存知",
    ],
    "lab_silo_and_ai_urgency": [
        "属研究室化",
        "研究室単位",
        "重複実装",
        "重複開発",
        "AIが開発を加速",
        "AIによる開発加速",
        "AIにより開発が加速",
    ],
    "jpmars_github_governance": [
        "JP-MARs",
        "GitHub",
        "issue",
        "pull request",
        "CI",
        "ライセンス",
    ],
    "executable_public_outputs": [
        "問題仕様",
        "参照実装",
        "ベンチマーク",
        "教材",
        "mcp-server",
        "ワークショップ",
    ],
    "reuse_and_upstream_first": [
        "既存OSS",
        "再利用",
        "上流貢献",
        "上流で改良",
        "meshio",
        "機能重複",
    ],
    "scientific_quality_gate": [
        "試験データ",
        "期待値",
        "許容誤差",
        "既知制約",
        "参照実装",
        "別機関",
        "CI",
    ],
    "ai_machine_reexecution": [
        "AI",
        "機械実行",
        "再実行",
        "Python/MCP",
        "環境構築支援",
        "依存関係",
    ],
    "tacit_knowledge_conversion": [
        "暗黙知",
        "解析条件",
        "境界条件",
        "メッシュ",
        "収束",
        "失敗例",
    ],
    "domestic_international_collaboration": [
        "国内外",
        "国際",
        "共同実装",
        "コードレビュー",
        "TU Wien",
        "TU Graz",
        "IGTE",
    ],
    "education_and_capacity": [
        "学生",
        "若手",
        "博士",
        "教育",
        "技術力",
        "実装力",
        "検証力",
    ],
    "environment_portability": [
        "mdx",
        "Windows",
        "Linux",
        "スパコン",
        "GPU",
        "ARM",
        "アーキテクチャ",
        "可搬",
        "移植",
    ],
    "maintenance_governance": [
        "運営",
        "保守",
        "CI",
        "リリース",
        "ドキュメント",
        "再現性",
        "貢献者",
    ],
}

_BUDGET_POLICY = (
    "申請予算は遠慮して小さく見せず、助成上限に近い額まで必要な計画として組む。"
    "ただし、上限近くでも不自然に見えないよう、単価・数量・月数/回数・年度配分・"
    "見積根拠を具体的に積算し、検証ループと社会実装に直結する経費として説明する。"
    "外部料金は公式料金表または見積書に基づき、提供者、URL、料金年度・改定日、参照日、"
    "税込/税抜、最低購入単位、有効期限、通貨・為替、端数処理を記録する。申請年度の料金が"
    "未公表なら現行料金による暫定積算と再確認時期を示し、公式単価、価格変動への予備幅、"
    "使用上限を区別する。費目別・年度別の総額を再計算し、直接経費上限、研究期間、間接経費"
    "の扱いと照合する。旅費は開催地の公表状況を区別し、人数、泊数、運賃、登録費、宿泊費、"
    "日当・現地交通へ分解する。"
    "最大費目を明示して中心的な研究行為へ対応付け、その他費目はサービス別に分解し、"
    "サーバ、CI、保存等の二重計上を避ける。AIと計算資源を併用する場合は、AI推論、"
    "常時稼働の低費用検証、短期集中計算を分け、各層の単位、期間、成果物を示す。"
    "GPU等の機種名を記す場合は、アクセラレータとホストCPUを公式仕様で区別し、"
    "実在する構成だけを予算化する。高速化の成功でなく、異機種間で保存すべき物理量・"
    "設計判断と不一致時の扱いを示す。"
    "科研費の基盤系種目では、採択時に申請額の約7割程度へ減額されて内定する場合が"
    "多い(充足率)。減額後も検証ループが成立する経費の優先順位を用意する。"
    "挑戦的研究は原則満額支給である。"
)

_BUDGET_AXIS_COMMENTS = {
    "ai_agent_costs": (
        "AI agent / LLM costs are thin. Tie Claude, Codex, Fable, or related tools "
        "to proposal drafting, design automation, validation, or public deliverables."
    ),
    "compute_resources": (
        "Compute costs are thin. Tie servers, cloud, GPU, or MDX-like resources to "
        "the verification loop and expected run volume."
    ),
    "poc_experiment": (
        "PoC / experiment costs are thin. Tie boards, parts, consumables, measurement, "
        "or evaluation work to concrete implementation evidence."
    ),
    "dissemination": (
        "Dissemination costs are thin. Tie travel, presentations, workshops, or reports "
        "to handoff and social implementation."
    ),
    "near_ceiling_strategy": (
        "予算を遠慮して小さく見せるのではなく、助成上限額に近い申請額が必要である方針を明記する。"
    ),
    "itemized_calculation": (
        "上限近くでも不自然に見えないよう、単価、数量、月数/回数、年度配分、見積根拠を具体的に積算する。"
    ),
    "pricing_provenance": (
        "外部料金は公式料金表または見積書へ遡れるよう、URL、料金年度・参照日、税込/税抜、"
        "最低購入単位、有効期限、通貨・為替、端数処理を記録する。"
    ),
}

_POWER_ELECTRONICS_FOCUS_TRIGGERS = [
    "パワーエレクトロニクス",
    "パワエレ",
    "厚銅",
    "LTspice",
    "Radia",
    "NGSolve",
    "CAE-AI",
]

_POWER_ELECTRONICS_FOCUS_AXES = {
    "main_theme_specificity": [
        "パワーエレクトロニクス基板",
        "パワエレ基板",
        "電源基板",
        "厚銅基板",
        "CAE-AI",
    ],
    "board_physics": [
        "寄生",
        "EMC",
        "熱",
        "電流集中",
        "放熱",
        "浮遊容量",
        "インダクタンス",
    ],
    "ai_mcp_loop": ["MCP", "AI", "LLM", "生成AI", "自律駆動", "自然言語"],
    "tool_chain": ["LTspice", "SPICE", "Radia", "NGSolve", "PEEC"],
    "social_target": ["1000人以下", "中小企業", "地域製造業", "製造業"],
    "poc_handoff": ["PoC", "試作", "技術プレゼン", "試験導入", "導入候補", "MotorAI"],
    "commercial_positioning": ["商用CAE", "置換", "入口", "届かない", "習得コスト"],
    "llm_native_advantage": ["Python-native", "現代的", "API", "ツール呼び出し"],
}


_KAKEN_REVIEW_CRITERIA_AXES = {
    "academic_importance": [
        "学術的重要性",
        "学術的意義",
        "学術的問い",
        "独自性",
        "独創性",
        "波及",
    ],
    "method_validity": [
        "研究方法",
        "研究計画",
        "妥当",
        "検証",
        "実証",
        "評価方法",
        "達成指標",
    ],
    "feasibility_environment": [
        "遂行能力",
        "研究環境",
        "研究実績",
        "準備状況",
        "実施体制",
        "研究設備",
        "予備",
    ],
    "internationality": [
        "国際性",
        "世界の研究",
        "国際共同",
        "国内外",
        "我が国独自",
        "日本独自",
    ],
}

_KAKEN_BRIEFING_NOTES = [
    (
        "研究計画の評定要素は3つ: (1)研究課題の学術的重要性、(2)研究方法の妥当性、"
        "(3)研究遂行能力及び研究環境の適切性。これとは別に国際性も評定される。"
    ),
    (
        "研究経費は、研究方法の妥当性の中で計画との整合性を見られ、さらに配分額の"
        "判断に用いる別枠の経費妥当性評価を受ける。独立した同格の総合評点軸ではない。"
    ),
    (
        "審査委員は約1ヶ月で多い場合100件程度の計画調書を審査する。"
        "専門外の読者でも読みやすい調書が圧倒的に採択されやすい。"
    ),
    (
        "カラーの図・写真は審査時に白黒印刷される種目がある。"
        "色の違いだけで系列を区別しない。"
    ),
    (
        "審査ではresearchmapが研究者番号で参照される。"
        "応募前に更新と研究者番号の登録を確認する。"
    ),
    (
        "「人権の保護及び法令等の遵守への対応」欄は例年審査委員からの指摘が"
        "非常に多い。該当なしの場合も判断根拠を一文添える。"
    ),
    (
        "基盤系種目は申請額の約7割程度への減額内定が多い(充足率)。"
        "挑戦的研究は原則満額支給だが採択率が低く、基盤研究との重複応募を検討する。"
    ),
]

_SUPPORTED_PROGRAMS = frozenset({
    "generic", "kaken_generic", "kaken_oss", "kaken_oss_platform",
    "kddi_digital",
})


def _validate_program(program: str) -> None:
    if program not in _SUPPORTED_PROGRAMS:
        choices = ", ".join(sorted(_SUPPORTED_PROGRAMS))
        raise ValueError(f"unknown grant program {program!r}; choose one of: {choices}")


def _section_axes_for_program(program: str) -> dict[str, list[str]]:
    _validate_program(program)
    if program == "kddi_digital":
        return _KDDI_DIGITAL_AXES
    if program == "kaken_generic":
        return _KAKEN_REVIEW_CRITERIA_AXES
    if program in {"kaken_oss", "kaken_oss_platform"}:
        return _KAKEN_OSS_PLATFORM_AXES
    return _GENERIC_AXES


def grant_writing_section_presence(text: str, program: str = "generic") -> dict:
    """Check whether a proposal draft contains the expected review axes."""
    text = _prose_for_lint(_read_text_if_path(text))
    low = text.lower()
    axes = _section_axes_for_program(program)
    axis_results = {}
    missing = []
    for axis, keywords in axes.items():
        matches = _contains_any(low, keywords)
        ok = bool(matches)
        axis_results[axis] = {
            "ok": ok,
            "matches": matches[:8],
            "keywords": keywords,
        }
        if not ok:
            missing.append(axis)
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "program": program,
        "score": score,
        "axes_total": len(axes),
        "axes_present": len(axes) - len(missing),
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "target": "all core review axes present at least once",
    }


def grant_writing_kddi_digital_check(text: str) -> dict:
    """KDDI Foundation Digital Innovation / social implementation check."""
    text = _read_text_if_path(text)
    presence = grant_writing_section_presence(text, program="kddi_digital")
    comments = []
    for axis in presence["missing_axes"]:
        comments.append(f"Missing KDDI social-implementation axis: {axis}")
    score = presence["score"]
    if "PoC" not in text and "実証" not in text and "試作" not in text:
        comments.append("Add a concrete PoC / field validation output.")
        score = max(0.0, score - 1.0)
    if "公開" not in text and "OSS" not in text and "レポジトリ" not in text:
        comments.append("Add a shareable output: repository, demo, report, or workshop material.")
        score = max(0.0, score - 1.0)
    return {
        "score": round(score, 1),
        "missing_required_count": len(comments),
        "comments": comments,
        "axis_results": presence["axis_results"],
        "source": "KDDI Digital Innovation social-implementation proposal axes",
    }


def grant_writing_kddi_power_electronics_focus_check(text: str) -> dict:
    """Check the current KDDI power-electronics-board CAE-AI framing.

    This is intentionally domain-specific.  Use it when the KDDI Digital
    Innovation proposal is about LTspice/Radia/NGSolve/MCP for power
    electronics boards.  It checks that "companies below 1000 employees" is
    kept as the implementation target, not as the main theme.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    axis_results = {}
    missing = []
    for axis, keywords in _POWER_ELECTRONICS_FOCUS_AXES.items():
        matches = _contains_any(low, keywords)
        axis_results[axis] = {
            "ok": bool(matches),
            "matches": matches[:8],
            "keywords": keywords,
        }
        if not matches:
            missing.append(axis)

    score = round(
        10.0 * (len(_POWER_ELECTRONICS_FOCUS_AXES) - len(missing))
        / len(_POWER_ELECTRONICS_FOCUS_AXES),
        1,
    )
    comments = [f"Power-electronics focus axis missing or thin: {axis}" for axis in missing]

    sentences = [s.strip() for s in re.split(r"[。．!?！？]", text) if s.strip()]
    theme_sentences = [
        s for s in sentences
        if any(token in s for token in ("主題", "目的", "テーマ案", "テーマ", "本研究"))
    ][:5]
    if theme_sentences:
        theme_blob = "。".join(theme_sentences)
        if (
            "パワーエレクトロニクス基板" not in theme_blob
            and "パワエレ基板" not in theme_blob
            and "CAE-AI" not in theme_blob
        ):
            comments.append(
                "The opening theme sentences do not clearly state the power-electronics-board CAE-AI subject."
            )
            score = max(0.0, score - 1.0)

    risky_phrases = []
    for phrase in ("一般的なCAE導入", "すべての製造業CAE", "商用CAEの代替"):
        if phrase in text:
            risky_phrases.append(phrase)
    if "商用CAE" in text and "置換する" in text and "直ちに置換するものではない" not in text:
        risky_phrases.append("商用CAEを置換する")
    if risky_phrases:
        comments.append(
            "Avoid broad or adversarial positioning; keep commercial CAE as powerful but hard to access: "
            + ", ".join(risky_phrases)
        )
        score = max(0.0, score - 1.0)

    recommendations = [
        "Make the subject: power-electronics-board circuit/EM/thermal CAE-AI environment.",
        "Keep companies below 1000 employees as the first user and implementation field.",
        "Frame commercial CAE as powerful but expensive/difficult; the proposal creates an AI/MCP entry point.",
        "Tie the PoC to a board-level outcome: parasitics, EMC risk, heat, measurement, and handoff.",
    ]
    return {
        "score": round(score, 1),
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "comments": comments,
        "theme_sentence_examples": theme_sentences,
        "recommendations": recommendations,
        "target": "specific power-electronics-board CAE-AI theme; 1000-person firms as implementation target",
        "source": "KDDI power-electronics-board CAE-AI framing check",
    }


def grant_writing_kaken_oss_platform_check(text: str) -> dict:
    """Check KAKENHI framing for an AI-era OSS research platform proposal."""
    text = _read_text_if_path(text)
    low = text.lower()
    presence = grant_writing_section_presence(text, program="kaken_oss")
    comments = [
        f"KAKEN OSS platform axis missing or thin: {axis}"
        for axis in presence["missing_axes"]
    ]
    score = presence["score"]

    public_outputs = [
        "公開リポジトリ",
        "問題仕様",
        "参照実装",
        "ベンチマーク",
        "software.html",
        "elemag/index.php",
        "geometry.php",
        "mcp-server",
    ]
    output_hits = _contains_any(low, public_outputs)
    if len(output_hits) < 4:
        comments.append(
            "Name at least four concrete public outputs: repository, problem specs, "
            "reference implementations, benchmarks, teaching pages, or mcp-server."
        )
        score = max(0.0, score - 1.0)

    if "技術報告" not in text and "報告書" not in text:
        comments.append(
            "Treat existing Japanese technical reports as a research knowledge source, "
            "then connect them to executable assets and review history."
        )
        score = max(0.0, score - 0.5)

    silo_hits = _contains_any(
        low,
        ["属研究室化", "研究室単位", "研究室内に閉じ", "重複実装", "重複開発"],
    )
    ai_urgency_hits = _contains_any(
        low,
        [
            "AIが開発を加速",
            "AIによる開発加速",
            "AIにより開発が加速",
            "AI時代",
            "高速に再生産",
            "増幅",
        ],
    )
    if not silo_hits or not ai_urgency_hits:
        comments.append(
            "Make the why-now explicit: AI can accelerate duplicated, lab-siloed "
            "implementations unless shared specifications, tests, and review convert "
            "them into collaborative assets."
        )
        score = max(0.0, score - 0.5)

    reuse_hits = _contains_any(
        low,
        ["既存OSS", "再利用", "上流貢献", "上流で改良", "upstream", "meshio", "機能重複"],
    )
    if not reuse_hits:
        comments.append(
            "Add an upstream-first gate: survey existing OSS, reuse or contribute "
            "upstream, and justify new code with a tested domain-specific gap."
        )
        score = max(0.0, score - 1.0)

    quality_groups = {
        "test_corpus": ["試験データ", "テストデータ", "検証データ"],
        "expected_tolerance": ["期待値", "許容誤差", "参照値"],
        "automated_validation": ["CI", "自動テスト", "再実行"],
        "scope_and_provenance": ["既知制約", "適用範囲", "由来", "出典"],
        "independent_reproduction": ["別機関", "第三者", "独立検証", "再現確認"],
    }
    quality_results = {
        name: _contains_any(low, keywords)
        for name, keywords in quality_groups.items()
    }
    if sum(bool(hits) for hits in quality_results.values()) < 4:
        comments.append(
            "Do not treat GitHub as a repository dump. Require test provenance, "
            "expected values/tolerances, CI, documented limits, and independent re-execution."
        )
        score = max(0.0, score - 1.0)

    responsibility_phrases = ["利用者の自己責任", "使う人の自己責任"]
    responsibility_negations = [
        "自己責任を品質保証の代替にしない",
        "自己責任を理由にreleaseしない",
        "利用者への責任転嫁ではなく",
    ]
    if (
        _contains_any(low, responsibility_phrases)
        and not _contains_any(low, responsibility_negations)
    ):
        comments.append(
            "A no-warranty license does not replace scientific validation; do not "
            "shift supported-scope verification responsibility to users."
        )
        score = max(0.0, score - 1.0)

    domestic_evidence = ["伊田", "HACApK", "CEFC", "ADVENTURE", "国内共同"]
    overseas_evidence = [
        "Hollaus",
        "TU Wien",
        "TU Graz",
        "IGTE",
        "openCFS",
        "NGSolve",
    ]
    domestic_hits = _contains_any(low, domestic_evidence)
    overseas_hits = _contains_any(low, overseas_evidence)
    if not domestic_hits:
        comments.append(
            "Add named domestic preliminary evidence such as the Ida/HACApK/CEFC "
            "collaboration or an ADVENTURE connection."
        )
        score = max(0.0, score - 0.5)
    if not overseas_hits:
        comments.append(
            "Add named overseas preliminary evidence such as Hollaus/TU Wien/IGTE, "
            "TU Graz/openCFS, or the NGSolve community."
        )
        score = max(0.0, score - 0.5)

    catch_up_negations = [
        "追いつくことではない",
        "追いつく計画ではない",
        "追いつく計画ではなく",
        "追いつくことを目的にしない",
        "欧州を一方向の到達目標",
    ]
    if "追いつく" in text and not any(phrase in text for phrase in catch_up_negations):
        comments.append(
            "Avoid a catch-up frame; state reciprocal collaboration and "
            "complementarity instead of Europe as a one-way benchmark."
        )
        score = max(0.0, score - 0.5)

    environment_groups = {
        "enterprise_os": ["Windows", "Linux", "企業"],
        "shared_compute": ["mdx", "クラウド", "仮想環境"],
        "hpc": ["スパコン", "スーパーコンピュータ", "HPC"],
        "accelerator_generation": ["GPU", "A100", "H200", "ICCG"],
        "architecture_change": ["ARM", "arm64", "アーキテクチャ", "可搬", "移植"],
    }
    environment_results = {
        name: _contains_any(low, keywords)
        for name, keywords in environment_groups.items()
    }
    if sum(bool(hits) for hits in environment_results.values()) < 3:
        comments.append(
            "Frame compute as a portable execution matrix spanning enterprise "
            "Windows/Linux, mdx, HPC, GPU generations, and future CPU architectures."
        )
        score = max(0.0, score - 1.0)

    hardware_purchase_terms = ["計算機を購入", "マシンを購入", "GPUを購入"]
    portability_terms = ["可搬", "移植", "再実行", "複数環境", "アーキテクチャ"]
    if _contains_any(low, hardware_purchase_terms) and not _contains_any(low, portability_terms):
        comments.append(
            "Do not make hardware acquisition the research objective; explain how "
            "specifications, containers, CI, and AI-assisted setup survive platform change."
        )
        score = max(0.0, score - 1.0)

    platform_mentions = len(re.findall(r"JP[-‐‑–—]?MARs", text, flags=re.IGNORECASE))
    radia_mentions = len(
        re.findall(r"(?<![A-Za-z])Radia(?:/radia-mcp|-mcp)?", text, flags=re.IGNORECASE)
    )
    platform_focus = {
        "ok": platform_mentions >= max(1, radia_mentions),
        "jpmars_mentions": platform_mentions,
        "radia_mentions": radia_mentions,
    }
    if not platform_focus["ok"]:
        comments.append(
            "Keep JP-MARs as the governing research platform; position Radia and "
            "other named software as preliminary evidence or demonstrators."
        )
        score = max(0.0, score - 1.0)

    radia_integration = {
        "mentioned": bool(radia_mentions),
        "positioned_as_integration_evidence": bool(
            radia_mentions
            and _contains_any(low, ["統合", "融合", "接続", "先行実証", "インタフェース"])
        ),
    }
    if radia_integration["mentioned"] and not radia_integration["positioned_as_integration_evidence"]:
        comments.append(
            "When Radia is named, use it as preliminary evidence that specialized "
            "software and methods can be integrated, not as the proposal's main product."
        )
        score = max(0.0, score - 0.5)

    recommendations = [
        "Open with how technical reports become reviewable, executable research assets in the AI era.",
        "Use GitHub/JP-MARs to preserve issue, implementation, validation, review, and contributor history.",
        "Define lab-level siloing and explain why AI otherwise accelerates duplicated and unverified development.",
        "Survey existing OSS first; reuse it or contribute upstream before authoring a new implementation.",
        "Separate license disclaimers from scientific responsibility and require test data, tolerances, limits, CI, and independent re-execution.",
        "Show both domestic and overseas preliminary collaborations with names and outputs.",
        "Treat AI as an accelerator for search, translation, environment setup, and re-execution; researchers retain physical validation.",
        "Define a portable execution matrix across Windows/Linux, mdx, HPC, GPU generations, and future architectures such as ARM.",
        "Split platform adaptation into issues and CI jobs so institutions and companies can maintain different environments through pull requests.",
        "Keep hardware acquisition subordinate to reproducible specifications, CI, containers, and long-term governance.",
    ]
    return {
        "score": round(score, 1),
        "missing_count": presence["missing_count"],
        "missing_axes": presence["missing_axes"],
        "axis_results": presence["axis_results"],
        "public_output_hits": output_hits,
        "domestic_evidence_hits": domestic_hits,
        "overseas_evidence_hits": overseas_hits,
        "environment_results": environment_results,
        "silo_hits": silo_hits,
        "ai_urgency_hits": ai_urgency_hits,
        "reuse_hits": reuse_hits,
        "quality_results": quality_results,
        "platform_focus": platform_focus,
        "radia_integration": radia_integration,
        "comments": comments,
        "recommendations": recommendations,
        "target": "KAKENHI OSS platform: overcome lab-level siloing through upstream-first reuse, scientific quality gates, reciprocal collaboration, and AI-assisted portability",
        "source": "KAKENHI AI-era OSS research-platform framing check",
    }


def grant_writing_internal_evidence_to_external_scale_check(text: str) -> dict:
    """Check whether an internal success is evidence for external transfer.

    The check is optional: it becomes applicable only when a draft claims an
    existing pilot, operation, or preliminary result. It then asks what is
    transferred, to whom, through which route, and how independent success is
    verified. This keeps the rule useful across grant programs and domains.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    axes = {
        "internal_evidence": [
            "予備成果",
            "運用実績",
            "実証済",
            "試行済",
            "既に実現",
            "すでに実現",
            "既に運用",
            "すでに運用",
            "機能している",
            "利用している",
            "活用している",
            "導入済",
        ],
        "observed_users_or_context": [
            "利用者",
            "ユーザー",
            "現場",
            "研究室",
            "学内",
            "社内",
            "組織内",
            "教員",
            "学生",
            "担当者",
            "企業",
        ],
        "transferable_unit": [
            "問題仕様",
            "参照実装",
            "入力データ",
            "試験データ",
            "検証結果",
            "実装",
            "コード",
            "手順",
            "API",
            "教材",
            "ドキュメント",
            "ベンチマーク",
            "知識基盤",
            "データセット",
        ],
        "external_route": [
            "他機関",
            "別機関",
            "研究室間",
            "複数機関",
            "国内外",
            "国際",
            "共同研究",
            "外部展開",
            "公開",
            "社会実装",
            "技術移転",
            "普及",
            "上流貢献",
        ],
        "external_actor": [
            "第三者",
            "他機関",
            "別機関",
            "外部利用者",
            "共同研究先",
            "企業",
            "海外研究者",
        ],
        "independent_validation": [
            "再実行",
            "再現",
            "比較",
            "評価",
            "検証",
            "レビュー",
            "受入",
            "フィードバック",
        ],
    }
    axis_results = {}
    for name, keywords in axes.items():
        hits = _contains_any(low, keywords)
        axis_results[name] = {
            "ok": bool(hits),
            "matches": hits[:8],
            "keywords": keywords,
        }

    applicable = axis_results["internal_evidence"]["ok"]
    if not applicable:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": axis_results,
            "comments": [],
            "target": (
                "when internal evidence is claimed, connect it to a transferable "
                "unit, external route and actor, and independent validation"
            ),
            "source": "generic internal-evidence-to-external-scale check",
        }

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "observed_users_or_context": "内部実証を誰がどの環境で利用したかを示す。",
        "transferable_unit": (
            "外部へ渡す単位を、仕様、実装、データ、手順等として明示する。"
        ),
        "external_route": (
            "内部実証を他機関・共同研究・社会実装へ接続する経路を示す。"
        ),
        "external_actor": "内部関係者ではない利用者・検証者を明示する。",
        "independent_validation": (
            "外部での再実行、比較、評価等により移転の成否を判定する。"
        ),
    }
    comments = [comments_by_axis[name] for name in missing if name in comments_by_axis]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "comments": comments,
        "target": (
            "internal success is feasibility evidence, not external validity; state "
            "what travels, to whom, by which route, and how success is verified"
        ),
        "source": "generic internal-evidence-to-external-scale check",
    }


def grant_writing_domain_outcome_chain_check(text: str) -> dict:
    """Check that a platform or tool proposal ends in domain knowledge.

    Infrastructure, OSS, APIs, repositories, and AI interfaces can make a
    study feasible, but they are rarely the academic outcome by themselves.
    This optional check becomes applicable when such enabling technology is
    prominent and asks for the target object, measurable quantity, conditional
    knowledge product, changed decision, falsification gate, and an explicit
    statement that the technology is subordinate to the research question.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    tool_terms = [
        "プラットフォーム",
        "研究基盤",
        "ソフトウェア",
        "OSS",
        "github",
        "gitlab",
        "MCP",
        "API",
        "リポジトリ",
        "モジュール",
        "インターフェース",
        "インタフェース",
    ]
    tool_hits = _contains_any(low, tool_terms)
    if not tool_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "tool_hits": [],
            "comments": [],
            "target": "enabling technology must terminate in a field-specific, falsifiable knowledge product",
            "source": "generic tool-to-domain-outcome chain check",
        }

    axes = {
        "domain_object": [
            "機器",
            "装置",
            "材料",
            "回路",
            "磁気浮上",
            "電動機",
            "患者",
            "診断",
            "製造工程",
            "観測対象",
            "制御対象",
        ],
        "measurable_domain_quantity": [
            "設計量",
            "目的量",
            "性能",
            "損失",
            "力",
            "温度",
            "平衡",
            "剛性",
            "感度",
            "精度",
            "収率",
            "寿命",
            "制約",
        ],
        "conditional_knowledge_product": [
            "選択則",
            "設計則",
            "成立条件",
            "成立域",
            "適用範囲",
            "適用境界",
            "不能条件",
            "選択域",
            "保存条件",
            "変化する条件",
            "閾値",
        ],
        "decision_consequence": [
            "設計判断",
            "設計順位",
            "順位を確定",
            "採否",
            "選択する",
            "切り替",
            "進むべき",
            "停止条件",
            "追加解析",
            "高忠実度",
        ],
        "falsifiable_gate": [
            "合格条件",
            "不合格",
            "反証",
            "許容差",
            "包含",
            "棄却",
            "失敗条件",
            "適用不能",
            "満たさない場合",
            "主張しない",
        ],
        "tool_subordination": [
            "問いではなく",
            "手段である",
            "手段とし",
            "実装手段",
            "検証手段",
            "道具である",
            "を用いる",
            "が支援する",
        ],
    }
    axis_results = {}
    for name, keywords in axes.items():
        hits = _contains_any(low, keywords)
        axis_results[name] = {
            "ok": bool(hits),
            "matches": hits[:8],
            "keywords": keywords,
        }

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "domain_object": "基盤を用いて研究する具体的な対象物・現象を明示する。",
        "measurable_domain_quantity": "対象分野で測る設計量・性能量・制約を明示する。",
        "conditional_knowledge_product": (
            "構築物でなく、成立条件、境界、選択則等の新しい知見を成果にする。"
        ),
        "decision_consequence": (
            "得られる知見が、どの設計・診断・採否を変えるかまで書く。"
        ),
        "falsifiable_gate": "知見を棄却または限定する事前の判定条件を置く。",
        "tool_subordination": (
            "OSS、AI、API、MCP、GitHub等は問いでなく検証手段だと明記する。"
        ),
    }
    comments = [comments_by_axis[name] for name in missing]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "tool_hits": tool_hits[:12],
        "comments": comments,
        "target": (
            "show a complete chain from enabling technology to a measured "
            "domain quantity, conditional knowledge, changed decision, and "
            "falsification gate"
        ),
        "source": "generic tool-to-domain-outcome chain check",
    }


def grant_writing_derived_metric_validation_check(text: str) -> dict:
    """Check calibration and falsification of a proposal-specific metric.

    A newly named interval, score, index, or composite criterion can be useful,
    but it must not be tuned and judged on the same examples. This optional
    check asks for an operational definition, observable components, a bounded
    calibration set, pre-test freezing, held-out validation, an acceptance
    threshold, and a consequence for failure.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    metric_terms = [
        "判定区間",
        "許容区間",
        "信頼区間",
        "評価指標",
        "判定指標",
        "複合指標",
        "合成指標",
        "リスクスコア",
        "評価スコア",
        "選択指数",
        "安全係数",
    ]
    metric_hits = _contains_any(low, metric_terms)
    if not metric_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "metric_hits": [],
            "comments": [],
            "target": "define, calibrate, freeze, independently validate, and falsify a proposal-specific metric",
            "source": "generic derived-metric validation check",
        }

    axes = {
        "operational_definition": [
            "定義する",
            "定義した",
            "算出式",
            "計算式",
            "合成則",
            "\\Delta",
            "I_q",
            "重み付き",
        ],
        "observable_components": [
            "対計算差",
            "成分",
            "残差",
            "測定値",
            "求解差",
            "離散化差",
            "連成差",
            "入力変数",
            "説明変数",
        ],
        "bounded_calibration_set": [
            "校正集合",
            "校正ケース",
            "校正データ",
            "参照データ",
            "基準問題",
            "訓練集合",
            "学習データ",
        ],
        "pretest_freeze": [
            "事前登録",
            "事前に固定",
            "事前固定",
            "凍結",
            "検証前に固定",
            "候補比較前に固定",
            "最終候補を得る前",
        ],
        "heldout_validation": [
            "保留データ",
            "独立検証",
            "検証集合",
            "テスト集合",
            "未使用データ",
            "外部検証",
            "holdout",
            "held-out",
        ],
        "acceptance_threshold": [
            "合格条件",
            "全点包含",
            "包含率",
            "許容差以内",
            "以上を合格",
            "以下を合格",
            "閾値",
            "判定基準",
        ],
        "failure_consequence": [
            "再校正せず",
            "主張しない",
            "不合格",
            "棄却",
            "適用境界",
            "反例",
            "失敗条件",
            "満たさない場合",
        ],
    }
    axis_results = {}
    for name, keywords in axes.items():
        hits = _contains_any(low, keywords)
        axis_results[name] = {
            "ok": bool(hits),
            "matches": hits[:8],
            "keywords": keywords,
        }

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "operational_definition": "新しい区間・指数・スコアの算出式または手順を定義する。",
        "observable_components": "指標を構成する観測可能な成分を列挙する。",
        "bounded_calibration_set": "校正に使うケース、範囲、件数を固定する。",
        "pretest_freeze": "独立検証を見る前に式・係数・閾値を凍結する。",
        "heldout_validation": "校正に使わない保留データまたは外部データで検証する。",
        "acceptance_threshold": "包含率、許容差、件数等の合否閾値を事前に置く。",
        "failure_consequence": "不合格時に再調整で救済せず、棄却・限定・代替経路を定める。",
    }
    comments = [comments_by_axis[name] for name in missing]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "metric_hits": metric_hits[:12],
        "comments": comments,
        "target": (
            "a proposal-specific metric must have an operational definition, "
            "separate calibration and held-out validation, a frozen acceptance "
            "rule, and a failure consequence"
        ),
        "source": "generic derived-metric validation check",
    }


def grant_writing_cross_organization_pilot_check(text: str) -> dict:
    """Check whether preliminary evidence crosses an organization boundary.

    User counts, repository publication, links, and internal operation are weak
    evidence for an inter-organizational research workflow. Around a claimed
    pilot or preparation result, this check looks for a named cross-boundary
    actor, a transferred artifact, a bounded task, an observed outcome, an
    independent action, and an explicit statement of what remains unproven.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    low = text.lower()
    prior_terms = [
        "予備実証",
        "予備試験",
        "予備成果",
        "準備状況",
        "既往実績",
        "実証済",
        "再実行した",
        "統合済",
    ]
    prior_hits = _contains_any(low, prior_terms)
    if not prior_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "prior_hits": [],
            "comments": [],
            "target": "a claimed collaboration pilot must show an artifact crossing an organization boundary and producing an observed result",
            "source": "generic cross-organization preliminary-evidence check",
        }

    sentences = [
        sentence.strip()
        for sentence in re.split(r"(?<=[。．!?！？])|\n+", text)
        if sentence.strip()
    ]
    marked = [
        i
        for i, sentence in enumerate(sentences)
        if _contains_any(sentence.lower(), prior_terms)
    ]
    selected: set[int] = set()
    for index in marked:
        selected.update(range(max(0, index - 1), min(len(sentences), index + 4)))
    evidence_text = " ".join(sentences[i] for i in sorted(selected)) or text
    evidence_low = evidence_text.lower()

    axes = {
        "cross_organization_actor": [
            "大学間",
            "他機関",
            "別機関",
            "複数機関",
            "三機関",
            "四機関",
            "共同研究先",
            "開発元",
            "提供者",
        ],
        "transferred_artifact": [
            "コード",
            "実装",
            "モデル",
            "入力データ",
            "行列",
            "試験",
            "ベンチマーク",
            "変更",
            "pull request",
            "プルリクエスト",
        ],
        "bounded_task": [
            "同一問題",
            "再実行",
            "比較した",
            "統合した",
            "差し替",
            "取り込",
            "再現した",
            "検証した",
        ],
        "observed_outcome": [
            "収束",
            "反復",
            "残差",
            "誤差",
            "一致",
            "不一致",
            "採択",
            "棄却",
            "満たさなかった",
            "取り込まれ",
        ],
        "independent_action": [
            "別機関レビュー",
            "独立レビュー",
            "第三者",
            "非実装担当",
            "別機関が",
            "他機関が",
            "提供の",
            "提供した",
        ],
        "remaining_gap": [
            "未達",
            "未検証",
            "未実施",
            "主張しない",
            "に留まる",
            "残る",
            "出発点",
            "一方",
            "本研究で",
        ],
    }
    axis_results = {}
    for name, keywords in axes.items():
        hits = _contains_any(evidence_low, keywords)
        axis_results[name] = {
            "ok": bool(hits),
            "matches": hits[:8],
            "keywords": keywords,
        }

    missing = [name for name, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "cross_organization_actor": "誰から誰へ資産が渡ったか、組織境界を明示する。",
        "transferred_artifact": "実際に渡したコード、モデル、データ、試験等を特定する。",
        "bounded_task": "同一問題での再実行、比較、変更等の限定した作業を示す。",
        "observed_outcome": "残差、誤差、採否、失敗等の観測結果を示す。",
        "independent_action": "資産作成者とは別の機関が行った実行・変更・レビューを示す。",
        "remaining_gap": "予備実証でまだ証明していない範囲を明記する。",
    }
    comments = [comments_by_axis[name] for name in missing]
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    reviewed = bool(re.search(
        r"(?:別機関|他機関|第三者|非実装担当)[^。．]{0,50}"
        r"(?:レビューした|査読した|採択した|棄却した|変更した)",
        evidence_text,
    ))
    if axis_results["observed_outcome"]["ok"] and reviewed:
        evidence_level = "L3: independently reviewed or adopted scientific result"
    elif axis_results["observed_outcome"]["ok"]:
        evidence_level = "L2: cross-organization re-execution with an observed outcome"
    elif axis_results["transferred_artifact"]["ok"]:
        evidence_level = "L1: artifact transfer or build workflow only"
    else:
        evidence_level = "L0: publication, link, or internal-use claim only"
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "prior_hits": prior_hits[:12],
        "evidence_excerpt": evidence_text[:1000],
        "evidence_level": evidence_level,
        "independent_review_or_adoption": reviewed,
        "comments": comments,
        "target": (
            "name the cross-organization actor and artifact, show a bounded "
            "independent task and observed result, and state what remains unproven"
        ),
        "source": "generic cross-organization preliminary-evidence check",
    }


def grant_writing_named_software_abstraction_check(
    text: str,
    software_names: str = "",
) -> dict:
    """Check whether named software is used at the right proposal level.

    Titles, summaries, academic questions, aims, novelty, and impact should
    normally name the technical category (for example, OSS, a high-order
    finite-element OSS, or an OSS coupling platform), not make one software
    package the research concept.  Named software remains appropriate when it
    makes methods, preliminary evidence, collaboration, rights, or costs
    reproducible.

    Args:
        text: Proposal text or an existing .md/.tex/.txt path.
        software_names: Optional comma-separated names added to the built-in
            software list.
    """
    text = _read_text_if_path(text)
    default_names = [
        "NGSolve",
        "ONELAB",
        "openCFS",
        "FreeFEM++",
        "Gmsh",
        "GetDP",
        "preCICE",
        "OpenMDAO",
        "COMSOL",
        "JMAG",
        "ANSYS",
        "Radia",
        "MATLAB",
        "Simulink",
    ]
    extra_names = [
        name.strip() for name in software_names.split(",") if name.strip()
    ]
    names = list(dict.fromkeys(default_names + extra_names))

    core_markers = [
        "概要",
        "研究背景",
        "学術的問い",
        "中心の問い",
        "研究目的",
        "本研究の目的",
        "研究課題",
        "主題",
        "独創性",
        "新規性",
        "波及効果",
        "国際性",
        "研究動向",
        "位置付け",
        "意義",
        "タイトル",
    ]
    implementation_markers = [
        "研究方法",
        "達成指標",
        "実証",
        "検証",
        "年度計画",
        "研究計画",
        "準備状況",
        "遂行能力",
        "研究実績",
        "既往実績",
        "予備成果",
        "予備試験",
        "共同実績",
        "共同研究",
        "発表",
        "掲載",
        "リンク",
        "権利",
        "ライセンス",
        "予算",
        "経費",
        "費用",
        "リスク",
        "代替策",
        "ベンチマーク",
        "比較する",
        "用いる",
        "使用する",
        "実装する",
        "接続した",
        "統合した",
        "経験した",
        "共著",
    ]
    heading_pattern = re.compile(
        r"\\(?:sub)*section\*?\{([^{}]+)\}|^\s*#{1,6}\s+(.+?)\s*$"
    )

    def marker_hits(fragment: str, markers: list[str]) -> list[str]:
        return [marker for marker in markers if marker in fragment]

    current_section = ""
    risks = []
    allowed_mentions = []
    unclassified_mentions = []
    for line_number, line in enumerate(text.splitlines(), 1):
        heading_match = heading_pattern.search(line)
        if heading_match:
            current_section = heading_match.group(1) or heading_match.group(2) or ""

        for name in names:
            name_pattern = re.compile(
                rf"(?<![A-Za-z0-9]){re.escape(name)}(?![A-Za-z0-9])",
                re.IGNORECASE,
            )
            for match in name_pattern.finditer(line):
                sentence_start = line.rfind("。", 0, match.start()) + 1
                sentence_stop = line.find("。", match.end())
                if sentence_stop < 0:
                    sentence_stop = len(line)
                else:
                    sentence_stop += 1
                excerpt = line[sentence_start:sentence_stop].strip()
                local_core = marker_hits(excerpt, core_markers)
                section_core = marker_hits(current_section, core_markers)
                local_implementation = marker_hits(excerpt, implementation_markers)
                section_implementation = marker_hits(
                    current_section, implementation_markers
                )
                mention = {
                    "software": name,
                    "line": line_number,
                    "section": current_section,
                    "excerpt": excerpt[:320],
                }
                if local_implementation or section_implementation:
                    mention["allowed_by"] = (
                        local_implementation or section_implementation
                    )
                    allowed_mentions.append(mention)
                elif local_core or section_core:
                    mention.update({
                        "severity": "MEDIUM",
                        "framing_hits": local_core or section_core,
                        "comment": (
                            f"{name} is a named implementation, not the research "
                            "concept; state the OSS or technical category here."
                        ),
                    })
                    risks.append(mention)
                else:
                    unclassified_mentions.append(mention)

    applicable = bool(risks or allowed_mentions or unclassified_mentions)
    score = None if not applicable else max(0.0, 10.0 - 2.5 * len(risks))
    risky_names = sorted({risk["software"] for risk in risks})
    comments = [
        "問い・背景・目的・独創性・波及効果では、個別ソフト名 "
        + "、".join(risky_names)
        + " をOSSまたは技術カテゴリへ抽象化する。"
    ] if risky_names else []
    recommendations = [
        "中心概念は「OSS」「高次有限要素OSS」「既存OSS連成基盤」等の役割・技術カテゴリで記述する。",
        "固有ソフト名は、研究方法、予備成果、共同実績、権利・ライセンス、費用根拠で再現性を担保するときに限る。",
        "固有名詞を残す場合は、先に一般カテゴリと研究上の役割を示し、その一実装として名称を置く。",
    ]
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "allowed_mentions": allowed_mentions,
        "unclassified_mentions": unclassified_mentions,
        "comments": comments,
        "recommendations": recommendations if risks else [],
        "target": (
            "category-level framing in the title, summary, question, aims, "
            "novelty, and impact; named software only where it supplies "
            "reproducible implementation or feasibility evidence"
        ),
        "source": "generic named-software abstraction-level check",
    }


def grant_writing_reviewer_vocabulary_check(text: str) -> dict:
    """Check whether proposal vocabulary is accessible to the likely reviewer.

    A domain reviewer can be expected to know the field, but not every software
    acronym, foreign-university abbreviation, laboratory shorthand, or named
    benchmark. This check also keeps benchmarks in a verification role.
    """
    text = _read_text_if_path(text)
    risks: list[dict] = []
    term_results: dict[str, dict] = {}

    def context_for(match: re.Match, radius: int = 180) -> str:
        start = max(0, match.start() - radius)
        stop = min(len(text), match.end() + radius)
        return re.sub(r"\s+", " ", text[start:stop]).strip()

    def add_risk(
        risk_type: str,
        term: str,
        match: re.Match,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
    ) -> None:
        risks.append({
            "type": risk_type,
            "term": term,
            "line": text.count("\n", 0, match.start()) + 1,
            "severity": severity,
            "excerpt": context_for(match)[:360],
            "comment": comment,
            "recommendation": recommendation,
        })

    acronym_specs = {
        "OSS": {
            "pattern": r"(?<![A-Za-z0-9_])OSS(?![A-Za-z0-9_])",
            "explanations": [
                "オープンソースソフトウェア",
                "オープンソースソフトウエア",
            ],
            "recommendation": "初出を「オープンソースソフトウェア（OSS）」とする。",
        },
        "LLM": {
            "pattern": r"(?<![A-Za-z0-9_])LLM(?![A-Za-z0-9_])",
            "explanations": ["大規模言語モデル"],
            "recommendation": "略語を避けて「大規模言語モデル」と記す。",
        },
        "MCP": {
            "pattern": r"(?<![A-Za-z0-9_])MCP(?![A-Za-z0-9_])",
            "explanations": [
                "知識・実行インターフェース",
                "知識・実行インタフェース",
                "知識・実行窓口",
                "知識と実行のインターフェース",
                "知識と実行のインタフェース",
            ],
            "recommendation": (
                "初出で「AIが利用する知識・実行インターフェース」等の"
                "研究上の役割を日本語で説明する。英語名の展開だけにしない。"
            ),
        },
    }
    for term, spec in acronym_specs.items():
        matches = list(re.finditer(spec["pattern"], text))
        if not matches:
            term_results[term] = {"present": False, "count": 0, "explained": True}
            continue
        first = matches[0]
        local = context_for(first)
        explained = any(word in local for word in spec["explanations"])
        term_results[term] = {
            "present": True,
            "count": len(matches),
            "explained": explained,
            "first_line": text.count("\n", 0, first.start()) + 1,
        }
        if not explained:
            add_risk(
                "unexplained_acronym",
                term,
                first,
                f"{term}の初出に、審査者が読める日本語の意味・役割がない。",
                spec["recommendation"],
            )

    mcp_storage_pattern = re.compile(
        r"(?<![A-Za-z0-9_])MCP(?![A-Za-z0-9_])(?:そのもの)?(?:に|へ)"
        r"[^。！？\n]{0,80}(?:蓄積|保存|格納|収録)(?:する|し|した|している)?"
    )
    for match in mcp_storage_pattern.finditer(text):
        add_risk(
            "mcp_described_as_storage",
            "MCP",
            match,
            (
                "MCPを情報の保存場所として説明している。MCPは知識や実行機能を"
                "AIへ提示する規約・インターフェースである。"
            ),
            (
                "保存先をリポジトリ、文書又はデータベースとして明示し、"
                "『実装判断・検証手順をMCPサーバーから利用可能にする』等と役割を分ける。"
            ),
            severity="HIGH",
        )

    institution_aliases = {
        "TU Graz": "グラーツ工科大学",
        "TU Wien": "ウィーン工科大学",
    }
    for alias, japanese_name in institution_aliases.items():
        matches = list(re.finditer(re.escape(alias), text, flags=re.IGNORECASE))
        term_results[alias] = {
            "present": bool(matches),
            "count": len(matches),
            "preferred": japanese_name,
        }
        if matches:
            add_risk(
                "foreign_institution_alias",
                alias,
                matches[0],
                f"{alias}は分野外の審査者には大学名として即読できない。",
                f"本文では「{japanese_name}」と記し、必要な箇所だけ原語を併記する。",
            )

    shorthand_specs = {
        "MMM": {
            "pattern": r"(?<![A-Za-z0-9_])MMM(?![A-Za-z0-9_])",
            "preferred": "磁気モーメント法",
            "comment": "研究室内略称のMMMは電磁気分野でも共有を仮定できない。",
        },
        "H(curl)": {
            "pattern": (
                r"H\s*\(\s*(?:\\(?:mathrm|operatorname)\s*\{\s*)?"
                r"curl\s*\}?\s*\)"
            ),
            "preferred": "辺要素",
            "comment": "関数空間記号より、審査者が設計手法として読める用語を先に置く。",
        },
    }
    for term, spec in shorthand_specs.items():
        matches = list(re.finditer(spec["pattern"], text, flags=re.IGNORECASE))
        term_results[term] = {
            "present": bool(matches),
            "count": len(matches),
            "preferred": spec["preferred"],
        }
        if matches:
            add_risk(
                "domain_shorthand",
                term,
                matches[0],
                spec["comment"],
                (
                    f"概要・目的・意義では「{spec['preferred']}」を用い、"
                    "数理記号は必要なら方法節で補足する。"
                ),
            )

    benchmark_pattern = r"TEAM\s*(?:Problem\s*)?28"
    benchmark_matches = list(re.finditer(benchmark_pattern, text, flags=re.IGNORECASE))
    benchmark_result = {
        "present": bool(benchmark_matches),
        "count": len(benchmark_matches),
        "framed_as_verification": True,
        "limitation_stated": True,
        "engineering_value_distinguished": True,
        "negative_disclaimer_present": False,
        "used_as_significance": False,
    }
    if benchmark_matches:
        benchmark_context = " ".join(
            context_for(match, radius=260) for match in benchmark_matches
        )
        verification_terms = [
            "単純化",
            "公開基準問題",
            "初期検証",
            "整合性",
            "校正",
            "基礎検証",
        ]
        limitation_terms = [
            "工学的有用性の根拠には用いない",
            "工学的有用性の根拠にしない",
            "工学的有用性を示すものではない",
            "有用性の根拠には用いない",
            "有用性の根拠にしない",
            "有用性を示すものではない",
        ]
        real_design_terms = [
            "主実証",
            "実設計",
            "実機寸法",
            "制約付き機器設計",
            "設計判断",
            "設計量",
            "加速器電磁石",
        ]
        significance_terms = [
            "中心の問い",
            "学術的問い",
            "本研究の目的",
            "独創性",
            "新規性",
            "波及効果",
            "主成果",
        ]
        framed = any(term in benchmark_context for term in verification_terms)
        negative_disclaimer = any(term in text for term in limitation_terms)
        positive_design_bridge = framed and any(
            term in benchmark_context for term in real_design_terms
        )
        limited = positive_design_bridge
        significance = any(
            term in context_for(match, radius=120)
            for match in benchmark_matches
            for term in significance_terms
        )
        benchmark_result.update({
            "framed_as_verification": framed,
            "limitation_stated": limited,
            "engineering_value_distinguished": limited,
            "negative_disclaimer_present": negative_disclaimer,
            "used_as_significance": significance,
        })
        first = benchmark_matches[0]
        if not framed:
            add_risk(
                "benchmark_unframed",
                "TEAM Problem 28",
                first,
                "固有ベンチマークの研究計画上の役割が明示されていない。",
                "単純化した公開基準問題を初期検証または校正に用いる、と位置付ける。",
            )
        if not limited:
            add_risk(
                "benchmark_without_limit",
                "TEAM Problem 28",
                first,
                "基準問題への一致と工学的有用性が区別されていない。",
                (
                    "基準問題の役割を初期検証とし、工学的価値を示す制約付き"
                    "実設計課題と設計判断を続けて記す。予備成果を否定する"
                    "免責文は置かない。"
                ),
                severity="HIGH",
            )
        if negative_disclaimer:
            disclaimer_match = next(
                match
                for term in limitation_terms
                if (match := re.search(re.escape(term), text)) is not None
            )
            add_risk(
                "benchmark_self_negating_disclaimer",
                "TEAM Problem 28",
                disclaimer_match,
                "予備検証の成果を、直後の免責文が自ら打ち消している。",
                (
                    "基準問題で何を確認できたかを肯定形で述べ、続けて主実証の"
                    "制約付き実設計と設計判断を示す。"
                ),
            )
        if significance:
            add_risk(
                "benchmark_as_significance",
                "TEAM Problem 28",
                first,
                "固有ベンチマークが問い・目的・独創性・波及効果の代わりになっている。",
                "固有名は方法・予備結果へ下げ、工学的意義を実設計量と意思決定で示す。",
                severity="HIGH",
            )
        if len(benchmark_matches) > 2:
            add_risk(
                "benchmark_overexposure",
                "TEAM Problem 28",
                benchmark_matches[2],
                "固有ベンチマーク名の反復が、研究の主役であるように見せている。",
                "固有名は方法節で一度定義し、以後は「公開基準問題」等と記す。",
            )

    term_results["TEAM Problem 28"] = benchmark_result
    applicable = any(result.get("present") for result in term_results.values())
    deductions = sum(2.0 if risk["severity"] == "HIGH" else 1.0 for risk in risks)
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(
        dict.fromkeys(risk["recommendation"] for risk in risks)
    )
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "term_results": term_results,
        "benchmark_result": benchmark_result,
        "target": (
            "assume field expertise, not familiarity with every software/AI acronym, "
            "foreign-institution alias, laboratory shorthand, or named benchmark; "
            "use benchmarks for verification and real design tasks for engineering value"
        ),
        "source": "generic reviewer-vocabulary and benchmark-role check",
    }


def grant_writing_persuasion_quality_check(text: str) -> dict:
    """Check reviewer-facing hierarchy, equations, and defensive prose.

    The check targets a common failure mode in technically careful proposals:
    evidence is immediately negated, equations arrive before their symbols and
    decision role, or exceptions crowd out the main claim. It is a heuristic
    quality gate, not a mathematical parser.
    """
    text = _read_text_if_path(text)
    risks: list[dict] = []

    def add_risk(
        risk_type: str,
        start: int,
        excerpt: str,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
        **details,
    ) -> None:
        item = {
            "type": risk_type,
            "line": text.count("\n", 0, start) + 1,
            "severity": severity,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:360],
            "comment": comment,
            "recommendation": recommendation,
        }
        item.update(details)
        risks.append(item)

    self_negating_patterns = [
        re.compile(
            r"(?:この|その)?(?:一致|結果|成果|実績|予備実証)"
            r"[^。！？\n]{0,60}(?:根拠|有用性|意義|価値)"
            r"[^。！？\n]{0,35}(?:しない|用いない|ではない|示さない)"
        ),
        re.compile(
            r"(?:工学的|学術的)?(?:有用性|意義|価値)[^。！？\n]{0,40}"
            r"(?:示すものではない|根拠(?:に|と)(?:は)?"
            r"(?:しない|ならない|用いない))"
        ),
        re.compile(
            r"(?:結果|成果|実績)[^。！？\n]{0,50}"
            r"(?:だけを示す|にすぎない)"
        ),
    ]
    seen_self_negating: set[tuple[int, int]] = set()
    for pattern in self_negating_patterns:
        for match in pattern.finditer(text):
            span = match.span()
            if span in seen_self_negating:
                continue
            seen_self_negating.add(span)
            add_risk(
                "self_negating_evidence",
                match.start(),
                match.group(0),
                "提示した成果を同じ段落の免責文が打ち消している。",
                (
                    "成果が確認した範囲を肯定形で述べ、次の文で主実証が検証する"
                    "設計量・判断・適用範囲へ接続する。"
                ),
                severity="HIGH",
            )

    equation_pattern = re.compile(
        r"(?<!\\)\\\[(?P<bracket>.*?)(?<!\\)\\\]"
        r"|\\begin\{(?P<env>equation\*?|align\*?|gather\*?)\}"
        r"(?P<environment>.*?)\\end\{(?P=env)\}"
        r"|\$\$(?P<dollar>.*?)\$\$",
        flags=re.DOTALL,
    )
    symbol_pattern = re.compile(
        r"(?:\\[A-Za-z]+|[A-Za-z])(?:_\{[^{}\n]+\}|_[A-Za-z0-9]+)"
    )
    on_ramp_terms = [
        "定義",
        "表す",
        "求める",
        "算出",
        "評価",
        "判断",
        "区間",
        "指標",
        "関係",
        "変化",
        "比較",
        "幅",
    ]
    interpretation_terms = [
        "と表す",
        "意味",
        "示す",
        "判断",
        "順位",
        "区間",
        "係数",
        "用いる",
        "分離",
        "重な",
    ]

    def compact_latex(value: str) -> str:
        value = re.sub(r"\\(?:left|right|,|;|!|quad|qquad)", "", value)
        return re.sub(r"[\s{}$]", "", value)

    equation_count = 0
    for match in equation_pattern.finditer(text):
        equation_count += 1
        body = next(
            group
            for group in (
                match.group("bracket"),
                match.group("environment"),
                match.group("dollar"),
            )
            if group is not None
        )
        before = text[max(0, match.start() - 500):match.start()]
        after = text[match.end():min(len(text), match.end() + 500)]
        before_prose = _prose_for_lint(before)
        after_prose = _prose_for_lint(after)
        equation_excerpt = match.group(0)

        if not any(term in before_prose for term in on_ramp_terms):
            add_risk(
                "equation_without_on_ramp",
                match.start(),
                equation_excerpt,
                "数式の前に、何を判断するための式かが示されていない。",
                (
                    "式の前に対象量、比較目的、式が答える判断を一文で述べてから"
                    "数式を置く。"
                ),
                severity="HIGH",
            )

        symbols = sorted(set(symbol_pattern.findall(body)))
        surrounding = compact_latex(before + after)
        missing_symbols = [
            symbol for symbol in symbols if compact_latex(symbol) not in surrounding
        ]
        if missing_symbols:
            add_risk(
                "equation_symbols_not_introduced",
                match.start(),
                equation_excerpt,
                "数式中の記号が周辺の文章で導入されていない。",
                (
                    "各記号を物理的意味と計算操作で先に定義する。添字だけで"
                    "解析条件を推測させない。"
                ),
                severity="HIGH",
                missing_symbols=missing_symbols,
            )

        if not any(term in after_prose for term in interpretation_terms):
            add_risk(
                "equation_without_interpretation",
                match.start(),
                equation_excerpt,
                "数式の後に、値をどう研究判断へ使うかが説明されていない。",
                (
                    "式の直後に、大小・分離・閾値がどの設計判断を変えるかを"
                    "平文で述べる。"
                ),
            )

    inline_math_pattern = re.compile(
        r"(?<!\$)\$(?!\$)(?P<body>[^$\n]{1,120})(?<!\$)\$(?!\$)"
    )
    inline_condition_pattern = re.compile(
        r"^\s*(?P<symbol>(?:\\[A-Za-z]+|[A-Za-z])"
        r"(?:_\{[^{}\n]+\}|_[A-Za-z0-9]+)?)\s*"
        r"(?:=|\\leq?|\\geq?|<|>)\s*"
        r"(?:[-+]?(?:\d+(?:\.\d*)?|\.\d+)|\\infty)\s*$"
    )
    unexplained_inline_condition_count = 0
    for match in inline_math_pattern.finditer(text):
        condition = inline_condition_pattern.match(match.group("body"))
        if condition is None:
            continue
        symbol = condition.group("symbol")
        base_match = re.match(r"(?:\\[A-Za-z]+|[A-Za-z])", symbol)
        if base_match is None:
            continue
        base_symbol = base_match.group(0)
        before = text[max(0, match.start() - 400):match.start()]
        prior_symbol_pattern = re.compile(
            r"(?<!\$)\$(?!\$)[^$\n]{0,80}"
            + re.escape(base_symbol)
            + r"[^$\n]{0,80}(?<!\$)\$(?!\$)"
        )
        if prior_symbol_pattern.search(before):
            continue
        unexplained_inline_condition_count += 1
        add_risk(
            "inline_condition_before_physical_meaning",
            match.start(),
            match.group(0),
            "記号だけの条件式が、その記号の工学的意味より先に現れている。",
            (
                "現象と設計上の効果を平文で述べ、量の名称と記号を独立に"
                "定義してから条件式を置く。概要・図・成果名には工学的効果を使う。"
            ),
            severity="HIGH",
            symbol=symbol,
        )

    paragraph_negative_terms = [
        "必達成果には含めない",
        "目的としない",
        "前提にしない",
        "依存しない",
        "対象外",
        "今後の課題",
        "失敗すれば",
        "含めない",
        "ではない",
        "しない",
        "不能",
        "未達",
    ]
    defensive_exemptions = [
        "リスク",
        "研究遂行上の課題",
        "代替策",
        "失敗時",
        "不成立時",
        "倫理",
        "法令",
        "安全",
        "個人情報",
        "知的財産",
        "権利",
    ]
    negative_pattern = re.compile(
        "|".join(re.escape(term) for term in paragraph_negative_terms)
    )
    cursor = 0
    defensive_paragraph_count = 0
    for raw_paragraph in re.split(r"\n\s*\n+", text):
        start = text.find(raw_paragraph, cursor)
        cursor = max(cursor, start + len(raw_paragraph))
        paragraph = _prose_for_lint(raw_paragraph)
        if len(paragraph) < 35 or any(
            term in paragraph for term in defensive_exemptions
        ):
            continue
        negative_hits = negative_pattern.findall(paragraph)
        if len(negative_hits) >= 3:
            defensive_paragraph_count += 1
            add_risk(
                "defensive_paragraph",
                max(0, start),
                paragraph,
                "否定・除外・未達の説明が一段落に集中し、主張が見えにくい。",
                (
                    "本文では成立させる主張と判定法を先に書く。例外と代替策は"
                    "研究遂行上のリスク段落へまとめる。"
                ),
                negative_hits=negative_hits,
            )

    optional_pattern = re.compile(
        r"条件付き追加検証|必達成果には含めない|余力があれば|条件が整えば|"
        r"可能であれば[^。！？\n]{0,40}(?:追加|検証|実施)"
    )
    optional_branch_count = 0
    for match in optional_pattern.finditer(text):
        optional_branch_count += 1
        add_risk(
            "optional_branch_in_core_plan",
            match.start(),
            match.group(0),
            "任意の追加課題が中心計画に入り、必達成果の輪郭を弱めている。",
            (
                "中核計画は必達の問い・実証・成果に絞り、条件付き拡張は"
                "リスク対応または将来展開へ移す。"
            ),
        )

    memo_specs = {
        "lettered_experiment": re.compile(r"実証[AB](?![A-Za-z0-9])"),
        "letter_pair": re.compile(
            r"(?<![A-Za-z0-9])A(?:・|/)B"
            r"(?=(?:で|の|を|は|へ|と|仕様|論文|成果|実証))"
        ),
        "coded_stage": re.compile(r"(?<![A-Za-z0-9])L[1-4](?![A-Za-z0-9])"),
        "parenthetical_milestone": re.compile(
            r"[（(]\s*(?:年度末|到達点|ゴール)\s*[:：]"
        ),
        "draft_marker": re.compile(
            r"(?i:\b(?:TODO|TBD|FIXME)\b)|要確認|暫定案|仮置き"
        ),
    }
    memo_shorthand_count = 0
    for memo_type, pattern in memo_specs.items():
        matches = list(pattern.finditer(text))
        if not matches:
            continue
        memo_shorthand_count += len(matches)
        first = matches[0]
        add_risk(
            "internal_memo_shorthand",
            first.start(),
            first.group(0),
            "内部管理用の略号・進捗メモが申請本文に残っている。",
            (
                "実証は科学的内容で呼び、年度到達点は括弧メモでなく"
                "完結した文で記す。下書き状態の印は本文から除く。"
            ),
            memo_type=memo_type,
            occurrence_count=len(matches),
        )

    prose = _prose_for_lint(text)
    acronym_pattern = re.compile(
        r"(?<![A-Za-z0-9])(?:[A-Z]{2,}[A-Za-z0-9]*"
        r"(?:[-/][A-Za-z0-9]+)*|H\([a-z]+\))(?![A-Za-z0-9])"
    )
    acronym_pile_count = 0
    # An inventory is not a sentence that hides its meaning behind acronyms.
    # 「Adventure, CST Studio, ELF/Magic, Elmer, EMCoS, EMSolution, ...」 in an
    # adopted proposal's 研究環境 is a list of the software the lab owns, and
    # naming them is the whole point.
    for sentence in re.split(r"(?<=[。．!?！？])|\n", prose):
        acronyms = sorted(set(acronym_pattern.findall(sentence)))
        if len(acronyms) < 6:
            continue
        if sentence.count(",") + sentence.count("、") >= 6:
            continue
        acronym_pile_count += 1
        start = text.find(sentence[:24])
        add_risk(
            "acronym_pile",
            max(0, start),
            sentence,
            "一文に固有略語が集中し、研究上の役割より名称の列挙が先に見える。",
            (
                "まず手法カテゴリ、比較軸、変わる設計判断を述べ、固有名は"
                "方法・実行可能性の箇所に分けて置く。"
            ),
            acronyms=acronyms,
        )

    deductions = sum(2.0 if risk["severity"] == "HIGH" else 1.0 for risk in risks)
    applicable = bool(text.strip())
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(
        dict.fromkeys(risk["recommendation"] for risk in risks)
    )
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "counts": {
            "display_equations": equation_count,
            "self_negating_evidence": len(seen_self_negating),
            "defensive_paragraphs": defensive_paragraph_count,
            "optional_branches": optional_branch_count,
            "memo_shorthand": memo_shorthand_count,
            "acronym_piles": acronym_pile_count,
            "unexplained_inline_conditions": unexplained_inline_condition_count,
        },
        "target": (
            "lead with a positive claim and reviewer decision; introduce every "
            "equation and inline condition by its physical meaning and symbols, "
            "interpret it immediately, collect exceptions in one risk section, "
            "and keep optional branches out of the core plan"
        ),
        "source": "generic reviewer-facing persuasion and equation-on-ramp check",
    }


_CLAIM_MARKERS = (
    "本研究は次を問う",
    "を問う",
    "中心の問い",
    "学術的な問い",
    "学術的問い",
    "研究目的は",
    "本研究の目的",
    "目的は",
)

# The noun that names what the answer will BE. Parallel statements of one
# question should use the same noun for the same semantic role. Different
# nouns may coexist when the prose distinguishes, for example, a scientific
# condition, an operational criterion, and an application limit.
_CLAIM_OUTCOME_NOUNS = (
    "条件", "境界", "範囲", "領域", "基準", "指標", "手順", "限界", "閾値",
    "選択則", "設計則", "法則",
)

# The operation the researcher performs on it.
_CLAIM_OPERATION_NOUNS = (
    "定量化", "記述", "検証", "同定", "決定", "確定", "評価", "予測", "測定",
)

_CLAIM_TERM = re.compile(r"[一-龥]{2,}|[ァ-ヴー]{3,}")
_CLAIM_TERM_STOPWORDS = frozenset({
    "本研究", "中心", "場合", "以下", "以上", "今回", "一方", "同様", "本文",
    "概要", "研究",
})


def _claim_terms(fragment: str) -> set[str]:
    """Technical-noun proxy: kanji runs and katakana runs, minus filler.

    Japanese has no spaces, so a naive run of any kana/kanji swallows whole
    clauses. Kanji and katakana runs approximate the content nouns well
    enough to compare two statements of the same claim.
    """
    return {
        term
        for term in _CLAIM_TERM.findall(fragment)
        if term not in _CLAIM_TERM_STOPWORDS
    }


def _claim_statements(text: str) -> list[dict]:
    """Locate statements of the central question / aim.

    A claim often spans several sentences (「中心の問いは次である。…。…か。」),
    so each marker takes its own sentence plus the following ones up to a
    sentence that closes the claim.
    """
    sentences = [s for s in re.split(r"(?<=[。．!?！？])", text) if s.strip()]

    def marker_of(fragment: str) -> str | None:
        return next((m for m in _CLAIM_MARKERS if m in fragment), None)

    # An opener defers the claim to what follows (「中心の問いは次である。」).
    opener = re.compile(r"(?:次である|次を問う|次のとおり|以下である)[。．]\s*$")
    closer = re.compile(r"(?:か|である|ことである|問う|明らかにする)[。．]\s*$")

    statements: list[dict] = []
    consumed: set[int] = set()
    for index, sentence in enumerate(sentences):
        if index in consumed:
            continue
        marker = marker_of(sentence)
        if marker is None:
            continue
        chunk = [sentence]
        if opener.search(sentence.strip()) or not closer.search(sentence.strip()):
            for offset in range(1, 4):
                nxt = index + offset
                if nxt >= len(sentences):
                    break
                # Never swallow the next claim: it is a separate statement.
                if marker_of(sentences[nxt]):
                    break
                chunk.append(sentences[nxt])
                consumed.add(nxt)
                if closer.search(sentences[nxt].strip()):
                    break
        body = "".join(chunk).strip()
        terms = _claim_terms(body)
        # A heading fragment or a passing mention is not a claim statement.
        if len(terms) < 4 or not closer.search(body):
            continue
        statements.append({
            "marker": marker,
            "sentence_index": index + 1,
            "text": body,
            "terms": terms,
        })
    return statements


def grant_writing_central_claim_consistency_check(text: str) -> dict:
    """Check that one central claim is not stated as two different claims.

    A proposal commonly states its question in the summary and again in the
    body. When parallel statements use different decisive nouns for what
    appears to be the same role -- one promises a 「境界」, the other a
    「条件」 -- a reviewer cannot tell whether the proposal has one question
    or two. This is not a universal ban on either noun: distinct roles may
    legitimately use distinct terms. Keyword-coverage checks score the
    ambiguous draft perfectly because every required word is present
    somewhere; the defect is that the words disagree with each other.

    The check is optional: it needs at least two claim statements.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    statements = _claim_statements(text)
    if len(statements) < 2:
        return {
            "applicable": False,
            "score": None,
            "statement_count": len(statements),
            "statements": [
                {k: v for k, v in s.items() if k != "terms"} for s in statements
            ],
            "risks": [],
            "comments": [],
            "target": (
                "define one central-question semantic contract; parallel "
                "restatements must preserve its decisive roles"
            ),
            "source": "central-claim consistency check",
        }

    risks: list[dict] = []
    pairs: list[dict] = []
    for i in range(len(statements)):
        for j in range(i + 1, len(statements)):
            a, b = statements[i], statements[j]
            shared = a["terms"] & b["terms"]
            union = a["terms"] | b["terms"]
            similarity = round(len(shared) / max(1, len(union)), 2)
            a_outcomes = {n for n in _CLAIM_OUTCOME_NOUNS if n in a["text"]}
            b_outcomes = {n for n in _CLAIM_OUTCOME_NOUNS if n in b["text"]}
            a_ops = {n for n in _CLAIM_OPERATION_NOUNS if n in a["text"]}
            b_ops = {n for n in _CLAIM_OPERATION_NOUNS if n in b["text"]}
            pair = {
                "statements": [a["sentence_index"], b["sentence_index"]],
                "markers": [a["marker"], b["marker"]],
                "similarity": similarity,
                "shared_terms": sorted(shared),
                "only_in_first": sorted(a["terms"] - b["terms"]),
                "only_in_second": sorted(b["terms"] - a["terms"]),
                "outcome_nouns": [sorted(a_outcomes), sorted(b_outcomes)],
                "operation_nouns": [sorted(a_ops), sorted(b_ops)],
            }
            pairs.append(pair)

            # Same topic (they share anchors) but the promised answer differs.
            if len(shared) >= 2 and a_outcomes and b_outcomes and not (
                a_outcomes & b_outcomes
            ):
                risks.append({
                    "type": "outcome_noun_divergence",
                    "severity": "HIGH",
                    "statements": pair["statements"],
                    "comment": (
                        "同じ問いの言い直しだが、答えの形を表す名詞が異なる: "
                        + "／".join(sorted(a_outcomes))
                        + " と "
                        + "／".join(sorted(b_outcomes))
                    ),
                    "recommendation": (
                        "同じ答えの役割なら同じ名詞を使う。別の役割を意図する"
                        "場合は、条件、判断基準、適用限界等の関係を定義する。"
                    ),
                    "excerpts": [a["text"][:160], b["text"][:160]],
                })
            elif len(shared) >= 2 and a_ops and b_ops and not (a_ops & b_ops):
                risks.append({
                    "type": "operation_noun_divergence",
                    "severity": "MEDIUM",
                    "statements": pair["statements"],
                    "comment": (
                        "問いに対して行う操作の語が異なる: "
                        + "／".join(sorted(a_ops))
                        + " と "
                        + "／".join(sorted(b_ops))
                    ),
                    "recommendation": (
                        "中心となる操作を固定し、記述や検証を併記する場合は"
                        "主操作との関係を示す。"
                    ),
                    "excerpts": [a["text"][:160], b["text"][:160]],
                })

            if similarity >= 0.8:
                risks.append({
                    "type": "verbatim_restatement",
                    "severity": "LOW",
                    "statements": pair["statements"],
                    "comment": "概要と本文がほぼ同一文になっている。",
                    "recommendation": (
                        "概要は全体の要約、本文は問いを導く論証と定義、と"
                        "役割を分ける。"
                    ),
                    "excerpts": [a["text"][:160], b["text"][:160]],
                })

    deductions = sum(
        3.0 if r["severity"] == "HIGH" else 1.5 if r["severity"] == "MEDIUM" else 0.5
        for r in risks
    )
    score = max(0.0, round(10.0 - deductions, 1))
    return {
        "applicable": True,
        "score": score,
        "statement_count": len(statements),
        "statements": [
            {k: v for k, v in s.items() if k != "terms"} for s in statements
        ],
        "pairs": pairs,
        "risk_count": len(risks),
        "risks": risks,
        "comments": list(dict.fromkeys(r["comment"] for r in risks)),
        "recommendations": list(dict.fromkeys(r["recommendation"] for r in risks)),
        "target": (
            "one central-question semantic contract, preserved across "
            "summary and body"
        ),
        "source": "central-claim consistency check",
    }


# Verbs that promise a result without naming the operation that produces it.
# Derived from an experienced PI's rewrite of a 2026 Power Academy proposal
# (2026-08-18): 「統合し」 went 3 -> 0 while 「双方向に連成」 went 0 -> 5 and
# 「連成」 2 -> 8. The reviewer question these verbs leave unanswered is
# always the same: integrated *how*?
_VAGUE_CLAIM_VERBS = (
    "統合する", "統合し", "統合を行う",
    "連携する", "連携し",
    "活用する", "活用し",
    "融合する", "融合し",
    "高度化する", "高度化し",
    "推進する", "推進し",
    "橋渡しする", "橋渡しし",
)

# Naming any of these makes the sentence answer "how".
_MECHANISM_MARKERS = (
    "連成", "接続", "射影", "写像", "写す", "縮約", "変換",
    "入力", "出力", "反映", "介して", "を通じて", "経由",
    "双方向", "電圧", "電流", "を渡す", "受け渡", "同じ",
    "同一", "共通", "基底", "境界条件", "パラメータ",
)


# The review criterion this check implements, quoted from a 2025 KAKENHI
# review disclosure: 「研究課題の核心をなす学術的『問い』は明確であり、学術的
# 独自性や創造性が認められるか」. Three of five reviewers marked that single
# item down, and it drove the academic-importance score to 1.60 against an
# adopted average of 2.83 -- the largest gap on the sheet.
_ORIGINALITY_MARKERS = (
    "独自", "独創", "新規", "初めて", "初の", "従来にない", "既存手法にない",
    "本研究に固有", "他に例をみない",
)

# A gap statement: prior work did X, but Y is not established. Without one,
# an originality word is an assertion rather than a position.
# 「従来」 must stand alone. The funded proposal writes 「従来提案された…は…
# 考慮できず」 and 「従来提案されている…をそのまま適用できない」, which a list
# of 従来手法／従来法／従来の misses entirely.
_PRIOR_WORK_MARKERS = (
    "既往研究", "既存研究", "先行研究", "従来", "これまで", "既報",
    "既存の", "現状の", "現行の",
)
# Measured against five real submitted proposals (one funded, four not).
# None of them phrase the limit as 「確立していない」. They negate a capability:
# 「考慮できず」「そのまま適用できない」「実現できなかった」「制限される」.
# The first vocabulary here was written from assumption and matched zero of
# the five, including the funded one.
_GAP_MARKERS = (
    # capability negation -- how real proposals actually write it
    "できない", "できず", "できなかった", "えない", "困難",
    "適用できない", "考慮できず", "制限される", "制約される",
    "限界がある", "対応できない", "十分に扱えない",
    # explicit absence
    "確立していない", "確立されていない", "体系化されていない",
    "明らかでない", "明らかにされていない", "得られていない",
    "十分でない", "十分ではない", "扱われていない",
    "できていない", "限られている", "残されている",
    "未解決", "至っていない", "難しい", "存在しない", "知られていない",
)


# A proposal can be thoroughly international without ever writing 「国際」;
# it names a region or a foreign institution instead. Triggering only on the
# explicit word would skip exactly the drafts that do this well.
# Words that assert an international dimension wherever they appear. The
# region names below are weaker: they also name a conference venue.
_INTERNATIONAL_STRONG_TRIGGERS = (
    "国際", "海外", "国外", "世界", "グローバル", "外国", "諸外国",
)
_INTERNATIONAL_TRIGGERS = (
    "国際", "海外", "国外", "世界", "グローバル", "外国", "諸外国",
    "欧州", "欧米", "米国", "アジア", "アメリカ", "ヨーロッパ",
    "ドイツ", "オーストリア", "フランス", "英国", "イタリア", "中国",
    "韓国", "台湾", "スイス", "オランダ", "北欧",
)

# What makes an international claim checkable rather than aspirational.
# An institution or a person one could collaborate with. A conference is not
# one of those: asking why COMPUMAG rather than a domestic substitute is not a
# question anybody can answer.
_NAMED_PARTNER = re.compile(
    r"[ァ-ヴー]{2,}(?:工科)?大学|"
    # 氏 must be the honorific, not the 氏 of 氏名, and it must sit right after
    # the name. A form's フリガナ field puts スガハラ ケンゴ on one line and 氏名
    # on the next, which was read as a foreign counterpart named ケンゴ氏 and
    # then judged for international reciprocity and irreplaceability.
    r"[ァ-ヴー]{3,}[ 　]{0,2}(?:氏(?!名)|教授|博士)|"
    r"(?:TU|ETH|MIT|EPFL)\s*[A-Za-z]*|"
    r"[A-Z][a-z]+\s+(?:University|Institute)"
)
# International venues and societies. Publishing there is evidence of
# international activity, but they are not counterparts.
_NAMED_INTERNATIONAL_VENUE = re.compile(r"(?:IGTE|CEFC|COMPUMAG|ICEM|IEEE)")
# A claimed relationship with someone abroad. These are the words that make
# 「相手先を名指ししていない」 a fair thing to say.
_INTERNATIONAL_RELATION_MARKERS = (
    "共同研究", "共著", "招へい", "招聘", "招請", "受入", "受け入れ",
    "派遣", "訪問", "滞在", "連携", "留学", "分担", "共同開発",
    # NOT bare 交流: in an electrical proposal it is alternating current, and
    # a glossary row for 誘導加熱 was read as an international exchange.
    "国際交流", "学術交流", "人的交流",
)
# Words that make a sentence about the applicant's own international activity,
# as opposed to the worldwide importance of a problem or the size of a market.
_INTERNATIONAL_ACTIVITY_MARKERS = _INTERNATIONAL_RELATION_MARKERS + (
    "発表", "参加", "投稿", "採択", "登壇", "議論", "レビュー",
)
_RECIPROCAL_MARKERS = (
    "相互", "双方向", "還流", "共同研究", "共著", "招へい", "招聘", "派遣",
    "受入", "交流", "分担", "共同開発", "共同実装",
)
_ONE_WAY_MARKERS = (
    "追いつく", "追随", "導入する", "学ぶ", "取り入れる", "輸入",
    "後追い", "キャッチアップ",
)
_INTERNATIONAL_OUTPUT_MARKERS = (
    "共著", "国際会議", "国際学会", "国際誌", "国際共同", "英文",
    "国際ベンチマーク", "国際レビュー", "査読",
)

# Whether an international output already exists decides how much it is
# worth. A record a reviewer can look up outranks an intention, and mixing
# the two lets a plan read as an achievement.
_ACHIEVED_MARKERS = (
    "採択", "掲載", "受理", "出版", "発表した", "共著論文",
    "実施した", "完了した", "招へいした", "訪問した", "得た", "行った",
)
_PLANNED_MARKERS = (
    "予定", "目指す", "したい", "見込み", "構想", "を計画",
)
_NATIONAL_VALUE_MARKERS = (
    "日本発", "我が国独自", "国内で発展", "国内発", "本邦",
    "日本独自", "国内の知見",
)


def grant_writing_international_standing_check(text: str) -> dict:
    """Check that an international claim is evidenced, not aspirational.

    Many programmes weigh an international dimension, and several state it
    as an explicit criterion. The KAKENHI wording decomposes into three
    things a proposal can actually show: leading the field, contributing
    through collaboration, and creating value distinct to one's own country.
    None of them are demonstrated by saying 「国際的に展開する」.

    The check is optional and fires only when the draft raises the subject.
    A 2025 review disclosure scored this axis 1.60 against an adopted 2.70,
    which is why aspiration without a named partner is treated as the
    defect it is.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    # A country name inside a travel or cost line is a venue, not a claim of
    # international standing. 「Conference（2026/5/17~22, フランス）：50万円」
    # in a budget table was read as one, and the proposal was then told to name
    # the counterpart institution it had never claimed to have.
    logistics = re.compile(
        r"円|旅費|参加費|宿泊|渡航|出張|"
        r"\d{4}\s*[/年]\s*\d{1,2}|\d{1,2}\s*[/月]\s*\d{1,2}"
    )
    # 「世界的な社会課題である」 says the problem matters everywhere, not that
    # the applicant works with anyone abroad, and the check went on to demand
    # the counterpart institution behind a claim never made. A trigger counts
    # only in a sentence that also describes activity or a relationship.
    trigger_lines = [s for s in re.split(r"[。．\n]", text) if s.strip()]
    trigger_hits = [
        t
        for t in _INTERNATIONAL_TRIGGERS
        if t in text
        and any(
            t in line
            and not logistics.search(line)
            and any(a in line for a in _INTERNATIONAL_ACTIVITY_MARKERS)
            for line in trigger_lines
        )
    ]
    if _NAMED_INTERNATIONAL_VENUE.search(text):
        trigger_hits = trigger_hits + [_NAMED_INTERNATIONAL_VENUE.search(text).group(0)]
    # Naming a foreign institution is itself the subject being raised. A draft
    # that says ウィーン工科大学 is international whether or not it also says
    # 国際 or a region name, and requiring the word skipped exactly those.
    named_trigger = _NAMED_PARTNER.search(text)
    if named_trigger:
        trigger_hits = trigger_hits + [named_trigger.group(0).strip()]

    # Scope by prose segment, not by sentence. Text extracted from a PDF table
    # or diagram carries no full stops, so a whole page becomes one "sentence"
    # and any two words in it appear to co-occur: a business-model figure put
    # 連携 beside an unrelated 海外 that way.
    segments = [
        s.strip()
        for s in re.split(r"(?<=[。．!?！？])|\n", text)
        if s and s.strip() and len(s.strip()) <= 200
    ]
    # The relationship must be an international one. A domestic 共同開発 with a
    # partner company, in a proposal that separately mentions 海外市場, is not a
    # foreign collaboration missing its counterpart's name.
    relations = [
        m
        for m in _INTERNATIONAL_RELATION_MARKERS
        if any(
            m in s and any(t in s for t in _INTERNATIONAL_TRIGGERS)
            for s in segments
        )
    ]
    # Citing foreign prior work is not a claim of international standing. A
    # domestic 基盤 proposal that surveyed 「フランス・Clenet教授らによる」 was
    # told it showed no international output, which is true and irrelevant.
    # The check opens on a claimed relationship, a named counterpart or venue,
    # or a claimed output -- not on a region name sitting near a verb.
    if not (
        relations
        or named_trigger
        or _NAMED_INTERNATIONAL_VENUE.search(text)
        or any(m in text for m in _INTERNATIONAL_OUTPUT_MARKERS)
    ):
        trigger_hits = []
    if not trigger_hits:
        return {
            "applicable": False,
            "score": None,
            "risks": [],
            "comments": [],
            "target": (
                "an international claim names its partners, its exchanged "
                "artefacts, and what flows back"
            ),
            "source": "international-standing check",
        }

    sentences = [s for s in re.split(r"(?<=[。．!?！？])", text) if s.strip()]
    partners = sorted({m.group(0).strip() for m in _NAMED_PARTNER.finditer(text)})
    reciprocal = [m for m in _RECIPROCAL_MARKERS if m in text]
    one_way = [m for m in _ONE_WAY_MARKERS if m in text]
    outputs = [m for m in _INTERNATIONAL_OUTPUT_MARKERS if m in text]
    # Publishing in IEEE Transactions or presenting at CEFC is international
    # output. Counting those venues as evidence of international activity while
    # reporting that the document shows none was a contradiction the adopted
    # 基盤 proposal exposed.
    outputs += sorted({m.group(0) for m in _NAMED_INTERNATIONAL_VENUE.finditer(text)})
    national = [m for m in _NATIONAL_VALUE_MARKERS if m in text]

    # Presenting abroad is international output; it has no counterpart to name.
    # Only a claimed relationship does. A proposal whose international content
    # was 「国際競争力強化に貢献する」 and 「想定する国内、海外市場」 was told
    # to name the partner institution behind a collaboration it never claimed.
    risks: list[dict] = []
    if relations and not partners:
        risks.append({
            "type": "no_named_counterpart",
            "severity": "HIGH",
            "comment": "国際連携に触れているが、相手先の機関名・研究者名がない。",
            "recommendation": (
                "「海外の研究者と連携する」ではなく、機関名と個人名を書く。"
                "審査者が確認できない連携は、意図の表明にとどまる。"
            ),
        })
    # Sentences that carry an international output, split by whether the
    # output exists yet. A plan and a record must not read alike.
    output_sentences = [
        s for s in sentences if any(o in s for o in _INTERNATIONAL_OUTPUT_MARKERS)
    ]
    achieved = [
        s for s in output_sentences if any(a in s for a in _ACHIEVED_MARKERS)
    ]
    planned = [
        s for s in output_sentences
        if any(p in s for p in _PLANNED_MARKERS)
        and not any(a in s for a in _ACHIEVED_MARKERS)
    ]

    if not outputs:
        risks.append({
            "type": "no_international_output",
            "severity": "MEDIUM",
            "comment": "国際的な成果物（共著、国際会議、国際誌等）が示されていない。",
            "recommendation": (
                "既にある共著・国際会議発表・国際レビューを挙げ、本計画で"
                "何を追加するかを書く。"
            ),
        })
    elif planned and not achieved:
        risks.append({
            "type": "international_output_all_planned",
            "severity": "MEDIUM",
            "comment": (
                "国際的な成果がすべて予定であり、既に成立したものがない。"
            ),
            "recommendation": (
                "採択済み・掲載済みのものがあれば、その状態を明記して分ける。"
                "無い場合は形成途上であることを正直に書き、本研究で到達する"
                "段階を示す。実績と予定を混ぜて書くと、予定が実績に読める。"
            ),
        })
    if not reciprocal:
        risks.append({
            "type": "no_reciprocity",
            "severity": "MEDIUM",
            "comment": "やり取りの双方向性（何を渡し、何が返るか）が書かれていない。",
            "recommendation": (
                "招へい・派遣・共同実装・共著など、双方向の往来を具体化する。"
            ),
        })
    # A catch-up frame concedes that the value flows one way, which is the
    # opposite of what an international-standing criterion asks for.
    if one_way and not national:
        risks.append({
            "type": "one_way_catch_up_frame",
            "severity": "MEDIUM",
            "comment": (
                "海外に追随する枠組みで書かれており、自国発の価値が示されていない。"
            ),
            "recommendation": (
                "「追いつく」ではなく、国内で発展した手法と海外の手法を相互検証し、"
                "双方へ還流する構図にする。"
            ),
            "one_way_markers": one_way[:5],
        })

    deductions = sum(
        3.0 if r["severity"] == "HIGH" else 1.5 for r in risks
    )
    score = max(0.0, round(10.0 - deductions, 1))
    return {
        "applicable": True,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "named_counterparts": partners[:12],
        "reciprocity_markers": reciprocal[:8],
        "one_way_markers": one_way[:5],
        "international_outputs": outputs[:8],
        "achieved_output_sentences": [
            re.sub(r"\s+", " ", x).strip()[:160] for x in achieved[:4]
        ],
        "planned_output_sentences": [
            re.sub(r"\s+", " ", x).strip()[:160] for x in planned[:4]
        ],
        "national_value_markers": national[:5],
        "comments": [r["comment"] for r in risks],
        "recommendations": [r["recommendation"] for r in risks],
        "target": (
            "named counterparts, real international outputs, two-way exchange, "
            "and a value that originates here rather than a catch-up plan"
        ),
        "source": "international-standing check",
    }


# What the applicant's side supplies that the partner cannot obtain locally.
# Derived from a 2026 bilateral case: an Austrian counterpart was asked by
# his own funder why Japan, and told that Europe already has CERN. A facility
# was never what he needed -- he needed a method line developed elsewhere to
# cross-validate against. The answer to "why this partner" is always a named
# asset, not a shared intention.
_TRANSFERABLE_ASSET_MARKERS = (
    "日本発", "国内で発展", "我が国独自", "国内発",
    "本研究室が開発", "代表者が開発", "が開発した",
    "提案者の一人", "開発者", "原著者",
)
_PARTNER_DEMAND_MARKERS = (
    "求められ", "要請", "招請", "招へい", "招聘", "打診",
    "関心を示", "議論したい", "共同研究の申し出", "取り上げられ",
    "採用され", "導入され", "参照され", "問い合わせ",
)
_SUBSTITUTE_QUESTION_MARKERS = (
    "他に代え", "代替できない", "他では得られない", "唯一",
    "独立に発展", "別系統", "異なる系譜", "相互検証",
)


def grant_writing_collaboration_irreplaceability_check(text: str) -> dict:
    """Check that a named collaboration says why that partner, both ways.

    A funder on either side asks the same thing: why this counterpart rather
    than someone closer to home. The answer is never a shared intention to
    cooperate. It is a named asset one side holds and the other cannot obtain
    locally, plus evidence that the other side actually wants it.

    The check is optional and applies only when a foreign counterpart is
    named. It does not judge whether the collaboration is good -- it asks
    whether the proposal answers the question a reviewer will ask.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    partner = _NAMED_PARTNER.search(text)
    if partner is None:
        return {
            "applicable": False,
            "score": None,
            "risks": [],
            "comments": [],
            "target": (
                "a named collaboration states the asset only this side holds "
                "and the evidence that the other side wants it"
            ),
            "source": "collaboration-irreplaceability check",
        }

    assets = [m for m in _TRANSFERABLE_ASSET_MARKERS if m in text]
    demand = [m for m in _PARTNER_DEMAND_MARKERS if m in text]
    substitute = [m for m in _SUBSTITUTE_QUESTION_MARKERS if m in text]

    risks: list[dict] = []
    if not assets:
        risks.append({
            "type": "no_asset_this_side_holds",
            "severity": "HIGH",
            "comment": (
                "相手が自国で得られない、こちら側固有の資産が示されていない。"
            ),
            "recommendation": (
                "手法名・ライブラリ名・原著者を挙げ、それが国内で発展した"
                "ものであることを書く。「連携する」だけでは、なぜその相手か"
                "にも、なぜこちらかにも答えていない。"
            ),
        })
    if not demand:
        risks.append({
            "type": "no_evidence_partner_wants_it",
            "severity": "MEDIUM",
            "comment": "相手側がそれを求めている証拠が書かれていない。",
            "recommendation": (
                "招請、共同研究の申し出、相手が自国で取り上げた事実など、"
                "相手側から来た動きを具体的に書く。こちらの意欲ではなく、"
                "相手の需要が連携の必然性を示す。"
            ),
        })
    if not substitute:
        risks.append({
            "type": "no_substitution_argument",
            "severity": "LOW",
            "comment": (
                "近場で代替できない理由が書かれていない。"
            ),
            "recommendation": (
                "独立に発展した別系統との相互検証である、など、"
                "同一国内や近隣機関では代えられない理由を一文で置く。"
            ),
        })

    deductions = sum(
        3.0 if r["severity"] == "HIGH" else 1.5 if r["severity"] == "MEDIUM" else 0.5
        for r in risks
    )
    score = max(0.0, round(10.0 - deductions, 1))
    return {
        "applicable": True,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "named_partner": partner.group(0).strip(),
        "asset_markers": assets[:8],
        "demand_markers": demand[:8],
        "substitution_markers": substitute[:8],
        "comments": [r["comment"] for r in risks],
        "recommendations": [r["recommendation"] for r in risks],
        "target": (
            "name the asset only this side holds, show the other side asking "
            "for it, and say why a closer substitute will not do"
        ),
        "source": "collaboration-irreplaceability check",
    }


_FORM_OVERFLOW_NOTICE = re.compile(
    r"「(?P<field>[^」]{2,40})」は(?P<limit>\d+)ページ以内で"
)
_TEX_FIELD_LIMIT = re.compile(
    r"\\section\{(?P<title>[^{}]+)\}(?P<tail>(?:[^\n]*\n){0,5})",
)
_TEX_LIMIT_VALUE = re.compile(r"[＜<]{2}\s*最大\s*(?P<pages>\d+)\s*ページ\s*[＞>]{2}")


def _flatten(text: str) -> str:
    return re.sub(r"\s+", "", text)


def grant_writing_page_limit_check(pdf_path: str, tex_dir: str = "") -> dict:
    """Check each field of a compiled proposal against its page allowance.

    A page limit is the one rule a funder enforces before anyone reads the
    science: a field that runs past its allowance can be returned unexamined.
    The limit is also a target. A field that leaves a page unused has thrown
    away space the applicant was given to argue in, which is the same defect
    seen from the other side.

    Two independent signals are used. Japanese form templates print their own
    notice (「<欄名>」はNページ以内で書いてください) onto the overflow page, and
    that string in the compiled PDF is proof on its own. Separately, the field
    spans measured from the PDF are compared with the allowances declared in
    the LaTeX source (＜＜最大　Nページ＞＞), which catches a form that stays
    silent.
    """
    pdf = pathlib.Path(pdf_path)
    if not pdf.is_file():
        raise FileNotFoundError(f"compiled proposal not found: {pdf_path}")

    import fitz  # PyMuPDF; imported here so text-only checks never need it.

    doc = fitz.open(str(pdf))
    pages: list[dict] = []
    for index in range(doc.page_count):
        page = doc.load_page(index)
        blocks = [b for b in page.get_text("blocks") if b[4].strip()]
        pages.append({
            "number": index + 1,
            "flat": _flatten(page.get_text()),
            "bottom": max((b[3] for b in blocks), default=0.0),
        })
    doc.close()
    if not pages:
        raise ValueError(f"compiled proposal has no pages: {pdf_path}")

    text_bottom = max(p["bottom"] for p in pages)

    declared: list[tuple[str, int]] = []
    source_dir = pathlib.Path(tex_dir) if tex_dir else pdf.parent
    for tex in sorted(source_dir.glob("*.tex")):
        body = tex.read_text(encoding="utf-8", errors="replace")
        for match in _TEX_FIELD_LIMIT.finditer(body):
            limit = _TEX_LIMIT_VALUE.search(match.group("tail"))
            if limit:
                declared.append((_flatten(match.group("title")), int(limit.group("pages"))))

    spans: list[dict] = []
    starts: list[tuple[int, str, int]] = []
    for title, limit in declared:
        start = next((p["number"] for p in pages if title in p["flat"]), None)
        if start is not None:
            starts.append((start, title, limit))
    starts.sort()
    for position, (start, title, limit) in enumerate(starts):
        end = starts[position + 1][0] - 1 if position + 1 < len(starts) else len(pages)
        last = pages[end - 1]
        spans.append({
            "field": title,
            "declared_max_pages": limit,
            "first_page": start,
            "last_page": end,
            "used_pages": end - start + 1,
            "last_page_fill": round(last["bottom"] / text_bottom, 2) if text_bottom else 0.0,
        })

    risks: list[dict] = []
    for page in pages:
        notice = _FORM_OVERFLOW_NOTICE.search(page["flat"])
        if notice:
            risks.append({
                "severity": "CRITICAL",
                "location": f"PDF p{page['number']}",
                "comment": (
                    f"様式が超過を印字している: 「{notice.group('field')}」は"
                    f"{notice.group('limit')}ページ以内。"
                ),
                "recommendation": "溢れた分の文をまるごと落とすか、別欄へ移す。",
            })

    for span in spans:
        if span["used_pages"] > span["declared_max_pages"]:
            risks.append({
                "severity": "CRITICAL",
                "location": f"PDF p{span['first_page']}-p{span['last_page']}",
                "comment": (
                    f"「{span['field']}」は{span['declared_max_pages']}ページ指定に対し"
                    f"{span['used_pages']}ページ占めている。"
                ),
                "recommendation": "文を圧縮せず、文・段落単位で落とすか別欄へ移す。",
            })
        elif span["used_pages"] < span["declared_max_pages"] or (
            span["last_page_fill"] < 0.6 and span["declared_max_pages"] >= 2
        ):
            # A single-page administrative field is often short because the
            # honest answer is short, so only an unused whole page or a slack
            # multi-page allowance is reported.
            risks.append({
                "severity": "MEDIUM",
                "location": f"PDF p{span['last_page']}",
                "comment": (
                    f"「{span['field']}」は{span['declared_max_pages']}ページ許容のうち"
                    f"{span['used_pages']}ページ、最終ページの充填率"
                    f"{span['last_page_fill']:.0%}。"
                ),
                "recommendation": "許容ページは埋める対象。証拠・数値・図で残りを使う。",
            })

    if not declared and not risks:
        return {
            "applicable": False,
            "score": None,
            "risk_count": 0,
            "risks": [],
            "page_count": len(pages),
            "fields": [],
            "comments": [],
            "recommendations": [],
            "reason": (
                "ページ上限の宣言（＜＜最大　Nページ＞＞）が見つからず、様式の超過印字もない。"
            ),
            "source": "page-limit check",
        }

    deductions = sum(
        5.0 if r["severity"] == "CRITICAL" else 1.5 if r["severity"] == "MEDIUM" else 0.5
        for r in risks
    )
    return {
        "applicable": True,
        "score": max(0.0, round(10.0 - deductions, 1)),
        "risk_count": len(risks),
        "risks": risks,
        "page_count": len(pages),
        "fields": spans,
        "comments": [r["comment"] for r in risks],
        "recommendations": [r["recommendation"] for r in risks],
        "target": "each field inside its allowance, and filling it",
        "source": "page-limit check",
    }


# Someone who carries work without a responsibility or budget share. A funder
# reads 研究分担者 as accountable and these as not.
_NON_MEMBER_ROLE = re.compile(r"連携研究者|研究協力者|協力者|アドバイザー|オブザーバ")
# A line that hands a named person a job, not prose that happens to use the
# role word: 「連携研究者　浅川伸一：機械学習に関する専門知識の供与」 assigns,
# while 「有能な研究協力者を有する」 describes. Only the former is checkable.
_ROLE_ASSIGNMENT = re.compile(
    r"(?:^|[\s　、，])(?P<role>連携研究者|研究協力者|アドバイザー|オブザーバ)"
    r"[\s　]{0,4}(?P<name>[^\s　：:、，。]{2,12})[：:]\s*(?P<desc>.+)$"
)
# The form explains these roles in its own instruction text, which describes
# nobody and therefore carries no capability.
_FORM_INSTRUCTION_HINT = re.compile(
    r"本欄には|記入して|記述して|してください|下さい|参照|場合には|必要に応じて"
)
# Vocabulary every proposal uses regardless of field: uncovered here means
# nothing.
_CAPABILITY_STOPWORDS = frozenset({
    "本研究", "研究", "開発", "手法", "方法", "解析", "評価", "検証", "計算",
    "設計", "技術", "問題", "課題", "目的", "計画", "実施", "実装", "適用",
    "利用", "使用", "対象", "結果", "効果", "内容", "以下", "以上", "場合",
    "今回", "一方", "同様", "現在", "近年", "従来", "既存", "各種", "全体",
    "期間", "年度", "本欄", "記述", "検討", "提案", "構築", "向上", "実現",
    "必要", "可能", "重要", "特徴", "状況", "分野", "国内", "国外",
    # Role-description filler: these say what the person does for the
    # project, not which capability they hold.
    "専門知識", "供与", "助言", "提供", "全般", "担当", "統括", "指導",
    "協力", "支援", "情報", "経験", "知見", "専門家",
})


def grant_writing_capability_responsibility_check(text: str) -> dict:
    """Check who carries the capability the novelty rests on.

    A proposal usually joins a field the applicant knows to one they do not,
    and reviewers read the capability criterion by asking whether the team can
    do the part that is new. The answer is a matter of roles: 研究代表者 and
    研究分担者 are accountable and funded, while 連携研究者, 研究協力者 and
    アドバイザー are not counted the same way.

    Measured on a rejected 基盤C whose novelty was machine learning applied to
    topology optimisation: the applicant's 23 listed items were patents and
    papers on 電磁界解析 and accelerators with no machine-learning entry, and
    the line that supplied the missing capability read 「連携研究者　浅川伸一：
    機械学習に関する専門知識の供与」. The same document's adopted counterpart
    named five people and gave every one of them a 分担 role.

    An earlier version of this check compared the novelty vocabulary against
    the words in the evidence list. It fired hardest on the *adopted*
    proposal, whose novelty is a compound it coined (マルチフィジクスモデル
    縮約) that no paper title could contain, and stayed silent on the rejected
    one. Lexical overlap does not measure capability, so only the role
    attribution -- which is mechanical and locatable -- is tested here.
    """
    text = _prose_for_lint(_read_text_if_path(text))

    segments = [
        s.strip()
        for s in re.split(r"(?<=[。．!?！？])|\n", text)
        if s and s.strip()
    ]
    assignments = [
        (s, m)
        for s in segments
        for m in [_ROLE_ASSIGNMENT.search(s)]
        if m and not _FORM_INSTRUCTION_HINT.search(s)
    ]
    if not assignments:
        return {
            "applicable": False,
            "score": None,
            "risk_count": 0,
            "risks": [],
            "carried_terms": [],
            "role_lines": [],
            "comments": [],
            "recommendations": [],
            "reason": (
                "連携研究者・研究協力者・アドバイザーへの担当割り当て行がないため"
                "判定しない。"
            ),
            "source": "capability-responsibility check",
        }

    risks: list[dict] = []
    carried: list[str] = []
    for line, match in assignments:
        role = match.group("role")
        # Only the description after the name states a capability; the name
        # itself and the role word are not ones.
        terms = sorted({
            term
            for term in _CLAIM_TERM.findall(match.group("desc"))
            if len(term) >= 3
            and term not in _CAPABILITY_STOPWORDS
            and not _NON_MEMBER_ROLE.fullmatch(term)
            and text.count(term) >= 3
        })
        if not terms:
            continue
        carried.extend(terms)
        risks.append({
            "type": "novelty_capability_on_non_member",
            "severity": "HIGH",
            "role": role,
            "terms": terms,
            "excerpt": line[:120],
            "comment": (
                f"「{'、'.join(terms)}」の能力が{role}に置かれている。"
            ),
            "recommendation": (
                "予算と責任を持つ研究分担者にする。連携研究者・協力者は"
                "遂行体制として数えにくく、その能力を要する主張は"
                "裏付けのない主張として読まれる。"
            ),
        })

    deductions = sum(3.0 for _ in risks)
    return {
        "applicable": True,
        "score": max(0.0, round(10.0 - deductions, 1)),
        "risk_count": len(risks),
        "risks": risks,
        "carried_terms": sorted(set(carried)),
        "role_lines": [line[:120] for line, _ in assignments[:4]],
        "comments": [r["comment"] for r in risks],
        "recommendations": [r["recommendation"] for r in risks],
        "target": (
            "every capability the claim depends on sits with someone who has "
            "a responsibility share"
        ),
        "source": "capability-responsibility check",
    }


def grant_writing_question_originality_check(text: str) -> dict:
    """Check that the central question carries an originality position.

    The criterion has two halves and this suite already covers the first:
    the question must be clear (see the central-claim consistency check).
    This one covers the second half -- that originality and creativity are
    recognisable -- which needs three things present and connected: a stated
    question, a claim of what is new, and a gap in prior work that the claim
    stands against. An originality adjective with no gap behind it is an
    assertion; a gap with no question attached belongs to somebody else's
    proposal.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    statements = _claim_statements(text)
    if not statements:
        return {
            "applicable": False,
            "score": None,
            "risks": [],
            "comments": [],
            "target": (
                "the central question states what is new and what prior work "
                "leaves unresolved"
            ),
            "source": "question-originality check (2025 review disclosure)",
        }

    sentences = [s for s in re.split(r"(?<=[。．!?！？])", text) if s.strip()]
    originality_hits = [m for m in _ORIGINALITY_MARKERS if m in text]
    prior_hits = [m for m in _PRIOR_WORK_MARKERS if m in text]
    gap_hits = [m for m in _GAP_MARKERS if m in text]

    # A gap statement is a position only when prior work and the gap are
    # joined. Either in one sentence -- 「既往研究は…進めてきたが、…は確立
    # していない」 -- or across two that a contrastive connective ties
    # together: 「既往研究では…進められてきた。一方、…は体系化されていない」.
    # Both are ordinary Japanese; requiring the single-sentence form only
    # would fail correct prose.
    contrastive = ("一方", "しかし", "だが", "ところが", "他方", "これに対し")
    gap_sentences = []
    for i, sentence in enumerate(sentences):
        has_prior = any(p in sentence for p in _PRIOR_WORK_MARKERS)
        has_gap = any(g in sentence for g in _GAP_MARKERS)
        if has_prior and has_gap:
            gap_sentences.append({
                "sentence_index": i + 1,
                "form": "single_sentence",
                "excerpt": re.sub(r"\s+", " ", sentence).strip()[:220],
            })
            continue
        if not has_prior:
            continue
        nxt = sentences[i + 1] if i + 1 < len(sentences) else ""
        if any(g in nxt for g in _GAP_MARKERS) and any(
            c in nxt for c in contrastive
        ):
            gap_sentences.append({
                "sentence_index": i + 1,
                "form": "contrastive_pair",
                "excerpt": re.sub(r"\s+", " ", sentence + nxt).strip()[:260],
            })

    risks: list[dict] = []
    if not originality_hits:
        risks.append({
            "type": "no_originality_claim",
            "severity": "HIGH",
            "comment": "何が新しいのかを述べた語がない。",
            "recommendation": (
                "独自性・新規性を一語で名指しする。審査項目は「学術的独自性や"
                "創造性が認められるか」であり、読み取れなければ低評価になる。"
            ),
        })
    if not gap_sentences:
        risks.append({
            # MEDIUM, not HIGH. Measured on ten submitted proposals (three
            # funded): this rule fires on funded work too, so it cannot carry
            # the severity of a defect. It marks a contrast worth writing,
            # not a reason the proposal will fail.
            "type": "no_gap_against_prior_work",
            "severity": "MEDIUM",
            "comment": (
                "既往研究の限界を一文で述べていない。"
                if not (prior_hits and gap_hits)
                else "既往研究への言及と未解決の指摘が別々の文に散っている。"
            ),
            "recommendation": (
                "「既往研究は〜を進めてきたが、〜は確立していない」の形で、"
                "先行研究と未解決点を同一文に置く。独自性はこの対比の上に立つ。"
            ),
        })

    deductions = sum(
        3.0 if r["severity"] == "HIGH" else 1.5 for r in risks
    )
    score = max(0.0, round(10.0 - deductions, 1))
    return {
        "applicable": True,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "statement_count": len(statements),
        "originality_markers": originality_hits[:8],
        "prior_work_markers": prior_hits[:8],
        "gap_statements": gap_sentences[:5],
        "comments": [r["comment"] for r in risks],
        "recommendations": [r["recommendation"] for r in risks],
        "target": (
            "the central question states what is new and what prior work "
            "leaves unresolved, in one contrast"
        ),
        "source": "question-originality check (2025 review disclosure)",
    }


_AMOUNT_PATTERN = re.compile(r"\d[\d,]*\s*(?:千円|万円|億円|円)")
_NECESSITY_MARKERS = (
    "必要性", "計上", "見積", "そのため", "必要である", "必要がある",
    "購入する", "使用する",
)
# Money the applicant expects to receive, not to spend.
_REVENUE_MARKERS = (
    "顧客", "販売", "売上", "収益", "価格", "ライセンス", "採算", "利益",
    "事業化", "市場規模", "単価設定", "課金",
)
_TRAVEL_MARKERS = ("旅費", "出張", "渡航")
_DISSEMINATION_MARKERS = (
    "学会", "国際会議", "研究会", "発表", "シンポジウム", "講演", "報告会",
)
_CALCULATION_BASIS_MARKERS = (
    "×", " x ", " X ", "単価", "月額", "年額", "人泊", "見積書",
    "公式料金", "料金表", "旅費規程", "利用料", "数量", "月数", "回数",
)
_QUANTITY_BASIS_PATTERN = re.compile(r"\d+\s*(?:名|回|件|台|個|月|日|泊)")


def grant_writing_budget_narrative_check(text: str) -> dict:
    """Check the necessity narrative that sits beside a budget table.

    The FY2027 JSPS Web entry guide requires both necessity and a calculation
    basis. A 2019 editor also correctly warned that copying bare totals into
    prose creates two values to maintain. Therefore an amount in the
    narrative is acceptable when it forms a recomputable basis (unit price,
    quantity, duration, official tariff, or quotation); only a bare amount is
    reported.

    The check is optional: it applies only where a necessity narrative and
    money both appear.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    sentences = [s for s in re.split(r"(?<=[。．!?！？])", text) if s.strip()]

    narrative = [
        (i, s) for i, s in enumerate(sentences)
        if any(m in s for m in _NECESSITY_MARKERS)
    ]
    if not narrative:
        return {
            "applicable": False,
            "score": None,
            "risks": [],
            "comments": [],
            "target": (
                "the narrative explains necessity and gives a recomputable "
                "calculation basis where an amount is stated"
            ),
            "source": (
                "FY2027 JSPS Web entry guide plus 2019 editor review of a "
                "funded proposal"
            ),
        }

    risks: list[dict] = []
    seen_excerpts: set[str] = set()
    for index, sentence in narrative:
        # A price charged is not a cost incurred. 「ライセンスビジネスの権利
        # 付与型は…1件あたり5,000千円に設定」 is a business model, and 計上 in
        # it refers to the company's own accounting, not to a budget line.
        if any(m in sentence for m in _REVENUE_MARKERS):
            continue
        if sentence in seen_excerpts:
            continue
        seen_excerpts.add(sentence)
        amounts = _AMOUNT_PATTERN.findall(sentence)
        has_calculation_basis = any(
            marker in sentence for marker in _CALCULATION_BASIS_MARKERS
        ) or bool(_QUANTITY_BASIS_PATTERN.search(sentence))
        if amounts and not has_calculation_basis:
            risks.append({
                "type": "amount_without_calculation_basis",
                "severity": "MEDIUM",
                "sentence_index": index + 1,
                "amounts": amounts[:4],
                "excerpt": re.sub(r"\s+", " ", sentence).strip()[:200],
                "comment": (
                    "必要性の説明に積算根拠のない金額がある: "
                    + "、".join(amounts[:3])
                ),
                "recommendation": (
                    "金額を残すなら、単価×数量×月数/回数、見積書、公式料金表等の"
                    "積算根拠を同じ記述に置く。合計額だけなら積算表へ集約する。"
                ),
            })

    travel = [s for s in sentences if any(t in s for t in _TRAVEL_MARKERS)]
    if travel and not any(
        any(d in s for d in _DISSEMINATION_MARKERS) for s in sentences
    ):
        risks.append({
            "type": "travel_without_dissemination_plan",
            "severity": "LOW",
            "comment": "旅費を計上しているが、発表・参加する場が書かれていない。",
            "recommendation": (
                "出張の行き先だけでなく、学会・国際会議など何のための移動かを"
                "書く。成果発表の予定がないなら、その旨を明示する。"
            ),
        })

    deductions = sum(
        1.5 if r["severity"] == "MEDIUM" else 0.5 for r in risks
    )
    score = max(0.0, round(10.0 - deductions, 1))
    return {
        "applicable": True,
        "score": score,
        "risk_count": len(risks),
        "risks": risks[:20],
        "necessity_sentence_count": len(narrative),
        "comments": list(dict.fromkeys(r["comment"] for r in risks)),
        "recommendations": list(dict.fromkeys(r["recommendation"] for r in risks)),
        "target": (
            "the narrative says what the money buys and why, amounts carry a "
            "recomputable basis, and travel names its purpose"
        ),
        "source": (
            "FY2027 JSPS Web entry guide plus 2019 editor review of a funded "
            "proposal"
        ),
    }


def grant_writing_template_residue_check(text: str) -> dict:
    """Find unfilled placeholders and leftover form instructions.

    Evidence from two real 2026 applications. A Power Academy draft still
    carried 「小計：○○○○千円（税込）」 and the form's own 「〜記入して
    ください」 sentences when it reached the co-investigator, who deleted
    them. A JSPS application was sent back by the office before any reviewer
    saw it, purely on form compliance. Both defect classes are mechanical,
    locatable, and fatal in a way no argument about research quality is.
    """
    text = _read_text_if_path(text)
    lines = text.splitlines()
    risks: list[dict] = []
    instructions: list[dict] = []

    placeholder = re.compile(
        r"[○◯]{2,}|[×✕]{3,}|＿{2,}|_{4,}|"
        r"[XxＸｘ]{3,}(?![A-Za-z0-9])|"
        r"【\s*(?:記入|入力|ここに|要記入)[^】]*】|"
        # A prefix match on 「（入力…）」 is not enough. 入力 is ordinary
        # technical vocabulary, and 「（入力したエネルギーに対するビーム強度）」
        # -- a term gloss in a submitted proposal -- was reported as an
        # unfilled field. A placeholder parenthetical holds the instruction
        # word and nothing else.
        r"（\s*(?:記入|入力)(?:して\s*(?:ください|下さい)|欄|例|箇所|事項)?\s*）|"
        r"(?:氏名|所属|役職|研究者番号|研究期間|申請金額|課題名)"
        r"[^。\n]{0,20}未定|"
        r"\bTBD\b|\bTODO\b"
    )
    # Instruction sentences are counted but NOT reported as defects. Measured
    # on the 2026 Power Academy rewrite: the co-investigator deleted 6 of them
    # and deliberately kept 13, and the two groups are indistinguishable in
    # flat text -- 「研究マップを…選択してください」 went, 「性別を…選択して
    # ください」 stayed, because one labels a prose box the applicant fills and
    # the other is part of the form's own answer structure. A rule that cannot
    # tell them apart cannot be a detector, so this becomes a count the author
    # judges, not a finding.
    instruction = re.compile(
        r"(?:記入|記載|入力|選択|要約|確認|参照)して\s*ください|"
        r"ご記入|ご覧いただき|お書きください"
    )

    for number, line in enumerate(lines, 1):
        stripped = line.strip()
        if not stripped:
            continue
        for match in placeholder.finditer(stripped):
            risks.append({
                "type": "unfilled_placeholder",
                "severity": "HIGH",
                "line": number,
                "match": match.group(0)[:40],
                "excerpt": stripped[:200],
                "comment": "未記入のプレースホルダが残っている。",
                "recommendation": (
                    "提出前に実際の値へ置き換える。金額欄の ○○○○ は"
                    "そのまま提出されると事務差し戻しの原因になる。"
                ),
            })
            break
        if instruction.search(stripped):
            instructions.append({"line": number, "excerpt": stripped[:200]})

    applicable = bool(text.strip())
    deductions = sum(2.0 for r in risks)
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    question = (
        f"様式の記入説明文が {len(instructions)} 件ある。"
        "自由記述欄の説明文なら提出前に削除し、選択・回答欄の様式文なら残す。"
        "文面だけでは区別できないため、欄ごとに確認する。"
    ) if instructions else ""
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks[:40],
        "instruction_sentence_count": len(instructions),
        "instruction_sentences": instructions[:20],
        "questions": [question] if question else [],
        "comments": list(dict.fromkeys(r["comment"] for r in risks)),
        "recommendations": list(dict.fromkeys(r["recommendation"] for r in risks)),
        "target": (
            "no unfilled placeholder survives into the submitted document; "
            "instruction sentences are counted for the author to judge"
        ),
        "source": "template-residue check (2026 Power Academy / JSPS evidence)",
    }


def grant_writing_vague_claim_verb_check(text: str) -> dict:
    """Flag 統合/連携/活用 that never say how.

    A proposal that "integrates three technologies" has told the reviewer
    nothing; one that "couples the two models bidirectionally through the
    winding current and induced voltage" has. This check finds the first
    shape and asks for the second. It is the edit an experienced PI makes
    first, and no keyword-coverage axis sees it, because the vague verb and
    the required nouns are all present -- separately.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    # A newline ends a sentence here too: a person's name on its own line
    # merged with the paragraph under it, and the merged run was reported
    # eighteen times across six real proposals.
    sentences = [s for s in re.split(r"(?<=[。．!?！？])|\n", text) if s.strip()]

    risks: list[dict] = []
    concrete: list[dict] = []
    # Work already done is a report, not a promise. 「両者を辺要素OSSへ統合し、
    # …再実行した」 names an artefact and states an outcome; asking it to say
    # how is asking it to re-describe finished work.
    completed = re.compile(
        r"(?:した|した。|してきた|実施した|再実行した|完了した|得た|"
        r"確認した|発表した|示した)"
    )
    # The same applies to a record written in the ongoing form. 「…知識を活用
    # して先端技術開発を行っている（S1,2）」 in a これまでの研究活動 field is a
    # record with a citation attached, not a promise about the project.
    record_ending = re.compile(r"(?:た|ている|ており|てきた|ています)$")
    trailing_note = re.compile(r"[（(][^（(）)]*[）)]\s*$")
    for index, sentence in enumerate(sentences):
        verb = next((v for v in _VAGUE_CLAIM_VERBS if v in sentence), None)
        if verb is None:
            continue
        if completed.search(sentence):
            continue
        stem = trailing_note.sub("", sentence.strip().rstrip("。．"))
        if record_ending.search(stem):
            continue
        # 「誘導加熱技術を活用する幅広い産業分野」 modifies a noun; it describes
        # who uses the technology, not what this project will do. A claim verb
        # followed immediately by a noun is adnominal, not the predicate.
        after = sentence[sentence.index(verb) + len(verb):]
        if after and re.match(r"[一-龥ァ-ヴー]", after):
            continue
        found = [m for m in _MECHANISM_MARKERS if m in sentence]
        entry = {
            "sentence_index": index + 1,
            "verb": verb,
            "mechanism_markers": found[:6],
            "excerpt": re.sub(r"\s+", " ", sentence).strip()[:220],
        }
        if found:
            concrete.append(entry)
            continue
        entry.update({
            "type": "claim_verb_without_mechanism",
            "severity": "MEDIUM",
            "comment": (
                f"「{verb}」が、何をどう渡すのかを書かずに使われている。"
            ),
            "recommendation": (
                "動詞を操作に置き換える。何と何を、どの物理量を介して、"
                "どちら向きに渡すのかを書く。"
                "例: 「三者の技術を統合し」→「両モデルを巻線電流と誘起電圧を"
                "介して双方向に連成し」。"
            ),
        })
        risks.append(entry)

    applicable = bool(risks or concrete)
    score = None if not applicable else max(0.0, round(10.0 - 1.5 * len(risks), 1))
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "concrete_uses": concrete,
        "comments": list(dict.fromkeys(r["comment"] for r in risks)),
        "recommendations": list(dict.fromkeys(r["recommendation"] for r in risks)),
        "target": (
            "every 統合/連携/活用 names the operation, the quantity exchanged, "
            "and the direction"
        ),
        "source": "vague-claim-verb check (Power Academy 2026 rewrite evidence)",
    }


def grant_writing_kaken_review_format_check(text: str) -> dict:
    """Check KAKENHI reviewer-format realities on a proposal draft.

    Combines the official B/C research-plan elements and internationality
    rating with the in-house R9/FY2027 call briefing: reviewers may read up to
    ~100 proposals in about a month; figures may be printed in monochrome for
    some categories; publication records are read through researchmap; the
    human-rights/legal-compliance box draws the most reviewer remarks; and the
    funding-overlap box has a fixed format. Fragment-level triggers gate each
    sub-check, so short excerpts stay clean; full-draft heuristics apply above
    ~1500 chars.
    """
    raw = _read_text_if_path(text)
    prose = _prose_for_lint(raw)
    risks: list[dict] = []

    def add_risk(
        risk_type: str,
        start: int,
        excerpt: str,
        comment: str,
        recommendation: str,
        severity: str = "MEDIUM",
        **details,
    ) -> None:
        item = {
            "type": risk_type,
            "line": raw.count("\n", 0, max(0, start)) + 1,
            "severity": severity,
            "excerpt": re.sub(r"\s+", " ", excerpt).strip()[:360],
            "comment": comment,
            "recommendation": recommendation,
        }
        item.update(details)
        risks.append(item)

    color_figure_pattern = re.compile(
        r"(?:赤|青|緑|黄|橙|紫|桃)(?:色|い)?(?:の)?"
        r"(?:実線|破線|点線|一点鎖線|線|棒|丸印|丸|矢印|領域|塗り|"
        r"プロット|曲線|マーカー|文字|枠)"
        r"|色分け|色で(?:区別|示|表)|カラーで(?:区別|示|表)"
    )
    mono_terms = [
        "白黒",
        "モノクロ",
        "グレースケール",
        "線種",
        "濃淡",
        "ハッチング",
        "マーカー形状",
    ]
    color_matches = list(color_figure_pattern.finditer(raw))
    if color_matches and not any(term in raw for term in mono_terms):
        first = color_matches[0]
        add_risk(
            "color_dependent_figure",
            first.start(),
            first.group(0),
            "色の違いだけで図の系列・領域を区別している。審査時に白黒印刷される種目がある。",
            "線種・マーカー・直接ラベル・濃淡で区別し、白黒でも判別できる図にする。",
            severity="HIGH",
            occurrence_count=len(color_matches),
        )

    subject_terms = [
        "アンケート",
        "質問紙",
        "被験者",
        "調査対象者",
        "インタビュー",
        "動物実験",
        "実験動物",
        "ヒト由来",
        "臨床",
        "個人情報",
    ]
    # NOTE: 「遵守」 alone is NOT a safeguard: the box heading itself
    # (「人権の保護及び法令等の遵守への対応」) contains it, so it would
    # suppress the check on every draft that quotes the heading.
    safeguard_terms = [
        "倫理審査",
        "倫理委員会",
        "動物実験委員会",
        "同意",
        "インフォームド",
        "承認",
        "匿名化",
        "個人情報保護",
        "適切に管理",
    ]
    subject_hits = [term for term in subject_terms if term in prose]
    if subject_hits and not any(term in prose for term in safeguard_terms):
        add_risk(
            "human_subjects_without_safeguard",
            max(0, raw.find(subject_hits[0])),
            subject_hits[0],
            "人・動物・個人情報を扱う記述があるのに、講じる対策・措置の記載が見当たらない。",
            "「人権の保護及び法令等の遵守への対応」欄に、倫理審査、同意取得、"
            "匿名化等の具体的な対策を記載する。例年、審査委員からの指摘が最も多い欄である。",
            severity="HIGH",
            subject_hits=subject_hits,
        )

    na_pattern = re.compile(r"該当\s*(?:は|事項は)?\s*(?:なし|ない|無し|ありません)")
    rationale_terms = [
        "ため",
        "ので",
        "対象としない",
        "対象とせず",
        "行わない",
        "用いない",
        "使用しない",
        "含まない",
        "扱わない",
        "のみであり",
        "のみで",
        "理由",
    ]
    ethics_context_terms = [
        "人権",
        "法令",
        "倫理",
        "安全対策",
        "個人情報",
        "被験者",
        "動物",
        "ヒト",
        "アンケート",
        "インタビュー",
    ]
    bare_na_matches = []
    # Inspect applicant prose, not raw LaTeX. A final-year field commonly
    # defines several ``\newcommand`` values as 「該当なし」; those form values
    # are not the human-rights/legal rationale this check is about.
    for match in na_pattern.finditer(prose):
        local_context = prose[max(0, match.start() - 500):match.end()]
        if not any(term in local_context for term in ethics_context_terms):
            continue
        sentence_start = max(prose.rfind("。", 0, match.start()) + 1, 0)
        sentence_stop = prose.find("。", match.end())
        if sentence_stop < 0:
            sentence_stop = len(prose)
        sentence = prose[sentence_start:sentence_stop]
        if not any(term in sentence for term in rationale_terms):
            bare_na_matches.append((match, sentence))
    if bare_na_matches:
        first, sentence = bare_na_matches[0]
        raw_start = raw.find(first.group(0))
        add_risk(
            "not_applicable_without_rationale",
            max(0, raw_start),
            sentence,
            "「該当なし」とだけ書かれ、そう判断した根拠がない。",
            "人を対象としない数値解析のみである等、該当なしと判断した根拠を"
            "一文添える。この欄は例年審査委員からの指摘が非常に多い。",
            occurrence_count=len(bare_na_matches),
        )

    # The final-year-early-application box follows a different rule from the
    # human-rights/legal box. If the applicant is not eligible for this box,
    # the official instruction says to retain the page and leave every field
    # blank. Writing 「該当なし」 or an explanatory sentence is itself a form
    # violation; adding a rationale does not cure it.
    final_year_heading = "研究計画最終年度前年度応募を行う場合の記述事項"
    final_year_start = raw.find(final_year_heading)
    final_year_rule_checked = final_year_start >= 0
    final_year_blank_violations: list[dict] = []
    if final_year_rule_checked:
        next_section = re.search(
            r"\\section\*?\{",
            raw[final_year_start + len(final_year_heading):],
        )
        final_year_end = (
            final_year_start
            + len(final_year_heading)
            + next_section.start()
            if next_section is not None
            else len(raw)
        )
        final_year_raw = raw[final_year_start:final_year_end]
        final_year_raw = re.sub(r"(?m)^\s*%.*$", " ", final_year_raw)
        final_year_na_pattern = re.compile(
            r"該当\s*(?:は|事項は)?\s*"
            r"(?:なし|無し|ない|ありません|しない|していない|せず)|対象外"
        )
        final_year_macro_pattern = re.compile(
            r"\\(?:re)?newcommand\*?\{\\(?P<name>最終年度[^{}]*)\}"
            r"(?:\[[^\]]*\])?\{(?P<value>[^{}]*)\}"
        )
        for match in final_year_macro_pattern.finditer(final_year_raw):
            value = match.group("value").strip()
            if not value:
                continue
            if final_year_na_pattern.search(value) or re.fullmatch(
                r"[-‐‑‒–—―ー－]+", value
            ):
                final_year_blank_violations.append({
                    "kind": "field_value",
                    "field": match.group("name"),
                    "value": value,
                })

        # _prose_for_lint removes TeX definitions and the form's own sentence
        # 「該当しない場合は…空欄のまま」, leaving only applicant prose.
        final_year_prose = _prose_for_lint(final_year_raw)
        for match in final_year_na_pattern.finditer(final_year_prose):
            sentence_start = max(
                final_year_prose.rfind("。", 0, match.start()) + 1, 0
            )
            sentence_stop = final_year_prose.find("。", match.end())
            if sentence_stop < 0:
                sentence_stop = len(final_year_prose)
            final_year_blank_violations.append({
                "kind": "explanatory_text",
                "value": final_year_prose[sentence_start:sentence_stop].strip(),
            })

        if final_year_blank_violations:
            first = final_year_blank_violations[0]
            first_value = first["value"]
            add_risk(
                "final_year_non_applicant_field_not_blank",
                max(0, raw.find(first_value, final_year_start)),
                first_value,
                (
                    "研究計画最終年度前年度応募に該当しない場合の欄へ、"
                    "「該当なし」等の値又は説明文が記入されている。"
                ),
                (
                    "第4欄のページ、表、見出しは削除せず、研究種目名、課題番号、"
                    "課題名、研究期間、「当初研究計画及び研究成果」、"
                    "「前年度応募する理由」の全記述欄を空欄にする。"
                ),
                severity="HIGH",
                occurrence_count=len(final_year_blank_violations),
                violations=final_year_blank_violations[:8],
                official_rule=(
                    "該当しない場合は記述欄を削除することなく、空欄のまま提出すること。"
                ),
            )

    pub_triggers = ["研究遂行能力", "研究業績", "主要論文", "代表論文", "発表論文", "業績"]
    pub_hit = next((term for term in pub_triggers if term in prose), None)
    identifier_pattern = re.compile(
        r"(?:19|20)\d{2}|vol\.?\s*\d|no\.?\s*\d|pp\.?\s*\d|doi|DOI|"
        r"\d+\s*巻|\d+\s*号|第\d+巻",
        flags=re.IGNORECASE,
    )
    if pub_hit is not None and not identifier_pattern.search(raw):
        add_risk(
            "publication_not_identifiable",
            max(0, raw.find(pub_hit)),
            pub_hit,
            "業績への言及があるが、業績を特定できる情報(誌名・年・巻号等)がない。",
            "審査はresearchmapを研究者番号で参照する。調書に業績を書く場合は、"
            "特定するための十分な情報(著者、誌名、年等)を添える。",
        )

    overlap_triggers = [
        "応募中の研究費",
        "受入予定の研究費",
        "応募・受入",
        "応募中及び受入",
        "応募中および受入",
    ]
    overlap_hit = next((term for term in overlap_triggers if term in raw), None)
    if overlap_hit is not None:
        role_pattern = re.compile(
            r"(?:大学|研究所|機構|高等専門学校|高専)[^。\n]{0,15}"
            r"(?:教授|准教授|講師|助教|研究員)"
        )
        missing_parts = []
        if "相違" not in raw:
            missing_parts.append("本応募課題との相違点")
        if not re.search(r"応募(?:する)?理由", raw):
            missing_parts.append("応募する理由")
        if not role_pattern.search(raw):
            missing_parts.append("所属組織・役職(例: ○○大学教授)")
        if missing_parts:
            add_risk(
                "funding_overlap_format",
                max(0, raw.find(overlap_hit)),
                overlap_hit,
                "応募・受入状況欄に必要な記載要素が欠けている: "
                + "、".join(missing_parts),
                "国外資金・民間財団助成・受託/共同研究費も全て記載し、2件目以降は"
                "研究内容の相違点と応募する理由、所属組織・役職を添える。"
                "代表課題は分担者を含む金額、分担課題は自身の経費のみを書く。",
                missing_parts=missing_parts,
            )

    low = prose.lower()
    axis_hits = {
        axis: _contains_any(low, keywords)
        for axis, keywords in _KAKEN_REVIEW_CRITERIA_AXES.items()
    }
    # Length alone cannot tell a proposal body from a compact application form.
    # A 1,715-character 住友財団 form -- one 要旨 box, then keywords, amounts and
    # a funding plan -- matched none of the review vocabularies, and reporting a
    # missing 研究遂行能力 axis there is a finding its author would argue with,
    # because the form gave them nowhere to write it. A document that already
    # speaks two review vocabularies is a proposal body, and absent axes then
    # warrant review.
    full_draft = len(prose) >= 1500 and sum(bool(h) for h in axis_hits.values()) >= 2
    criteria_axis_results: dict[str, dict] = {}
    if full_draft:
        for axis, keywords in _KAKEN_REVIEW_CRITERIA_AXES.items():
            hits = axis_hits[axis]
            criteria_axis_results[axis] = {
                "ok": bool(hits),
                "matches": hits[:8],
                "keywords": keywords,
            }
        missing_criteria = [
            axis
            for axis, result in criteria_axis_results.items()
            if not result["ok"]
        ]
        if missing_criteria:
            add_risk(
                "review_criteria_axis_missing",
                0,
                "、".join(missing_criteria),
                "研究計画3要素(学術的重要性・方法の妥当性・遂行能力/環境)と"
                "別評定の国際性のうち、読み取れない軸がある: "
                + "、".join(missing_criteria),
                "各セクションがどの評定要素で読まれるかを意識し、研究計画3要素と"
                "国際性に対応する記述を置く。",
                missing_axes=missing_criteria,
            )
        emphasis_pattern = re.compile(
            r"下線|太字|ゴシック|アンダーライン|\\underline|\\textbf|"
            r"[図表]\s*\d|Fig\.?\s*\d"
        )
        if not emphasis_pattern.search(raw):
            add_risk(
                "no_emphasis_or_figures",
                0,
                prose[:80],
                "長い本文に、強調(下線・太字・ゴシック)や図表参照が見当たらない。",
                "審査委員は約1ヶ月で最大100件程度を読む。要点への下線・太字と"
                "図表で、一読で主張が追える構成にする。",
            )

    deductions = sum(2.0 if risk["severity"] == "HIGH" else 1.0 for risk in risks)
    applicable = bool(raw.strip())
    score = None if not applicable else max(0.0, round(10.0 - deductions, 1))
    comments = list(dict.fromkeys(risk["comment"] for risk in risks))
    recommendations = list(
        dict.fromkeys(risk["recommendation"] for risk in risks)
    )
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "comments": comments,
        "recommendations": recommendations,
        "full_draft_heuristics_applied": full_draft,
        "criteria_axis_results": criteria_axis_results,
        "final_year_blank_rule_checked": final_year_rule_checked,
        "final_year_blank_violations": final_year_blank_violations,
        "briefing_notes": list(_KAKEN_BRIEFING_NOTES),
        "target": (
            "a proposal a reviewer can judge on the three plan axes plus "
            "internationality at "
            "~100-proposals-per-month reading speed: monochrome-safe figures, "
            "identifiable publications, an explicit human-rights/legal box, "
            "and a complete funding-overlap box"
        ),
        "source": (
            "official KAKENHI B/C criteria plus in-house R9/FY2027 "
            "review-format briefing"
        ),
    }


def grant_writing_literature_gap_evidence_check(text: str) -> dict:
    """Check whether literature-survey evidence supports the claimed gap.

    A missing keyword or method name in a bounded corpus can help select
    examples, but it does not by itself establish field-wide non-adoption, a
    causal implementation barrier, or an academic gap.  This check looks for
    those inference jumps in a local sentence window.  It is not a systematic
    review or citation-quality validator.
    """
    text = _prose_for_lint(_read_text_if_path(text))
    sentences = [
        sentence.strip()
        for sentence in re.split(r"(?<=[。．!?！？])|\n+", text)
        if sentence.strip()
    ]

    keyword_groups = {
        "search_method": [
            "全文検索",
            "キーワード検索",
            "検索語",
            "同一語彙",
            "検出頁",
            "目視確認",
            "出現頻度",
            "ヒット件数",
            "full-text search",
            "keyword search",
        ],
        "corpus_scope": [
            "調査対象",
            "対象文献",
            "対象資料",
            "技術報告",
            "選定した文献",
            "限定コーパス",
            "bounded corpus",
            "selected literature",
        ],
        # An existential claim, as opposed to a report of having searched.
        # Derived from the corpus: every absence claim in eight real proposals
        # is phrased this way, and none is phrased as a search report, so the
        # non_detection vocabulary below matched nothing at all.
        "bare_absence": [
            "存在しない",
            "存在せず",
            "他に類を見な",
            "類を見ない",
            "前例がな",
            "報告例がな",
            "研究例がな",
            "事例がな",
            "皆無",
            "知られていない",
            "does not exist",
            "no such",
        ],
        "non_detection": [
            "確認できなかった",
            "確認されなかった",
            "見当たらなかった",
            "検出されなかった",
            "現れなかった",
            "ヒットしなかった",
            "記載がない",
            "言及がない",
            "採用例がない",
            "not found",
            "no occurrence",
            "was not identified",
        ],
        "field_scope": [
            "国内では",
            "日本では",
            "我が国では",
            "国内全",
            "当該分野",
            "分野全体",
            "コミュニティ全体",
            "海外では",
            "欧州では",
            "in japan",
            "field-wide",
            "community-wide",
        ],
        "adoption_claim": [
            "普及していない",
            "使われていない",
            "利用されていない",
            "採用されていない",
            "浸透していない",
            "未利用",
            "デファクト",
            "候補から外れ",
            "取り込まれない",
            "遅れている",
            "ガラパゴス",
            "not widely used",
            "not adopted",
        ],
        "causal_barrier": [
            "知識不足",
            "実装障壁",
            "言語障壁",
            "情報不足",
            "不足により",
            "不足から",
            "起因",
            "阻害",
            "妨げ",
            "implementation barrier",
            "language barrier",
        ],
        "academic_gap": [
            "研究障壁",
            "学術的空白",
            "研究上の空白",
            "未解決問題",
            "未解決の課題",
            "academic gap",
            "research barrier",
        ],
        "inference_link": [
            "示す",
            "裏付ける",
            "根拠",
            "意味する",
            "したがって",
            "ゆえに",
            "demonstrates",
            "therefore",
            "evidence of",
        ],
        "limited_role": [
            "背景補助",
            "補助証拠",
            "実証対象の選定",
            "実証対象選定",
            "事例選定",
            "予備調査",
            "本調査範囲",
            "網羅調査ではなく",
            "網羅的ではなく",
            "主根拠としない",
            "scope examples",
            "case selection",
        ],
    }
    count_pattern = re.compile(
        r"(?:全|計)?\s*\d+\s*(?:冊|報|件|編|誌|論文|文献|資料|報告)"
    )

    def keyword_hits(fragment: str, group: str) -> list[str]:
        low = fragment.lower()
        return [word for word in keyword_groups[group] if word.lower() in low]

    absence_indices = [
        index
        for index, sentence in enumerate(sentences)
        if keyword_hits(sentence, "non_detection")
        or keyword_hits(sentence, "bare_absence")
    ]
    candidate_windows = []
    unbacked: list[dict] = []
    for index in absence_indices:
        start = max(0, index - 4)
        stop = min(len(sentences), index + 5)
        window = "".join(sentences[start:stop])
        search_hits = keyword_hits(window, "search_method")
        corpus_hits = keyword_hits(window, "corpus_scope")
        count_hits = count_pattern.findall(window)
        if not search_hits and not corpus_hits and not count_hits:
            # An absence asserted with no account of how the applicant looked.
            # Four of eight real proposals do this -- 「統合的なマルチスケール
            # モデル縮約法が存在しない」, 「直接的な競合製品は存在しない」 --
            # and a reviewer who knows one counterexample loses the sentence
            # and some of the trust around it.
            bare = keyword_hits(sentences[index], "bare_absence")
            if bare:
                unbacked.append({
                    "sentence": index + 1,
                    "excerpt": sentences[index].strip()[:200],
                    "absence_hits": bare,
                })
            continue
        candidate_windows.append({
            "sentence": index + 1,
            "excerpt": window[:360],
            "search_hits": search_hits,
            "corpus_hits": corpus_hits + count_hits,
            "non_detection_hits": keyword_hits(window, "non_detection"),
            "field_scope_hits": keyword_hits(window, "field_scope"),
            "adoption_hits": keyword_hits(window, "adoption_claim"),
            "causal_hits": keyword_hits(window, "causal_barrier"),
            "academic_gap_hits": keyword_hits(window, "academic_gap"),
            "inference_hits": keyword_hits(window, "inference_link"),
            "limited_role_hits": keyword_hits(window, "limited_role"),
        })

    risk_specs = {
        "field_generalization": (
            "限定した文献群での非検出から、国内・分野全体の未普及や採用状況へ一般化している。",
            lambda item: bool(item["field_scope_hits"] and item["adoption_hits"]),
        ),
        "unsupported_causal_inference": (
            "語の非検出から、知識・言語・実装上の障壁という原因へ飛躍している。",
            lambda item: bool(item["causal_hits"] and item["inference_hits"]),
        ),
        "absence_as_academic_gap": (
            "限定コーパスでの不在を、未解決の研究障壁・学術的空白の主根拠にしている。",
            lambda item: bool(
                item["academic_gap_hits"] and item["inference_hits"]
            ),
        ),
    }
    risks = []
    for risk_type, (comment, predicate) in risk_specs.items():
        evidence = [item for item in candidate_windows if predicate(item)]
        if evidence:
            risks.append({
                "type": risk_type,
                "severity": "HIGH",
                "comment": comment,
                "evidence": evidence[:3],
            })

    if unbacked:
        risks.append({
            "type": "absence_claimed_without_search",
            "severity": "MEDIUM",
            "comment": (
                "「存在しない」と断定しているが、どう調べたかが書かれていない。"
            ),
            "evidence": unbacked[:3],
        })

    applicable = bool(candidate_windows or unbacked)
    score = None if not applicable else max(0.0, 10.0 - 3.0 * len(risks))
    rewrite_strategy = [
        (
            "断定した不在には、調べた範囲（検索語、対象、年）を一文添えるか、"
            "「知る限り」に落とす。審査者は反例を一つ知っていれば足りる。"
        ),
        "限定コーパスで語を確認できないことを、分野全体の未普及・優劣・研究障壁の主根拠にしない。",
        "文献調査は背景補助、候補選定、予備調査に位置づけ、対象範囲と限界を明示する。",
        "学術的空白は、既往研究が解いていない条件、理論、比較可能性、または検証可能な仮説として独立に述べる。",
        "普及実態を主張する場合は、検索式、選定基準、対象範囲、代替語、再現可能な集計を別途示す。",
    ]
    return {
        "applicable": applicable,
        "score": score,
        "risk_count": len(risks),
        "risks": risks,
        "candidate_windows": candidate_windows[:5],
        "unbacked_absence_claims": unbacked[:5],
        "comments": [risk["comment"] for risk in risks],
        "rewrite_strategy": rewrite_strategy if risks else [],
        "target": (
            "use bounded literature searches for background or case selection; "
            "establish the academic gap through an unresolved scientific "
            "condition, theory, comparison, or testable hypothesis"
        ),
        "source": "generic literature-gap evidence-scope check",
    }


def grant_writing_collaborative_integration_risk_check(text: str) -> dict:
    """Check recurring risks in collaborative software-integration proposals.

    The check is domain-neutral and applies only when a draft proposes coupling,
    integration, interoperability, or cross-organization reuse. It distinguishes
    the academic question from the implementation mechanism and checks costs on
    both sides of the interface, ecosystem positioning, scope control, negative
    results, team readiness, evaluation ethics, and asset provenance.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    applicable_hits = _contains_any(
        low,
        [
            "結合",
            "統合",
            "連携",
            "連成",
            "相互運用",
            "interoperability",
            "integration",
            "coupling",
        ],
    )

    axis_groups = {
        "academic_question_vs_mechanism": [
            ["学術的問い", "中心の問い", "研究上の問い", "仮説", "research question"],
            ["実装手段", "検証手段", "研究手段", "method", "mechanism"],
        ],
        "provider_and_reuse_cost": [
            ["所有者側", "提供者側", "開発元", "provider", "maintainer"],
            ["初期整備", "初期費用", "保守", "maintenance", "setup cost"],
            ["総負担", "再利用", "損益分岐", "reuse", "total cost"],
        ],
        "existing_ecosystem_boundary": [
            ["既存基盤", "既存規格", "既存oss", "existing framework", "existing standard"],
            ["置き換えず", "置換しない", "再利用", "補完", "boundary", "reuse"],
        ],
        "core_vs_optional_scope": [
            ["中核", "必達", "成立条件", "core", "required"],
            [
                "独立課題", "別課題", "発展候補", "条件付き", "optional",
                "exploratory",
                # A proposal that says "NVH等は波及効果とする" has made exactly
                # this split; the axis must not miss it on vocabulary alone.
                "波及効果", "対象外", "今後の展開", "将来展開",
            ],
        ],
        "negative_result_value": [
            ["結合不能", "適用境界", "不能理由", "反例", "不成立", "negative result"],
            ["成果", "判定", "同定", "明らか", "result", "outcome"],
        ],
        "team_readiness": [
            ["共著", "既往成果", "既発表", "共同遂行", "prior work", "track record"],
            ["担当", "役割", "責任", "role", "responsibility"],
            ["着手", "準備", "基礎", "既に", "ready", "readiness"],
        ],
        "evaluation_unit_and_ethics": [
            ["評価単位", "分析単位", "課題単位", "unit of analysis"],
            ["個人を評価", "個人の能力", "個人情報", "individual productivity"],
            ["倫理", "同意", "該当性", "ethics", "consent"],
        ],
        "asset_provenance_and_fallback": [
            ["権利", "保守主体", "所有者", "provenance", "maintainer"],
            ["参照実装", "公開ベンチマーク", "代替", "fallback", "reference implementation"],
        ],
    }

    axis_results = {}
    for axis, groups in axis_groups.items():
        group_hits = [_contains_any(low, group) for group in groups]
        axis_results[axis] = {
            "ok": all(bool(hits) for hits in group_hits),
            "group_matches": [hits[:8] for hits in group_hits],
            "groups": groups,
        }

    # Naming people is not the trigger; MEASURING them is. A sentence like
    # 「教員・学生が利用している」 says who uses a tool, and demanding an ethics
    # determination for it is a false positive. Terms that are inherently about
    # human-subject measurement stand alone; merely naming a person category
    # counts only next to a measurement verb in the same sentence.
    _MEASURED_PEOPLE_TERMS = [
        "工程時間",
        "手作業時間",
        "生産性",
        "被験者",
        "アンケート",
    ]
    _PERSON_CATEGORY_TERMS = ["学生", "若手", "参加者", "教員"]
    _MEASUREMENT_TERMS = [
        "評価する",
        "評価を",
        "計測",
        "測定",
        "記録し",
        "記録する",
        "比較する",
        "分析単位",
        "調査",
        "収集",
    ]
    people_process_hits = _contains_any(low, _MEASURED_PEOPLE_TERMS)
    if not people_process_hits:
        for sentence in re.split(r"(?<=[。．!?！？])", text):
            s_low = sentence.lower()
            if _contains_any(s_low, _PERSON_CATEGORY_TERMS) and _contains_any(
                s_low, _MEASUREMENT_TERMS
            ):
                people_process_hits = _contains_any(s_low, _PERSON_CATEGORY_TERMS)
                break
    if not people_process_hits:
        axis_results["evaluation_unit_and_ethics"].update(
            {"ok": True, "not_applicable": True}
        )

    if not applicable_hits:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": axis_results,
            "comments": [],
            "target": (
                "for collaborative integration proposals, separate the research "
                "question from tooling and test lifecycle cost, scope, negative "
                "results, ethics, and provenance"
            ),
            "source": "generic collaborative-integration risk check",
        }

    missing = [axis for axis, result in axis_results.items() if not result["ok"]]
    comments_by_axis = {
        "academic_question_vs_mechanism": (
            "製品・プロトコル名を問いにせず、学術的問いと実装手段を分ける。"
        ),
        "provider_and_reuse_cost": (
            "利用側だけでなく、提供者側の初期整備・保守を含む総負担を再利用回数に対して評価する。"
        ),
        "existing_ecosystem_boundary": (
            "既存規格・連携基盤が担う範囲を認め、置換せず再利用する範囲と研究上の空白を示す。"
        ),
        "core_vs_optional_scope": (
            "中核実証と条件付きの発展候補を分け、発展候補を必達成果の成立条件にしない。"
        ),
        "negative_result_value": (
            "順位不変、結合不能、反例等も、条件・原因・適用境界を同定できれば成果となる成功条件を置く。"
        ),
        "team_readiness": (
            "各担当者について、既往成果、利用可能資産、役割、着手可能性を一続きで示す。"
        ),
        "evaluation_unit_and_ethics": (
            "工程を比較する場合は個人でなく課題・成果物を評価単位とし、倫理該当性と同意手続を確認する。"
        ),
        "asset_provenance_and_fallback": (
            "固有資産の権利・保守主体を確認し、利用不能時の公開ベンチマークまたは参照実装を用意する。"
        ),
    }
    comments = [comments_by_axis[axis] for axis in missing]
    score = round(10.0 * (len(axis_groups) - len(missing)) / len(axis_groups), 1)
    return {
        "applicable": True,
        "score": score,
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": axis_results,
        "people_process_hits": people_process_hits,
        "comments": comments,
        "target": (
            "a collaborative integration proposal with an academic question distinct "
            "from tooling, full lifecycle cost, ecosystem boundaries, controlled scope, "
            "valuable negative results, team readiness, ethical evaluation, and fallback assets"
        ),
        "source": "generic collaborative-integration risk check",
    }


_BUDGET_COST_TOKENS = (
    "円", "費", "単価", "積算", "計上", "予算", "経費", "見積", "金額", "内訳",
)


def _mentions_cost_nearby(text: str, keyword: str) -> bool:
    """True when a keyword is costed in the sentence that mentions it.

    Bare presence is not evidence of a budget rationale: a methods section
    says 評価 and AI constantly without costing anything. The sentence is the
    right window -- a cost token one sentence away belongs to a different
    claim -- so this is what separates "mentions it" from "budgets for it".
    """
    needle = keyword.lower()
    for sentence in re.split(r"(?<=[。．!?！？\n])", text):
        if needle in sentence.lower() and any(
            token in sentence for token in _BUDGET_COST_TOKENS
        ):
            return True
    return False


def grant_writing_budget_alignment_check(text: str) -> dict:
    """Check that budget items are tied to verification and implementation.

    The check is optional. A proposal section that carries no budget content
    at all -- a research-plan or feasibility section, for example -- is not a
    thin budget; it is the wrong document for this question, so the check
    reports ``applicable: False`` rather than a low score. When it does apply,
    a cost keyword only counts if a money token sits next to it.
    """
    text = _read_text_if_path(text)
    low = text.lower()
    axes = {
        "ai_agent_costs": ["claude", "codex", "fable", "生成ai", "llm", "ai"],
        "compute_resources": ["mdx", "計算資源", "gpu", "クラウド", "サーバ"],
        "poc_experiment": ["試作", "基板", "部品", "消耗品", "計測", "評価"],
        "dissemination": ["旅費", "発表", "技術プレゼン", "ワークショップ", "報告"],
        "near_ceiling_strategy": [
            "上限",
            "助成上限",
            "上限額",
            "限度額",
            "満額",
            "ほぼ上限",
            "上限いっぱい",
            "上限近く",
        ],
        "itemized_calculation": [
            "内訳",
            "単価",
            "数量",
            "月数",
            "回数",
            "年度配分",
            "見積",
            "積算",
            "算出",
            "根拠",
            "税込",
            "税抜",
        ],
        "pricing_provenance": [
            "公式料金",
            "料金表",
            "見積書",
            "参照日",
            "改定日",
            "料金年度",
            "最低購入",
            "有効期限",
            "為替",
            "端数処理",
        ],
    }
    budget_markers = [
        "予算", "経費", "費目", "直接経費", "設備備品", "消耗品費", "旅費",
        "人件費", "謝金", "その他", "千円", "万円", "単価", "積算", "計上",
    ]
    money_pattern = re.compile(r"\d[\d,\.]*\s*(?:千円|万円|円|kJPY|JPY)")
    marker_hits = _contains_any(low, budget_markers)
    has_money = bool(money_pattern.search(text))
    if not has_money and len(marker_hits) < 2:
        return {
            "applicable": False,
            "score": None,
            "missing_count": 0,
            "missing_axes": [],
            "axis_results": {},
            "comments": [],
            "budget_marker_hits": marker_hits,
            "budget_policy": _BUDGET_POLICY,
            "target": (
                "budget rationale is judged only where budget content exists; "
                "a research-plan or feasibility section carries none by design"
            ),
        }

    results = {}
    missing = []
    for axis, keywords in axes.items():
        if axis == "pricing_provenance":
            required_groups = {
                "source": ["公式料金", "料金表", "見積書", "https://", "url"],
                "vintage": ["参照日", "改定日", "料金年度", "年度料金", "2025年度", "2026年度"],
                "accounting": [
                    "税込",
                    "税抜",
                    "最低購入",
                    "有効期限",
                    "為替",
                    "端数処理",
                ],
            }
            group_matches = {
                group: _contains_any(low, group_keywords)
                for group, group_keywords in required_groups.items()
            }
            matches = [match for values in group_matches.values() for match in values]
            ok = all(group_matches.values())
            results[axis] = {
                "ok": ok,
                "matches": matches,
                "keywords": keywords,
                "required_groups": required_groups,
                "group_matches": group_matches,
            }
        elif axis in {
            "ai_agent_costs",
            "compute_resources",
            "poc_experiment",
            "dissemination",
        }:
            # A resource keyword is budget evidence only next to a money token.
            # Otherwise a methods section "passes" on 評価 / AI it never costs.
            matches = [
                kw for kw in keywords
                if kw.lower() in low and _mentions_cost_nearby(text, kw)
            ]
            ok = bool(matches)
            results[axis] = {
                "ok": ok,
                "matches": matches,
                "keywords": keywords,
                "requires_cost_context": True,
            }
        else:
            matches = _contains_any(low, keywords)
            ok = bool(matches)
            results[axis] = {"ok": ok, "matches": matches, "keywords": keywords}
        if not ok:
            missing.append(axis)
    score = round(10.0 * (len(axes) - len(missing)) / len(axes), 1)
    comments = [
        _BUDGET_AXIS_COMMENTS.get(axis, f"Budget rationale missing or thin: {axis}")
        for axis in missing
    ]
    if "効率化" in text and "検証" not in text and "PoC" not in text:
        comments.append("AI費用が一般的な効率化に見える。検証ループの実行経費として説明する。")
        score = max(0.0, score - 1.0)
    return {
        "applicable": True,
        "score": round(score, 1),
        "missing_count": len(missing),
        "missing_axes": missing,
        "axis_results": results,
        "comments": comments,
        "budget_marker_hits": marker_hits,
        "budget_policy": _BUDGET_POLICY,
        "target": (
            "every major cost maps to AI/tool execution, compute, PoC, or dissemination; "
            "the requested amount may be close to the ceiling when itemized, justified, "
            "and traceable to dated official prices or quotations"
        ),
    }


_ARGUMENT_EVIDENCE_ROLES = {
    "central_question": {
        "description": "解明対象となる中心の問いまたは主張",
        "terms": ["中心の問い", "学術的問い", "研究上の問い", "何を明らか", "問う"],
    },
    "prior_gap": {
        "description": "既往研究との限定された対比と未解決点",
        "terms": [
            "未解決", "明らかでない", "体系化されていない", "検証されていない",
            "限界", "研究障壁", "一方", "これに対し",
        ],
    },
    "method_operation": {
        "description": "明示した入力・モデルに対して行う操作",
        "terms": [
            "比較する", "構成する", "結合する", "射影", "算出", "測定する",
            "検証する", "導入する", "評価する", "同定する", "実装する",
        ],
    },
    "decision_rule": {
        "description": "主張を支持または棄却できる観測量・判定規則",
        "terms": [
            "判定", "許容差", "閾値", "順位", "一致", "不一致", "採否",
            "反例", "成立条件", "適用境界", "信頼区間",
        ],
    },
    "knowledge_output": {
        "description": "ソフトウェア成果物を越えて得る分野知",
        "terms": [
            "設計則", "選択則", "適用条件", "成立条件", "適用境界", "知見",
            "体系化", "指針", "条件を明らか",
        ],
    },
    "preliminary_evidence": {
        "description": "提案研究の実行可能性を支える完了済みの根拠",
        "terms": [
            "実装した", "完了した", "確認した", "再現した", "発表した",
            "発表予定", "採択", "共著", "予備実証", "既往成果",
        ],
    },
    "preparation_plan_link": {
        "description": "完了済みの準備実績を開始可能な研究項目へ結ぶ文",
        "terms": [
            "実績により", "成果により", "これにより研究", "したがって本研究",
            "から着手できる", "を開始できる", "実行できること", "遂行できること",
            "を担保する", "準備が整っている", "基盤が既に整っている",
        ],
    },
    "responsibility": {
        "description": "研究の成立に不可欠な能力を担う構成員",
        "terms": ["担当", "担う", "責任", "研究代表者", "研究分担者", "役割"],
    },
    "negative_result": {
        "description": "結合や仮説が成立しない場合にも得られる知識",
        "terms": [
            "不成立", "結合不能", "不能理由", "反例", "適用境界", "失敗例",
            "成立しない", "否定された場合",
        ],
    },
}


def _argument_segments(text: str) -> list[str]:
    prose = _prose_for_lint(_read_text_if_path(text))
    segments = [
        re.sub(r"\s+", " ", item).strip()
        for item in re.split(r"(?<=[。！？!?])|\n", prose)
    ]
    return [item for item in segments if 12 <= len(item) <= 500]


def grant_writing_argument_evidence_map(text: str) -> dict:
    """Map argument roles to excerpts without scoring scientific validity.

    This bridges deterministic lint and an LLM/human close read. Lexical hits
    only locate candidate evidence; they do not prove that the question is
    original, the method is valid, or the evidence supports the claim.
    """
    segments = _argument_segments(text)
    evidence_map = {}
    for role, spec in _ARGUMENT_EVIDENCE_ROLES.items():
        matched = [
            segment for segment in segments
            if any(term.lower() in segment.lower() for term in spec["terms"])
        ]
        hits = sorted({
            term
            for term in spec["terms"]
            if any(term.lower() in segment.lower() for segment in matched)
        })
        evidence_map[role] = {
            "description": spec["description"],
            "candidate_count": len(matched),
            "terms_hit": hits,
            "excerpts": matched[:3],
        }

    untraced = [
        role for role, result in evidence_map.items()
        if result["candidate_count"] == 0
    ]
    prompts = [
        (
            "概要、目的、方法、年度計画で、中心の問いを決める名詞が同じか。"
        ),
        (
            "各問いについて、何にどの操作を行い、どの観測量の変化を見るか。"
        ),
        (
            "判定規則は主張を棄却できるか。どの結果でも成功になる設計ではないか。"
        ),
        (
            "最終成果はリポジトリやソフトだけでなく、条件・境界・設計則などの分野知か。"
        ),
        (
            "完了済みの根拠と担当者は、特に新規部分を含む必須能力を全て覆うか。"
        ),
        (
            "各準備実績は、どの研究項目を直ちに開始・遂行できる根拠かを明示しているか。"
        ),
    ]
    if (
        evidence_map["preliminary_evidence"]["candidate_count"] > 0
        and evidence_map["preparation_plan_link"]["candidate_count"] == 0
    ):
        prompts.insert(
            0,
            (
                "準備実績の候補はあるが、研究項目への橋渡し文を確認できない。"
                "実績ごとに、何を開始・遂行できる根拠かを対応付ける。"
            ),
        )
    if untraced:
        prompts.insert(
            0,
            "語彙上の候補を確認できない役割: " + ", ".join(untraced)
            + "。該当節を通読する。ここでの未検出は欠陥を意味しない。",
        )

    return {
        "applicable": bool(segments),
        "evidence_map": evidence_map,
        "untraced_roles": untraced,
        "manual_review_prompts": prompts,
        "warning": (
            "候補文は通読の索引であり、論理的妥当性の証拠でも点数でもない。"
        ),
        "source": "反復的な申請書レビューから構成した論証追跡マップ",
    }


def grant_writing_recommendation_letter_template(
    program: str = "kddi_digital",
    applicant: str = "菅原賢悟准教授",
) -> str:
    """Return a one-page recommendation-letter draft template."""
    if program != "kddi_digital":
        return (
            f"本学所属の{applicant}による助成申請を推薦いたします。\n\n"
            "同氏は、研究課題の遂行に必要な専門性、研究実績、研究環境を備えており、"
            "本助成によって得られる成果は当該分野の発展に資するものです。"
            "本学としても、研究実施に必要な環境整備および研究活動を支援いたします。\n\n"
            "以上の理由により、本申請を推薦いたします。"
        )
    return (
        f"本学理工学部 電気電子通信工学科 {applicant}による、"
        "KDDI財団デジタルイノベーション社会実装助成への申請を推薦いたします。\n\n"
        "同氏は、三菱電機におけるPCB EMCを含む電磁気設計・CAE実務、"
        "本学における電磁界解析・熱解析研究、ならびにRadia、NGSolve、LTspice等を"
        "MCPサーバを介してAI駆動で扱う研究開発環境の構築実績を有しています。"
        "本申請は、地域製造業における回路・電磁界・熱協調設計の属人性を低減し、"
        "PoC基板評価と公開可能な設計支援基盤を通じて社会実装へ接続するものです。\n\n"
        "本学は、同氏の研究室が有する実験・評価環境、学生を含む研究実施体制、"
        "および外部連携活動を支援し、本助成期間中の計画遂行を後押しします。"
        "以上の理由により、本申請を推薦いたします。"
    )


# The test that decides which list a check belongs in: can it point at a
# place in the text and say what is wrong there, such that the author agrees
# without argument? 「この文は91字」 and 「概要は境界、本文は条件」 pass.
# 「公開成果が3つで4つ未満」 does not -- that is an opinion wearing the
# clothes of a measurement.
_DETECTOR_TOOLS = frozenset({
    "sentence",
    "bedrock",
    "weak",
    "claim",
    "vague",
    "format",
    "persuasion",
    "vocabulary",
    "abstraction",
    "literature",
    "residue",
    "narrative",
    "originality",
    "international",
    "irreplaceable",
    "pages",
    "capability",
})

_DETECTOR_RESULT_KEYS = frozenset({
    "sentence",
    "bedrock",
    "weak",
    "central_claim_consistency",
    "vague_claim_verb",
    "kaken_review_format",
    "persuasion_quality",
    "reviewer_vocabulary",
    "named_software_abstraction",
    "literature_gap_evidence",
    "template_residue",
    "budget_narrative",
    "question_originality",
    "international_standing",
    "collaboration_irreplaceability",
    "page_limit",
    "capability_responsibility",
})


def _find_compiled_pdf(text_or_path: str, pdf: str) -> pathlib.Path | None:
    """Locate the compiled proposal whose page allowances should be checked.

    An explicit path always wins. Otherwise the PDF is only inferred when the
    source directory holds exactly one of them, because guessing which of
    several PDFs is the submission would report page counts for the wrong
    document.
    """
    if pdf:
        candidate = pathlib.Path(pdf)
        if not candidate.is_file():
            raise FileNotFoundError(f"compiled proposal not found: {pdf}")
        return candidate

    # Proposal text must never reach the filesystem: statting it once per call
    # is a network round trip on a NAS-hosted tree, and the suite calls this on
    # every report.
    if len(text_or_path) > 260 or "\n" in text_or_path:
        return None
    source = pathlib.Path(text_or_path)
    if source.suffix.lower() not in {".md", ".tex", ".txt"} or not source.is_file():
        return None
    candidates = sorted(source.parent.glob("*.pdf"))
    return candidates[0] if len(candidates) == 1 else None


def grant_writing_health_report(
    text_or_path: str,
    program: str = "generic",
    skip: str = "",
    pdf: str = "",
) -> dict:
    """Integrated grant-writing health report.

    Args:
        text_or_path: Proposal text or an existing .md/.tex/.txt path.
        program: ``generic``, ``kaken_generic``, ``kddi_digital``, or
            ``kaken_oss``. Use ``kaken_generic`` for ordinary KAKENHI drafts;
            ``kaken_oss`` adds checks specific to the current OSS-platform
            proposal.
        skip: comma-separated tool ids to skip, e.g. ``sentence,literature``.
        pdf: compiled proposal to check page allowances against. When omitted
            and the source is a path, a single sibling PDF is used.
    """
    _validate_program(program)
    text = _read_text_if_path(text_or_path)
    skip_set = {s.strip().lower() for s in skip.split(",") if s.strip()}
    valid_skip_ids = {
        "abstraction", "argument_map", "bedrock", "budget", "capability",
        "claim", "domain", "focus", "format", "integration", "international",
        "irreplaceable", "kaken", "kddi", "literature", "metric", "narrative",
        "originality", "pages", "persuasion", "pilot", "residue", "scale",
        "japanese", "readability", "momentum", "sections", "sentence", "vague",
        "vocabulary", "weak",
    }
    unknown_skip_ids = sorted(skip_set - valid_skip_ids)
    if unknown_skip_ids:
        raise ValueError(
            "unknown grant-writing skip id(s): " + ", ".join(unknown_skip_ids)
        )

    detailed_results: dict[str, dict] = {}
    detailed_scores: dict[str, float] = {}
    priority_issues: list[dict] = []

    if "sections" not in skip_set:
        sections = grant_writing_section_presence(text, program=program)
        detailed_results["sections"] = sections
        detailed_scores["sections"] = sections["score"]
        for axis in sections["missing_axes"]:
            priority_issues.append({
                "tool": "sections",
                "name": "section_presence",
                "severity": "HIGH",
                "score": sections["score"],
                "comments": [f"Add reviewer-visible content for: {axis}"],
            })

    if program == "kddi_digital" and "kddi" not in skip_set:
        kddi = grant_writing_kddi_digital_check(text)
        detailed_results["kddi_digital"] = kddi
        detailed_scores["kddi_digital"] = kddi["score"]
        if kddi["comments"]:
            priority_issues.append({
                "tool": "kddi",
                "name": "kddi_digital_check",
                "severity": _severity_from_score(kddi["score"]),
                "score": kddi["score"],
                "comments": kddi["comments"][:5],
            })

    if (
        program == "kddi_digital"
        and "focus" not in skip_set
        and _contains_any(text.lower(), _POWER_ELECTRONICS_FOCUS_TRIGGERS)
    ):
        focus = grant_writing_kddi_power_electronics_focus_check(text)
        detailed_results["power_electronics_focus"] = focus
        detailed_scores["power_electronics_focus"] = focus["score"]
        if focus["comments"]:
            priority_issues.append({
                "tool": "focus",
                "name": "kddi_power_electronics_focus_check",
                "severity": _severity_from_score(focus["score"]),
                "score": focus["score"],
                "comments": focus["comments"][:5],
            })

    if program in {"kaken_oss", "kaken_oss_platform"} and "kaken" not in skip_set:
        kaken = grant_writing_kaken_oss_platform_check(text)
        detailed_results["kaken_oss_platform"] = kaken
        detailed_scores["kaken_oss_platform"] = kaken["score"]
        if kaken["comments"]:
            priority_issues.append({
                "tool": "kaken",
                "name": "kaken_oss_platform_check",
                "severity": _severity_from_score(kaken["score"]),
                "score": kaken["score"],
                "comments": kaken["comments"][:5],
            })

    if (
        program in {"kaken_oss", "kaken_oss_platform"}
        and "abstraction" not in skip_set
    ):
        abstraction = grant_writing_named_software_abstraction_check(text)
        detailed_results["named_software_abstraction"] = abstraction
        if abstraction["applicable"]:
            detailed_scores["named_software_abstraction"] = abstraction["score"]
            if abstraction["risks"]:
                priority_issues.append({
                    "tool": "abstraction",
                    "name": "named_software_abstraction_check",
                    "severity": "MEDIUM",
                    "score": abstraction["score"],
                    "comments": abstraction["comments"][:5],
                })

    if (
        program in {"kaken_oss", "kaken_oss_platform"}
        and "vocabulary" not in skip_set
    ):
        vocabulary = grant_writing_reviewer_vocabulary_check(text)
        detailed_results["reviewer_vocabulary"] = vocabulary
        if vocabulary["applicable"]:
            detailed_scores["reviewer_vocabulary"] = vocabulary["score"]
            if vocabulary["risks"]:
                priority_issues.append({
                    "tool": "vocabulary",
                    "name": "reviewer_vocabulary_check",
                    "severity": "MEDIUM",
                    "score": vocabulary["score"],
                    "comments": vocabulary["comments"][:5],
                })

    if "persuasion" not in skip_set:
        persuasion = grant_writing_persuasion_quality_check(text)
        detailed_results["persuasion_quality"] = persuasion
        if persuasion["applicable"]:
            detailed_scores["persuasion_quality"] = persuasion["score"]
            if persuasion["risks"]:
                priority_issues.append({
                    "tool": "persuasion",
                    "name": "persuasion_quality_check",
                    "severity": _severity_from_score(persuasion["score"]),
                    "score": persuasion["score"],
                    "comments": persuasion["comments"][:5],
                })

    if "narrative" not in skip_set:
        narrative = grant_writing_budget_narrative_check(text)
        detailed_results["budget_narrative"] = narrative
        if narrative["applicable"]:
            detailed_scores["budget_narrative"] = narrative["score"]
            if narrative["risks"]:
                priority_issues.append({
                    "tool": "narrative",
                    "name": "budget_narrative_check",
                    "severity": max(
                        (r["severity"] for r in narrative["risks"]),
                        key=lambda x: {"HIGH": 2, "MEDIUM": 1, "LOW": 0}[x],
                    ),
                    "score": narrative["score"],
                    "comments": narrative["comments"][:5],
                })

    if "residue" not in skip_set:
        residue = grant_writing_template_residue_check(text)
        detailed_results["template_residue"] = residue
        if residue["applicable"]:
            detailed_scores["template_residue"] = residue["score"]
            if residue["risks"]:
                priority_issues.append({
                    "tool": "residue",
                    "name": "template_residue_check",
                    "severity": max(
                        (r["severity"] for r in residue["risks"]),
                        key=lambda s: {"HIGH": 2, "MEDIUM": 1, "LOW": 0}[s],
                    ),
                    "score": residue["score"],
                    "comments": residue["comments"][:5],
                })

    if "vague" not in skip_set:
        vague = grant_writing_vague_claim_verb_check(text)
        detailed_results["vague_claim_verb"] = vague
        if vague["applicable"]:
            detailed_scores["vague_claim_verb"] = vague["score"]
            if vague["risks"]:
                priority_issues.append({
                    "tool": "vague",
                    "name": "vague_claim_verb_check",
                    "severity": "MEDIUM",
                    "score": vague["score"],
                    "comments": vague["comments"][:5],
                })

    if "pages" not in skip_set:
        compiled = _find_compiled_pdf(text_or_path, pdf)
        if compiled is not None:
            limits = grant_writing_page_limit_check(str(compiled))
            detailed_results["page_limit"] = limits
            if limits["applicable"]:
                detailed_scores["page_limit"] = limits["score"]
                if limits["risks"]:
                    priority_issues.append({
                        "tool": "pages",
                        "name": "page_limit_check",
                        "severity": max(
                            (r["severity"] for r in limits["risks"]),
                            key=lambda x: {
                                "CRITICAL": 3, "HIGH": 2, "MEDIUM": 1, "LOW": 0,
                            }[x],
                        ),
                        "score": limits["score"],
                        "comments": limits["comments"][:5],
                    })

    if "capability" not in skip_set:
        capability = grant_writing_capability_responsibility_check(text)
        detailed_results["capability_responsibility"] = capability
        if capability["applicable"]:
            detailed_scores["capability_responsibility"] = capability["score"]
            if capability["risks"]:
                priority_issues.append({
                    "tool": "capability",
                    "name": "capability_responsibility_check",
                    "severity": "HIGH",
                    "score": capability["score"],
                    "comments": capability["comments"][:5],
                })

    if "irreplaceable" not in skip_set:
        irrep = grant_writing_collaboration_irreplaceability_check(text)
        detailed_results["collaboration_irreplaceability"] = irrep
        if irrep["applicable"]:
            detailed_scores["collaboration_irreplaceability"] = irrep["score"]
            if irrep["risks"]:
                priority_issues.append({
                    "tool": "irreplaceable",
                    "name": "collaboration_irreplaceability_check",
                    "severity": max(
                        (r["severity"] for r in irrep["risks"]),
                        key=lambda x: {"HIGH": 2, "MEDIUM": 1, "LOW": 0}[x],
                    ),
                    "score": irrep["score"],
                    "comments": irrep["comments"][:5],
                })

    if "international" not in skip_set:
        intl = grant_writing_international_standing_check(text)
        detailed_results["international_standing"] = intl
        if intl["applicable"]:
            detailed_scores["international_standing"] = intl["score"]
            if intl["risks"]:
                priority_issues.append({
                    "tool": "international",
                    "name": "international_standing_check",
                    "severity": max(
                        (r["severity"] for r in intl["risks"]),
                        key=lambda x: {"HIGH": 2, "MEDIUM": 1, "LOW": 0}[x],
                    ),
                    "score": intl["score"],
                    "comments": intl["comments"][:5],
                })

    if "originality" not in skip_set:
        originality = grant_writing_question_originality_check(text)
        detailed_results["question_originality"] = originality
        if originality["applicable"]:
            detailed_scores["question_originality"] = originality["score"]
            if originality["risks"]:
                priority_issues.append({
                    "tool": "originality",
                    "name": "question_originality_check",
                    "severity": max(
                        (r["severity"] for r in originality["risks"]),
                        key=lambda x: {"HIGH": 2, "MEDIUM": 1, "LOW": 0}[x],
                    ),
                    "score": originality["score"],
                    "comments": originality["comments"][:5],
                })

    if "claim" not in skip_set:
        claim = grant_writing_central_claim_consistency_check(text)
        detailed_results["central_claim_consistency"] = claim
        if claim["applicable"]:
            detailed_scores["central_claim_consistency"] = claim["score"]
            if claim["risks"]:
                priority_issues.append({
                    "tool": "claim",
                    "name": "central_claim_consistency_check",
                    "severity": max(
                        (r["severity"] for r in claim["risks"]),
                        key=lambda s: {"HIGH": 2, "MEDIUM": 1, "LOW": 0}[s],
                    ),
                    "score": claim["score"],
                    "comments": claim["comments"][:5],
                })

    if "format" not in skip_set:
        review_format = grant_writing_kaken_review_format_check(text)
        detailed_results["kaken_review_format"] = review_format
        if review_format["applicable"]:
            detailed_scores["kaken_review_format"] = review_format["score"]
            if review_format["risks"]:
                priority_issues.append({
                    "tool": "format",
                    "name": "kaken_review_format_check",
                    "severity": _severity_from_score(review_format["score"]),
                    "score": review_format["score"],
                    "comments": review_format["comments"][:5],
                })

    if "scale" not in skip_set:
        scale = grant_writing_internal_evidence_to_external_scale_check(text)
        detailed_results["internal_to_external_scale"] = scale
        if scale["applicable"]:
            detailed_scores["internal_to_external_scale"] = scale["score"]
            if scale["comments"]:
                priority_issues.append({
                    "tool": "scale",
                    "name": "internal_evidence_to_external_scale_check",
                    "severity": _severity_from_score(scale["score"]),
                    "score": scale["score"],
                    "comments": scale["comments"][:5],
                })

    if "domain" not in skip_set:
        domain = grant_writing_domain_outcome_chain_check(text)
        detailed_results["domain_outcome_chain"] = domain
        if domain["applicable"]:
            detailed_scores["domain_outcome_chain"] = domain["score"]
            if domain["comments"]:
                priority_issues.append({
                    "tool": "domain",
                    "name": "domain_outcome_chain_check",
                    "severity": _severity_from_score(domain["score"]),
                    "score": domain["score"],
                    "comments": domain["comments"][:5],
                })

    if "metric" not in skip_set:
        metric = grant_writing_derived_metric_validation_check(text)
        detailed_results["derived_metric_validation"] = metric
        if metric["applicable"]:
            detailed_scores["derived_metric_validation"] = metric["score"]
            if metric["comments"]:
                priority_issues.append({
                    "tool": "metric",
                    "name": "derived_metric_validation_check",
                    "severity": _severity_from_score(metric["score"]),
                    "score": metric["score"],
                    "comments": metric["comments"][:5],
                })

    if "pilot" not in skip_set:
        pilot = grant_writing_cross_organization_pilot_check(text)
        detailed_results["cross_organization_pilot"] = pilot
        if pilot["applicable"]:
            detailed_scores["cross_organization_pilot"] = pilot["score"]
            if pilot["comments"]:
                priority_issues.append({
                    "tool": "pilot",
                    "name": "cross_organization_pilot_check",
                    "severity": _severity_from_score(pilot["score"]),
                    "score": pilot["score"],
                    "comments": pilot["comments"][:5],
                })

    if "literature" not in skip_set:
        literature = grant_writing_literature_gap_evidence_check(text)
        detailed_results["literature_gap_evidence"] = literature
        if literature["applicable"]:
            detailed_scores["literature_gap_evidence"] = literature["score"]
            if literature["risks"]:
                priority_issues.append({
                    "tool": "literature",
                    "name": "literature_gap_evidence_check",
                    "severity": "HIGH",
                    "score": literature["score"],
                    "comments": literature["comments"][:5],
                })

    if "integration" not in skip_set:
        integration = grant_writing_collaborative_integration_risk_check(text)
        detailed_results["collaborative_integration_risk"] = integration
        if integration["applicable"]:
            detailed_scores["collaborative_integration_risk"] = integration["score"]
            if integration["comments"]:
                priority_issues.append({
                    "tool": "integration",
                    "name": "collaborative_integration_risk_check",
                    "severity": _severity_from_score(integration["score"]),
                    "score": integration["score"],
                    "comments": integration["comments"][:5],
                })

    if "sentence" not in skip_set:
        sent = grant_writing_analyze_sentences(text)
        detailed_results["sentence"] = sent
        if "error" not in sent:
            over = sent["over_threshold_count"]
            score = max(0.0, round(10.0 - min(over, 8) * 1.0, 1))
            detailed_scores["sentence"] = score
            if over:
                priority_issues.append({
                    "tool": "sentence",
                    "name": "analyze_sentences",
                    "severity": "MEDIUM" if over < 4 else "HIGH",
                    "score": score,
                    "comments": [f"{over} sentence(s) exceed the threshold."],
                })

    if "readability" not in skip_set:
        readability = grant_writing_adjacent_reviewer_readability_check(text)
        detailed_results["adjacent_reviewer_readability"] = readability

    if "japanese" not in skip_set:
        japanese = grant_writing_japanese_readability_score(
            text,
            document_type="grant_proposal",
        )
        detailed_results["japanese_readability"] = japanese
        if japanese["applicable"]:
            detailed_scores["japanese_readability"] = japanese["score"] / 10

    if "momentum" not in skip_set:
        momentum = grant_writing_reviewer_momentum_check(text)
        detailed_results["reviewer_momentum"] = momentum

    if "weak" not in skip_set:
        weak = grant_writing_count_weak_expressions(text)
        detailed_results["weak"] = weak
        score = max(0.0, round(10.0 - min(weak["total_weak_expressions"], 10), 1))
        detailed_scores["weak"] = score
        if weak["total_weak_expressions"]:
            priority_issues.append({
                "tool": "weak",
                "name": "count_weak_expressions",
                "severity": "MEDIUM" if weak["total_weak_expressions"] < 5 else "HIGH",
                "score": score,
                "comments": [f"Weak / non-committal expressions: {weak['by_pattern']}"],
            })

    if "bedrock" not in skip_set:
        bedrock = grant_writing_lint_bedrock(_prose_for_lint(text))
        detailed_results["bedrock"] = bedrock
        issue_count = bedrock.get("issue_count", 0)
        score = max(0.0, round(10.0 - issue_count * 1.5, 1))
        detailed_scores["bedrock"] = score
        if issue_count:
            priority_issues.append({
                "tool": "bedrock",
                "name": "lint_bedrock",
                "severity": _severity_from_score(score),
                "score": score,
                "comments": [i.get("rule", "issue") for i in bedrock.get("issues", [])[:5]],
            })

    if "budget" not in skip_set:
        budget = grant_writing_budget_alignment_check(text)
        detailed_results["budget"] = budget
        # A section with no budget content is not a thin budget. Scoring it
        # would drag the overall score down and raise a HIGH issue about
        # itemization the section cannot carry.
        if budget.get("applicable", True):
            detailed_scores["budget"] = budget["score"]
            if budget["comments"]:
                priority_issues.append({
                    "tool": "budget",
                    "name": "budget_alignment_check",
                    "severity": _severity_from_score(budget["score"]),
                    "score": budget["score"],
                    "comments": budget["comments"][:5],
                })

    if "argument_map" not in skip_set:
        argument_map = grant_writing_argument_evidence_map(text)
        detailed_results["argument_evidence_map"] = argument_map

    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2, "UNKNOWN": 3, "LOW": 4}

    # Split what the suite produces into the two kinds it actually contains.
    # A DETECTOR points at a place in the text and says what is wrong there;
    # the author can check it and fix it without argument. A QUESTION asks
    # whether the proposal covers a topic; keyword presence cannot answer
    # that, so it is surfaced as a prompt and never scored. Averaging the two
    # into one number was what made the old overall_score meaningless: a
    # draft with a fatal inconsistency could score 10 while a clean section
    # scored 8.6 on question noise.
    findings = [i for i in priority_issues if i["tool"] in _DETECTOR_TOOLS]
    questions = [i for i in priority_issues if i["tool"] not in _DETECTOR_TOOLS]
    findings.sort(
        key=lambda x: (sev_rank.get(x["severity"], 99), -1 * (x.get("score") or 0))
    )
    questions.sort(key=lambda x: x["name"])
    for q in questions:
        q["kind"] = "question"
        q.pop("severity", None)
        q.pop("score", None)

    defect_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0, "CRITICAL": 0}
    for item in findings:
        defect_counts[item["severity"]] = defect_counts.get(item["severity"], 0) + 1
    defect_counts["total"] = len(findings)

    detector_scores = {
        name: value
        for name, value in detailed_scores.items()
        if name in _DETECTOR_RESULT_KEYS
    }
    defect_score = (
        round(sum(detector_scores.values()) / len(detector_scores), 1)
        if detector_scores
        else 10.0
    )

    if defect_counts["total"] == 0:
        summary = (
            "No located defects. This says the mechanics are clean; it does "
            "not say the argument holds."
        )
    else:
        summary = (
            f"{defect_counts['total']} located defect(s): "
            f"{defect_counts['CRITICAL']} critical, {defect_counts['HIGH']} high, "
            f"{defect_counts['MEDIUM']} medium, {defect_counts['LOW']} low. "
            "Fix the findings; the questions are prompts, not defects."
        )

    japanese = detailed_results.get("japanese_readability", {})
    return {
        "defect_counts": defect_counts,
        "findings": findings,
        "questions": questions,
        "defect_score": defect_score,
        "score_max": 10,
        "japanese_readability_score": japanese.get("score"),
        "japanese_readability_status": japanese.get("status", "skipped"),
        "summary_comment": summary,
        "program": program,
        "detailed_scores": detailed_scores,
        "detailed_results": detailed_results,
        "tools_run": sorted(detailed_results),
        "tools_skipped": sorted(skip_set),
        "manual_review_prompts": (
            detailed_results.get("argument_evidence_map", {})
            .get("manual_review_prompts", [])
        ),
        "hint": (
            "findings locate defects and are worth fixing; questions cannot be "
            "answered by keyword presence and are for the author to judge. "
            "defect_score measures located mechanical defects only -- it is not "
            "a judgement of the research, and editing to raise it is wasted work."
        ),
        "source": "radia_mcp.grant_writing public document server",
    }
