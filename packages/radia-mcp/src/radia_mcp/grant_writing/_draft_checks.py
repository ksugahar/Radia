"""Read-only draft planning diagnostics; never automatic quality scores."""

from __future__ import annotations

import math
import re


def _draft_prose(text: str) -> str:
    # Deferred import keeps public registration in tools without a cycle.
    from .tools import _prose_for_lint, _read_text_if_path

    return _prose_for_lint(_read_text_if_path(text))


# A proposal may legitimately pursue more than one question. What a reviewer
# cannot do is read two of them as one. The defect named here is the mismatch
# between an announcement of singularity and the plurality that follows it --
# 「本研究の問いは明快である。」 followed by two questions joined by 「さらに」.
# The consistency check above cannot see this: it compares two statements of
# the SAME question and reports when their decisive nouns disagree. Here each
# statement is internally consistent; they are simply different questions.
_SINGULAR_ANNOUNCEMENTS = (
    "明快である",
    "明確である",
    "単純である",
    "一つである",
    "ひとつである",
    "一点である",
    "次である",
    "次のとおりである",
    "以下である",
    "次を問う",
)

# Accept both explicit question marks and the Japanese interrogative ending.
_QUESTION_ENDING = re.compile(r"(?:か|かどうか)[。．?？]?\s*$")

# 「さらに、…できるか。」 is the specific shape a second central question takes
# when it is appended to a first one rather than subordinated to it.
_ADDITIVE_OPENERS = (
    "さらに",
    "また",
    "加えて",
    "第二に",
    "もう一つ",
    "もうひとつ",
    "同時に",
    "あわせて",
    "併せて",
)

# 「成果は二つある」 in one section is a plan. Repeated across the aims, the
# outcomes and the funding-overlap sections, it is the shape of the whole
# document, and the reviewer must carry two of everything to the end.
_ENUMERATION_DECLARATION = re.compile(
    r"(?:二つ|三つ|ふたつ|みっつ|２つ|2つ|３つ|3つ|二点|三点)"
    r"(?:ある|あります|である|に分か|からなる|を挙げ|示す)"
)


def grant_writing_central_question_singularity_check(text: str) -> dict:
    """Check that a proposal announcing one central question states one.

    A reviewer reading 「本研究の問いは明快である。」 commits to holding one
    question. When the next sentences deliver two, the promise of clarity is
    spent before the science is read, and every later section inherits the
    ambiguity: which of the two does the schedule serve, which does the budget
    buy, which one fails if the other succeeds?

    This is deliberately not a ban on multiple questions. Lexical candidates
    require a close read: this check cannot determine whether the questions
    are scientifically distinct or one is subordinate to another. No quality
    score or submission verdict is produced.

    The check is optional: it needs an announced question or aim.
    """
    from .tools import _CLAIM_MARKERS, _prose_for_lint, _read_text_if_path

    # Navigation labels are not applicant assertions. Strip Markdown headings
    # before prose normalization, which can erase their structural markers.
    # Limit this to the claim detector: other diagnostics need field headings.
    raw = _read_text_if_path(text).replace("\r\n", "\n").replace("\r", "\n")
    raw = re.sub(r"(?m)^[ \t]{0,3}#{1,6}(?:[ \t]+[^\r\n]*|[ \t]*)$", "", raw)
    raw = re.sub(
        r"(?m)^[ \t]{0,3}[^\s\r\n][^\r\n]*\r?\n[ \t]{0,3}(?:=+|-+)[ \t]*$",
        "",
        raw,
    )
    text = _prose_for_lint(raw)
    sentences = [
        fragment
        for line in text.split("\n")
        for fragment in re.split(r"(?<=[。．!?！？])", line)
        if fragment.strip()
    ]

    def is_question(fragment: str) -> bool:
        return bool(_QUESTION_ENDING.search(fragment.strip()))

    announcements: list[dict] = []
    for index, sentence in enumerate(sentences):
        stripped = sentence.strip()
        marker = next((m for m in _CLAIM_MARKERS if m in stripped), None)
        if marker is None:
            continue
        singular = next((a for a in _SINGULAR_ANNOUNCEMENTS if a in stripped), None)
        announcements.append(
            {
                "sentence_index": index + 1,
                "marker": marker,
                "singular_phrase": singular,
                "text": stripped[:160],
            }
        )

    risks: list[dict] = []
    # The window is the claim zone: a central question is stated close to its
    # announcement, not a page later.
    window = 6
    for announcement in announcements:
        start = announcement["sentence_index"]
        zone = sentences[start - 1 : start + window]
        # Do not borrow questions from the next explicitly announced claim.
        for offset, fragment in enumerate(zone[1:], start=1):
            if any(marker in fragment for marker in _CLAIM_MARKERS) and any(
                phrase in fragment for phrase in _SINGULAR_ANNOUNCEMENTS
            ):
                zone = zone[:offset]
                break
        questions = [q.strip() for q in zone if is_question(q)]
        additive = [
            q.strip()
            for q in zone
            if is_question(q) and any(q.strip().startswith(o) for o in _ADDITIVE_OPENERS)
        ]
        if announcement["singular_phrase"] and len(questions) >= 2:
            risks.append(
                {
                    "type": "announced_singular_but_multiple",
                    "severity": "HIGH",
                    "sentence_index": announcement["sentence_index"],
                    "question_count": len(questions),
                    "comment": (
                        "「"
                        + announcement["singular_phrase"]
                        + "」と宣言した直後に問いが"
                        + str(len(questions))
                        + "つ並んでいる。"
                    ),
                    "recommendation": (
                        "中心の問いを一つに決める。残りは、その問いを解くための"
                        "実装条件、検証条件、適用範囲として従属させる。"
                    ),
                    "excerpts": [q[:160] for q in questions[:3]],
                }
            )
        elif additive:
            risks.append(
                {
                    "type": "additive_second_question",
                    "severity": "MEDIUM",
                    "sentence_index": announcement["sentence_index"],
                    "question_count": len(questions),
                    "comment": ("中心の問いの直後に、接続語で追加された問いがある。"),
                    "recommendation": ("追加の問いを主たる問いへ従属させるか、主従を明示する。"),
                    "excerpts": [q[:160] for q in additive[:2]],
                }
            )

    enumerations = [
        {"sentence_index": index + 1, "text": sentence.strip()[:120]}
        for index, sentence in enumerate(sentences)
        if _ENUMERATION_DECLARATION.search(sentence)
    ]
    return {
        "applicable": bool(announcements),
        "score": None,
        "automatic_score_prohibited": True,
        "status": (
            "review_required" if risks else "no_candidates" if announcements else "not_applicable"
        ),
        "announcement_count": len(announcements),
        "announcements": announcements,
        "enumeration_declarations": enumerations,
        "risk_count": len(risks),
        "risks": risks,
        "comments": list(dict.fromkeys(r["comment"] for r in risks)),
        "recommendations": list(dict.fromkeys(r["recommendation"] for r in risks)),
        "target": (
            "one announced central question delivers one question; "
            "further questions are subordinated to it"
        ),
        "source": "central-question singularity check",
    }


# Configurable planning assumption, not a validated typesetting guarantee.
_CHARS_PER_PAGE_DEFAULT = 1600

# A figure does not merely displace its own area: it takes the column with it.
# A half-page figure is worth about half a page of prose.
_FIGURE_PAGE_COST_DEFAULT = 0.5


def grant_writing_draft_length_budget_check(
    text: str,
    page_limit: float = 0.0,
    chars_per_page: int = _CHARS_PER_PAGE_DEFAULT,
    reserved_pages: float = 0.0,
    figure_count: int = 0,
    figure_page_cost: float = _FIGURE_PAGE_COST_DEFAULT,
) -> dict:
    """Estimate whether a plain-text draft can fit its form's page limit.

    ``grant_writing_page_limit_check`` is the authority, but it needs a
    compiled PDF, so it can only speak once the draft has been poured into the
    form -- days before a deadline, when cutting a third of the prose means
    cutting the argument rather than the wording. This check runs on the
    Markdown or text draft from the first day, and it exists to be wrong early
    rather than right late.

    ``reserved_pages`` is the part of the allowance the prose will never see:
    a 研究歴 and publication list, a mandatory budget table, a signature block.
    Those pages are consumed by content that is not written yet and cannot be
    compressed, so counting them as available is how a draft "fits" on paper
    and overflows in the form.

    The estimate is deliberately coarse. It answers "is a 2,000-character cut
    coming?", not "will line 43 wrap".
    """
    for name, value in (
        ("page_limit", page_limit),
        ("chars_per_page", chars_per_page),
        ("reserved_pages", reserved_pages),
        ("figure_count", figure_count),
        ("figure_page_cost", figure_page_cost),
    ):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError(f"{name} must be finite and nonnegative")
    if chars_per_page <= 0 or int(chars_per_page) != chars_per_page:
        raise ValueError("chars_per_page must be a positive integer")
    if int(figure_count) != figure_count:
        raise ValueError("figure_count must be a nonnegative integer")
    prose = _draft_prose(text)

    # Headings, list bullets and table pipes survive into the form as layout,
    # not as prose, so they are counted out of the character budget.
    body_lines = [line for line in prose.split("\n") if not line.strip().startswith("#")]
    body = "\n".join(body_lines)
    stripped = re.sub(r"[*_`|>\-\s]", "", body)
    japanese = sum(1 for c in stripped if "\u3040" <= c <= "\u30ff" or "\u4e00" <= c <= "\u9fff")
    counted = len(stripped)

    if page_limit <= 0:
        return {
            "applicable": False,
            "score": None,
            "automatic_score_prohibited": True,
            "status": "not_applicable",
            "risks": [],
            "risk_count": 0,
            "counted_characters": counted,
            "japanese_characters": japanese,
            "comments": [],
            "recommendations": [],
            "reason": ("page_limit が未指定のため判定しない。様式のページ数を渡す。"),
            "target": (
                "a draft whose character count fits the prose pages its form " "actually leaves it"
            ),
            "source": "draft-length budget check",
        }

    figure_pages = figure_count * figure_page_cost
    prose_pages = max(0.0, page_limit - reserved_pages - figure_pages)
    allowance = int(prose_pages * chars_per_page)
    overflow = counted - allowance
    ratio = round(counted / allowance, 2) if allowance > 0 else None

    risks: list[dict] = []
    if allowance <= 0:
        risks.append(
            {
                "type": "no_prose_pages_left",
                "severity": "HIGH",
                "comment": ("確保済みページと図で上限を使い切っており、本文の余地がない。"),
                "recommendation": ("reserved_pages と図の数を見直すか、様式のページ配分を変える。"),
            }
        )
    elif overflow > 0:
        severity = "HIGH" if overflow > allowance * 0.2 else "MEDIUM"
        risks.append(
            {
                "type": "draft_exceeds_prose_allowance",
                "severity": severity,
                "overflow_characters": overflow,
                "comment": (
                    "本文が"
                    + str(counted)
                    + "字あり、様式に収まる推定"
                    + str(allowance)
                    + "字を"
                    + str(overflow)
                    + "字超えている。"
                ),
                "recommendation": (
                    "章立てを様式の記入欄へ対応付けたうえで削る。"
                    "字数合わせの圧縮ではなく、二次的な結果と重複説明を落とす。"
                ),
            }
        )
    return {
        "applicable": True,
        "score": None,
        "automatic_score_prohibited": True,
        "status": (
            "no_prose_capacity"
            if allowance <= 0
            else "exceeds_estimate" if overflow > 0 else "within_estimate"
        ),
        "counted_characters": counted,
        "japanese_characters": japanese,
        "page_limit": page_limit,
        "reserved_pages": reserved_pages,
        "figure_pages": figure_pages,
        "prose_pages_available": prose_pages,
        "chars_per_page": chars_per_page,
        "estimated_allowance": allowance,
        "overflow_characters": max(0, overflow),
        "fill_ratio": ratio,
        "risk_count": len(risks),
        "risks": risks,
        "comments": [r["comment"] for r in risks],
        "recommendations": [r["recommendation"] for r in risks],
        "estimate_basis": (
            f"{chars_per_page} chars/page x {prose_pages} prose pages "
            f"(limit {page_limit} - reserved {reserved_pages} "
            f"- figures {figure_pages})"
        ),
        "warning": ("概算である。確定判定は組版後の " "grant_writing_page_limit_check で行う。"),
        "target": (
            "a draft whose character count fits the prose pages its form " "actually leaves it"
        ),
        "source": "draft-length budget check",
    }


# A funder's form is a fixed set of boxes, and a proposal is judged box by
# box. The generic section check above asks whether a proposal contains a
# purpose, a method and a schedule at all; it cannot know that 立石 asks for
# 「目指す姿」 to be tied to its own theme, or that 研究歴 and a publication
# list must fit inside the same three pages as the science. A draft can score
# a perfect 10 on the generic axes and still have no text for half the boxes.
#
# Each entry maps a printed field name to the words a draft would use if it
# had written that field. Synonyms are deliberate: a proposal rarely repeats
# the form's own heading.
_FORM_FIELD_PRESETS: dict[str, tuple[tuple[str, tuple[str, ...]], ...]] = {
    # 公益財団法人立石科学技術振興財団 研究助成(A)(B)(C) 申請書(研究課題).
    # Fields transcribed from the 2027 form: 1.研究目的 (1)-(4),
    # 2.研究計画 (1)-(2), 3.その他 (1)-(2).
    "tateisi_research": (
        ("1-(1) 目指す姿", ("目指す姿", "長期的", "将来像", "人間と機械", "調和")),
        ("1-(2) 提供する価値", ("提供する価値", "提供価値", "社会", "価値")),
        ("1-(3) 目標と解決すべき課題", ("目標", "技術課題", "解決すべき", "従来")),
        ("1-(4) 独創性", ("独創", "特色", "独自", "新規性")),
        ("2-(1) マイルストーンと方策", ("マイルストーン", "方策", "年度", "月")),
        ("2-(2) 共同研究者", ("共同研究者", "連携", "所属", "専門分野")),
        ("3-(1) 研究歴", ("研究歴", "これまでの研究", "経歴", "従事")),
        (
            "3-(2) 発表論文・著書・招待講演・受賞歴",
            ("論文", "著書", "招待講演", "受賞", "査読"),
        ),
    ),
}


def grant_writing_form_field_coverage_check(
    text: str,
    fields: str = "",
    preset: str = "",
) -> dict:
    """Locate lexical candidates for form fields, never certify coverage.

    Supply comma-separated field names or a named preset, not both. Matching
    neither proves presence nor proves absence: synonymous scientific prose
    may have no keyword hits, and an empty template may contain all of them.
    The caller must compare each candidate with the actual current form.
    """
    if preset and fields.strip():
        raise ValueError("supply fields or preset, not both")
    if preset and preset not in _FORM_FIELD_PRESETS:
        raise ValueError(f"unknown preset: {preset}")
    spec = _FORM_FIELD_PRESETS.get(preset, ())
    if not preset and fields.strip():
        names = list(dict.fromkeys(n.strip() for n in fields.split(",") if n.strip()))
        if not names:
            raise ValueError("fields must contain a nonempty field name")
        spec = tuple((name, (name,)) for name in names)
    prose = _draft_prose(text)
    candidates, unmatched = [], []
    for order, (name, keywords) in enumerate(spec, start=1):
        hits = [key for key in keywords if key in prose]
        item = {"field": name, "form_order": order, "matched": hits}
        if hits:
            item["first_position_ratio"] = round(
                min(prose.find(key) for key in hits) / max(1, len(prose)), 2
            )
            candidates.append(item)
        else:
            unmatched.append(item)
    return {
        "applicable": bool(spec),
        "score": None,
        "automatic_score_prohibited": True,
        "status": "manual_review_required" if spec else "not_applicable",
        "preset": preset or None,
        "available_presets": sorted(_FORM_FIELD_PRESETS),
        "field_count": len(spec),
        "candidate_count": len(candidates),
        "unmatched_count": len(unmatched),
        "candidate_fields": candidates,
        "unmatched_fields": unmatched,
        "risk_count": 0,
        "risks": [],
        "comments": [],
        "recommendations": (
            ["Compare every field with the actual form and read its answer."] if spec else []
        ),
        "warning": "Keyword hits are candidates, not coverage; no hits do not establish missing content.",
        "source": "form-field lexical candidate check",
    }


# Japanese application forms ask for 査読の有無 entry by entry, and for most
# venues the answer is fixed by the venue rather than by the paper. 電気学会
#研究会資料 carries no review at all; an applicant who writes 査読あり beside one
# has made a factual error that a reviewer in the same society notices at once.
# A conference proceedings is the genuinely ambiguous case and is reported as
# such rather than guessed.
_VENUE_REVIEW_CONVENTIONS = (
    (
        re.compile(
            r"電気学会研究会資料|研究会資料|電学研資|信学技報|技術研究報告",
            re.IGNORECASE,
        ),
        "not_reviewed",
        "研究会資料・技術研究報告は査読を行わない。",
    ),
    (
        re.compile(r"arXiv|preprint|プレプリント", re.IGNORECASE),
        "not_reviewed",
        "プレプリントは査読前の公開である。",
    ),
    (
        re.compile(
            r"IEEE\s+Trans(?:actions)?\b|IEICE\s+Trans(?:actions)?\b"
            r"|IEEJ\s+Trans(?:actions)?\b|電気学会論文誌",
            re.IGNORECASE,
        ),
        "reviewed",
        "学会論文誌は査読あり。",
    ),
    (
        re.compile(
            r"J(?:ournal)?\.?\s*(?:of\s+)?Magn(?:etics)?\.?\s*Soc"
            r"|日本磁気学会論文誌|Journal of Magnetic Resonance",
            re.IGNORECASE,
        ),
        "reviewed",
        "学術誌は査読あり。",
    ),
    (
        re.compile(
            r"COMPUMAG|CEFC|OIPE|IGTE|ICEAA|PIERS|LDIA|ISEF"
            r"|Int(?:ernational)?\.?\s+Conf|Proc\.|proceedings",
            re.IGNORECASE,
        ),
        "venue_dependent",
        (
            "国際会議。採否審査はあるが、査読の有無と単位（digest か full paper か）"
            "は会議規定による。投稿時の案内で確認する。"
        ),
    ),
)

_PUBLICATION_VENUE_CUE = re.compile(
    r"IEEE|IEICE|IEEJ|Trans(?:actions)?\b|Journal\b|Proc\.|proceedings|"
    r"研究会資料|技術研究報告|論文誌|学会誌|doi\s*:",
    re.IGNORECASE,
)
_NUMBERED_LIST_ENTRY = re.compile(r"^\s*(?:[-*]\s*)?(?:\d+[.．)、]|\[\d+\])")
_QUOTED_TITLE = re.compile(r"「[^」]+」|“[^”]+”|\"[^\"]+\"")


def _looks_like_publication_entry(line: str) -> bool:
    """Distinguish bibliography rows from numbered proposal prose."""
    if _PUBLICATION_VENUE_CUE.search(line):
        return True
    if not _NUMBERED_LIST_ENTRY.search(line):
        return False
    return bool(_QUOTED_TITLE.search(line) or line.count(",") >= 2)


def _publication_candidate_lines(raw: str) -> list[tuple[int, str]]:
    """Prefer an explicit numbered bibliography when a full proposal is given.

    A collaborator biography may mention an IEEE Transactions paper without
    being an achievement-list entry.  If the input contains a numbered
    bibliography, that stronger structure owns the audit; a standalone
    unnumbered citation remains supported when no numbered list is present.
    """
    lines = [
        (number, line.strip())
        for number, line in enumerate(raw.splitlines(), 1)
        if line.strip()
    ]
    candidates = [(n, line) for n, line in lines if _looks_like_publication_entry(line)]
    numbered = [(n, line) for n, line in candidates if _NUMBERED_LIST_ENTRY.search(line)]
    return numbered or candidates


def grant_writing_peer_review_convention_hints(text: str) -> dict:
    """List the peer-review convention of each venue named in an achievement list.

    The form asks the applicant to state 査読の有無 per entry, and getting one
    wrong is not a matter of emphasis: 研究会資料 is not reviewed, and claiming
    otherwise is a false statement about a specific paper.

    This reports what each named venue's convention is. It decides nothing.
    A proceedings can be either, and only the author knows what a given paper
    actually went through -- so entries are returned for the author to answer,
    never as a verdict.
    """
    raw = text
    from .tools import _read_text_if_path

    raw = _read_text_if_path(raw)
    entries: list[dict] = []
    for number, stripped in _publication_candidate_lines(raw):
        for pattern, convention, note in _VENUE_REVIEW_CONVENTIONS:
            match = pattern.search(stripped)
            if match:
                entries.append({
                    "line": number,
                    "venue_signal": match.group(0),
                    "convention": convention,
                    "note": note,
                    "excerpt": stripped[:160],
                })
                break

    counts: dict[str, int] = {}
    for entry in entries:
        counts[entry["convention"]] = counts.get(entry["convention"], 0) + 1
    return {
        "applicable": bool(entries),
        "score": None,
        "automatic_judgment_prohibited": True,
        "status": "manual_review_required",
        "entry_count": len(entries),
        "by_convention": counts,
        "entries": entries,
        "recommendations": [
            "Write 査読の有無 on every entry; the form asks for it explicitly.",
            "venue_dependent の行は、投稿時の案内か採録通知で確認して確定する。",
        ],
        "warning": (
            "Venue conventions only. A matched venue does not establish what a "
            "particular paper went through, and an unmatched line is not "
            "evidence that its venue is unreviewed."
        ),
        "source": "venue peer-review convention hints",
    }


# An achievement list is read as a record of what exists. An entry dated after
# the application, or numbered SA-27-0xx, is a plan wearing the clothes of a
# result; left unmarked among published work it invites the reviewer to
# discount the entries that are real.
_PUBLICATION_PLACEHOLDER = re.compile(
    r"0xx|x{2,}|＿＿|_{3,}|未定|未確定|TBD|\?\?+", re.IGNORECASE
)
_YEAR_IN_LINE = re.compile(r"(?<!\d)(20\d{2})(?!\d)")


def grant_writing_future_dated_publication_check(
    text: str,
    application_year: int = 0,
) -> dict:
    """List achievement entries dated ahead of the application, or unfilled.

    Two things make an entry unverifiable, and they usually arrive together:
    a year later than the application, and a number the venue has not issued
    yet. Neither is misconduct -- a paper genuinely in press belongs on the
    list -- but an entry that states neither 発表予定 nor 投稿中 reads as a
    published result the reviewer cannot find.

    ``application_year`` defaults to the current year when omitted.
    """
    from datetime import datetime

    from .tools import _read_text_if_path

    raw = _read_text_if_path(text)
    year_now = application_year or datetime.now().astimezone().year

    # An entry that already says what it is does not need the reviewer to ask.
    disclosed = re.compile(
        r"発表予定|投稿中|投稿予定|採録決定|採択済|in press|accepted|submitted",
        re.IGNORECASE,
    )

    entries: list[dict] = []
    for number, stripped in _publication_candidate_lines(raw):
        # The caller may pass a complete proposal, not just its bibliography.
        # A future grant year or a Markdown placeholder is not a publication.
        # Require either list-entry syntax or a bibliographic venue signal
        # before interpreting dates and placeholders as publication metadata.
        years = [int(y) for y in _YEAR_IN_LINE.findall(stripped)]
        future = [y for y in years if y > year_now]
        placeholder = _PUBLICATION_PLACEHOLDER.findall(stripped)
        if not future and not placeholder:
            continue
        entries.append({
            "line": number,
            "future_years": future,
            "placeholders": sorted(set(placeholder)),
            "status_disclosed": bool(disclosed.search(stripped)),
            "excerpt": stripped[:160],
        })

    undisclosed = [e for e in entries if not e["status_disclosed"]]
    return {
        "applicable": bool(entries),
        "score": None,
        "automatic_judgment_prohibited": True,
        "status": "manual_review_required",
        "application_year": year_now,
        "entry_count": len(entries),
        "undisclosed_count": len(undisclosed),
        "entries": entries,
        "recommendations": [
            "各行に発表状況（発表予定・投稿中・採録決定など）を明記する。",
            "確定していない書誌は、順位を下げるか、確定するまで載せない。",
        ],
        "warning": (
            "A future year is not a defect by itself: a paper in press belongs "
            "on the list. What this locates is an entry a reviewer cannot "
            "verify and that does not say so."
        ),
        "source": "future-dated / placeholder publication audit",
    }


_SOFTWARE_CATEGORIES = re.compile(
    r"オープンソース|OSS|ソフトウェア|解析基盤|計算基盤|研究基盤|"
    r"ソルバ|ライブラリ|パッケージ|プラットフォーム|ツール|コード|"
    r"open[- ]source|software|solver|library|package|platform|toolkit",
    re.IGNORECASE,
)
_SOFTWARE_OPERATIONS = re.compile(
    r"(?:磁場|電磁界|形状|電流|軌道|データ|モデル|方程式|解析|計算|設計|"
    r"最適化|可視化|連成|制御|予測|評価|生成|保存|管理|変換).{0,28}"
    r"(?:計算|解析|求め|解[くき]|扱[うい]|設計|最適化|可視化|連成|制御|"
    r"予測|評価|生成|保存|管理|変換|実行)|"
    r"(?:compute|solve|analyse|analyze|design|optimi[sz]e|simulate|visuali[sz]e|"
    r"couple|control|predict|evaluate|generate|store|manage|convert).{0,80}",
    re.IGNORECASE,
)
_IMPLEMENTATION_LABEL = re.compile(
    r"(?:Python|MATLAB|C\+\+|Java|Julia)[-‐‑–— ]?native|"
    r"(?:Python|MATLAB|C\+\+|Java|Julia)ネイティブ",
    re.IGNORECASE,
)


def _named_term_pattern(name: str) -> re.Pattern[str]:
    boundary_left = r"(?<![A-Za-z0-9])" if name[:1].isascii() else ""
    boundary_right = r"(?![A-Za-z0-9])" if name[-1:].isascii() else ""
    return re.compile(boundary_left + re.escape(name) + boundary_right, re.IGNORECASE)


def _named_diagnostic_sentences(raw: str) -> list[dict]:
    sentences: list[dict] = []
    for line_number, line in enumerate(raw.splitlines(), 1):
        for fragment in re.split(r"(?<=[。．.!?！？])", line):
            text = re.sub(r"\s+", " ", fragment).strip()
            if text:
                sentences.append({"line": line_number, "text": text})
    return sentences


def grant_writing_named_software_first_use_check(
    text: str,
    software_names: str = "",
) -> dict:
    """Show whether a named tool explains what it does at first use.

    A label such as ``公開研究基盤Radia`` tells the reviewer how the applicant
    values the software, not what the software accepts, computes, or returns.
    This diagnostic therefore looks for both a generic category and a concrete
    operation in the sentence that first names each tool.  It does not verify
    that the description is technically true.  Architecture labels such as
    ``Python-native`` are surfaced for source verification rather than accepted
    as a functional explanation.

    ``software_names`` may add comma-separated project-specific tools or
    methods to the built-in list.
    """
    from .tools import _read_text_if_path

    raw = _read_text_if_path(text)
    defaults = [
        "Radia", "NGSolve", "ONELAB", "openCFS", "FreeFEM++", "Gmsh",
        "GetDP", "preCICE", "OpenMDAO", "COMSOL", "JMAG", "ANSYS",
        "MATLAB", "Simulink",
    ]
    extras = [name.strip() for name in software_names.split(",") if name.strip()]
    names = list(dict.fromkeys(defaults + extras))
    sentences = _named_diagnostic_sentences(raw)

    entries: list[dict] = []
    for name in names:
        pattern = _named_term_pattern(name)
        first = next((item for item in sentences if pattern.search(item["text"])), None)
        if first is None:
            continue
        category_hits = [m.group(0) for m in _SOFTWARE_CATEGORIES.finditer(first["text"])]
        operation_hits = [m.group(0) for m in _SOFTWARE_OPERATIONS.finditer(first["text"])]
        implementation_claims = [
            m.group(0) for m in _IMPLEMENTATION_LABEL.finditer(first["text"])
        ]
        entries.append({
            "software": name,
            "line": first["line"],
            "first_use": first["text"][:320],
            "category_hits": category_hits,
            "operation_hits": operation_hits,
            "functional_identity_present": bool(category_hits and operation_hits),
            "implementation_claims_requiring_source_check": implementation_claims,
        })

    missing = [entry for entry in entries if not entry["functional_identity_present"]]
    source_checks = [
        entry for entry in entries
        if entry["implementation_claims_requiring_source_check"]
    ]
    return {
        "applicable": bool(entries),
        "score": None,
        "automatic_judgment_prohibited": True,
        "status": "manual_review_required" if entries else "not_applicable",
        "entry_count": len(entries),
        "missing_functional_identity_count": len(missing),
        "source_check_count": len(source_checks),
        "entries": entries,
        "manual_review_prompts": [
            "初出の一文だけで、専門外の審査者が入力・処理・出力の少なくとも一つを説明できるか。",
            "言語や実装方式のラベルは、機能説明の代わりにせず、リポジトリ等の一次情報と一致するか。",
        ],
        "recommendations": [
            "『公開研究基盤』だけで終えず、『磁石の形状と電流から三次元磁場を計算するソフトウェア』のように機能を先に書く。",
            "固有名は機能説明の後に置くか、『Radiaは、…するソフトウェアである』と同じ文で定義する。",
        ] if missing or source_checks else [],
        "warning": (
            "This is a first-use readability map, not a technical-fact checker. "
            "A detected category and operation can still describe the software incorrectly."
        ),
        "source": "named-software first-use functional-identity audit",
    }


_CURRENT_STATUS = re.compile(
    r"完成しつつある|実装済み|検証済み|開発済み|整備済み|接続済み|"
    r"公開済み|既に|すでに|現時点|現在利用|利用可能|予備検討|"
    r"動作している|動いている|利用している|継続して開発|"
    r"already|implemented|validated|available|current(?:ly)?",
    re.IGNORECASE,
)
_FUTURE_STATUS = re.compile(
    r"本研究(?:で|では)|助成期間|研究期間|今後|これから|新たに|着手|"
    r"構築する|実現する|"
    r"開発する|統合する|検証する|評価する|拡張する|目指す|予定|"
    r"機能を加える|つなぐ|"
    r"will|to be (?:developed|implemented|validated|integrated)|planned|proposed",
    re.IGNORECASE,
)


def grant_writing_capability_status_map(
    text: str,
    capability_names: str = "",
) -> dict:
    """Map named capabilities to current evidence and proposed work.

    The check prevents a mature foundation, a planned optimization, and a
    planned validation from collapsing into one present-tense product claim.
    It only lists sentences and lexical status signals.  Whether the claimed
    baseline is true, or whether the future work is sufficiently ambitious,
    remains a source-based human judgment.

    Pass comma-separated names such as ``Radia,HDiv-MMM,EnergyStop``.
    """
    from .tools import _read_text_if_path

    raw = _read_text_if_path(text)
    names = [name.strip() for name in capability_names.split(",") if name.strip()]
    if not names:
        return {
            "applicable": False,
            "score": None,
            "automatic_judgment_prohibited": True,
            "status": "not_applicable",
            "capabilities": [],
            "manual_review_prompts": [
                "capability_names に、現在地と助成期間中の到達点を分けたい名称を指定する。"
            ],
            "source": "current/future capability status map",
        }

    sentences = _named_diagnostic_sentences(raw)
    capabilities: list[dict] = []
    for name in names:
        pattern = _named_term_pattern(name)
        statements = []
        for sentence in sentences:
            if not pattern.search(sentence["text"]):
                continue
            current_hits = [m.group(0) for m in _CURRENT_STATUS.finditer(sentence["text"])]
            future_hits = [m.group(0) for m in _FUTURE_STATUS.finditer(sentence["text"])]
            status = (
                "current_and_future" if current_hits and future_hits else
                "current" if current_hits else
                "future" if future_hits else
                "unstated"
            )
            statements.append({
                "line": sentence["line"],
                "status": status,
                "current_signals": current_hits,
                "future_signals": future_hits,
                "excerpt": sentence["text"][:320],
            })
        capabilities.append({
            "capability": name,
            "mention_count": len(statements),
            "has_current_statement": any(
                s["status"] in {"current", "current_and_future"} for s in statements
            ),
            "has_future_statement": any(
                s["status"] in {"future", "current_and_future"} for s in statements
            ),
            "unstated_status_count": sum(s["status"] == "unstated" for s in statements),
            "statements": statements,
        })

    return {
        "applicable": True,
        "score": None,
        "automatic_judgment_prohibited": True,
        "status": "manual_review_required",
        "capabilities": capabilities,
        "manual_review_prompts": [
            "完成済み・完成しつつある基盤と、助成期間中に初めて行う最適化・統合・検証を別文で特定できるか。",
            "現在地の根拠は実装・試験・公開物に結び付き、将来到達点は研究計画と評価指標に結び付いているか。",
            "現在形の『求める』『探索する』『扱う』が、未実施の成果を完成済みに見せていないか。",
        ],
        "warning": (
            "Lexical status signals only. This map does not verify repository "
            "state, experimental completion, or scientific readiness."
        ),
        "source": "current/future capability status map",
    }
