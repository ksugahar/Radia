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
    from .tools import _CLAIM_MARKERS

    text = _draft_prose(text)
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
