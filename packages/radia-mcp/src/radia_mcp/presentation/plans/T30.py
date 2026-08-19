"""Presentation T30: Script vs slide coverage (verbal side 強化)
(presentation_script_vs_slide_coverage).

speaker_note (台本) が slide 内容を網羅しているかを診断。
台本に言及されない slide content は「発表中言い忘れる」リスク。
逆に、slide に無いが台本にある内容は補足説明として OK。

検出ロジック:
  - 各 slide の body text tokens と speaker_note tokens を抽出
  - per slide: coverage = speaker_note でカバーされた body token の比率
  - missing_from_script: slide にあるが note で言及されない tokens
  - extra_in_script: note にあるが slide に無い tokens (補足、基本 OK)
"""

from __future__ import annotations

import pathlib
import re


def _iter_shapes(shape_collection):
    try:
        from radia_mcp.presentation.tools import _walk_shapes
        yield from _walk_shapes(shape_collection)
    except Exception:
        for s in shape_collection:
            yield s


# A Latin token bounded by ASCII, NOT by \b.  Python's \w matches CJK, so
# "のHACApK" has no word boundary before the H and \b[A-Z] never fires --
# which made every acronym spoken inside a Japanese sentence invisible.
_LATIN = re.compile(r"(?<![A-Za-z0-9])[A-Z][A-Za-z0-9]{2,}(?![A-Za-z0-9])")
_KANJI = re.compile(r"[一-鿿]{2,}")
_KATAKANA = re.compile(r"[ァ-ヴー]{3,}")


def _extract_tokens(text: str) -> set[str]:
    tokens: set[str] = set()
    for pattern in (_KANJI, _KATAKANA, _LATIN):
        for m in pattern.finditer(text):
            tokens.add(m.group(0))
    return tokens


def _get_slide_body(slide) -> str:
    parts: list[str] = []
    try:
        from radia_mcp.presentation.tools import _slide_title
        title = _slide_title(slide)
    except Exception:
        title = ""
    for shape in _iter_shapes(slide.shapes):
        try:
            if shape.has_text_frame:
                txt = shape.text_frame.text or ""
                if txt.strip() != title:
                    parts.append(txt)
        except Exception:
            continue
    return "\n".join(parts)


def _get_notes(slide) -> str:
    try:
        if slide.has_notes_slide:
            return slide.notes_slide.notes_text_frame.text or ""
    except Exception:
        pass
    return ""


def presentation_script_vs_slide_coverage(pptx_path: str) -> dict:
    """台本 (speaker_note) が slide 内容を網羅しているかを per slide 診断。

    Returns:
        dict: per_slide / avg_coverage / low_coverage_slides /
              uncovered_tokens / comments
    """
    try:
        import pptx as _pptx
    except ImportError:
        return {"error": "python-pptx not installed."}
    p = pathlib.Path(pptx_path)
    if not p.exists():
        return {"error": f"file not found: {pptx_path}"}
    prs = _pptx.Presentation(str(p))

    per_slide: list[dict] = []
    low_coverage_slides: list[dict] = []

    # Deck furniture -- the affiliation in the corner, a running section label
    # -- appears on nearly every slide and is not spoken.  Scoring it against
    # the script marks every slide down for the same non-omission.
    slides = list(prs.slides)
    seen: dict[str, int] = {}
    for slide in slides:
        for token in _extract_tokens(_get_slide_body(slide)):
            seen[token] = seen.get(token, 0) + 1
    boilerplate = {t for t, n in seen.items()
                   if len(slides) >= 4 and n >= 0.8 * len(slides)}

    for i, slide in enumerate(slides, 1):
        if i == 1:
            # The cover carries the formal title, the English subtitle and the
            # author list.  skill.md says the speaker need not read the formal
            # title aloud, so scoring it as unspoken content is a false alarm --
            # the title checkers skip slide 1 for the same reason.
            per_slide.append({
                "slide_no": i, "coverage": None,
                "skip_reason": "cover_slide",
            })
            continue
        body = _get_slide_body(slide).strip()
        notes = _get_notes(slide).strip()
        body_tokens = _extract_tokens(body)
        notes_tokens = _extract_tokens(notes)

        # Skip if either is empty
        if not body_tokens and not notes_tokens:
            per_slide.append({
                "slide_no": i, "coverage": None,
                "skip_reason": "empty_both",
            })
            continue
        if not body_tokens:
            per_slide.append({
                "slide_no": i, "coverage": None,
                "skip_reason": "title_only_slide",
                "notes_chars": len(notes),
            })
            continue
        if not notes_tokens:
            per_slide.append({
                "slide_no": i,
                "coverage": 0.0,
                "body_tokens": len(body_tokens),
                "notes_tokens": 0,
                "missing_from_script": sorted(body_tokens)[:10],
                "warning": "speaker_note 空で slide に内容あり",
            })
            low_coverage_slides.append({
                "slide_no": i,
                "coverage": 0.0,
                "reason": "speaker_note 未記入",
                "missing_tokens": sorted(body_tokens)[:5],
            })
            continue

        body_tokens = body_tokens - boilerplate
        if not body_tokens:
            per_slide.append({
                "slide_no": i, "coverage": None,
                "skip_reason": "only_deck_furniture",
            })
            continue

        # A token counts as spoken if the script contains it at all -- saying
        # "Gram二重積分" is saying "二重積分".  Set intersection alone required
        # the two segmentations to agree, so a longer phrase read as a gap.
        covered = {t for t in body_tokens if t in notes_tokens or t in notes}
        missing = body_tokens - covered
        extra = notes_tokens - body_tokens
        coverage = len(covered) / len(body_tokens)

        per_slide.append({
            "slide_no": i,
            "coverage": round(coverage, 3),
            "body_tokens": len(body_tokens),
            "notes_tokens": len(notes_tokens),
            "covered_tokens": sorted(covered)[:5],
            "missing_from_script": sorted(missing)[:5],
            "extra_in_script": sorted(extra)[:5],
        })

        if coverage < 0.5:
            low_coverage_slides.append({
                "slide_no": i,
                "coverage": round(coverage, 3),
                "reason": f"slide 内容の {coverage*100:.0f}% しか言及されない",
                "missing_tokens": sorted(missing)[:5],
            })

    # Overall stats
    valid_coverages = [s["coverage"] for s in per_slide
                       if s.get("coverage") is not None]
    avg_coverage = (
        sum(valid_coverages) / len(valid_coverages) if valid_coverages else 0
    )
    n_low = len(low_coverage_slides)

    # Score
    score = round(avg_coverage * 10, 1)
    if n_low > len(valid_coverages) * 0.3 and valid_coverages:
        score -= 2
    score = max(0.0, round(score, 1))

    # Comments
    comments: list[str] = []
    comments.append(
        f"平均 coverage: {avg_coverage*100:.0f}% "
        f"({len(valid_coverages)} slides 評価対象) | "
        f"low-coverage (<50%): {n_low} slides"
    )

    if n_low >= 1:
        comments.append(
            f"⚠ {n_low} slides で台本 coverage <50%。"
            "発表中言い忘れるリスク。speaker_note に言及追加推奨。"
        )
        for lc in low_coverage_slides[:5]:
            comments.append(
                f"  • slide {lc['slide_no']} ({lc['coverage']*100:.0f}%): "
                f"missing {lc['missing_tokens']}"
            )

    if avg_coverage >= 0.7 and n_low == 0:
        comments.append(
            f"✓ 平均 coverage {avg_coverage*100:.0f}%、全 slide で balance OK"
        )

    hint = (
        "body token と note token の overlap で coverage を推定。"
        "missing は『slide に書いてあるが台本に無い → 発表中飛ばすリスク』。"
        "extra は『台本にあるが slide に無い → 補足説明、基本 OK』。"
        "coverage 50% 未満は警戒、70%+ で良好。"
    )

    return {
        "score": score,
        "score_max": 10,
        "avg_coverage_pct": round(avg_coverage * 100, 1),
        "n_slides_evaluated": len(valid_coverages),
        "n_low_coverage_slides": n_low,
        "deck_boilerplate_ignored": sorted(boilerplate),
        "low_coverage_slides": low_coverage_slides,
        "per_slide": per_slide[:30],
        "comments": comments,
        "hint": hint,
        "source": (
            "presentation リハーサル支援。既存 speaker_note_ratio (T7) は "
            "note の有無、本 tool は note と slide の内容一致度を見る。"
            "(presentation 2026-04 採用、v0.26.0)。"
        ),
    }
