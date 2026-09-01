"""Does the talk turn? A 起承転結 check that reads the bottom line, not the title.

The existing structure checks look for phase words in slide titles -- "Results",
"Method", "Background". The lab's own message-hierarchy rule forbids exactly
that: a title is a short noun phrase and the claim goes in the bottom takeaway.
So a deck can score 10/10 on hierarchy and near zero on structure, and neither
number says anything about whether the talk has an arc.

This reads the takeaway line instead, where the lab's format puts the assertion,
and looks for the four moves rather than for vocabulary:

  起  the object and the existing approach, stated without complaint
  承  the complaint: what the existing approach cannot do
  転  the turn: a proposal that answers that complaint
  結  evidence for the proposal, its limits, and the summary

What is measured is the SHAPE -- a complaint before the turn, a turn that
answers it, evidence after it, a close that echoes it -- not whether the words
match a list. A deck whose turn comes on the last slide, or which never states a
problem, or whose summary has drifted from what was proposed, fails here while
passing every keyword test.
"""
from __future__ import annotations

import re

# The complaint: something the current approach cannot do, or does badly.
_COMPLAINT = re.compile(
    r"(?i)\bcannot\b|\bcan not\b|\bfails?\b|\bno[t]? \w+ed\b|\bstops? at\b"
    r"|\bby hand\b|\bceiling\b|\blimit(?:ation)?s?\b|\bprice\b|\bwrong\b"
    r"|\bmisses\b|\bbreaks?\b|\bnot the\b|できない|限界|欠陥|課題|問題")
# The turn: an instruction or a resolution, often answering the complaint.
_TURN = re.compile(
    r"(?i)^(?:do not|use|take|add|replace|instead)\b|\buse both\b"
    r"|\bcomes out of\b|\bremoves?\b|\bnothing is (?:left to be )?fit"
    r"|\bnothing is fitted\b|\binstead\b|\bwe propose\b|提案|代わりに")
# Evidence: a takeaway that carries a measured number.
_NUMBER = re.compile(r"\d+(?:\.\d+)?\s*(?:%|x\b|×|倍)")


def _moves(takeaway: str) -> set[str]:
    t = (takeaway or "").strip()
    out = set()
    if _COMPLAINT.search(t):
        out.add("complaint")
    if _TURN.search(t):
        out.add("turn")
    if _NUMBER.search(t):
        out.add("evidence")
    return out


def _overlap(a: str, b: str) -> float:
    """Content-word overlap, as a fraction of the shorter side."""
    stop = {"the", "a", "an", "is", "are", "of", "to", "and", "in", "on",
            "that", "it", "this", "with", "for", "not", "no", "one", "at",
            "from", "by", "its", "has", "have", "be", "so", "but", "which"}
    wa = {w for w in re.findall(r"[a-z]+", a.lower()) if w not in stop}
    wb = {w for w in re.findall(r"[a-z]+", b.lower()) if w not in stop}
    if not wa or not wb:
        return 0.0
    return len(wa & wb) / min(len(wa), len(wb))


def presentation_kishotenketsu_check(pptx_path: str,
                                     backup_title: str = "Backup") -> dict:
    """Report whether the deck has a 起承転結 arc, and where it turns.

    backup_title : slides from this one onward are treated as backup and left
                   out; they are answers to questions, not part of the arc.
    """
    try:
        from pptx import Presentation
    except ImportError:
        return {"error": "python-pptx not installed."}

    prs = Presentation(pptx_path)
    slides = []
    for i, s in enumerate(prs.slides, 1):
        title = s.shapes.title.text.strip() if s.shapes.title else ""
        texts = [sh.text_frame.text.strip() for sh in s.shapes
                 if sh.has_text_frame and sh.text_frame.text.strip()]
        body = [t for t in texts if t != title]
        # the takeaway is the lab's bottom banner: the last text block
        takeaway = body[-1].splitlines()[-1].strip() if body else ""
        slides.append({"slide": i, "title": title, "takeaway": takeaway})

    cut = next((k for k, s in enumerate(slides)
                if s["title"].strip() == backup_title), len(slides))
    main = slides[1:cut] if cut > 1 else slides       # drop the title card
    if len(main) < 4:
        return {"error": f"only {len(main)} content slides; too few to have an arc."}

    for s in main:
        s["moves"] = sorted(_moves(s["takeaway"]))

    first_complaint = next((k for k, s in enumerate(main)
                            if "complaint" in s["moves"]), None)
    turn = None
    if first_complaint is not None:
        turn = next((k for k, s in enumerate(main)
                     if k > first_complaint and "turn" in s["moves"]), None)
    evidence = [k for k, s in enumerate(main)
                if "evidence" in s["moves"] and (turn is None or k > turn)]

    # 起 is whatever precedes the first complaint; 結 begins at the first
    # evidence slide after the turn.
    arc = {}
    if first_complaint is not None:
        arc["ki"] = [main[k]["slide"] for k in range(first_complaint)]
        end_sho = turn if turn is not None else len(main)
        arc["sho"] = [main[k]["slide"] for k in range(first_complaint, end_sho)]
    if turn is not None:
        end_ten = evidence[0] if evidence else len(main)
        arc["ten"] = [main[k]["slide"] for k in range(turn, end_ten)]
        arc["ketsu"] = [main[k]["slide"] for k in range(end_ten, len(main))]

    # does the turn answer the complaint, and does the close echo the turn?
    # A deck may raise several complaints and answer them on different slides.
    # Pairing only the first of each reports "the turn answers nothing" for a
    # talk that answers everything, just not in that pair.
    answers, unanswered = 0.0, []
    complaints = [k for k, sl in enumerate(main) if "complaint" in sl["moves"]]
    proposals = [k for k, sl in enumerate(main)
                 if "turn" in sl["moves"] and (turn is None or k >= turn)]
    for c in complaints:
        if turn is not None and c > turn:
            continue                      # a limitation stated after the turn
        best = max((_overlap(main[c]["takeaway"], main[t]["takeaway"])
                    for t in proposals), default=0.0)
        answers = max(answers, best)
        if best < 0.15:
            unanswered.append({"slide": main[c]["slide"],
                               "takeaway": main[c]["takeaway"],
                               "best_overlap": round(best, 2)})
    echo = 0.0
    if turn is not None:
        echo = max((_overlap(main[t]["takeaway"], main[-1]["takeaway"])
                    for t in proposals), default=0.0)

    checks = {
        "起 exists (setup before any complaint)":
            first_complaint is not None and first_complaint > 0,
        "承 exists (a stated complaint)": first_complaint is not None,
        "転 exists (a proposal after the complaint)": turn is not None,
        "結 exists (measured evidence after the turn)": bool(evidence),
        "転 answers 承 (shared content words)": answers >= 0.15,
        "結 echoes 転 (the close returns to the proposal)": echo >= 0.15,
        "転 lands in the middle (not the last quarter)":
            turn is not None and turn < len(main) * 0.75,
    }
    score = round(10.0 * sum(checks.values()) / len(checks), 1)

    comments = [f"{'OK  ' if v else 'FAIL'} {k}" for k, v in checks.items()]
    if first_complaint is None:
        comments.append(
            "承 が無い。従来法の何ができないかを下端文で言い切ると、転が効く。")
    elif turn is None:
        comments.append(
            "転が見つからない。提案を下端文に命令形か『〜が決まる』の形で置く。")
    for u in unanswered:
        comments.append(
            f"スライド {u['slide']} の欠陥に答える提案が見当たらない"
            f"（最大の語重なり {u['best_overlap']:.2f}）: 「{u['takeaway'][:52]}」")
    if turn is not None and echo < 0.15:
        comments.append(
            f"まとめが転に戻っていない（語の重なり {echo:.2f}）。"
            "最後の下端文は提案の言い直しにする。")

    outline = [
        {"phase": ph, "slides": arc.get(ph, []),
         "takeaways": [s["takeaway"] for s in main if s["slide"] in arc.get(ph, [])]}
        for ph in ("ki", "sho", "ten", "ketsu")
    ]
    return {
        "score": score,
        "score_max": 10,
        "content_slides": len(main),
        "backup_from_slide": slides[cut]["slide"] if cut < len(slides) else None,
        "arc": {k: v for k, v in arc.items()},
        "turn_slide": main[turn]["slide"] if turn is not None else None,
        "turn_position": round((turn + 1) / len(main), 2) if turn is not None else None,
        "turn_answers_complaint": round(answers, 2),
        "unanswered_complaints": unanswered,
        "close_echoes_turn": round(echo, 2),
        "checks": checks,
        "comments": comments,
        "outline": outline,
        "hint": "タイトルではなく下端 takeaway で判定する。研究室の規約は"
                "タイトルを名詞句に限るので、phase キーワード方式は原理的に"
                "使えない。承→転→結 の形そのものを見る。",
        "source": "起承転結（四段構成）を学会発表に適用。"
                  "presentation_check_slide_message_hierarchy と対で使う: "
                  "あちらは 1 枚の中の役割分担、こちらは 枚をまたぐ筋。",
    }
