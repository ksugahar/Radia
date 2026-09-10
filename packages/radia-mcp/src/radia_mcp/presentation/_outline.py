"""途中の目次: 背景がどこで終わり、提案がどこから始まるかを枚で示しているか。

松尾先生 (2026-09-10, IGTE'26 のデッキ査読):

    よくあるように途中に目次/アウトラインを挟むと分かりやすいと思いました。
    例えば、どこまでのスライドが現状の問題(あるいは背景)で、どこからが今回の
    解決法の話になるのか、明示する方が聴衆には分かりやすいと思いました。

冒頭に一度だけ目次を出す形は、この指摘には答えていない。聞き手が「今どこに
いるのか」を見失うのは冒頭ではなく中盤である ―― 従来法の話が三枚続いたあと、
提案が始まる境目で、いま聞いているのが「まだ問題の話」なのか「もう解決法の
話」なのかが分からなくなる。だからこの検査は、目次が**在るか**ではなく
**転の直前に在るか**を見る。

境目は 起承転結 の 転 から取る（``presentation_kishotenketsu_check`` と同じ
定義を使い、``read_deck`` で同じ読み方をする）。転はタイトルの語ではなく下端
の主張文から決まるので、研究室の「タイトルは名詞句・主張は下端帯」の様式と
衝突しない ―― 「Method」「Results」といった章題を並べる流儀を前提にした検査
なら、この様式のデッキでは常に不合格になってしまう。

目次の枚は二通りの見つけ方をする。題が目次語（Outline / 目次 / 本日の流れ …）
であるか、あるいは本文の行がこのデッキ自身の他の枚の題を並べているか。後者が
要るのは、区切りの枚に「Where we are going」のような固有の題を付けることが
あるからで、題だけを見ると見落とす。

短い発表では枚を一枚使う余裕がないこともある。その場合の答えは「検査を無視
する」ではなく、**境目の枚の下端文で章が変わったと言い切る**か、置き場所を
`suggested_outline` が指す位置に絞ることである。何秒かかるかは
``presentation_estimate_per_slide_time`` で測れる。
"""
from __future__ import annotations

import re

from ._kishotenketsu import _overlap, presentation_kishotenketsu_check, read_deck

# 題が目次だと名乗っている枚。
_OUTLINE_TITLE = re.compile(
    r"(?i)\boutline\b|\bagenda\b|\bcontents\b|\broad ?map\b|\boverview\b"
    r"|\bwhat follows\b|\bwhere we are going\b|\bthe plan of the talk\b"
    r"|目次|構成|本日の流れ|発表の流れ|アウトライン|お話しする順|全体像")
# 片側が「現状・背景・問題」だと言っているか。
_PROBLEM_SIDE = re.compile(
    r"(?i)\bbackground\b|\bproblem\b|\bmotivation\b|\bstate of the art\b"
    r"|\bprior work\b|\bwhat is known\b|\bwhere it fails\b|\bthe difficulty\b"
    r"|背景|課題|問題|現状|従来|これまで")
# もう片側が「提案・解決法・結果」だと言っているか。
_SOLUTION_SIDE = re.compile(
    r"(?i)\bmethod\b|\bapproach\b|\bproposal\b|\bwe propose\b|\bour \w+\b"
    r"|\bthis talk\b|\bsolution\b|\bresults?\b|\bevidence\b|\bwhat we do\b"
    r"|提案|手法|解決|本研究|本発表|結果|検証")


def _looks_like_outline(slide: dict, other_titles: list[str]) -> tuple[bool, str]:
    """Is this slide an outline / agenda / section divider?

    Two signals, because a divider need not be titled "Outline": the title
    says so, or the body lists the deck's own slide titles back to it.
    """
    if _OUTLINE_TITLE.search(slide["title"] or ""):
        return True, "title"
    lines = [ln.strip() for ln in (slide["text"] or "").splitlines()
             if len(ln.strip()) > 3]
    echoed = sum(1 for ln in lines
                 if max((_overlap(ln, t) for t in other_titles), default=0.0) >= 0.5)
    if echoed >= 2:
        return True, f"lists {echoed} of the deck's own titles"
    return False, ""


def presentation_check_outline_slide(pptx_path: str,
                                     backup_title: str = "Backup",
                                     max_slides_before_turn: int = 1) -> dict:
    """Report whether a mid-deck outline marks the problem/solution boundary.

    max_slides_before_turn : how far ahead of the turn a divider still counts
        as marking it. 1 means the slide immediately before the turn (the
        usual placement); 0 demands the divider be the turn slide itself.
    """
    try:
        from pptx import Presentation  # noqa: F401
    except ImportError:
        return {"error": "python-pptx not installed."}

    slides, cut, main = read_deck(pptx_path, backup_title)
    if len(main) < 4:
        return {"error": f"only {len(main)} content slides; too few to need "
                         "an outline."}

    arc = presentation_kishotenketsu_check(pptx_path, backup_title)
    turn = arc.get("turn_slide")

    titles = [s["title"] for s in main if s["title"]]
    found = []
    for s in main:
        others = [t for t in titles if t != s["title"]]
        ok, why = _looks_like_outline(s, others)
        if ok:
            text = f"{s['title']}\n{s['text']}"
            found.append({
                "slide": s["slide"], "title": s["title"], "detected_by": why,
                "names_the_problem_side": bool(_PROBLEM_SIDE.search(text)),
                "names_the_solution_side": bool(_SOLUTION_SIDE.search(text)),
            })

    # 「途中」= 転の直前。冒頭にしか無い目次は、境目では思い出されない。
    at_boundary = [f for f in found
                   if turn is not None
                   and 0 <= turn - f["slide"] <= max_slides_before_turn]
    marks_split = [f for f in found
                   if f["names_the_problem_side"] and f["names_the_solution_side"]]

    problem = (arc.get("arc", {}).get("ki", []) + arc.get("arc", {}).get("sho", []))
    solution = (arc.get("arc", {}).get("ten", []) + arc.get("arc", {}).get("ketsu", []))
    by_no = {s["slide"]: s["title"] for s in main}
    suggested = {
        "insert_before_slide": turn,
        "problem": [{"slide": n, "title": by_no.get(n, "")} for n in problem],
        "solution": [{"slide": n, "title": by_no.get(n, "")} for n in solution],
    }

    checks = {
        "目次・区切りの枚がある": bool(found),
        "その一枚が転の直前にある（冒頭だけではない）": bool(at_boundary),
        "問題の側と解決法の側を両方名指ししている": bool(marks_split),
    }
    score = round(10.0 * sum(checks.values()) / len(checks), 1)
    comments = [f"{'OK  ' if v else 'FAIL'} {k}" for k, v in checks.items()]

    if turn is None:
        comments.append(
            "転が見つからないので境目を決められない。"
            "presentation_kishotenketsu_check を先に通すこと ―― "
            "どこから解決法かを枚で示す前に、話が転じている必要がある。")
    elif not found:
        comments.append(
            f"目次の枚が無い。{turn} 枚目から解決法の話に変わるので、その直前に"
            f"一枚入れ、{problem and problem[0]}–{turn - 1} 枚目までが問題、"
            f"{turn} 枚目からが提案、と明示する。"
            "`suggested_outline` に並べる項目を入れてある。")
    elif not at_boundary:
        where = ", ".join(str(f["slide"]) for f in found)
        comments.append(
            f"目次はある（{where} 枚目）が、境目（{turn} 枚目の直前）に無い。"
            "冒頭で一度読み上げた目次は、三枚あとの境目では思い出されない。"
            "同じ枚を、今どこにいるかを変えてもう一度出すのが普通のやり方で、"
            "発話は各回十数秒で済む。")
    if found and not marks_split:
        comments.append(
            "目次が章の名前を並べているだけで、どこまでが問題でどこからが"
            "解決法かを言っていない。二つの側にそれぞれ見出しを付けること"
            "（松尾先生 2026-09-10）。")

    return {
        "score": score,
        "score_max": 10,
        "content_slides": len(main),
        "turn_slide": turn,
        "outline_slides": found,
        "outline_at_boundary": [f["slide"] for f in at_boundary],
        "suggested_outline": suggested,
        "checks": checks,
        "comments": comments,
        "hint": "境目は転（起承転結）から取る。タイトルの語ではなく下端の"
                "主張文で決まるので、「タイトルは名詞句・主張は下端帯」という"
                "研究室の様式のままで使える。枚を足す余裕が無いときは、"
                "境目の枚の下端文で章が変わったと言い切るのでもよい。",
        "source": "松尾先生のデッキ査読 2026-09-10。"
                  "presentation_kishotenketsu_check と対で使う: あちらは"
                  "筋が転じているか、こちらはその転じ目が聴衆に見えているか。",
    }
