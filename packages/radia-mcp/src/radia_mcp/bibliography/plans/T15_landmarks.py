"""Landmark papers of this field, reachable by description rather than by key.

Someone asks for "the Italian's original play-model paper" or "the H-matrix
proposal", not for bobbio1997play or hackbusch1999sparse. A literal term search
cannot bridge that: the entry says nothing about being Italian, or about being
the first.

So the connection is written down. Each landmark carries the words a person
would actually use -- in Japanese and English, including the author's name, the
method's acronym, and what the paper is FOR -- and bibliography_search consults
this alongside the entry text.

This is a cache, not an authority. It covers the works a Radia manuscript cites
as the origin of a method it uses; anything else is found by ordinary search, or
supplied and added. Every key here is checked against the bibliography at import
time, so a rename cannot leave the table pointing at nothing.
"""
from __future__ import annotations

# key -> the phrases someone would search by. Terms are matched case-folded as
# substrings, so "プレイモデル" also catches "プレイモデルの原著".
LANDMARKS: dict[str, list[str]] = {
    # --- hysteresis -------------------------------------------------------
    "bobbio1997play": [
        "プレイモデル 原著", "play model original", "play stop hysteron",
        "プレイヒステロン", "ストップモデル", "イタリア", "italian",
        "bobbio", "napoli", "naples",
    ],
    "krasnoselskii1989systems": [
        "プレイ作用素 起源", "play operator origin", "hysteresis operator",
        "krasnoselskii", "pokrovskii", "ヒステリシス作用素",
    ],
    "visintin1994differential": [
        "ヒステリシス 数学", "differential models of hysteresis", "visintin",
    ],
    "brokate1996hysteresis": [
        "ヒステリシス 相転移", "hysteresis phase transitions", "brokate",
        "sprekels",
    ],
    "preisach1935magnetische": [
        "preisach 原著", "プライザッハ 原著", "preisach original",
        "magnetische nachwirkung",
    ],
    "mayergoyz1986preisach": [
        "preisach 数学モデル", "mayergoyz", "preisach mathematical",
    ],
    "jiles1986theory": ["jiles atherton", "ja モデル", "jiles"],
    "stoner1948mechanism": ["stoner wohlfarth", "単磁区", "coherent rotation"],
    "matsuo2005representation": [
        "ベクトルプレイ", "vector play", "stop model 表現定理",
        "representation theorem",
    ],
    # --- benchmarks -------------------------------------------------------
    "karl1996description": [
        "team 28 原著", "team28 original", "team workshop problem 28",
        "電磁浮上 ベンチマーク", "electrodynamic levitation device",
    ],
    # --- finite elements and differential forms ---------------------------
    "nedelec1980mixed": [
        "辺要素 原著", "edge element original", "nedelec", "ネデレック",
        "mixed finite elements in r3",
    ],
    "bossavit1988whitney": [
        "whitney form", "ホイットニー要素", "bossavit whitney",
    ],
    "raviart1977mixed": ["raviart thomas", "rt 要素", "面要素 原著"],
    "arnold2006finite": [
        "feec", "有限要素外微分", "finite element exterior calculus", "arnold",
    ],
    # --- matrix compression and acceleration ------------------------------
    "hackbusch1999sparse": [
        "h行列 提案", "階層行列 原著", "h-matrix proposal", "hackbusch",
        "hierarchical matrix original",
    ],
    "bebendorf2000approximation": [
        "aca 原著", "adaptive cross approximation", "bebendorf", "交差近似",
    ],
    "greengard1987fast": [
        "fmm 原著", "高速多重極", "fast multipole", "greengard",
    ],
    # --- circuits and model order reduction -------------------------------
    "ruehli1974equivalent": ["peec 原著", "peec original", "ruehli"],
    "odabasioglu1998prima": ["prima", "passive reduced order interconnect"],
    "feldmann1995efficient": ["pvl", "lanczos mor", "pade via lanczos"],
    "kameari2018cauer": [
        "cln 原著", "cln original", "cauer ladder 原著",
        "cauer ladder network original", "亀有",
    ],
    "kuriyama2019cauer": [
        "cln 多展開点", "multiple expansion points", "栗山",
    ],
    "henneron2015model": ["pod 縮約", "pod model order reduction", "henneron"],
    # --- surface impedance and multiscale ---------------------------------
    "senior1962impedance": ["sibc senior", "表面インピーダンス 原著", "senior"],
    "mitzner1967integral": ["sibc mitzner", "mitzner"],
    "yuferev2009surface": ["sibc 教科書", "surface impedance textbook", "yuferev"],
    "hollaus2018some": [
        "hollaus マルチスケール", "msfem", "multiscale finite element",
        "積層鉄心 マルチスケール",
    ],
    "hollaus2026nonlinear": [
        "hollaus 実効表面インピーダンス", "esim hollaus",
        "nonlinear effective surface impedance",
    ],
    "hiruma2023extended": [
        "xfem 比留間", "em-xfem", "extended finite element eddy current",
    ],
    # --- open boundary and sources ----------------------------------------
    "berenger1994perfectly": ["pml 原著", "berenger", "perfectly matched layer"],
    "quarteroni1999domain": ["領域分割", "domain decomposition", "steklov"],
    "urankar1980vector": [
        "円弧 線電流", "arc segment vector potential", "urankar",
        "ビオサバール 解析解",
    ],
    "warburg1899ueber": ["warburg", "ワールブルク", "定位相要素 起源"],
    # --- software ----------------------------------------------------------
    "chubar1998three": ["radia 原著", "radia original", "chubar"],
    "schoberl2014implementation": ["ngsolve", "netgen"],
    "kamon1994fasthenry": ["fasthenry", "インダクタンス抽出"],
    # --- numerics ----------------------------------------------------------
    "walker2011anderson": ["anderson 加速", "anderson acceleration"],
    "henrotte1993new": [
        "henrotte 軸対称", "軸対称 要素", "axisymmetric element",
        "henrotte axisymmetric",
    ],
}


def landmark_keys_for(terms: list[str]) -> set[str]:
    """Keys whose description phrases contain every one of the given terms."""
    out = set()
    for key, phrases in LANDMARKS.items():
        joined = " ".join(phrases).lower()
        if all(t in joined for t in terms):
            out.add(key)
    return out


def bibliography_landmarks(topic: str = "") -> str:
    """List the landmark papers, or those matching a topic word.

    These are the works a Radia manuscript cites as the origin of a method it
    uses. The list is a convenience, not a boundary: anything absent is found by
    ordinary search or supplied and added.
    """
    from .._bibparse import read_bib_file, _strip_latex
    from .T14_canonical import CANONICAL

    entries = {e.key: e for e in read_bib_file(CANONICAL)}
    t = topic.strip().lower()
    rows, missing = [], []
    for key, phrases in LANDMARKS.items():
        if t and t not in " ".join(phrases).lower() and t not in key.lower():
            continue
        e = entries.get(key)
        if e is None:
            missing.append(key)
            continue
        rows.append((key, _strip_latex(e.fields.get("title", ""))[:60],
                     e.fields.get("year", "?"), phrases))

    out = [f"bibliography_landmarks: {len(rows)} 件"
           + (f"（topic={topic!r}）" if topic else "")]
    for key, title, year, phrases in sorted(rows, key=lambda r: r[0]):
        out.append(f"  {key}  ({year})")
        out.append(f"      {title}")
        out.append(f"      引ける語: {', '.join(phrases[:5])}")
    if missing:
        out.append(f"\n  表が指す先が書誌に無い: {', '.join(missing)}")
        out.append("  改名か削除で壊れている。表を直すこと。")
    return "\n".join(out)
