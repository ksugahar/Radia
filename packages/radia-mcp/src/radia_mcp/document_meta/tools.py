"""Tool functions for the document_meta sub-module.

Cross-cutting helpers that do not belong to a specific document type
(grant / paper / presentation): deadline tracking, version diffing,
template scaffolding, and full-tree linting.
"""
from __future__ import annotations

import datetime
import difflib
import json
import pathlib
import re


def document_meta_deadline_countdown(deadline_iso: str,
                                       context: str = "") -> dict:
    """任意の締切までの日数と推奨アクションを返す。

    Args:
        deadline_iso: ISO-8601 date ("2026-05-15") or datetime ("2026-05-15T17:00").
        context: "jsps_s2" | "ieee_tmag" | "ieej_kaishi" | "aps_prb" | "" (generic)

    Returns:
        {deadline, today, days_remaining, weeks_remaining, urgency,
         phase, recommended_action}
    """
    try:
        if "T" in deadline_iso:
            dl = datetime.datetime.fromisoformat(deadline_iso)
        else:
            dl = datetime.datetime.fromisoformat(deadline_iso + "T23:59:59")
    except ValueError as exc:
        return {"error": f"bad deadline_iso {deadline_iso!r}: {exc}"}
    now = datetime.datetime.now()
    delta = dl - now
    days = delta.days
    weeks = days / 7.0

    if days < 0:
        urgency = "past_due"
        phase = "expired"
        action = "Deadline passed. Reassess scope or target next round."
    elif days < 3:
        urgency = "critical"
        phase = "final_polish"
        action = "Freeze scope. Proof-read + validate_pdf_pages + count_underlines only."
    elif days < 7:
        urgency = "high"
        phase = "polish"
        action = "Stop feature work. Run full lint: count_underlines, analyze_sentences, check_overfull_hbox, find_undefined_acronyms."
    elif days < 21:
        urgency = "medium"
        phase = "revision"
        action = "Second round of content revision. Co-author feedback integration + figure finalization."
    elif days < 60:
        urgency = "low"
        phase = "drafting"
        action = "First full draft. Establish §structure, figure kamishibai, bib skeleton."
    else:
        urgency = "planning"
        phase = "ideation"
        action = "Outline + literature survey + figure placeholders."

    ctx_hints = {
        "jsps_s2": "招へい事業: 5 pages, underline density 3-6/page, acronym first-use definition.",
        "ieee_tmag": "IEEE TMAG: 4 pages (letter) / open (regular), IEEEtran.cls, Abstract <=250w.",
        "ieej_kaishi": "電気学会論文誌: 6-8 pages typical, jsarticle + yohei class, 和文 main.",
        "aps_prb": "APS PRB: RevTeX4-2, no page limit, Abstract <=600w.",
    }
    return {
        "deadline": dl.isoformat(timespec="minutes"),
        "today": now.isoformat(timespec="minutes"),
        "days_remaining": days,
        "weeks_remaining": round(weeks, 1),
        "urgency": urgency,
        "phase": phase,
        "recommended_action": action,
        "context_hint": ctx_hints.get(context, ""),
    }


def document_meta_diff_versions(path_a: str, path_b: str,
                                  max_hunks: int = 20) -> dict:
    """2 つのテキスト file の unified diff を返す (作文 version 比較)。

    Use case: v5 → v6 で section 1 の句読点と下線をどう直したかを確認。

    Args:
        path_a: 旧版 (tex / md / txt)。
        path_b: 新版。
        max_hunks: 返却する diff hunk の上限。

    Returns:
        {hunks, added_lines, removed_lines, chars_a, chars_b,
         underline_count_a, underline_count_b, sentence_count_a/b}
    """
    pa = pathlib.Path(path_a).expanduser()
    pb = pathlib.Path(path_b).expanduser()
    if not pa.exists():
        return {"error": f"path_a not found: {pa}"}
    if not pb.exists():
        return {"error": f"path_b not found: {pb}"}
    try:
        text_a = pa.read_text(encoding="utf-8")
        text_b = pb.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        return {"error": f"decode error: {exc}"}

    ln_a = text_a.splitlines()
    ln_b = text_b.splitlines()
    diff_iter = difflib.unified_diff(ln_a, ln_b,
                                       fromfile=pa.name, tofile=pb.name,
                                       lineterm="", n=3)
    hunks: list[str] = []
    current: list[str] = []
    added = removed = 0
    for line in diff_iter:
        if line.startswith("@@"):
            if current:
                hunks.append("\n".join(current))
                if len(hunks) >= max_hunks:
                    break
            current = [line]
        elif line.startswith("+") and not line.startswith("+++"):
            added += 1
            current.append(line)
        elif line.startswith("-") and not line.startswith("---"):
            removed += 1
            current.append(line)
        else:
            current.append(line)
    if current and len(hunks) < max_hunks:
        hunks.append("\n".join(current))

    ul_re = re.compile(r"\\uline\{|\\underline\{|\\ul\{")
    sent_re = re.compile(r"[。．\.!?！？]")
    return {
        "chars_a": len(text_a),
        "chars_b": len(text_b),
        "lines_a": len(ln_a),
        "lines_b": len(ln_b),
        "added_lines": added,
        "removed_lines": removed,
        "underline_count_a": len(ul_re.findall(text_a)),
        "underline_count_b": len(ul_re.findall(text_b)),
        "sentence_count_a": len(sent_re.findall(text_a)),
        "sentence_count_b": len(sent_re.findall(text_b)),
        "hunks": hunks,
        "truncated": len(hunks) >= max_hunks,
    }


_TEMPLATES = {
    "jsps_s2": r"""% JSPS 外国人招へい事業 (S2) 申請書 skeleton — jsarticle + ulem
\documentclass[11pt,a4paper]{jsarticle}
\usepackage{graphicx}
\usepackage{ulem}
\usepackage{wrapfig}
\newcommand{\sectionbox}[2]{\par\noindent\fbox{\begin{minipage}{.98\linewidth}\textbf{#1}\\ \small #2\end{minipage}}\par\smallskip}

\begin{document}
\sectionbox{1. 研究目的・内容}{招へい研究者と申請者が共同で行う研究の目的、内容}
% Hook: 何が問題か / なぜ招へいが必要か
% Method: ESIM (Effective Surface Impedance Method) + POD + Multi-scale
% Why Karl: NGSolve core dev, Hollaus-SIBC の原典
% Validation platform value: FEM reference + target method + real data on one platform

\sectionbox{2. 研究計画}{期間中の活動計画、役割分担、成果目標}
% D1--D7: 環境 setup / D8--D21: 実装 / D22--D28: 検証 + writeup

\end{document}
""",
    "ieee_tmag": r"""% IEEE TMAG journal submission skeleton — IEEEtran.cls
\documentclass[journal,10pt]{IEEEtran}
\usepackage{graphicx,amsmath,amssymb,cite}
\begin{document}
\title{Title goes here}
\author{First Author\textsuperscript{1}, Second Author\textsuperscript{2}}
\maketitle
\begin{abstract}
[<=250 words. Background <=25\% / Question / Method / Key result / Implication]
\end{abstract}
\begin{IEEEkeywords} eddy currents, finite element methods \end{IEEEkeywords}
\section{Introduction} % ~15%
\section{Formulation} % ~25%
\section{Results} % ~30%
\section{Discussion} % ~20%
\section{Conclusion} % ~5%
\bibliographystyle{IEEEtran}
\bibliography{refs}
\end{document}
""",
    "ieej_kaishi": r"""% 電気学会論文誌 skeleton — yohei + jsarticle
\documentclass[twocolumn,10pt,a4paper]{jsarticle}
\usepackage{graphicx}
\title{論文題目}
\author{著者名\thanks{正員, 近畿大学}}
\begin{document}
\maketitle
\begin{abstract}
Abstract in English. <=250 words.
\end{abstract}
\section*{キーワード}
\section{はじめに} % 序論
\section{定式化}
\section{数値結果}
\section{考察}
\section{むすび}
\bibliographystyle{junsrt}
\bibliography{refs}
\end{document}
""",
    "aps_prb": r"""% APS PRB skeleton — RevTeX 4-2
\documentclass[aps,prb,twocolumn,showpacs,amsmath,amssymb]{revtex4-2}
\usepackage{graphicx}
\begin{document}
\title{Title}
\author{First Author} \affiliation{Kindai Univ.}
\date{\today}
\begin{abstract} [<=600 words] \end{abstract}
\maketitle
\section{Introduction}
\section{Model}
\section{Results}
\section{Discussion}
\section{Conclusions}
\bibliographystyle{apsrev4-2}
\bibliography{refs}
\end{document}
""",
    "ieej_sa_pres": r"""% IEEJ SA 研究会 presentation skeleton (Beamer)
\documentclass[11pt,aspectratio=169]{beamer}
\usetheme{default}
\usepackage{graphicx}
\title{発表題目}
\author{発表者名 (近畿大学)}
\date{IEEJ 静止器・回転機合同研究会}
\begin{document}
\frame{\titlepage}
\section{背景}
\begin{frame}{背景と目的}
% Hook: 1 slide。数値 punch を入れる
\end{frame}
\section{手法}
\begin{frame}{手法}
\end{frame}
\section{結果}
\begin{frame}{結果}
\end{frame}
\section{まとめ}
\begin{frame}{まとめ (Take-home)}
% 3 bullets + next step
\end{frame}
\end{document}
""",
}


def document_meta_template_loader(kind: str = "",
                                    out_path: str = "") -> dict:
    """学術 document の定型 skeleton を返す。

    Args:
        kind: "jsps_s2" | "ieee_tmag" | "ieej_kaishi" | "aps_prb" | "ieej_sa_pres"
              空文字列で一覧を返す。
        out_path: 保存先 (省略時は template 文字列を返却のみ)。

    Returns:
        {kind, template, out_path} / {available}
    """
    if not kind:
        return {
            "available": list(_TEMPLATES.keys()),
            "usage": "document_meta_template_loader(kind='jsps_s2', out_path='...')",
        }
    if kind not in _TEMPLATES:
        return {"error": f"unknown kind={kind!r}",
                "available": list(_TEMPLATES.keys())}
    tpl = _TEMPLATES[kind]
    if out_path:
        outp = pathlib.Path(out_path).expanduser()
        outp.parent.mkdir(parents=True, exist_ok=True)
        outp.write_text(tpl, encoding="utf-8")
        return {"kind": kind, "out_path": str(outp), "bytes": len(tpl.encode("utf-8"))}
    return {"kind": kind, "template": tpl, "chars": len(tpl)}


# ---------------------------------------------------------------------------
# lint_all -- redesigned 2026-06-02 for radia-mcp.
#
# The LAB-private original hard-imported grant_writing.tools for the "common"
# checks AND the "grant" domain.  grant_writing stays private, so radia-mcp
# must not depend on it.  Here lint_all is a DECLARATIVE REGISTRY of lint tools
# that all live inside radia-mcp (presentation / paper_writing / poster),
# dispatched by lazy import.  A "grant" file still gets the generic checks plus
# a pointer to the private server; radia-mcp never imports grant_writing.
# ---------------------------------------------------------------------------

# (module, function, input-kind) -- input-kind is "text" or "path".
_COMMON_LINTS = [
    ("radia_mcp.presentation.tools", "presentation_analyze_sentences", "text"),
    ("radia_mcp.presentation.tools", "presentation_count_weak_expressions", "text"),
    ("radia_mcp.presentation.tools", "presentation_find_undefined_acronyms", "text"),
    ("radia_mcp.presentation.tools", "presentation_check_notation_variants", "text"),
]
_DOMAIN_LINTS = {
    "paper": [
        ("radia_mcp.paper_writing.tools", "paper_writing_check_imrad_balance", "text"),
        ("radia_mcp.paper_writing.tools", "paper_writing_check_passive_voice_ratio", "text"),
        ("radia_mcp.paper_writing.tools", "paper_writing_check_paragraph_length", "text"),
    ],
    "presentation": [
        ("radia_mcp.presentation.tools", "presentation_count_underlines", "path"),
    ],
    "poster": [
        ("radia_mcp.poster.tools", "poster_lint", "path"),
    ],
    # "grant" intentionally absent: grant_writing lint is LAB-private.
}
# Domains whose dedicated linter ships only in the private mcp-server-document.
_PRIVATE_DOMAIN_NOTE = {
    "grant": ("grant_writing lint is LAB-private (mcp-server-document) and is "
              "not shipped in radia-mcp; only the generic/common checks above "
              "were run."),
}


def _run_one_lint(modpath: str, fname: str, kind: str,
                  text: str, path: "pathlib.Path"):
    """Lazily import modpath and call fname on the text or the path.

    Raises ImportError if the subpackage is unavailable (e.g. radia-mcp[document]
    not installed) and lets the caller record it as 'unavailable'.
    """
    import importlib
    mod = importlib.import_module(modpath)
    fn = getattr(mod, fname)
    if kind == "text":
        return fn(text)
    if kind == "path":
        return fn(str(path))
    raise ValueError(f"unknown input kind {kind!r}")


def document_meta_lint_all(path: str,
                             domain: str = "",
                             page_limit: int = 0) -> dict:
    """Run every applicable radia-mcp lint over one text / TeX file.

    Declarative registry (``_COMMON_LINTS`` + ``_DOMAIN_LINTS``) of lints that
    all live inside radia-mcp (presentation / paper_writing / poster),
    dispatched by lazy import -- so this has NO dependency on the LAB-private
    grant_writing.  A "grant" file still gets the generic checks plus a pointer
    to the private server.

    Args:
        path: .tex / .md / .txt file.
        domain: "paper" | "presentation" | "poster" | "grant" | ""
                ("" auto-detects from content).
        page_limit: if > 0 and a sibling .pdf exists, also report its page count.

    Returns:
        dict with domain_detected, one entry per lint that ran, lints_run,
        lints_unavailable, optional note (private domain), and a summary.
    """
    p = pathlib.Path(path).expanduser()
    if not p.exists():
        return {"error": f"file not found: {p}"}
    try:
        text = p.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        return {"error": f"decode error: {exc}"}

    # Auto-detect domain from content.
    if not domain:
        low = text.lower()
        if "\\documentclass{beamer}" in text or "\\begin{frame}" in text:
            domain = "presentation"
        elif "\\sectionbox" in text or "JSPS" in text or "申請" in text:
            domain = "grant"
        elif "tikzposter" in low or "a0poster" in low or "betterposter" in low:
            domain = "poster"
        elif "ieeetran" in low or "\\abstract" in low or "revtex" in low:
            domain = "paper"
        else:
            domain = "paper"  # sensible default

    findings: dict = {"domain_detected": domain}
    ran: list = []
    unavailable: list = []

    for modpath, fname, kind in _COMMON_LINTS + _DOMAIN_LINTS.get(domain, []):
        try:
            findings[fname] = _run_one_lint(modpath, fname, kind, text, p)
            ran.append(fname)
        except ImportError as exc:
            unavailable.append(f"{fname} ({modpath}: {exc})")
        except Exception as exc:  # the lint ran but errored on this input
            findings[f"error_{fname}"] = str(exc)

    if domain in _PRIVATE_DOMAIN_NOTE:
        findings["note"] = _PRIVATE_DOMAIN_NOTE[domain]

    # Optional PDF page-limit check -- self-contained via PyMuPDF.
    if page_limit > 0:
        pdf_path = p.with_suffix(".pdf")
        if pdf_path.exists():
            try:
                import fitz
                npages = fitz.open(str(pdf_path)).page_count
                findings["pdf_pages"] = {"pages": npages, "limit": page_limit,
                                          "ok": npages <= page_limit}
            except Exception as exc:
                findings["error_pdf_pages"] = str(exc)

    # Roll-up of problem categories.
    problems = 0
    for val in findings.values():
        if isinstance(val, dict):
            if val.get("undefined_count", 0) > 0:
                problems += 1
            if val.get("weak_count", 0) > 0:
                problems += 1
            if val.get("deviations"):
                problems += 1
            if val.get("total_findings", 0) > 0:
                problems += 1
    findings["lints_run"] = ran
    findings["lints_unavailable"] = unavailable
    findings["summary"] = {"problem_categories": problems,
                            "file": str(p),
                            "chars": len(text),
                            "n_lints_run": len(ran)}
    return findings
