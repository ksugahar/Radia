"""Tool functions for the document_meta sub-module.

Cross-cutting helpers that do not belong to a specific document type
(grant / paper / presentation): deadline tracking, version diffing,
template scaffolding, and full-tree linting.
"""
from __future__ import annotations

import ast
import datetime
import difflib
import hashlib
import json
import pathlib
import re
import subprocess


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
# checks AND the "grant" domain.  grant_writing is now a first-class public
# radia-mcp subpackage, so lint_all is a DECLARATIVE REGISTRY over the document
# family (grant_writing / presentation / paper_writing / poster), dispatched by
# lazy import.
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
    "grant": [
        ("radia_mcp.grant_writing.tools", "grant_writing_section_presence", "text"),
        ("radia_mcp.grant_writing.tools", "grant_writing_budget_alignment_check", "text"),
        ("radia_mcp.grant_writing.tools", "grant_writing_health_report", "text"),
    ],
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
    all live inside radia-mcp (grant_writing / presentation / paper_writing /
    poster), dispatched by lazy import.

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
            if val.get("total_weak_expressions", 0) > 0:
                problems += 1
            if val.get("missing_count", 0) > 0:
                problems += 1
            if val.get("missing_required_count", 0) > 0:
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


# ---------------------------------------------------------------------------
# Radia repository documentation/notebook audit helpers.
# ---------------------------------------------------------------------------

_AUDIT_SKIP_PARTS = {
    ".git", ".claude", ".pytest_cache", ".ipynb_checkpoints",
    "__pycache__", "build", "build-msvc",
}


def _resolve_repo_root(repo_root: str = "") -> pathlib.Path:
    """Resolve a Radia-style repo root from an explicit path, cwd, or this file."""
    if repo_root:
        return pathlib.Path(repo_root).expanduser().resolve()

    starts: list[pathlib.Path] = [pathlib.Path.cwd()]
    starts.extend(pathlib.Path(__file__).resolve().parents)

    seen: set[pathlib.Path] = set()
    for start in starts:
        try:
            p = start.resolve()
        except OSError:
            p = start
        for cand in (p, *p.parents):
            if cand in seen:
                continue
            seen.add(cand)
            if (
                (cand / "pyproject.toml").is_file()
                and (cand / "docs").is_dir()
                and (cand / "src" / "radia").is_dir()
            ):
                return cand
    return pathlib.Path(".").resolve()


def _rel(path: pathlib.Path, root: pathlib.Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except Exception:
        return path.as_posix()




def _is_git_ignored(path: pathlib.Path, repo_root: pathlib.Path) -> bool:
    """Return True when git ignores path.  False outside a git worktree."""
    if not (repo_root / ".git").exists():
        return False
    try:
        rel = _rel(path, repo_root)
        result = subprocess.run(
            ["git", "-C", str(repo_root), "check-ignore", "--quiet", "--", rel],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return result.returncode == 0
    except Exception:
        return False


def _git_tracked_paths(repo_root: pathlib.Path,
                       pathspecs: list[str]) -> set[str] | None:
    """Return tracked paths matching pathspecs, or None outside a git checkout."""
    if not (repo_root / ".git").exists():
        return None
    try:
        result = subprocess.run(
            ["git", "ls-files", "--", *pathspecs],
            cwd=str(repo_root),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    return {line.strip().replace("\\", "/") for line in result.stdout.splitlines() if line.strip()}


def _tracked_notebook_set(repo_root: pathlib.Path,
                          scan_root: pathlib.Path) -> set[str] | None:
    scan_rel = _rel(scan_root, repo_root).rstrip("/")
    if scan_rel in {"", "."}:
        pathspecs = ["*.ipynb", "**/*.ipynb"]
    else:
        pathspecs = [f"{scan_rel}/*.ipynb", f"{scan_rel}/**/*.ipynb"]
    return _git_tracked_paths(repo_root, pathspecs)


def _iter_notebooks(scan_root: pathlib.Path,
                    repo_root: pathlib.Path,
                    include_gitignored: bool,
                    tracked_only: bool = False):
    tracked = _tracked_notebook_set(repo_root, scan_root) if tracked_only else None
    for p in sorted(scan_root.rglob("*.ipynb")):
        if {part.lower() for part in p.parts} & _AUDIT_SKIP_PARTS:
            continue
        if tracked is not None and _rel(p, repo_root) not in tracked:
            continue
        if not include_gitignored and _is_git_ignored(p, repo_root):
            continue
        yield p


def _read_notebook_result_summary(path: pathlib.Path, repo_root: pathlib.Path) -> dict:
    """Summarise whether an ipynb has saved execution results."""
    try:
        nb = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {
            "path": _rel(path, repo_root),
            "error": f"notebook read/parse failed: {exc}",
            "result_saved": False,
        }

    metadata = nb.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    radia_metadata = metadata.get("radia", {})
    if not isinstance(radia_metadata, dict):
        radia_metadata = {}
    notebook_role = str(radia_metadata.get("notebook_role", "")).strip().lower()
    webgui_required = (
        notebook_role == "example"
        or radia_metadata.get("webgui_required") is True
    )
    webgui_field_required = radia_metadata.get("webgui_field_required") is True

    cells = nb.get("cells", [])
    code_cells = [c for c in cells if c.get("cell_type") == "code"]
    executed = [c for c in code_cells if c.get("execution_count") is not None]
    output_cells = [c for c in code_cells if c.get("outputs")]
    code_sources = []
    for cell in code_cells:
        source = cell.get("source", "")
        code_sources.append("".join(source) if isinstance(source, list) else str(source))

    webgui_import_pattern = re.compile(
        r"(?m)^\s*(?:from\s+(?:ngsolve|netgen)\.webgui\s+import\s+[^\n]*\bDraw\b"
        r"|import\s+(?:ngsolve|netgen)\.webgui(?:\s+as\s+\w+)?)"
    )
    webgui_draw_pattern = re.compile(
        r"\b(?:Draw|(?:ngsolve|netgen)\.webgui\.Draw|\w+\.Draw)\s*\("
    )
    webgui_import_present = any(
        webgui_import_pattern.search(source) for source in code_sources
    )
    webgui_draw_indices = [
        idx for idx, source in enumerate(code_sources)
        if webgui_draw_pattern.search(source)
    ]
    webgui_draw_present = webgui_import_present and bool(webgui_draw_indices)

    parameterized_field_draw_indices: list[int] = []
    for idx, source in enumerate(code_sources):
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            function = node.func
            is_draw = (
                isinstance(function, ast.Name) and function.id == "Draw"
            ) or (
                isinstance(function, ast.Attribute) and function.attr == "Draw"
            )
            if not is_draw:
                continue
            keyword_names = {kw.arg for kw in node.keywords if kw.arg is not None}
            if len(node.args) >= 2 and "name" in keyword_names:
                parameterized_field_draw_indices.append(idx)
                break
    parameterized_field_draw_indices = sorted(set(parameterized_field_draw_indices))
    parameterized_field_draw_present = (
        webgui_import_present and bool(parameterized_field_draw_indices)
    )

    def _has_webgui_output(cell: dict) -> bool:
        for output in cell.get("outputs", []):
            if output.get("output_type") not in {"display_data", "execute_result"}:
                continue
            data = output.get("data", {})
            if not isinstance(data, dict):
                continue
            widget = data.get("application/vnd.jupyter.widget-view+json")
            if isinstance(widget, dict) and widget.get("model_id"):
                return True
            html = data.get("text/html", "")
            html_text = "".join(html) if isinstance(html, list) else str(html)
            if "webgui" in html_text.lower() or "netgen" in html_text.lower():
                return True
        return False

    executed_webgui_draw_cells = [
        code_cells[idx] for idx in webgui_draw_indices
        if code_cells[idx].get("execution_count") is not None
    ]
    webgui_draw_cells_with_rich_outputs = [
        cell for cell in executed_webgui_draw_cells if _has_webgui_output(cell)
    ]
    executed_parameterized_field_draw_cells = [
        code_cells[idx] for idx in parameterized_field_draw_indices
        if code_cells[idx].get("execution_count") is not None
    ]
    parameterized_field_draw_cells_with_rich_outputs = [
        cell for cell in executed_parameterized_field_draw_cells
        if _has_webgui_output(cell)
    ]
    webgui_ready = (
        not webgui_required
        or (
            webgui_draw_present
            and bool(executed_webgui_draw_cells)
            and bool(webgui_draw_cells_with_rich_outputs)
        )
    )
    webgui_field_ready = (
        not webgui_field_required
        or (
            parameterized_field_draw_present
            and bool(executed_parameterized_field_draw_cells)
            and bool(parameterized_field_draw_cells_with_rich_outputs)
        )
    )
    error_outputs = []
    for idx, cell in enumerate(code_cells):
        for out in cell.get("outputs", []):
            if out.get("output_type") == "error":
                error_outputs.append({
                    "code_cell_index": idx,
                    "ename": out.get("ename", ""),
                    "evalue": out.get("evalue", ""),
                })

    result_saved = bool(code_cells and executed and output_cells and not error_outputs)
    return {
        "path": _rel(path, repo_root),
        "cells": len(cells),
        "code_cells": len(code_cells),
        "executed_code_cells": len(executed),
        "code_cells_with_outputs": len(output_cells),
        "output_count": sum(len(c.get("outputs", [])) for c in code_cells),
        "error_output_count": len(error_outputs),
        "error_outputs": error_outputs[:5],
        "result_saved": result_saved,
        "notebook_role": notebook_role,
        "webgui_required": webgui_required,
        "webgui_field_required": webgui_field_required,
        "webgui_import_present": webgui_import_present,
        "webgui_draw_present": webgui_draw_present,
        "webgui_draw_cell_count": len(webgui_draw_indices),
        "executed_webgui_draw_cell_count": len(executed_webgui_draw_cells),
        "webgui_draw_cells_with_rich_outputs": len(webgui_draw_cells_with_rich_outputs),
        "webgui_ready": webgui_ready,
        "parameterized_field_draw_present": parameterized_field_draw_present,
        "parameterized_field_draw_cell_count": len(parameterized_field_draw_indices),
        "executed_parameterized_field_draw_cell_count": len(
            executed_parameterized_field_draw_cells
        ),
        "parameterized_field_draw_cells_with_rich_outputs": len(
            parameterized_field_draw_cells_with_rich_outputs
        ),
        "webgui_field_ready": webgui_field_ready,
    }


def _collect_json_keys(obj) -> set[str]:
    keys: set[str] = set()
    if isinstance(obj, dict):
        for key, val in obj.items():
            keys.add(str(key))
            keys.update(_collect_json_keys(val))
    elif isinstance(obj, list):
        for val in obj[:100]:
            keys.update(_collect_json_keys(val))
    return keys


def _result_json_summary(path: pathlib.Path,
                         repo_root: pathlib.Path,
                         notebook_path: pathlib.Path | None = None) -> dict:
    date_keys = {
        "generated_at_utc", "generated_at", "executed_at_utc", "executed_at",
        "timestamp", "date", "created_at", "updated_at",
    }
    version_keys = {
        "version", "schema_version", "radia_version", "python_version",
        "runtime_version", "package_versions", "versions", "platform",
    }
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {
            "path": _rel(path, repo_root),
            "parse_ok": False,
            "error": str(exc),
            "has_date": False,
            "has_version": False,
            "has_version_date": False,
            "has_notebook_sha256": False,
            "notebook_sha256_matches": False,
        }
    keys = _collect_json_keys(data)
    has_date = bool(keys & date_keys)
    has_version = bool(keys & version_keys)
    notebook_sha = data.get("notebook_sha256") if isinstance(data, dict) else None
    has_notebook_sha = isinstance(notebook_sha, str) and bool(notebook_sha)
    notebook_sha_matches = False
    if has_notebook_sha and notebook_path is not None and notebook_path.exists():
        try:
            notebook_sha_matches = notebook_sha == _sha256_bytes(notebook_path.read_bytes())
        except Exception:
            notebook_sha_matches = False
    return {
        "path": _rel(path, repo_root),
        "parse_ok": True,
        "has_date": has_date,
        "has_version": has_version,
        "has_version_date": has_date and has_version,
        "has_notebook_sha256": has_notebook_sha,
        "notebook_sha256_matches": notebook_sha_matches,
        "matched_date_keys": sorted(keys & date_keys),
        "matched_version_keys": sorted(keys & version_keys),
    }


def _find_notebook_result_jsons(nb_path: pathlib.Path) -> list[pathlib.Path]:
    parent = nb_path.parent
    stem = nb_path.stem.lower()
    candidates: set[pathlib.Path] = set()
    for p in parent.glob("*.json"):
        name = p.stem.lower()
        if (
            name == stem
            or name.startswith(stem + "_")
            or name.startswith(stem + "-")
            or "result" in name
            or "summary" in name
            or "evidence" in name
            or "validation" in name
        ):
            candidates.add(p)
    return sorted(candidates)




def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()












def document_meta_notebook_result_audit(repo_root: str = "",
                                          notebook_root: str = "",
                                          require_json: bool = False,
                                          include_gitignored: bool = False,
                                          tracked_only: bool = True,
                                          max_items: int = 50) -> dict:
    """Audit docs notebooks for saved results and WebGUI scenes.

    Project policy:

    * every method/showcase ``docs/**/*.ipynb`` should be a result-saving
      notebook, meaning it has executed code cells with saved outputs and no
      saved error outputs;
    * docs notebooks are public demonstrations, not benchmark evidence. They
      do not require adjacent JSON or a runtime threshold. Numerical evidence
      and benchmark JSON belong under ``validation_test/``;
    * notebooks marked ``metadata.radia.notebook_role = "example"`` or
      ``metadata.radia.webgui_required = true`` must contain an executed
      ``ngsolve.webgui.Draw``/``netgen.webgui.Draw`` cell with saved rich
      output.
    * notebooks marked ``metadata.radia.webgui_field_required = true`` must
      contain an executed ``Draw(field, mesh, name=..., ...)`` call with saved
      rich output. A bare one-argument Draw call is not a field scene.

    Args:
        repo_root: Radia repository root.  Empty means auto-detect from cwd.
        notebook_root: Optional subdirectory to scan, absolute or repo-relative.
            Empty defaults to ``docs/`` when it exists.
        require_json: Optional compatibility diagnostic for an explicitly
            sidecar-managed corpus. False for the Radia docs policy.
        include_gitignored: Include notebooks ignored by git (LAB-local docs,
            work-in-progress notebook galleries).  Defaults false for public
            repo policy checks.
        tracked_only: In a git checkout, scan only tracked notebooks by default.
            Set false for a deliberate WIP-tree audit.
        max_items: Maximum detailed gap rows to include.

    Returns:
        A repository-wide notebook result-persistence report.
    """
    root = _resolve_repo_root(repo_root)
    if notebook_root:
        scan_root = pathlib.Path(notebook_root).expanduser()
        if not scan_root.is_absolute():
            scan_root = root / scan_root
    elif (root / "docs").is_dir():
        scan_root = root / "docs"
    else:
        scan_root = root
    if not scan_root.exists():
        return {"error": f"notebook_root not found: {scan_root}", "repo_root": str(root)}

    notebooks = list(_iter_notebooks(
        scan_root,
        root,
        include_gitignored,
        tracked_only=tracked_only,
    ))
    rows: list[dict] = []
    gaps: list[dict] = []
    for nb in sorted(notebooks):
        nb_summary = _read_notebook_result_summary(nb, root)
        sidecars = (
            [
                _result_json_summary(p, root, nb)
                for p in _find_notebook_result_jsons(nb)
            ]
            if require_json
            else []
        )
        good_json_count = sum(1 for s in sidecars if s.get("has_version_date"))
        synced_json_count = sum(
            1 for s in sidecars
            if s.get("has_version_date") and s.get("notebook_sha256_matches")
        )
        if nb_summary.get("error"):
            status = "notebook_parse_error"
        elif not nb_summary.get("code_cells"):
            status = "no_code_cells"
        elif not nb_summary.get("result_saved"):
            status = "needs_saved_outputs"
        elif (nb_summary.get("webgui_required")
              and not nb_summary.get("webgui_draw_present")):
            status = "needs_webgui_draw"
        elif (nb_summary.get("webgui_required")
              and not nb_summary.get("webgui_ready")):
            status = "needs_executed_webgui_output"
        elif (nb_summary.get("webgui_field_required")
              and not nb_summary.get("parameterized_field_draw_present")):
            status = "needs_parameterized_webgui_field_draw"
        elif (nb_summary.get("webgui_field_required")
              and not nb_summary.get("webgui_field_ready")):
            status = "needs_executed_parameterized_webgui_field_output"
        elif require_json and not sidecars:
            status = "needs_result_json_sidecar"
        elif require_json and good_json_count == 0:
            status = "result_json_missing_version_or_date"
        elif require_json and synced_json_count == 0:
            status = "result_json_not_synced_to_notebook"
        else:
            status = "ok_result_saved"

        row = {
            **nb_summary,
            "status": status,
            "result_json_count": len(sidecars),
            "result_json_with_version_date_count": good_json_count,
            "result_json_synced_to_notebook_count": synced_json_count,
            "result_jsons": sidecars[:10],
        }
        rows.append(row)
        if status != "ok_result_saved" and len(gaps) < max_items:
            gaps.append(row)

    summary = {
        "repo_root": str(root),
        "scan_root": _rel(scan_root, root),
        "notebooks_scanned": len(rows),
        "ok_result_saved": sum(1 for r in rows if r["status"] == "ok_result_saved"),
        "needs_saved_outputs": sum(1 for r in rows if r["status"] == "needs_saved_outputs"),
        "webgui_required": sum(1 for r in rows if r.get("webgui_required")),
        "webgui_ready": sum(
            1 for r in rows if r.get("webgui_required") and r.get("webgui_ready")
        ),
        "webgui_field_required": sum(
            1 for r in rows if r.get("webgui_field_required")
        ),
        "webgui_field_ready": sum(
            1 for r in rows
            if r.get("webgui_field_required") and r.get("webgui_field_ready")
        ),
        "needs_webgui_draw": sum(
            1 for r in rows if r["status"] == "needs_webgui_draw"
        ),
        "needs_executed_webgui_output": sum(
            1 for r in rows if r["status"] == "needs_executed_webgui_output"
        ),
        "needs_parameterized_webgui_field_draw": sum(
            1 for r in rows
            if r["status"] == "needs_parameterized_webgui_field_draw"
        ),
        "needs_executed_parameterized_webgui_field_output": sum(
            1 for r in rows
            if r["status"] == "needs_executed_parameterized_webgui_field_output"
        ),
        "needs_result_json_sidecar": sum(
            1 for r in rows if r["status"] == "needs_result_json_sidecar"
        ),
        "result_json_missing_version_or_date": sum(
            1 for r in rows if r["status"] == "result_json_missing_version_or_date"
        ),
        "result_json_not_synced_to_notebook": sum(
            1 for r in rows if r["status"] == "result_json_not_synced_to_notebook"
        ),
        "no_code_cells": sum(1 for r in rows if r["status"] == "no_code_cells"),
        "parse_errors": sum(1 for r in rows if r["status"] == "notebook_parse_error"),
    }
    return {
        "policy": {
            "scope": "tracked docs/**/*.ipynb method/showcase notebooks",
            "notebook": "docs ipynb files should store executed outputs",
            "example_webgui": "example notebooks must store executed ngsolve.webgui.Draw or netgen.webgui.Draw rich output",
            "field_webgui": "field notebooks must store Draw(field, mesh, name=..., ...) rich output with explicit view arguments",
            "json": "docs notebooks require no sidecar; validation_test owns machine-readable evidence JSON",
            "sync": "existing notebook-linked JSON is informational unless require_json is explicitly enabled",
            "tracked_only_default": "true; set tracked_only=false for explicit WIP-tree audits",
        },
        "summary": summary,
        "gaps": gaps,
        "gaps_truncated": len(gaps) < sum(1 for r in rows if r["status"] != "ok_result_saved"),
        "notebooks": rows[:max_items],
        "notebooks_truncated": len(rows) > max_items,
    }
