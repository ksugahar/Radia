"""Tool functions for the document_meta sub-module.

Cross-cutting helpers that do not belong to a specific document type
(grant / paper / presentation): deadline tracking, version diffing,
template scaffolding, and full-tree linting.
"""
from __future__ import annotations

import datetime
import difflib
import hashlib
import json
import pathlib
import platform
import re
import subprocess
import sys


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


# ---------------------------------------------------------------------------
# Radia repo documentation/notebook migration helpers.
#
# These are intentionally lightweight filesystem audits.  They do not execute
# notebooks or move files; they give agents a repeatable first pass before
# promoting examples into result-bearing docs notebooks or relocating panel
# surfaces.
# ---------------------------------------------------------------------------

_AUDIT_TEXT_SUFFIXES = {
    ".py", ".pyi", ".md", ".rst", ".txt", ".toml", ".json", ".ps1", ".yml", ".yaml"
}
_AUDIT_SKIP_PARTS = {
    ".git", ".claude", ".pytest_cache", ".ipynb_checkpoints",
    "__pycache__", "build", "build-msvc",
}


def _resolve_repo_root(repo_root: str = "") -> pathlib.Path:
    """Resolve a Radia-style repo root from an explicit path, cwd, or this file."""
    starts: list[pathlib.Path] = []
    if repo_root:
        starts.append(pathlib.Path(repo_root).expanduser())
    starts.append(pathlib.Path.cwd())
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
            if (cand / "examples").is_dir() and (cand / "docs").is_dir():
                return cand
    return pathlib.Path(repo_root or ".").expanduser().resolve()


def _rel(path: pathlib.Path, root: pathlib.Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except Exception:
        return path.as_posix()


def _iter_files(root: pathlib.Path, suffixes: set[str]):
    if not root.exists():
        return
    for p in root.rglob("*"):
        if p.is_file() and p.suffix.lower() in suffixes:
            parts = {part.lower() for part in p.parts}
            if parts & _AUDIT_SKIP_PARTS:
                continue
            yield p


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


def _iter_notebooks(scan_root: pathlib.Path,
                    repo_root: pathlib.Path,
                    include_gitignored: bool):
    for p in sorted(scan_root.rglob("*.ipynb")):
        if {part.lower() for part in p.parts} & _AUDIT_SKIP_PARTS:
            continue
        if not include_gitignored and _is_git_ignored(p, repo_root):
            continue
        yield p


def _count_tree_files(root: pathlib.Path, pattern: str) -> int:
    if not root.exists():
        return 0
    return sum(1 for p in root.rglob(pattern) if p.is_file())


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

    cells = nb.get("cells", [])
    code_cells = [c for c in cells if c.get("cell_type") == "code"]
    executed = [c for c in code_cells if c.get("execution_count") is not None]
    output_cells = [c for c in code_cells if c.get("outputs")]
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


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _coerce_output_text(value) -> str:
    if isinstance(value, list):
        return "".join(str(v) for v in value)
    return str(value)


def _summarize_notebook_output(output: dict,
                               include_output_text: bool,
                               max_output_chars: int) -> dict:
    out_type = output.get("output_type", "")
    summary = {"output_type": out_type}
    if out_type == "stream":
        text = _coerce_output_text(output.get("text", ""))
        summary.update({
            "name": output.get("name", ""),
            "text_chars": len(text),
            "text_sha256": _sha256_text(text),
        })
        if include_output_text:
            summary["text_preview"] = text[:max_output_chars]
            summary["text_truncated"] = len(text) > max_output_chars
        return summary

    if out_type in {"display_data", "execute_result"}:
        data = output.get("data", {})
        summary["data_keys"] = sorted(data.keys())
        if "text/plain" in data:
            text = _coerce_output_text(data["text/plain"])
            summary.update({
                "text_plain_chars": len(text),
                "text_plain_sha256": _sha256_text(text),
            })
            if include_output_text:
                summary["text_plain_preview"] = text[:max_output_chars]
                summary["text_plain_truncated"] = len(text) > max_output_chars
        image_hashes = {}
        for key in ("image/png", "image/jpeg", "image/svg+xml", "text/html"):
            if key in data:
                raw = _coerce_output_text(data[key])
                image_hashes[key] = {
                    "chars": len(raw),
                    "sha256": _sha256_text(raw),
                }
        if image_hashes:
            summary["rich_output_hashes"] = image_hashes
        return summary

    if out_type == "error":
        traceback_text = "\n".join(output.get("traceback", []))
        summary.update({
            "ename": output.get("ename", ""),
            "evalue": output.get("evalue", ""),
            "traceback_sha256": _sha256_text(traceback_text),
            "traceback_lines": len(output.get("traceback", [])),
        })
        return summary

    summary["raw_keys"] = sorted(output.keys())
    return summary


def _runtime_versions() -> dict:
    versions = {
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
    }
    try:
        import importlib.metadata as importlib_metadata
        versions["radia_mcp_version"] = importlib_metadata.version("radia-mcp")
    except Exception:
        versions["radia_mcp_version"] = None
    try:
        import radia
        versions["radia_version"] = getattr(radia, "__version__", None)
    except Exception:
        versions["radia_version"] = None
    return versions


def document_meta_write_notebook_result_json(notebook_path: str,
                                               output_path: str = "",
                                               overwrite: bool = False,
                                               include_output_text: bool = True,
                                               max_output_chars: int = 4000) -> dict:
    """Write a durable JSON sidecar summarising a saved-result notebook.

    This does **not** execute the notebook.  It records the notebook's current
    saved outputs, output hashes/previews, execution counts, and runtime/version
    metadata in adjacent JSON for future debugging.  Use it after executing a
    notebook with ``jupyter nbconvert --execute --inplace`` or an equivalent
    workflow.

    Args:
        notebook_path: Path to an executed ``.ipynb``.
        output_path: Optional JSON output path.  Defaults to
            ``<notebook stem>_result.json`` beside the notebook.
        overwrite: If false, fail rather than overwrite an existing sidecar.
        include_output_text: Store text/plain and stdout previews in the JSON.
            Rich media are hash-only.
        max_output_chars: Per-output text preview limit.

    Returns:
        Dict with the output path, summary counts, and policy status.
    """
    p = pathlib.Path(notebook_path).expanduser()
    if not p.exists():
        return {"error": f"notebook not found: {p}"}
    if p.suffix.lower() != ".ipynb":
        return {"error": f"not an ipynb: {p}"}
    try:
        raw = p.read_text(encoding="utf-8")
        nb = json.loads(raw)
    except Exception as exc:
        return {"error": f"notebook read/parse failed: {exc}", "path": str(p)}

    out = pathlib.Path(output_path).expanduser() if output_path else p.with_name(f"{p.stem}_result.json")
    if out.exists() and not overwrite:
        return {
            "error": f"output exists (pass overwrite=True): {out}",
            "output_path": str(out),
        }

    root = _resolve_repo_root(str(p.parent))
    nb_summary = _read_notebook_result_summary(p, root)
    cells = nb.get("cells", [])
    code_outputs: list[dict] = []
    for cell_index, cell in enumerate(cells):
        if cell.get("cell_type") != "code":
            continue
        cell_outputs = [
            _summarize_notebook_output(o, include_output_text, max_output_chars)
            for o in cell.get("outputs", [])
        ]
        code_outputs.append({
            "cell_index": cell_index,
            "execution_count": cell.get("execution_count"),
            "source_sha256": _sha256_text(_coerce_output_text(cell.get("source", ""))),
            "source_preview": _coerce_output_text(cell.get("source", ""))[:1000],
            "outputs": cell_outputs,
        })

    payload = {
        "schema": "radia.notebook_result.v1",
        "generated_at_utc": datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "notebook": _rel(p, root),
        "notebook_sha256": _sha256_bytes(p.read_bytes()),
        "versions": _runtime_versions(),
        "summary": {
            **{k: v for k, v in nb_summary.items() if k != "error_outputs"},
            "metadata_keys": sorted(nb.get("metadata", {}).keys()),
        },
        "outputs": code_outputs,
        "policy": {
            "result_saved": bool(nb_summary.get("result_saved")),
            "json_has_generated_at_utc": True,
            "json_has_versions": True,
        },
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    return {
        "output_path": str(out),
        "notebook": str(p),
        "result_saved": payload["policy"]["result_saved"],
        "code_cells": payload["summary"].get("code_cells", 0),
        "outputs_recorded": sum(len(c["outputs"]) for c in code_outputs),
        "generated_at_utc": payload["generated_at_utc"],
        "versions": payload["versions"],
    }


def document_meta_write_docs_notebook_result_jsons(repo_root: str = "",
                                                     notebook_root: str = "docs",
                                                     overwrite: bool = True,
                                                     dry_run: bool = False,
                                                     max_notebooks: int = 0,
                                                     include_gitignored: bool = False,
                                                     include_output_text: bool = True,
                                                     max_output_chars: int = 4000) -> dict:
    """Batch-write synchronized result JSON sidecars for executed docs notebooks.

    This is the mechanical first pass for the docs-notebook policy.  It scans
    ``docs/**/*.ipynb`` by default, skips notebooks that do not yet have saved
    execution outputs, and writes ``<notebook stem>_result.json`` sidecars with
    ``generated_at_utc``, runtime versions, output hashes/previews, and
    ``notebook_sha256``.  Domain-specific detailed JSON can still live beside
    this manifest; the manifest is the generic synchronization/debug record.

    Args:
        repo_root: Radia repository root.  Empty means auto-detect from cwd.
        notebook_root: Subdirectory to scan. Defaults to ``docs``.
        overwrite: Refresh existing generated sidecars.
        dry_run: Report planned writes without changing files.
        max_notebooks: Optional cap for incremental batches (0 = no cap).
        include_gitignored: Include notebooks ignored by git. Defaults false so
            public docs policy batches do not create ignored sidecars under
            LAB-local research trees.
        include_output_text: Store stdout/text/plain previews in sidecars.
        max_output_chars: Per-output text preview limit.

    Returns:
        Dict with written/skipped/error rows and an audit summary.
    """
    root = _resolve_repo_root(repo_root)
    scan_root = pathlib.Path(notebook_root).expanduser()
    if not scan_root.is_absolute():
        scan_root = root / scan_root
    if not scan_root.exists():
        return {"error": f"notebook_root not found: {scan_root}", "repo_root": str(root)}

    notebooks = list(_iter_notebooks(scan_root, root, include_gitignored))
    written: list[dict] = []
    skipped: list[dict] = []
    errors: list[dict] = []
    planned = 0
    for nb in notebooks:
        summary = _read_notebook_result_summary(nb, root)
        if not summary.get("result_saved"):
            skipped.append({
                "path": _rel(nb, root),
                "reason": "notebook_has_no_saved_outputs",
                "executed_code_cells": summary.get("executed_code_cells", 0),
                "code_cells_with_outputs": summary.get("code_cells_with_outputs", 0),
            })
            continue
        if max_notebooks and planned >= max_notebooks:
            skipped.append({"path": _rel(nb, root), "reason": "max_notebooks_reached"})
            continue
        planned += 1
        out = nb.with_name(f"{nb.stem}_result.json")
        if dry_run:
            written.append({
                "notebook": _rel(nb, root),
                "output_path": _rel(out, root),
                "dry_run": True,
            })
            continue
        result = document_meta_write_notebook_result_json(
            str(nb),
            output_path=str(out),
            overwrite=overwrite,
            include_output_text=include_output_text,
            max_output_chars=max_output_chars,
        )
        if "error" in result:
            errors.append({"notebook": _rel(nb, root), **result})
        else:
            written.append({
                **result,
                "notebook": _rel(nb, root),
                "output_path": _rel(pathlib.Path(result["output_path"]), root),
            })

    return {
        "policy": {
            "scope": "tracked docs/**/*.ipynb result-bearing notebooks",
            "sync": "each written JSON records notebook_sha256 for the current ipynb",
        },
        "repo_root": str(root),
        "scan_root": _rel(scan_root, root),
        "dry_run": dry_run,
        "summary": {
            "notebooks_scanned": len(notebooks),
            "written_or_planned": len(written),
            "skipped": len(skipped),
            "errors": len(errors),
        },
        "written": written,
        "skipped": skipped[:50],
        "skipped_truncated": len(skipped) > 50,
        "errors": errors,
    }


def document_meta_notebook_result_audit(repo_root: str = "",
                                          notebook_root: str = "",
                                          require_json: bool = True,
                                          include_gitignored: bool = False,
                                          max_items: int = 50) -> dict:
    """Audit docs notebooks for saved results and synchronized result JSON.

    Project policy:

    * every method/showcase ``docs/**/*.ipynb`` should be a result-saving
      notebook, meaning it has executed code cells with saved outputs and no
      saved error outputs;
    * if the docs notebook computes or verifies something, keep durable JSON
      beside it with date/version metadata for future debugging.  Preferred
      keys are ``generated_at_utc`` plus ``radia_version``/``python_version``
      or a ``versions``/``package_versions`` object;
    * the JSON and result-bearing notebook are a synchronized pair, so the JSON
      records ``notebook_sha256`` matching the current ``.ipynb``.

    Args:
        repo_root: Radia repository root.  Empty means auto-detect from cwd.
        notebook_root: Optional subdirectory to scan, absolute or repo-relative.
            Empty defaults to ``docs/`` when it exists.
        require_json: If true, result-saving notebooks without a matching result
            JSON sidecar are reported as gaps.
        include_gitignored: Include notebooks ignored by git (LAB-local docs,
            work-in-progress notebook galleries).  Defaults false for public
            repo policy checks.
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

    notebooks = list(_iter_notebooks(scan_root, root, include_gitignored))
    rows: list[dict] = []
    gaps: list[dict] = []
    for nb in sorted(notebooks):
        nb_summary = _read_notebook_result_summary(nb, root)
        sidecars = [
            _result_json_summary(p, root, nb)
            for p in _find_notebook_result_jsons(nb)
        ]
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
            "json": "computed docs results should have adjacent JSON with generated_at_utc and version/runtime metadata",
            "sync": "the JSON must include notebook_sha256 matching the current result-bearing ipynb",
        },
        "summary": summary,
        "gaps": gaps,
        "gaps_truncated": len(gaps) < sum(1 for r in rows if r["status"] != "ok_result_saved"),
        "notebooks": rows[:max_items],
        "notebooks_truncated": len(rows) > max_items,
    }


def _find_topic_notebooks(repo_root: pathlib.Path, topic: str) -> list[pathlib.Path]:
    docs = repo_root / "docs"
    matches: list[pathlib.Path] = []
    direct = docs / topic
    if direct.exists():
        matches.extend(p for p in direct.rglob("*.ipynb") if p.is_file())
    from_examples = docs / "notebooks" / "from_examples" / topic
    if from_examples.exists():
        matches.extend(p for p in from_examples.rglob("*.ipynb") if p.is_file())

    # Some historic notebooks use a broader docs directory while retaining the
    # topic token in the filename/path.  Keep this conservative to avoid
    # accidentally claiming unrelated notebooks.
    token = topic.lower().replace("_", "-")
    for p in docs.rglob("*.ipynb"):
        rel = _rel(p, repo_root).lower().replace("_", "-")
        if token in rel and p not in matches:
            matches.append(p)
    return sorted(matches)


def document_meta_examples_notebook_audit(repo_root: str = "",
                                            topic: str = "",
                                            max_items: int = 30) -> dict:
    """Audit examples -> docs/ipynb promotion state for the Radia repo.

    The audit enforces two project documentation rules:

    * notebooks promoted from examples should be result-saving artifacts
      (executed code cells with saved outputs and no saved error outputs);
    * Python tightly coupled to a notebook should live beside it under
      ``docs/<topic>/`` as a helper, while reusable behavior should be
      promoted into a ``src/`` API instead of remaining as a loose example.

    Args:
        repo_root: Radia repository root.  Empty means auto-detect from cwd.
        topic: Optional single ``examples/<topic>`` directory name.
        max_items: Maximum number of detailed problem rows to include.

    Returns:
        A dict with per-topic status, result-saving notebook summaries, and
        Python placement recommendations.
    """
    root = _resolve_repo_root(repo_root)
    examples = root / "examples"
    docs = root / "docs"
    if not examples.is_dir() or not docs.is_dir():
        return {
            "error": "repo root must contain examples/ and docs/",
            "repo_root": str(root),
        }

    if topic:
        topic_dirs = [examples / topic]
    else:
        topic_dirs = sorted(p for p in examples.iterdir() if p.is_dir())

    rows: list[dict] = []
    needs_notebook: list[dict] = []
    notebooks_without_results: list[dict] = []
    python_reviews: list[dict] = []

    for topic_dir in topic_dirs:
        if not topic_dir.is_dir():
            rows.append({
                "topic": topic_dir.name,
                "status": "missing_example_topic",
                "example_dir": _rel(topic_dir, root),
            })
            continue

        name = topic_dir.name
        py_files = sorted(
            p for p in topic_dir.rglob("*.py")
            if "__pycache__" not in {part.lower() for part in p.parts}
        )
        example_ipynbs = sorted(topic_dir.rglob("*.ipynb"))
        notebooks = _find_topic_notebooks(root, name)
        notebook_summaries = [
            _read_notebook_result_summary(p, root) for p in notebooks
        ]
        result_saved_count = sum(1 for s in notebook_summaries if s.get("result_saved"))
        docs_helper_dir = docs / name
        docs_py = sorted(docs_helper_dir.rglob("*.py")) if docs_helper_dir.exists() else []

        if not notebooks:
            status = "needs_result_saving_notebook" if py_files else "readme_or_data_only"
        elif result_saved_count == len(notebooks):
            status = "notebook_result_saved"
        elif result_saved_count:
            status = "partial_notebook_results"
        else:
            status = "notebook_without_saved_results"

        if py_files and notebooks:
            py_policy = "review_each_example_py_for_docs_helper_or_src_api"
        elif py_files:
            py_policy = "promote_to_result_saving_notebook_before_prune"
        elif docs_py:
            py_policy = "docs_local_helpers_present"
        else:
            py_policy = "no_topic_python"

        row = {
            "topic": name,
            "status": status,
            "example_dir": _rel(topic_dir, root),
            "has_readme": (topic_dir / "README.md").exists(),
            "example_py_count": len(py_files),
            "example_ipynb_count": len(example_ipynbs),
            "docs_notebook_count": len(notebooks),
            "result_saved_notebook_count": result_saved_count,
            "docs_helper_py_count": len(docs_py),
            "python_policy": py_policy,
            "sample_example_py": [_rel(p, root) for p in py_files[:5]],
            "sample_docs_helper_py": [_rel(p, root) for p in docs_py[:5]],
            "notebooks": notebook_summaries,
        }
        rows.append(row)

        if status == "needs_result_saving_notebook" and len(needs_notebook) < max_items:
            needs_notebook.append({
                "topic": name,
                "example_py_count": len(py_files),
                "sample_example_py": row["sample_example_py"],
            })
        if status in {"notebook_without_saved_results", "partial_notebook_results"}:
            if len(notebooks_without_results) < max_items:
                notebooks_without_results.append({
                    "topic": name,
                    "status": status,
                    "notebooks": notebook_summaries,
                })
        if py_policy in {
            "review_each_example_py_for_docs_helper_or_src_api",
            "promote_to_result_saving_notebook_before_prune",
        } and len(python_reviews) < max_items:
            python_reviews.append({
                "topic": name,
                "python_policy": py_policy,
                "example_py_count": len(py_files),
                "sample_example_py": row["sample_example_py"],
            })

    summary = {
        "repo_root": str(root),
        "topics_scanned": len(rows),
        "topics_with_docs_notebooks": sum(1 for r in rows if r.get("docs_notebook_count", 0) > 0),
        "topics_with_all_result_saved_notebooks": sum(
            1 for r in rows
            if r.get("docs_notebook_count", 0) > 0
            and r.get("docs_notebook_count") == r.get("result_saved_notebook_count")
        ),
        "topics_needing_result_saving_notebook": sum(
            1 for r in rows if r.get("status") == "needs_result_saving_notebook"
        ),
        "topics_with_notebook_result_gaps": sum(
            1 for r in rows
            if r.get("status") in {"notebook_without_saved_results", "partial_notebook_results"}
        ),
        "topics_needing_python_placement_review": sum(
            1 for r in rows
            if r.get("python_policy") in {
                "review_each_example_py_for_docs_helper_or_src_api",
                "promote_to_result_saving_notebook_before_prune",
            }
        ),
    }
    return {
        "policy": {
            "scope": "docs/<topic> result-bearing method/showcase notebooks",
            "notebook": "docs notebooks must be result-saving: executed cells with saved outputs, no saved error outputs",
            "json": "computed results should have adjacent JSON with generated_at_utc and version/runtime metadata",
            "sync": "the JSON sidecar and result-bearing ipynb are committed as a synchronized pair",
            "python": "notebook-coupled .py goes under docs/<topic>/; reusable behavior is promoted to src/ API",
        },
        "summary": summary,
        "needs_result_saving_notebook": needs_notebook,
        "notebooks_without_saved_results": notebooks_without_results,
        "python_placement_reviews": python_reviews,
        "topics": rows if topic else rows[:max_items],
        "topics_truncated": (not topic and len(rows) > max_items),
    }


def _literal_reference_hits(root: pathlib.Path,
                            needles: list[str],
                            max_hits: int = 50) -> list[dict]:
    hits: list[dict] = []
    scan_roots = [root / "pyproject.toml", root / "src", root / "packages",
                  root / "docs", root / "tests", root / "validation_test"]
    for base in scan_roots:
        paths = [base] if base.is_file() else list(_iter_files(base, _AUDIT_TEXT_SUFFIXES))
        for p in paths:
            try:
                text = p.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                try:
                    text = p.read_text(encoding="cp932")
                except Exception:
                    continue
            except Exception:
                continue
            for lineno, line in enumerate(text.splitlines(), start=1):
                if any(n in line for n in needles):
                    hits.append({
                        "path": _rel(p, root),
                        "line": lineno,
                        "text": line.strip()[:240],
                    })
                    if len(hits) >= max_hits:
                        return hits
    return hits


def document_meta_panel_layout_audit(repo_root: str = "",
                                       max_hits: int = 50) -> dict:
    """Audit impact of moving panel surfaces toward repo-root ``panels/``.

    This does not move files.  It reports the current panel tree, whether a
    root-level ``panels/`` tree already exists, and literal references that
    must be migrated or shimmed.

    Args:
        repo_root: Radia repository root.  Empty means auto-detect from cwd.
        max_hits: Maximum old-path reference hits to include.

    Returns:
        A migration-readiness report for the root-level panels policy.
    """
    root = _resolve_repo_root(repo_root)
    old_dir = root / "src" / "radia" / "panels"
    new_dir = root / "panels"
    old_refs = _literal_reference_hits(
        root,
        ["src/radia/panels", "src\\radia\\panels", "radia/panels"],
        max_hits=max_hits,
    )
    new_refs = _literal_reference_hits(
        root,
        ["panels/", "panels\\"],
        max_hits=max_hits,
    ) if new_dir.exists() else []

    old_py = _count_tree_files(old_dir, "*.py")
    old_calc = len(list(old_dir.glob("calc_*.py"))) if old_dir.exists() else 0
    old_samples = _count_tree_files(old_dir / "samples", "*")
    old_notebooks = _count_tree_files(old_dir / "notebooks", "*.ipynb")

    if not old_dir.exists():
        readiness = "no_src_panel_tree_found"
    elif old_refs:
        readiness = "needs_staged_migration_or_compat_shim"
    else:
        readiness = "ready_for_layout_move_review"

    return {
        "policy": {
            "desired_layout": "panel operating surface lives at repo-root panels/",
            "compatibility": "keep imports/package-data/tests green during migration; reusable calc logic belongs in src/ APIs",
            "notebook": "panel notebooks should also be result-saving artifacts when they demonstrate a method/result",
        },
        "repo_root": str(root),
        "current_panel_dir": _rel(old_dir, root),
        "target_panel_dir": _rel(new_dir, root),
        "target_panel_dir_exists": new_dir.exists(),
        "current_counts": {
            "py_files": old_py,
            "calc_scripts": old_calc,
            "sample_files": old_samples,
            "notebooks": old_notebooks,
        },
        "old_path_reference_count_reported": len(old_refs),
        "old_path_reference_hits": old_refs,
        "target_path_reference_count_reported": len(new_refs),
        "target_path_reference_hits": new_refs,
        "readiness": readiness,
        "recommended_sequence": [
            "classify each calc_*.py as panel-only wrapper or reusable src API candidate",
            "move user-facing panel assets/samples/notebooks toward panels/ with compatibility imports as needed",
            "update pyproject package-data, validation paths, and MCP/panel knowledge in the same change",
            "run panel CLI diff/golden checks before deleting any old-path shim",
        ],
    }
