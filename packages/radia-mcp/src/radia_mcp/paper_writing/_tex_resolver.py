"""Multi-file .tex resolver + abstract auto-extraction.

Added 2026-05-26 (radia-mcp v0.93.0).  Two helpers that close
known gaps in em_submission_gate / individual lint tools:

  1. ``resolve_input_chain(main_tex)`` -- recursively inline
     ``\\input{...}`` and ``\\include{...}`` directives so that
     downstream tools can scan the full document body.
  2. ``extract_abstract_from_tex(tex_source)`` -- pull the
     ``\\begin{abstract}...\\end{abstract}`` body for the
     abstract-only checks (validate_abstract_length,
     check_abstract_background_ratio, count_weak_expressions on
     the abstract block specifically).

Both are exposed as MCP tools so the orchestrator no longer needs
the user to manually pass abstract_text and so multi-chapter
theses / long IEEE papers using \\input get scanned in full.
"""

from __future__ import annotations

import os
import pathlib
import re
from typing import Optional


# Maximum recursion depth for \\input / \\include chains.  Caps a
# pathological cycle from blowing the stack.  Realistic papers
# never exceed depth 3 (book -> chapter -> section).
_MAX_INPUT_DEPTH = 8


_INPUT_PATTERN = re.compile(
    r"\\(?:input|include)\s*\{\s*([^{}\s][^{}]*?)\s*\}"
)


def resolve_input_chain(
    main_tex_path: str,
    max_depth: int = _MAX_INPUT_DEPTH,
    encoding: str = "utf-8",
) -> dict:
    """Recursively inline \\input{} and \\include{} to a single string.

    Walks the document tree starting from main_tex_path.  Each
    \\input{X} is resolved to <dir>/X.tex (or <dir>/X if it already
    ends in .tex).  Duplicate inputs (the same file referenced twice)
    are inlined ONCE then noted as duplicate in the report (LaTeX
    would re-typeset them but for static analysis a single pass is
    sufficient and avoids infinite-loop risk).

    Args:
        main_tex_path: path to the main .tex (the one with
            \\documentclass{...} and \\begin{document}).
        max_depth: recursion cap (default 8).  Raises if exceeded.
        encoding: file encoding (utf-8 default; falls back to
            cp932 for legacy Japanese files via 'replace' errors).

    Returns:
        dict with:
          "ok": bool
          "main": main_tex_path echoed
          "merged_tex": single-string body with every input inlined
          "files_resolved": list of {"path", "depth", "size_bytes"}
          "files_missing": list of {"path", "referenced_from"}
          "files_duplicate": list of {"path", "first_at",
                                       "duplicate_at"}
          "total_chars": int (size of merged output)

    Notes:
        * Comments after \\input{...} % comment are stripped.
        * \\input inside a comment line (leading %) is skipped.
        * Relative paths are resolved against the compilation working
          directory, represented here by the main file's directory. TeX does
          not change its working directory when it enters an input file.
    """
    main = os.path.abspath(main_tex_path)
    main_dir = os.path.dirname(main)
    if not os.path.exists(main):
        return {
            "ok": False,
            "error": f"main_tex_path not found: {main_tex_path}",
        }

    resolved: list[dict] = []
    missing: list[dict] = []
    duplicate: list[dict] = []
    seen_paths: dict[str, int] = {}   # abs_path -> first inlined depth

    def _read(path: str) -> str:
        raw = pathlib.Path(path).read_bytes()
        try:
            return raw.decode(encoding, errors="strict")
        except UnicodeDecodeError:
            if encoding.casefold().replace("-", "") != "utf8":
                raise
            return raw.decode("cp932", errors="strict")

    def _resolve(path: str, depth: int, parent: Optional[str] = None) -> str:
        abs_path = os.path.abspath(path)
        if abs_path in seen_paths:
            duplicate.append({
                "path": abs_path,
                "first_at": seen_paths[abs_path],
                "duplicate_at": depth,
                "referenced_from": parent or "(root)",
            })
            return f"\n% [resolve_input_chain: duplicate input {abs_path} skipped]\n"
        if depth > max_depth:
            raise RecursionError(
                f"resolve_input_chain: max_depth={max_depth} exceeded "
                f"at {abs_path} (cyclic \\input?)"
            )
        if not os.path.exists(abs_path):
            missing.append({
                "path": abs_path,
                "referenced_from": parent or "(root)",
            })
            return f"\n% [resolve_input_chain: missing file {abs_path}]\n"

        seen_paths[abs_path] = depth
        text = _read(abs_path)
        resolved.append({
            "path": abs_path,
            "depth": depth,
            "size_bytes": len(text.encode("utf-8", errors="replace")),
        })

        # Replace each \\input{X} with the recursive resolution.
        def _sub(m):
            # Skip if the \\input is inside a comment line: find
            # the start of the line and check for '%' before.
            line_start = text.rfind("\n", 0, m.start()) + 1
            line_before = text[line_start:m.start()]
            # If there's an unescaped %, this is a comment
            if re.search(r"(?<!\\)%", line_before):
                return m.group(0)
            inc = m.group(1).strip()
            if not os.path.splitext(inc)[1]:
                inc = inc + ".tex"
            child_path = inc if os.path.isabs(inc) else os.path.join(main_dir, inc)
            return _resolve(child_path, depth + 1, parent=abs_path)

        return _INPUT_PATTERN.sub(_sub, text)

    try:
        merged = _resolve(main, depth=0)
    except (RecursionError, OSError, UnicodeError) as e:
        return {
            "ok": False,
            "error": str(e),
            "files_resolved": resolved,
            "files_missing": missing,
        }

    result = {
        "ok": not missing,
        "main": main_tex_path,
        "merged_tex": merged,
        "files_resolved": resolved,
        "files_missing": missing,
        "files_duplicate": duplicate,
        "total_chars": len(merged),
    }
    if missing:
        result["error"] = f"{len(missing)} TeX input file(s) could not be resolved"
    return result


def paper_writing_resolve_input_chain(
    main_tex_path: str,
    max_depth: int = _MAX_INPUT_DEPTH,
    return_merged_text: bool = False,
    encoding: str = "utf-8",
) -> dict:
    """MCP tool wrapper for resolve_input_chain().

    Args:
        main_tex_path: main .tex file.
        max_depth: recursion cap (default 8).
        return_merged_text: True to include the full merged_tex in
            the response (default False -- can be 100s of KB; better
            to call this once, write the merged file, then pass that
            path to the per-tool checks).
        encoding: 'utf-8' default; legacy theses may need 'cp932'.

    Returns:
        dict with files_resolved / missing / duplicate metadata.
        merged_tex included only if return_merged_text=True.
    """
    r = resolve_input_chain(
        main_tex_path, max_depth=max_depth, encoding=encoding
    )
    if return_merged_text:
        return r
    # Strip merged_tex from response (just metadata)
    return {k: v for k, v in r.items() if k != "merged_tex"}


# ============================================================
# Abstract auto-extraction
# ============================================================


_ABSTRACT_PATTERNS = [
    # \begin{abstract} ... \end{abstract}
    re.compile(
        r"\\begin\{abstract\}(.*?)\\end\{abstract\}",
        re.DOTALL,
    ),
    # IEEEtran: \begin{IEEEabstract} ... \end{IEEEabstract}
    re.compile(
        r"\\begin\{IEEEabstract\}(.*?)\\end\{IEEEabstract\}",
        re.DOTALL,
    ),
    # IEEEkeywords precedes the body; abstract is in \abstract{...}
    re.compile(
        r"\\abstract\s*\{(.*?)\}\s*(?=\\(?:IEEEkeywords|maketitle|section))",
        re.DOTALL,
    ),
    # AAAI-style \begin{quote}\textbf{Abstract}\\ ... \end{quote}
    re.compile(
        r"\\textbf\{Abstract\}\s*\\\\\s*(.*?)\\end\{quote\}",
        re.DOTALL,
    ),
]


def extract_abstract_from_tex(tex_source: str) -> Optional[str]:
    """Extract the abstract body from a LaTeX source.

    Tries 4 patterns in order:
      1. \\begin{abstract}...\\end{abstract}  (standard)
      2. \\begin{IEEEabstract}...\\end{IEEEabstract}
      3. \\abstract{...}                       (IEEEtran legacy)
      4. \\textbf{Abstract}\\\\ ... \\end{quote}  (AAAI-style)

    Returns None if no abstract found.  Comments and leading/trailing
    whitespace are stripped from the returned text.
    """
    for pat in _ABSTRACT_PATTERNS:
        m = pat.search(tex_source)
        if m:
            body = m.group(1)
            # Strip line comments
            body = re.sub(r"(?<!\\)%[^\n]*", "", body)
            return body.strip()
    return None


def paper_writing_extract_abstract(
    tex_path: str,
    encoding: str = "utf-8",
    auto_resolve_inputs: bool = True,
) -> dict:
    """Extract the abstract body from a .tex file.

    Args:
        tex_path: main .tex file.
        encoding: 'utf-8' default; 'cp932' for legacy Japanese.
        auto_resolve_inputs: True (default) recursively inlines
            \\input{} so abstracts split into a separate file
            (e.g. ``\\input{abstract}``) are still found.

    Returns:
        dict with:
          "ok": bool
          "tex_path": echoed
          "abstract": str OR None
          "n_chars": int (length of abstract; 0 if None)
          "method": "direct" | "resolved_inputs" | "not_found"
          "advice": str next-step hint
    """
    if not os.path.exists(tex_path):
        return {
            "ok": False,
            "error": f"tex_path not found: {tex_path}",
            "abstract": None,
        }

    try:
        raw = pathlib.Path(tex_path).read_bytes()
        try:
            src = raw.decode(encoding, errors="strict")
        except UnicodeDecodeError:
            if encoding.casefold().replace("-", "") != "utf8":
                raise
            src = raw.decode("cp932", errors="strict")
    except (OSError, UnicodeError) as exc:
        return {
            "ok": False,
            "error": f"failed to read tex: {exc}",
            "abstract": None,
        }

    # Try direct first
    abs_text = extract_abstract_from_tex(src)
    if abs_text:
        return {
            "ok": True,
            "tex_path": tex_path,
            "abstract": abs_text,
            "n_chars": len(abs_text),
            "method": "direct",
            "advice": (f"Abstract extracted ({len(abs_text)} chars).  Pass "
                       f"directly to paper_writing_em_submission_gate via "
                       f"abstract_text=... parameter."),
        }

    # Fall back to resolve_input_chain
    if auto_resolve_inputs:
        r = resolve_input_chain(tex_path, encoding=encoding)
        if r.get("ok"):
            merged = r["merged_tex"]
            abs_text = extract_abstract_from_tex(merged)
            if abs_text:
                return {
                    "ok": True,
                    "tex_path": tex_path,
                    "abstract": abs_text,
                    "n_chars": len(abs_text),
                    "method": "resolved_inputs",
                    "files_resolved": len(r["files_resolved"]),
                    "advice": (f"Abstract found inside an \\input subfile "
                               f"after resolving {len(r['files_resolved'])} "
                               f"included files."),
                }

    return {
        "ok": False,
        "tex_path": tex_path,
        "abstract": None,
        "n_chars": 0,
        "method": "not_found",
        "advice": (
            "No abstract found.  Tried 4 patterns: \\begin{abstract}, "
            "\\begin{IEEEabstract}, \\abstract{...}, \\textbf{Abstract}.  "
            "If the abstract is in a non-standard form, extract manually "
            "and pass via abstract_text= directly."
        ),
    }
