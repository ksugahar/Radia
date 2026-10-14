"""Digest-specific lint tools for conference papers (IGTE / COMPUMAG /
CEFC / IEEE Transactions).

Added 2026-05-27.  Three companion checks that surface the most-frequent
reviewer comments on short-format papers (2-page digest, 4-6-page
extended abstract, full transactions):

  1. paper_writing_check_ref_label_consistency(tex_path)
     - \\label{} ↔ \\ref{} cross-reference integrity
     - dangling refs become "[??]" in the PDF (compile-time silent)
     - orphan labels (declared but never cited) flag unused figures

  2. paper_writing_check_ieee_keywords(tex_path)
     - \\begin{IEEEkeywords}...\\end{IEEEkeywords} exists + 3-7 entries
     - alternate form: \\keywords{} (Elsevier / Springer)
     - per-keyword length sanity check

  3. paper_writing_check_pdf_unresolved_markers(pdf_path)
     - post-compile scan of the rendered PDF for "[?]" / "[??]" markers
       (and the Hebrew "??" variant some templates produce).  Catches
       missing references that survived the LaTeX compile but render as
       broken in the final PDF.  Companion to the static
       paper_writing_check_citation_keys_exist + check_ref_label_consistency.

Why a single file:
  All three are digest-review lints called together by the digest
  submission gate; keeping them adjacent makes it obvious that they
  cover three layers of the same problem (label declaration → cite
  command → final PDF render).
"""
from __future__ import annotations

import os
import re


# ============================================================
# Tool 1: \ref{} ↔ \label{} consistency
# ============================================================

# All ref-family commands.  Captures the key inside {...}.
_REF_RE = re.compile(
    r"\\(?:ref|eqref|autoref|cref|Cref|pageref|nameref|vref|"
    r"figref|tabref|secref)"
    r"(?:\[[^\]]*\])?"
    r"\{([^\}]+)\}"
)
# \label{key}
_LABEL_RE = re.compile(r"\\label\{([^\}]+)\}")

# Patterns that should be stripped before scanning to avoid false matches
# inside comments / verbatim blocks.
_STRIP_BEFORE_SCAN = [
    re.compile(r"(?<!\\)%[^\n]*"),                                  # comments
    re.compile(r"\\begin\{verbatim\}.+?\\end\{verbatim\}", re.DOTALL),
    re.compile(r"\\begin\{lstlisting\}.+?\\end\{lstlisting\}", re.DOTALL),
    re.compile(r"\\verb\|[^|]*\|"),
]


def _strip_for_lint(src: str) -> str:
    out = src
    for pat in _STRIP_BEFORE_SCAN:
        out = pat.sub(lambda m: " " * (m.end() - m.start()), out)
    return out


def paper_writing_check_ref_label_consistency(
    tex_path: str,
    auto_resolve_inputs: bool = True,
    encoding: str = "utf-8",
    max_reported: int = 50,
) -> dict:
    """Verify that every ``\\ref{key}`` resolves to a ``\\label{key}`` and
    that every label is actually referenced.

    Args:
        tex_path: main .tex file.
        auto_resolve_inputs: True (default) inlines ``\\input{}`` first.
        encoding: file encoding (``utf-8`` / ``cp932``).
        max_reported: cap per-list report length.

    Returns:
        dict with::

            {
              "tex_path": echoed,
              "n_labels": int,
              "n_refs": int,
              "n_dangling_refs": int,    # \\ref{} without matching \\label{}
              "n_orphan_labels": int,    # \\label{} never \\ref'd
              "dangling_refs": [{"key", "n_occurrences", "first_offset",
                                 "first_context"}, ...],
              "orphan_labels": [str, ...],
              "all_dangling_refs": [...],
              "all_orphan_labels": [...],
              "status": "clean" | "warning" | "fail",
              "advice": str,
            }

    Status:
      * ``fail``: any dangling ref (renders as ``[??]`` in the PDF).
      * ``warning``: all refs resolve, but some labels are never cited
        — orphan figures/equations are space-wasters in a digest.
      * ``clean``: 1-to-1.
    """
    if not os.path.exists(tex_path):
        return {"error": f"tex_path not found: {tex_path}"}
    with open(tex_path, encoding=encoding, errors="replace") as fh:
        src = fh.read()

    if auto_resolve_inputs:
        try:
            from ._tex_resolver import resolve_input_chain
            r = resolve_input_chain(tex_path, encoding=encoding)
            if not r.get("ok"):
                return {"error": f"failed to resolve TeX inputs: {r.get('error')}"}
            src = r["merged_tex"]
        except Exception as exc:  # noqa: BLE001
            return {"error": f"failed to resolve TeX inputs: {exc}"}

    stripped = _strip_for_lint(src)

    labels: set[str] = set()
    for m in _LABEL_RE.finditer(stripped):
        labels.add(m.group(1).strip())

    ref_occurrences: dict[str, list[int]] = {}
    for m in _REF_RE.finditer(stripped):
        for key in (part.strip() for part in m.group(1).split(",")):
            if key:
                ref_occurrences.setdefault(key, []).append(m.start())

    ref_keys = set(ref_occurrences.keys())
    dangling = sorted(ref_keys - labels)
    orphan = sorted(labels - ref_keys)

    dangling_records = []
    for k in dangling[:max_reported]:
        first = ref_occurrences[k][0]
        a = max(0, first - 50)
        b = min(len(src), first + 60)
        dangling_records.append({
            "key": k,
            "n_occurrences": len(ref_occurrences[k]),
            "first_offset": first,
            "first_context": src[a:b].replace("\n", " "),
        })

    if dangling:
        status = "fail"
        advice = (
            f"{len(dangling)} ref(s) with no matching \\label{{}} — these "
            "render as '[??]' in the PDF.  Add the missing label OR fix "
            "the typo in the ref key.  Most common cause: copy-pasted "
            "\\ref{fig:foo} but the label is \\label{fig:bar}."
        )
    elif orphan:
        status = "warning"
        advice = (
            f"All {len(ref_keys)} ref(s) resolve, but {len(orphan)} "
            "label(s) are never cited.  In a digest, every figure / "
            "equation / table MUST be referenced in prose (Fig.~\\ref{} / "
            "Eq.~(\\ref{}) / Table~\\ref{}).  Orphan labels often mean: "
            "(a) the figure was added but not yet discussed, OR (b) the "
            "discussion was deleted but the figure remains — both waste "
            "digest page space."
        )
    else:
        status = "clean"
        advice = (
            f"All {len(ref_keys)} ref(s) resolve and all {len(labels)} "
            "label(s) are cited — cross-references are in sync."
        )

    return {
        "tex_path": tex_path,
        "n_labels": len(labels),
        "n_refs": len(ref_keys),
        "n_dangling_refs": len(dangling),
        "n_orphan_labels": len(orphan),
        "dangling_refs": dangling_records,
        "orphan_labels": orphan[:max_reported],
        "all_dangling_refs": dangling,
        "all_orphan_labels": orphan,
        "truncated_dangling": len(dangling) > max_reported,
        "truncated_orphan": len(orphan) > max_reported,
        "status": status,
        "advice": advice,
    }


# ============================================================
# Tool 2: IEEE Index Terms / Keywords check
# ============================================================

# IEEEtran: \begin{IEEEkeywords} ... \end{IEEEkeywords}
_IEEE_KW_BLOCK_RE = re.compile(
    r"\\begin\{IEEEkeywords\}(.+?)\\end\{IEEEkeywords\}",
    re.DOTALL,
)
# Elsevier / Springer: \keywords{ ... }  (also \PACS{} / Mathematics Subject Classification etc. — keep loose)
_ELSEVIER_KW_RE = re.compile(
    r"\\keywords\s*\{(.+?)\}",
    re.DOTALL,
)


def paper_writing_check_ieee_keywords(
    tex_path: str,
    auto_resolve_inputs: bool = True,
    encoding: str = "utf-8",
    min_keywords: int = 3,
    max_keywords: int = 7,
    min_kw_chars: int = 3,
    max_kw_chars: int = 50,
) -> dict:
    """Verify the paper has an Index Terms / Keywords block of acceptable
    size and per-entry length.

    Recognised templates:
      * IEEEtran: ``\\begin{IEEEkeywords} kw1, kw2, ... \\end{IEEEkeywords}``
      * Elsevier / Springer: ``\\keywords{kw1; kw2; ...}``

    Args:
        tex_path: main .tex file.
        auto_resolve_inputs: True (default) inlines ``\\input{}`` first.
        encoding: file encoding.
        min_keywords: minimum acceptable count (IEEE rule of thumb 3).
        max_keywords: maximum acceptable count (IEEE rule of thumb 7;
            COMPUMAG digest typically 5).
        min_kw_chars: per-entry minimum length (catches "X" stubs).
        max_kw_chars: per-entry maximum length (catches sentences-as-keywords).

    Returns:
        dict with::

            {
              "tex_path": echoed,
              "found_block": "IEEEkeywords" | "keywords" | None,
              "keywords": [str, ...],
              "n_keywords": int,
              "warnings": [str, ...],  # per-rule warnings
              "status": "missing" | "warning" | "clean",
              "advice": str,
            }
    """
    if not os.path.exists(tex_path):
        return {"error": f"tex_path not found: {tex_path}"}
    with open(tex_path, encoding=encoding, errors="replace") as fh:
        src = fh.read()

    if auto_resolve_inputs:
        try:
            from ._tex_resolver import resolve_input_chain
            r = resolve_input_chain(tex_path, encoding=encoding)
            if not r.get("ok"):
                return {"error": f"failed to resolve TeX inputs: {r.get('error')}"}
            src = r["merged_tex"]
        except Exception as exc:  # noqa: BLE001
            return {"error": f"failed to resolve TeX inputs: {exc}"}

    src = _strip_for_lint(src)

    block: str | None = None
    found_block: str | None = None
    m = _IEEE_KW_BLOCK_RE.search(src)
    if m:
        block = m.group(1)
        found_block = "IEEEkeywords"
    else:
        m = _ELSEVIER_KW_RE.search(src)
        if m:
            block = m.group(1)
            found_block = "keywords"

    if not block:
        return {
            "tex_path": tex_path,
            "found_block": None,
            "keywords": [],
            "n_keywords": 0,
            "warnings": [],
            "status": "missing",
            "advice": (
                "No Index Terms / Keywords block found.  Required by IEEE "
                "Transactions, IGTE / COMPUMAG / CEFC digest.  Add e.g.:\n"
                "  \\begin{IEEEkeywords}\n"
                "  finite element analysis, induction heating, "
                "magneto-quasi-static, surface impedance\n"
                "  \\end{IEEEkeywords}\n"
                "(3-7 entries, lowercase comma-separated)."
            ),
        }

    # Split on comma or semicolon; trim whitespace; drop empty
    raw_entries = re.split(r"[;,]", block)
    keywords = [e.strip().rstrip(".") for e in raw_entries]
    keywords = [k for k in keywords if k]
    n = len(keywords)

    warnings: list[str] = []
    if n < min_keywords:
        warnings.append(
            f"Only {n} keyword(s) — IEEE recommends >= {min_keywords}."
        )
    if n > max_keywords:
        warnings.append(
            f"{n} keywords — IEEE recommends <= {max_keywords} "
            "(too many dilutes indexing weight)."
        )
    for k in keywords:
        if len(k) < min_kw_chars:
            warnings.append(f"Keyword too short: {k!r} (< {min_kw_chars} chars)")
        elif len(k) > max_kw_chars:
            warnings.append(
                f"Keyword too long: {k[:60]!r}... (> {max_kw_chars} chars) "
                "— sentences-as-keywords get rejected by IEEE Xplore indexer."
            )

    if warnings:
        status = "warning"
        advice = (
            f"Keyword block present ({found_block}, {n} entries) but "
            f"{len(warnings)} issue(s) found.  See 'warnings' field for "
            "specifics.  Aim for 3-7 lowercase noun phrases of 3-50 chars each."
        )
    else:
        status = "clean"
        advice = (
            f"Keyword block OK: {found_block} with {n} entries, each "
            "within length bounds."
        )

    return {
        "tex_path": tex_path,
        "found_block": found_block,
        "keywords": keywords,
        "n_keywords": n,
        "warnings": warnings,
        "status": status,
        "advice": advice,
    }


# ============================================================
# Tool 3: PDF unresolved-marker scan
# ============================================================


# Markers that indicate unresolved \cite or \ref in the rendered PDF.
# LaTeX default: "[?]" (single ?) for unknown cite, "[??]" for unknown ref.
# Some templates render "(?)" or just "?" — but matching plain "?" causes
# huge false-positive rate, so we keep the patterns conservative.
_PDF_MARKER_PATTERNS = [
    (re.compile(r"\[\?\?\]"), "[??]", "unresolved \\ref or \\eqref"),
    (re.compile(r"\[\?\]"),   "[?]",  "unresolved \\cite"),
]


def paper_writing_check_pdf_unresolved_markers(
    pdf_path: str,
    max_reported: int = 50,
    context_chars: int = 60,
) -> dict:
    """Scan a rendered PDF for [?] / [??] markers (unresolved cite / ref).

    Companion to the static
    ``paper_writing_check_citation_keys_exist`` and
    ``paper_writing_check_ref_label_consistency`` tools — those check
    the .tex BEFORE compile.  This one checks the .pdf AFTER compile, so
    it catches cases where the .tex looks fine but bibtex/biber wasn't
    re-run, or where ``\\nocite`` failed.

    Implementation: extract per-page text via pymupdf and grep for the
    marker patterns.  Pure text scan; ignores text layout, so the marker
    must literally be in the rendered text.

    Args:
        pdf_path: rendered PDF (after pdflatex + bibtex + 2nd pdflatex).
        max_reported: cap the number of marker locations reported.
        context_chars: characters around each marker for the context
            snippet (default 60).

    Returns:
        dict with::

            {
              "pdf_path": echoed,
              "n_pages": int,
              "n_markers_total": int,
              "by_pattern": {"[?]": int, "[??]": int},
              "markers": [
                  {"pattern": "[??]",
                   "meaning": "unresolved \\ref or \\eqref",
                   "page": int (1-based),
                   "context": str},
                  ...
              ],
              "status": "clean" | "fail",
              "advice": str,
            }
    """
    if not os.path.exists(pdf_path):
        return {"error": f"pdf_path not found: {pdf_path}"}

    try:
        import pymupdf  # type: ignore
    except ImportError:
        return {
            "error": "pymupdf not installed.  pip install pymupdf",
            "pdf_path": pdf_path,
        }

    try:
        doc = pymupdf.open(pdf_path)
    except Exception as e:  # noqa: BLE001
        return {"error": f"failed to open PDF: {e}", "pdf_path": pdf_path}

    by_pattern: dict[str, int] = {marker: 0 for _, marker, _ in _PDF_MARKER_PATTERNS}
    markers: list[dict] = []

    try:
        n_pages = doc.page_count
        for i in range(n_pages):
            text = doc[i].get_text("text") or ""
            for pat, label, meaning in _PDF_MARKER_PATTERNS:
                for m in pat.finditer(text):
                    by_pattern[label] += 1
                    if len(markers) < max_reported:
                        a = max(0, m.start() - context_chars)
                        b = min(len(text), m.end() + context_chars)
                        ctx = text[a:b].replace("\n", " ")
                        markers.append({
                            "pattern": label,
                            "meaning": meaning,
                            "page": i + 1,
                            "context": ctx,
                        })
    finally:
        doc.close()

    n_total = sum(by_pattern.values())
    if n_total == 0:
        status = "clean"
        advice = (
            "No '[?]' / '[??]' markers found in the rendered PDF — every "
            "\\cite{} and \\ref{} resolved correctly at compile time."
        )
    else:
        status = "fail"
        parts = []
        if by_pattern["[?]"]:
            parts.append(
                f"{by_pattern['[?]']} unresolved \\cite{{}} marker(s) — "
                "likely cause: bibtex/biber not re-run after .bib edit, OR "
                "the cite key truly is missing.  Run "
                "'paper_writing_check_citation_keys_exist' on the .tex+.bib "
                "to identify the missing key, then re-compile with "
                "pdflatex -> bibtex -> pdflatex -> pdflatex."
            )
        if by_pattern["[??]"]:
            parts.append(
                f"{by_pattern['[??]']} unresolved \\ref{{}} marker(s) — "
                "likely cause: \\label{} missing or typo'd.  Run "
                "'paper_writing_check_ref_label_consistency' on the .tex to "
                "find the dangling refs."
            )
        advice = "  ".join(parts)

    return {
        "pdf_path": pdf_path,
        "n_pages": n_pages,
        "n_markers_total": n_total,
        "by_pattern": by_pattern,
        "markers": markers,
        "truncated": n_total > max_reported,
        "status": status,
        "advice": advice,
    }
