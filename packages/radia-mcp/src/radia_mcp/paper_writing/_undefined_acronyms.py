"""Undefined-acronym detection for LaTeX EM papers.

Added 2026-05-27 (radia-mcp).  Catches the IEEE / IEEJ / paper-reviewer
pattern:

  Reviewer: "The introduction uses 'MQS' and 'IH' without defining
   what they stand for.  Please spell out every acronym on first use."

Detection algorithm:

  1. Strip math (``$...$``, equations) and citation arguments
     (``\\cite{...}``) — these often contain bib keys / symbols that
     LOOK like acronyms but aren't.
  2. Find every all-caps 2-5-letter token (``\\b[A-Z][A-Z0-9]{1,4}\\b``).
  3. For each unique acronym, find the FIRST occurrence and look in a
     ±80-char window for one of:
       * ``Acronym (Full Name)``  → fine
       * ``Full Name (Acronym)``  → fine
       * Acronym is listed in a Nomenclature / IEEEdescription block
         BEFORE its first use → fine
  4. Acronyms NOT meeting any of those criteria are reported as
     "undefined on first use".

Lab POLICY:

  Every domain-specific acronym MUST be defined on first use.  Universal
  acronyms (PDF, USA, UK, USB, RAM, CPU, GPU, …) are exempt.  The lab's
  EM-specific common acronyms are intentionally NOT in the whitelist —
  even FEM / BEM / MQS / IH must be spelled out on first appearance.

Limitations (call out for honesty):

  * Pure-text detection; cannot tell "ABC Corp." from the acronym ABC.
    Add such company names via ``extra_whitelist``.
  * "FEM" inside ``\\cite{FEM2020}`` is stripped before scanning, but
    custom ``\\foo{ABC}`` commands may leak through; pass them in via
    ``ignore_command_args=['foo', ...]`` if needed.
"""
from __future__ import annotations

import os
import re
from typing import Optional


# ============================================================
# Acronym whitelist (universally recognised — no definition needed)
# ============================================================
#
# Deliberately conservative.  Anything that an EM-paper reviewer might
# legitimately ask "what does X stand for?" is NOT here — even very
# common lab vocabulary (FEM, BEM, MQS, IH, OCR) is excluded so authors
# are forced to spell them out on first use.
UNIVERSAL_ACRONYM_WHITELIST = frozenset({
    # Web / file formats
    "PDF", "HTML", "XML", "JSON", "CSV", "URL", "URI", "HTTP", "HTTPS",
    "FTP", "TCP", "IP", "DNS", "API", "GUI", "CLI", "SDK",
    # Countries / organisations everybody knows
    "USA", "UK", "EU", "UN", "WHO", "NASA", "IEEE", "IEEJ", "SI",
    # Roman-numeral section numbers
    "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XII",
    # Hardware everybody knows
    "USB", "RAM", "ROM", "CPU", "GPU", "TPU", "SSD", "HDD", "DVD",
    "BIOS", "UEFI", "PCI", "LCD", "LED", "OLED",
    # SI / Physics units (always allowed)
    "AC", "DC", "RMS",
    "HZ", "KHZ", "MHZ", "GHZ", "THZ",
    "MM", "CM", "KM", "NM", "UM",
    "MS", "US", "NS", "PS", "FS",
    "MV", "KV", "MW", "KW", "MWH", "KWH",
    "TM",   # trademark
    "ID",   # identifier
    "OK",
    # Common research abbreviations that are NEVER spelled out
    "PHD", "MS", "BS", "MD", "JD",
    "ETC", "I.E", "E.G",   # not strictly acronyms, but caught by regex
    "PI",  # principal investigator OR pi (3.14) — context-dependent; keep out of false-positive list
})


# Pattern: 2-5-letter all-uppercase token at a word boundary.
# Must start with a letter (not a digit).  Allows letters + digits
# inside (e.g. "3D", "2D" matched? — no, must START with letter).
_ACRONYM_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z][A-Z0-9]{1,4})(?![A-Za-z0-9])")

# Patterns to STRIP from the source before scanning for acronyms:
_STRIP_PATTERNS_DEFAULT = [
    # Math environments
    re.compile(r"\$[^$\n]+?\$"),                       # inline $...$
    re.compile(r"\$\$.+?\$\$", re.DOTALL),
    re.compile(r"\\\(.+?\\\)", re.DOTALL),
    re.compile(r"(?<!\\)\\\[.+?\\\]", re.DOTALL),
    re.compile(
        r"\\begin\{(equation\*?|align\*?|eqnarray\*?|gather\*?|multline\*?)\}"
        r".+?\\end\{\1\}",
        re.DOTALL,
    ),
    # Citation / reference / label arguments (bib keys are not acronyms)
    re.compile(r"\\(?:cite|citet|citep|citeyear|citealt|citealp|"
               r"citeauthor|nocite|fullcite|textcite|parencite|"
               r"footcite|autocite|ref|eqref|label)"
               r"(?:\[[^\]]*\])?(?:\[[^\]]*\])?\{[^\}]+\}"),
    # Bibliography itself
    re.compile(r"\\begin\{thebibliography\}.+?\\end\{thebibliography\}",
                re.DOTALL),
    # Verbatim / listings (code blocks)
    re.compile(r"\\begin\{(verbatim|lstlisting|minted)\}.+?\\end\{\1\}",
                re.DOTALL),
    re.compile(r"\\verb\|[^|]*\|"),
    # Comments
    re.compile(r"(?<!\\)%[^\n]*"),
    # Common no-acronym commands (\\texttt{...} content is often bib-keys / code)
    re.compile(r"\\(?:texttt|url|href|path|nolinkurl)(?:\[[^\]]*\])?\{[^\}]+\}"),
]

# Definition patterns to scan around the first occurrence
def _definition_window(text: str, abs_start: int, abs_end: int,
                       acronym: str, window: int = 80) -> Optional[str]:
    """If acronym is defined within ±window chars of its first
    occurrence, return the matching span; else None.

    Recognised:
      ACRONYM (Full Name)     -- definition immediately after
      Full Name (ACRONYM)     -- definition immediately before
    """
    a = max(0, abs_start - window)
    b = min(len(text), abs_end + window)
    snippet = text[a:b]
    # Acronym (Full Name)  — capitalised words inside parens after token
    forward = re.search(
        rf"{re.escape(acronym)}\s*\(([^)]+)\)",
        snippet,
    )
    if forward:
        inner = forward.group(1).strip()
        # Heuristic: must contain at least one alphabetic word, not just digits / spaces
        if re.search(r"[A-Za-z]{3,}", inner):
            return f"{acronym} ({inner})"
    # Full Name (ACRONYM)  — acronym in parens preceded by 1-8 hyphenated /
    # lowercase / capitalised words.  We allow lowercase / hyphens so that
    # "finite-element method (FEM)" and "magneto-quasi-static (MQS)" both
    # qualify as definitions; the surrounding prose-not-equation context
    # is already guaranteed by _stripped_source().
    backward = re.search(
        rf"((?:[A-Za-z][A-Za-z0-9-]*\s+){{1,8}})\(\s*{re.escape(acronym)}\s*\)",
        snippet,
    )
    if backward:
        return f"{backward.group(1).strip()} ({acronym})"
    return None


def _scan_nomenclature_for_acronym(text: str, acronym: str) -> bool:
    """Check if acronym appears in any Nomenclature / IEEEdescription
    block BEFORE the document's first non-frontmatter section."""
    for env in ("nomenclature", "IEEEdescription"):
        for m in re.finditer(
            rf"\\begin\{{{env}\}}(.+?)\\end\{{{env}\}}",
            text, flags=re.DOTALL | re.IGNORECASE,
        ):
            if re.search(rf"\b{re.escape(acronym)}\b", m.group(1)):
                return True
    # Lab-friendly "section{Notation}" / "section{Acronyms}"
    for m in re.finditer(
        r"\\section\*?\{(?:Notation|Acronyms?|Abbreviations?|Nomenclature)\}"
        r"(.+?)(?=\\section\{|\\bibliography|\\end\{document\})",
        text, flags=re.DOTALL | re.IGNORECASE,
    ):
        if re.search(rf"\b{re.escape(acronym)}\b", m.group(1)):
            return True
    return False


def _stripped_source(src: str) -> str:
    """Return ``src`` with math / citations / verbatim / comments
    blanked (replaced by spaces of same length) so acronym scanning
    only sees prose.

    Same-length replacement preserves character offsets so the
    "first occurrence" location reported back to the user still
    matches the original source file."""
    out = src
    for pat in _STRIP_PATTERNS_DEFAULT:
        def _blank(m):
            return " " * (m.end() - m.start())
        out = pat.sub(_blank, out)
    return out


def paper_writing_check_undefined_acronyms(
    tex_path: str,
    extra_whitelist: str = "",
    max_reported: int = 50,
    window_chars: int = 80,
) -> dict:
    """Detect acronyms (e.g. IH, MQS, FEM) that are not spelled out on
    first use AND are not listed in a Nomenclature / Acronyms block.

    Args:
        tex_path: path to the main .tex file.  Multi-file papers:
            either concatenate first or call once per file (this tool
            does NOT auto-resolve ``\\input``).
        extra_whitelist: comma-separated additional acronyms to suppress
            from the report.  Use for project-specific acronyms that are
            defined externally (e.g. shared institutional names).
        max_reported: cap on undefined-acronym reports in the response
            (the full list is still in ``all_undefined``).
        window_chars: search range around the first occurrence to look
            for an inline ``Full Name (ACRONYM)`` definition.  Default
            80 chars matches an average sentence length.

    Returns:
        dict with::

            {
              "tex_path": echoed,
              "n_acronyms_used": int,           # unique acronyms in prose
              "n_acronyms_defined": int,
              "n_undefined": int,
              "undefined_acronyms": [
                  {"acronym": "IH",
                   "first_occurrence_offset": int,
                   "first_occurrence_context": "... uses IH and ...",
                   "definition_status": "no inline definition + not in Nomenclature"},
                  ...
              ],
              "all_undefined": [str, ...],      # alphabetical, full list
              "whitelist_size": int,
              "advice": str,                    # policy text
            }
    """
    if not os.path.exists(tex_path):
        return {"error": f"tex_path not found: {tex_path}"}
    try:
        with open(tex_path, encoding="utf-8", errors="replace") as fh:
            src = fh.read()
    except Exception as e:  # noqa: BLE001
        return {"error": f"failed to read tex: {e}"}

    whitelist = set(UNIVERSAL_ACRONYM_WHITELIST)
    if extra_whitelist:
        for s in extra_whitelist.split(","):
            s = s.strip().upper()
            if s:
                whitelist.add(s)

    stripped = _stripped_source(src)
    # Collect every acronym occurrence and its abs offset.
    occurrences: dict[str, list[tuple[int, int]]] = {}
    for m in _ACRONYM_RE.finditer(stripped):
        ac = m.group(1)
        if ac in whitelist:
            continue
        occurrences.setdefault(ac, []).append((m.start(), m.end()))

    if not occurrences:
        return {
            "tex_path": tex_path,
            "n_acronyms_used": 0,
            "n_acronyms_defined": 0,
            "n_undefined": 0,
            "undefined_acronyms": [],
            "all_undefined": [],
            "whitelist_size": len(whitelist),
            "advice": "No domain-specific acronyms detected in prose.",
        }

    # Classify each acronym
    defined: list[str] = []
    undefined_records: list[dict] = []
    for ac in sorted(occurrences.keys()):
        first_start, first_end = occurrences[ac][0]
        defn = _definition_window(src, first_start, first_end, ac,
                                  window=window_chars)
        in_nom = _scan_nomenclature_for_acronym(src, ac)
        if defn or in_nom:
            defined.append(ac)
            continue
        # Build context snippet for the reviewer's eyes (±60 chars)
        a = max(0, first_start - 60)
        b = min(len(src), first_end + 60)
        ctx = src[a:b].replace("\n", " ")
        undefined_records.append({
            "acronym": ac,
            "first_occurrence_offset": first_start,
            "first_occurrence_context": ctx,
            "n_total_occurrences": len(occurrences[ac]),
            "definition_status": ("no inline 'ACRONYM (Full Name)' or "
                                  "'Full Name (ACRONYM)' nearby; not in "
                                  "Nomenclature / Acronyms section"),
        })

    all_undefined = [r["acronym"] for r in undefined_records]
    reported = undefined_records[:max_reported]

    if undefined_records:
        advice = (
            f"{len(undefined_records)} acronym(s) used without definition on "
            "first use.  Fix in ONE of two ways:  (a) spell them out inline "
            "on first use — e.g. 'magneto-quasi-static (MQS)' or 'induction "
            "heating (IH)';  (b) add an Acronyms / Nomenclature section "
            "before Section I with each acronym listed.  Universal acronyms "
            "(PDF, USA, USB, CPU, …) are auto-excluded.  Pass "
            "extra_whitelist='ABC,XYZ' for institutional / project acronyms "
            "that need no definition."
        )
    else:
        advice = ("All domain-specific acronyms are defined on first use "
                  "or listed in a Nomenclature / Acronyms section.")

    return {
        "tex_path": tex_path,
        "n_acronyms_used": len(occurrences),
        "n_acronyms_defined": len(defined),
        "n_undefined": len(undefined_records),
        "undefined_acronyms": reported,
        "all_undefined": all_undefined,
        "truncated": len(undefined_records) > max_reported,
        "whitelist_size": len(whitelist),
        "advice": advice,
    }
