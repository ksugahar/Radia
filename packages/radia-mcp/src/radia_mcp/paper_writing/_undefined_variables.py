"""Undefined-variable detection for LaTeX EM papers.

Added 2026-05-26 (radia-mcp v0.92.0).  Catches the IEEE / IEEJ
reviewer pattern:

  Reviewer: "Eq. (5) introduces the symbol $\\eta_s$ but I cannot find
   anywhere a definition of what this represents.  Please define
   every symbol that appears in an equation either via a Nomenclature
   section, a 'where ...' clause, or inline before first use."

Detection algorithm (regex + scope tracking):

  1. Extract all math symbols (single letters + multi-letter Greek)
     from $...$, $$...$$, \\(...\\), \\[...\\], equation,
     equation*, align, align*, eqnarray, gather environments.
  2. Build the DEFINED set from:
     * Nomenclature blocks (\\begin{nomenclature}...\\end{nomenclature}
       OR \\begin{IEEEdescription}...\\end{IEEEdescription}).
     * "where" clauses after equations ("where $X$ is the ...").
     * Inline definitions ("Let $X$ denote ..." or "$X$ denotes ...").
     * Plain-text "is the" patterns ("$\\sigma$ is the conductivity").
     * A built-in WHITELIST of universally-recognized symbols
       (\\pi, \\omega, \\mu_0, \\epsilon_0, \\hbar, e, j, i, etc.).
  3. Report the symbols that appear in math but never in any of
     the definition sources.

Limitations (call out for honesty):

  * The check is **syntactic, not semantic** -- if you write
    "where $X$ is" but mean a different X, that's not caught.
  * The check is **per-paper, not cross-paper** -- citation chains
    that define $X$ in a referenced paper are not detected.
  * Indices on subscripts (e.g. $H_n$ where $n$ is a generic index)
    can produce false positives; pass them via the `whitelist`
    parameter.

Web-research basis (2026-05-26):

  * **paper-linter** (misc0110/paper-linter, GitHub) -- general
    LaTeX lint, not symbol-specific.
  * AAAI Press LaTeX template -- "every symbol must be defined".
  * IEEE Editorial Style Manual -- nomenclature requirement.
"""

from __future__ import annotations

import os
import re
from typing import Optional


# ============================================================
# Symbol whitelist (universally recognized)
# ============================================================

# Symbols that are NEVER flagged as undefined regardless of context.
# These are conventionally understood across all EM / physics papers.
UNIVERSAL_WHITELIST = frozenset({
    # Greek constants / variables in universal use
    r"\pi", r"\Pi",
    r"\omega", r"\Omega",
    r"\epsilon", r"\varepsilon",
    r"\mu", r"\nu",
    r"\sigma", r"\Sigma",
    r"\rho",
    r"\tau",
    r"\lambda", r"\Lambda",
    r"\theta", r"\Theta", r"\phi", r"\Phi", r"\varphi", r"\psi", r"\Psi",
    r"\alpha", r"\beta", r"\gamma", r"\Gamma", r"\delta", r"\Delta",
    r"\eta", r"\zeta", r"\xi", r"\Xi", r"\chi", r"\kappa",
    # SI-style constants
    r"\mu_0", r"\epsilon_0", r"\varepsilon_0", r"\hbar",
    # Loop / index / counter (almost universally generic)
    "i", "j", "k", "l", "m", "n",
    # Time / frequency / position (universal in EM body text)
    "t", "f", "x", "y", "z", "r",
    # Math operators
    r"\sin", r"\cos", r"\tan", r"\sinh", r"\cosh", r"\tanh",
    r"\exp", r"\log", r"\ln", r"\arg",
    r"\det", r"\tr", r"\rank",
    r"\Re", r"\Im",
    # Calc differentials (universal)
    r"\mathrm{d}", "d",
    # Vector operators (universal)
    r"\nabla", r"\partial",
    # Boldface common vectors (universal in EM)
    r"\mathbf{E}", r"\mathbf{B}", r"\mathbf{H}", r"\mathbf{D}",
    r"\mathbf{J}", r"\mathbf{M}", r"\mathbf{A}", r"\mathbf{F}",
    r"\boldsymbol{E}", r"\boldsymbol{B}", r"\boldsymbol{H}",
    r"\boldsymbol{D}", r"\boldsymbol{J}", r"\boldsymbol{M}",
    r"\boldsymbol{A}", r"\boldsymbol{F}",
    # Plain unbolded forms (legacy but accepted)
    "E", "B", "H", "D", "J", "M", "A", "F",
})


# ============================================================
# Regex patterns
# ============================================================

# Math environments to extract
_MATH_ENV_PATTERNS = [
    # inline $...$ (single-dollar)
    re.compile(r"(?<!\\)(?<!\$)\$([^$\n]+?)(?<!\\)\$"),
    # inline \(...\)
    re.compile(r"\\\((.+?)\\\)", re.DOTALL),
    # displayed $$...$$
    re.compile(r"\$\$(.+?)\$\$", re.DOTALL),
    # displayed \[...\]
    re.compile(r"\\\[(.+?)\\\]", re.DOTALL),
    # equation / equation*
    re.compile(r"\\begin\{equation\*?\}(.+?)\\end\{equation\*?\}",
                re.DOTALL),
    # align / align*
    re.compile(r"\\begin\{align\*?\}(.+?)\\end\{align\*?\}",
                re.DOTALL),
    # eqnarray / eqnarray*
    re.compile(r"\\begin\{eqnarray\*?\}(.+?)\\end\{eqnarray\*?\}",
                re.DOTALL),
    # gather / gather*
    re.compile(r"\\begin\{gather\*?\}(.+?)\\end\{gather\*?\}",
                re.DOTALL),
]

# Definition-source patterns
_NOMENCLATURE_PATTERNS = [
    re.compile(r"\\begin\{nomenclature\}(.+?)\\end\{nomenclature\}",
                re.DOTALL | re.IGNORECASE),
    re.compile(r"\\begin\{IEEEdescription\}(.+?)\\end\{IEEEdescription\}",
                re.DOTALL),
    # Lab-friendly informal "Notation" / "Symbols" section
    re.compile(
        r"\\section\*?\{(?:Notation|Symbols|Nomenclature)\}"
        r"(.+?)(?=\\section\{|\\bibliography|\\end\{document\})",
        re.DOTALL | re.IGNORECASE,
    ),
]

# Where-clause patterns (after equation)
_WHERE_PATTERN = re.compile(
    r"\bwhere\s+(.+?)(?:\.|\\(?:section|subsection|paragraph)|"
    r"\\end\{|\\begin\{equation)",
    re.DOTALL | re.IGNORECASE,
)

# "Let $X$ be / denote ..." or "$X$ is the ..."
_INLINE_DEF_PATTERNS = [
    re.compile(r"Let\s+\$([^$]+)\$\s+(?:be|denote)", re.IGNORECASE),
    re.compile(
        r"\$([^$]+)\$\s+(?:denotes?|represents?|is\s+the|"
        r"is\s+a|is\s+defined|stands?\s+for)",
        re.IGNORECASE,
    ),
    # Reverse: "the conductivity $\sigma$"
    re.compile(r"(?:the|a|an)\s+\w+(?:\s+\w+){0,4}\s+\$([^$]+)\$"),
]


def _extract_math_blocks(tex_source: str) -> list[str]:
    """Pull out every math block from a .tex source."""
    blocks = []
    for pat in _MATH_ENV_PATTERNS:
        for m in pat.finditer(tex_source):
            blocks.append(m.group(1))
    return blocks


def _extract_symbols_from_math(math_block: str) -> set[str]:
    """Extract LaTeX symbols (commands + single letters) used in math.

    Returns a set of normalised symbol names.  Multi-letter Greek
    commands (\\sigma, \\omega) preserved verbatim; single ASCII
    letters returned as 1-char strings.
    """
    symbols = set()

    # 1. Backslash-prefixed commands (\\alpha, \\sigma, \\mathbf{E}, etc.)
    #    Capture the command name plus an optional {...} arg AND an
    #    optional _subscript (single letter / digit OR {...}).
    #    Subscripts matter because "\\eta_s" and "\\eta" are different
    #    physical quantities -- "\\eta" might be whitelisted but
    #    "\\eta_s" (surface impedance) must be defined.
    cmd_pat = re.compile(
        r"\\([A-Za-z]+)(?:\{([^}]+)\})?(_(?:\w|\{[^}]*\}))?"
    )
    for m in cmd_pat.finditer(math_block):
        cmd = m.group(1)
        arg = m.group(2)
        subscript = m.group(3) or ""   # e.g. "_s" or "_{12}"
        # Skip operators / structural commands that aren't symbols
        if cmd in {
            "frac", "sqrt", "cdot", "times", "div", "cdots",
            "ldots", "vdots", "ddots", "infty", "to", "rightarrow",
            "leftarrow", "Rightarrow", "Leftarrow",
            "leq", "geq", "neq", "ll", "gg", "approx", "sim",
            "equiv", "propto", "in", "notin",
            "sum", "int", "iint", "iiint", "oint", "prod",
            "lim", "max", "min", "sup", "inf",
            "left", "right", "big", "Big", "bigg", "Bigg",
            "quad", "qquad", "hspace", "vspace",
            "begin", "end", "label", "ref", "eqref",
            "text", "textrm", "textit", "textbf",
            "hat", "tilde", "bar", "vec", "dot", "ddot",
            "overline", "underline", "overbrace", "underbrace",
            "stackrel", "underset", "overset",
        }:
            continue
        if arg is not None:
            # \\mathbf{E} -> store as "\\mathbf{E}"
            # \\mathbf{E}_s -> store as "\\mathbf{E}_s"
            symbols.add(f"\\{cmd}{{{arg}}}{subscript}")
        else:
            symbols.add(f"\\{cmd}{subscript}")

    # 2. Single-letter Roman variables (italic by default in math mode)
    #    Catch a-z and A-Z when NOT preceded by a backslash.
    letter_pat = re.compile(
        r"(?<![\\A-Za-z\d])([A-Za-z])(?![A-Za-z])"
    )
    # Strip out parts inside backslash-commands first, otherwise
    # \\sigma -> we'd extract 's', 'i', 'g', 'm', 'a' which is wrong.
    cleaned = re.sub(r"\\[A-Za-z]+(?:\{[^}]*\})?", " ", math_block)
    for m in letter_pat.finditer(cleaned):
        symbols.add(m.group(1))

    return symbols


def _extract_defined_symbols(tex_source: str) -> set[str]:
    """Build set of symbols defined anywhere in the .tex source."""
    defined = set()

    # 1. Nomenclature / IEEEdescription blocks
    for pat in _NOMENCLATURE_PATTERNS:
        for m in pat.finditer(tex_source):
            block = m.group(1)
            # Each line typically: \\item[$X$] description
            for im in re.finditer(r"\\item\s*\[\s*\$([^$]+)\$\s*\]", block):
                # Add symbols extracted from $X$
                inner = im.group(1)
                defined |= _extract_symbols_from_math(inner)

    # 2. "where" clauses (limit per-clause scope to ~500 chars)
    for m in _WHERE_PATTERN.finditer(tex_source):
        clause = m.group(1)[:500]
        # Within the clause, every $...$ symbol is being defined
        for im in re.finditer(r"\$([^$\n]+?)\$", clause):
            defined |= _extract_symbols_from_math(im.group(1))

    # 3. Inline definitions
    for pat in _INLINE_DEF_PATTERNS:
        for m in pat.finditer(tex_source):
            defined |= _extract_symbols_from_math(m.group(1))

    return defined


# ============================================================
# Main MCP tool
# ============================================================


def paper_writing_check_undefined_variables(
    tex_path: str,
    extra_whitelist: str = "",
    report_first_use: bool = True,
    max_reported: int = 50,
) -> dict:
    """Detect math symbols that are used but never defined.

    Reviewer-pattern target: "Eq. (5) introduces $\\eta_s$ but I cannot
    find where it is defined.  Please define every symbol."

    Args:
        tex_path: main .tex file (single-file paper).  Multi-file
            papers: concatenate first or call once per file.
        extra_whitelist: comma-separated extra symbols to suppress
            from the report.  Useful for project-specific generic
            indices (e.g. "p,q" for finite-element basis order).
            Add to UNIVERSAL_WHITELIST.
        report_first_use: True (default) finds the FIRST occurrence
            of each undefined symbol and reports the surrounding
            sentence for context.
        max_reported: cap on the number of undefined symbols
            reported (sorted alphabetically; rest hidden).

    Returns:
        dict with:
          "tex_path": echoed
          "n_math_blocks": int (count of equations / inline math)
          "n_symbols_used": int (unique symbols in math)
          "n_symbols_defined": int (unique symbols defined)
          "n_undefined": int
          "undefined_symbols": list of {
              "symbol": str (LaTeX form, e.g. "\\eta_s" or "X"),
              "first_occurrence_context": str (~100 chars around
                  first occurrence in source, only if
                  report_first_use=True),
              "first_occurrence_char_offset": int,
          }
          "advice": text hint with the lab POLICY rule

    Lab POLICY:
        Every symbol that appears in an equation MUST be defined in
        ONE of:
          * Nomenclature / IEEEdescription block before Section I.
          * "where $X$ is ..." clause immediately after equation.
          * Inline before first use: "Let $X$ denote ..." or
            "The conductivity $\\sigma$ is ...".
        Symbols on the UNIVERSAL_WHITELIST (\\pi, \\omega, t, f, x,
        y, z, ...) are exempt; everything else triggers a reviewer
        comment.
    """
    if not os.path.exists(tex_path):
        return {"error": f"tex_path not found: {tex_path}"}
    try:
        with open(tex_path, encoding="utf-8", errors="replace") as fh:
            src = fh.read()
    except Exception as e:  # noqa: BLE001
        return {"error": f"failed to read tex: {e}"}

    # Build whitelist
    whitelist = set(UNIVERSAL_WHITELIST)
    if extra_whitelist:
        for s in extra_whitelist.split(","):
            s = s.strip()
            if s:
                whitelist.add(s)

    # Pull math blocks + symbols
    math_blocks = _extract_math_blocks(src)
    symbols_used = set()
    for blk in math_blocks:
        symbols_used |= _extract_symbols_from_math(blk)
    # Strip whitelist
    symbols_to_check = symbols_used - whitelist

    # Extract defined symbols
    symbols_defined = _extract_defined_symbols(src)

    # Find undefined
    undefined = sorted(symbols_to_check - symbols_defined)

    # Locate first occurrence for each undefined symbol
    reports = []
    for sym in undefined[:max_reported]:
        rec = {"symbol": sym}
        if report_first_use:
            # Find first occurrence in math context
            # For backslash commands, search for the literal form
            if sym.startswith("\\"):
                # Escape for regex
                search_pat = re.escape(sym)
            else:
                # Single letter -- find in math context: $...X...$
                # Conservative: look for "$" preceded letter "$"-ish
                search_pat = (r"\$[^$]*?(?<![A-Za-z\\])"
                              + re.escape(sym)
                              + r"(?![A-Za-z])[^$]*?\$")
            m = re.search(search_pat, src)
            if m:
                start = max(0, m.start() - 50)
                end = min(len(src), m.end() + 50)
                rec["first_occurrence_context"] = src[start:end].replace(
                    "\n", " ")
                rec["first_occurrence_char_offset"] = m.start()
            else:
                rec["first_occurrence_context"] = "(could not locate)"
                rec["first_occurrence_char_offset"] = -1
        reports.append(rec)

    advice = (
        "Every symbol used in an equation MUST be defined.  Add the "
        "missing symbols to a Nomenclature / IEEEdescription block "
        "BEFORE Section I, OR add a 'where $X$ is ...' clause "
        "immediately after the equation that introduces them.  Use "
        "extra_whitelist for project-specific generic indices that "
        "are conventionally generic (e.g. element-order p, q).  "
        "Failure to define symbols is the #11 reviewer comment in "
        "paper_writing_em_paper_style('reviewer_patterns')."
    )

    return {
        "tex_path": tex_path,
        "n_math_blocks": len(math_blocks),
        "n_symbols_used": len(symbols_used),
        "n_symbols_defined": len(symbols_defined),
        "n_undefined": len(undefined),
        "undefined_symbols": reports,
        "all_undefined": undefined if len(undefined) <= max_reported else
                          undefined[:max_reported] + ["...(truncated)"],
        "whitelist_size": len(whitelist),
        "truncated": len(undefined) > max_reported,
        "advice": advice,
    }
