"""Tier 7 — poster_from_paper_tex.

Extract Introduction / Methods / Results / Conclusion sections from an
existing paper .tex and emit a poster skeleton with each contentbox
prefilled. Aggressive abstraction: the user must hand-edit, but starts
from real text rather than placeholders.
"""
from __future__ import annotations

import pathlib
import re


_SECTION_RE = re.compile(
    r"\\section\*?\{([^}]+)\}", re.IGNORECASE
)
_KEEP_FIRST_N_PARAGRAPHS = 2
_MAX_CHARS_PER_BOX = 600

_SECTION_ALIASES = {
    "introduction": "intro",
    "background": "intro",
    "method": "method",
    "methods": "method",
    "materials": "method",
    "formulation": "method",
    "result": "result",
    "results": "result",
    "experiment": "result",
    "experiments": "result",
    "discussion": "discussion",
    "conclusion": "conclusion",
    "conclusions": "conclusion",
    "summary": "conclusion",
}


def _normalize(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def _escape(s: str) -> str:
    return s.replace("%", "\\%").replace("#", "\\#").replace("&", "\\&").replace("_", "\\_")


def poster_from_paper_tex(paper_tex: str,
                          out_tex_path: str | None = None,
                          title: str | None = None) -> str:
    """Convert a paper .tex into a Kelvin-style poster skeleton.

    Heuristics:
      - Detect ``\\section{Introduction}`` and aliases (背景, はじめに, etc.)
      - Keep first ~2 paragraphs per section, capped at 600 chars
      - Drop equations / tikz / tabular blocks
      - Preserve ``\\cite{...}`` markers for the user to relocate

    Args:
        paper_tex: Path to the source paper ``.tex``.
        out_tex_path: Optional output path. None returns the source.
        title: Optional override; default uses paper's ``\\title{}`` or
               fallback "Poster Title".

    Returns:
        Status string when ``out_tex_path`` is given; .tex otherwise.
    """
    p = pathlib.Path(paper_tex)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    src = p.read_text(encoding="utf-8", errors="replace")

    # Title
    if title is None:
        tm = re.search(r"\\title\{([^}]+)\}", src)
        title = tm.group(1) if tm else "Poster Title"

    # Locate sections
    matches = list(_SECTION_RE.finditer(src))
    sections: list[tuple[str, str]] = []
    for i, m in enumerate(matches):
        name = _normalize(m.group(1))
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(src)
        body = src[start:end]
        # Strip heavy env
        body = re.sub(
            r"\\begin\{(equation|align|figure|tabular|tikzpicture|verbatim|comment|table)\*?\}.*?\\end\{\1\*?\}",
            "", body, flags=re.DOTALL,
        )
        # Take first N paragraphs.
        paras = [_normalize(p_) for p_ in re.split(r"\n\s*\n", body) if _normalize(p_)]
        snippet = "\n\n".join(paras[:_KEEP_FIRST_N_PARAGRAPHS])
        if len(snippet) > _MAX_CHARS_PER_BOX:
            snippet = snippet[:_MAX_CHARS_PER_BOX].rsplit(" ", 1)[0] + "…"
        sections.append((name, snippet))

    # Filter to canonical poster sections.
    keep = []
    for name, body in sections:
        for alias, canon in _SECTION_ALIASES.items():
            if alias in name.lower():
                keep.append((canon, name, body))
                break

    if not keep:
        return f"Error: no recognized sections found in {p}"

    # Emit template body.
    parts = [
        "\\documentclass[a1paper,portrait]{bxjsarticle}",
        "\\usepackage[dvipdfmx]{graphicx,color}",
        "\\usepackage{amsmath,multicol,enumitem,bm,geometry,tikz,tcolorbox,caption,anyfontsize,xcolor}",
        "\\geometry{margin=2cm,top=1cm,bottom=1cm}",
        "\\tcbuselibrary{breakable,skins}",
        "\\usepackage{newtxtext,newtxmath}",
        "\\renewcommand{\\familydefault}{\\sfdefault}",
        "\\definecolor{headercolor}{RGB}{0,51,153}",
        "\\definecolor{sectionblue}{RGB}{51,102,204}",
        "\\tcbset{",
        "  contentbox/.style={colback=white,colframe=sectionblue,boxrule=2pt,arc=3mm,boxsep=10pt,left=10pt,right=10pt,top=10pt,bottom=10pt},",
        "  subsectionbox/.style={colback=sectionblue!15,coltext=sectionblue!80!black,boxrule=0pt,arc=3mm,left=12pt,right=12pt,top=6pt,bottom=6pt,fontupper=\\large\\bfseries},",
        "}",
        "\\begin{document}",
        "\\begin{flushleft}",
        f"  {{\\fontsize{{60}}{{72}}\\selectfont\\bfseries\\color{{headercolor}} {_escape(title)}}}",
        "\\end{flushleft}",
        "\\hrule",
        "\\vspace{10pt}",
        "\\begin{multicols}{2}",
    ]
    for canon, name, body in keep:
        parts.append("\\begin{tcolorbox}[contentbox,breakable]")
        parts.append("  \\begin{tcolorbox}[subsectionbox]")
        parts.append(f"    {{\\fontsize{{30pt}}{{36pt}}\\selectfont {_escape(name)}}}")
        parts.append("  \\end{tcolorbox}")
        parts.append("  {\\fontsize{24pt}{33pt}\\selectfont\\par")
        parts.append(f"  {_escape(body)}")
        parts.append("  \\par}")
        parts.append("\\end{tcolorbox}")
    parts.append("\\end{multicols}")
    parts.append("\\end{document}")

    src_out = "\n".join(parts)
    if out_tex_path is None:
        return src_out
    out = pathlib.Path(out_tex_path)
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(src_out, encoding="utf-8")
    return f"poster_from_paper_tex: wrote {out}\n  sections imported: {len(keep)} ({', '.join(s[0] for s in keep)})"
