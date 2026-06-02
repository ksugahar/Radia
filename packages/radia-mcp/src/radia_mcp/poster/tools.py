"""Tool functions for the poster sub-skill of mcp-server-document.

All functions are plain callables; they get wrapped by ``@mcp.tool()`` in
``poster/__init__.py``'s ``register()`` (called from the top-level
``server.py``). No FastMCP import here.

Origin: COMSOL Conference 2025 Tokyo poster ``Kelvin.tex`` (2025-12-05),
promoted into a reusable A1-portrait Japanese academic poster template +
lint + compile pipeline.
"""
from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
from typing import Any


_HERE = pathlib.Path(__file__).resolve().parent
_TEMPLATES = _HERE / "templates"


# Tier 1 — lint expansion ----------------------------------------------------
from .plans.T1_word_count import poster_word_count  # noqa: E402,F401
from .plans.T2_line_length import poster_line_length  # noqa: E402,F401
from .plans.T3_caption_self_contained import poster_caption_self_contained  # noqa: E402,F401
from .plans.T4_jp_font_check import poster_jp_font_check  # noqa: E402,F401
from .plans.T5_fontsize_by_distance import poster_fontsize_by_distance  # noqa: E402,F401
from .plans.T6_typography_lints import poster_typography_lints  # noqa: E402,F401

# Tier 2 — accessibility -----------------------------------------------------
from .plans.T7_color_contrast_wcag import poster_color_contrast_wcag  # noqa: E402,F401
from .plans.T8_colorblind_hint import poster_colorblind_hint  # noqa: E402,F401
from .plans.T9_color_count_321 import poster_color_count_321  # noqa: E402,F401

# Tier 3 — BetterPoster ------------------------------------------------------
from .plans.T10_template_betterposter import poster_template_betterposter  # noqa: E402,F401
from .plans.T11_betterposter_billboard_lint import poster_betterposter_billboard_lint  # noqa: E402,F401
from .plans.T12_zone_balance import poster_zone_balance_check  # noqa: E402,F401

# Tier 4 — QR ----------------------------------------------------------------
from .plans.T13_qr_inject import poster_qr_inject  # noqa: E402,F401
from .plans.T14_qr_audit import poster_qr_audit  # noqa: E402,F401

# Tier 5 — Print readiness ---------------------------------------------------
from .plans.T15_print_readiness_audit import poster_print_readiness_audit  # noqa: E402,F401
from .plans.T16_font_embed_check import poster_font_embed_check  # noqa: E402,F401

# Tier 6 — Intelligence Layer ------------------------------------------------
from .plans.T17_health_report import poster_health_report  # noqa: E402,F401
from .plans.T18_root_cause_diagnosis import poster_root_cause_diagnosis  # noqa: E402,F401
from .plans.T19_next_5_actions import poster_next_5_actions  # noqa: E402,F401
from .plans.T20_rewrite_suggest import poster_rewrite_suggest  # noqa: E402,F401
from .plans.T21_adaptive_health_report import poster_adaptive_health_report  # noqa: E402,F401
from .plans.T22_run_full_workflow import poster_run_full_workflow  # noqa: E402,F401

# Tier 7 — Entry width -------------------------------------------------------
from .plans.T23_from_pptx import poster_from_pptx  # noqa: E402,F401
from .plans.T24_from_paper_tex import poster_from_paper_tex  # noqa: E402,F401
from .plans.T25_elevator_pitch import poster_elevator_pitch_generate  # noqa: E402,F401
from .plans.T26_qa_anticipation import poster_qa_anticipation_list  # noqa: E402,F401


# ---------------------------------------------------------------------------
# 0. Skill doc accessor
# ---------------------------------------------------------------------------
def poster_skill_doc() -> str:
    """Return the poster sub-skill's ``skill.md`` documentation as text."""
    return (_HERE / "skill.md").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# 1. Template — Kelvin.tex (A1 portrait, bxjsarticle + tcolorbox + multicols)
# ---------------------------------------------------------------------------
# Single-line regexes that target the field positions in the Kelvin.tex
# template. Each pattern captures the leading sigil + spacing so the
# substitution preserves font-size / color directives verbatim.
_FIELD_PATTERNS: dict[str, re.Pattern[str]] = {
    # \fontsize{60}{72}\selectfont\bfseries\color{headercolor} <TITLE>
    "title": re.compile(
        r"(\{\\fontsize\{60\}\{72\}\\selectfont\\bfseries\\color\{headercolor\}\s+)([^}]+?)(\s*\})"
    ),
    # \fontsize{40}{48}\selectfont\bfseries\color{headercolor} \quad ~~<SUBTITLE>
    "subtitle": re.compile(
        r"(\{\\fontsize\{40\}\{48\}\\selectfont\\bfseries\\color\{headercolor\}\s+\\quad\s+~~)([^}]+?)(\s*\})"
    ),
    # \fontsize{35}{42}\selectfont <AUTHOR>
    "author": re.compile(
        r"(\{\\fontsize\{35\}\{42\}\\selectfont\s+)([^}]+?)(\s*\})"
    ),
    # Bottom footer banner inside the gradient tikz node.
    "footer": re.compile(
        r"(\\fontsize\{24\}\{40\}\\selectfont\s+)(Excerpt from proceedings of the[^}]+?)(\s*\})",
        re.DOTALL,
    ),
}


def poster_template_kelvin(out_path: str | None = None,
                           replace: dict[str, str] | None = None) -> str:
    """Return (or write) the A1-portrait Japanese poster template.

    The template is the verbatim ``Kelvin.tex`` from COMSOL Conference
    2025 Tokyo. Optional ``replace`` keys ({title, subtitle, author,
    abstract, footer}) perform conservative single-field substitutions
    that preserve the surrounding font-size / color directives.

    Args:
        out_path: Optional destination ``.tex`` path. If omitted, the
                  rendered template is returned as a string. If the
                  destination directory does not exist it is created.
        replace: Optional ``{field: replacement}`` dict. Unknown keys are
                 ignored with a warning at the head of the returned
                 status string. The ``abstract`` key replaces the entire
                 first paragraph after the ``はじめに`` ``sectionbox``.

    Returns:
        - When ``out_path`` is None: the rendered ``.tex`` source.
        - When ``out_path`` is given: a status string with the output
          path, number of substitutions made, and any warnings.
    """
    tex_src = (_TEMPLATES / "poster_a1_portrait.tex").read_text(encoding="utf-8")

    warnings: list[str] = []
    n_subs = 0

    if replace:
        unknown = set(replace) - (set(_FIELD_PATTERNS) | {"abstract"})
        if unknown:
            warnings.append(f"unknown replace keys ignored: {sorted(unknown)}")

        for field, pat in _FIELD_PATTERNS.items():
            if field not in replace:
                continue
            new_val = replace[field]
            new_src, n = pat.subn(
                lambda m, v=new_val: f"{m.group(1)}{v}{m.group(3)}",
                tex_src,
                count=1,
            )
            if n == 0:
                warnings.append(f"field {field!r}: pattern did not match (template edited?)")
            else:
                tex_src = new_src
                n_subs += 1

        if "abstract" in replace:
            # Replace the paragraph block between the はじめに sectionbox
            # and the next \vspace{5pt} ... \begin{multicols}{2}.
            abs_pat = re.compile(
                r"(\\end\{tcolorbox\}\s*\{\\fontsize\{24pt\}\{33pt\}\\selectfont\\par\s*\n)"
                r"(.*?)"
                r"(\n\s*\\par\s*\})",
                re.DOTALL,
            )
            tex_src, n = abs_pat.subn(
                lambda m, v=replace["abstract"]: f"{m.group(1)}{v}{m.group(3)}",
                tex_src,
                count=1,
            )
            if n == 0:
                warnings.append("field 'abstract': pattern did not match (template edited?)")
            else:
                n_subs += 1

    if out_path is None:
        head = ""
        if warnings:
            head = "% poster_template_kelvin warnings:\n" + "".join(
                f"%   - {w}\n" for w in warnings
            ) + "\n"
        return head + tex_src

    out = pathlib.Path(out_path)
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tex_src, encoding="utf-8")

    lines = [f"Wrote poster template to: {out}", f"  substitutions: {n_subs}"]
    if warnings:
        lines.append("  warnings:")
        lines.extend(f"    - {w}" for w in warnings)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 2. Lint — poster-specific (1m reading distance, figure dominance, flow)
# ---------------------------------------------------------------------------
_FONTSIZE_RE = re.compile(r"\\fontsize\{(\d+)(?:pt)?\}\{(\d+)(?:pt)?\}")
_INCLUDEGRAPHICS_RE = re.compile(r"\\includegraphics\b[^{]*\{([^}]+)\}")
_SECTIONBOX_RE = re.compile(r"\\begin\{tcolorbox\}\[sectionbox\]")
_CONTENTBOX_RE = re.compile(r"\\begin\{tcolorbox\}\[contentbox", re.IGNORECASE)
_MULTICOLS_RE = re.compile(r"\\begin\{multicols\}\{(\d+)\}")
_CITE_RE = re.compile(r"\\item\s+([^\n]+)")
_REFLIST_RE = re.compile(
    r"\\begin\{enumerate\}\[label=\{\[\\arabic\*\]\}.*?\\end\{enumerate\}",
    re.DOTALL,
)

# Severity thresholds.
_BODY_FONT_MIN = 24   # pt — body text at 1m reading distance
_REF_FONT_MIN = 20    # pt — references can be slightly smaller
_TITLE_FONT_MIN = 50  # pt — main title must be readable from 3m


def poster_lint(tex_path: str) -> str:
    """Lint a poster .tex against poster-specific (not slide) criteria.

    Checks five axes derived from Kelvin.tex experience:

    1. **Font size floors** — body text ≥24pt, references ≥20pt, title ≥50pt.
       Violations are HIGH severity (1m読みづらさ直結).
    2. **Figure dominance** — at least one ``\\includegraphics`` per
       ``contentbox``. Text-only contentbox is MEDIUM severity.
    3. **Column structure** — exactly one ``\\begin{multicols}{2}`` (or 3)
       so the reading flow is unambiguous.
    4. **Section count** — 3-5 ``sectionbox`` recommended (はじめに,
       本論 1-3, まとめ). <2 = thin, >6 = overcrowded.
    5. **Reference distribution** — references should appear inside
       contentbox (per-section), not gathered at the end. Counts
       ``\\item`` inside each contentbox.

    Args:
        tex_path: Path to the poster ``.tex`` source.

    Returns:
        Multi-line status string with PASS/FAIL per axis + severity tag
        + remediation hint. Designed to be readable in raw form, not JSON.
    """
    src_path = pathlib.Path(tex_path)
    if not src_path.is_absolute():
        src_path = pathlib.Path.cwd() / src_path
    if not src_path.exists():
        return f"Error: file not found: {src_path}"
    tex = src_path.read_text(encoding="utf-8")

    issues: list[tuple[str, str, str]] = []  # (severity, axis, message)

    # 1) Font size floors
    fontsizes = [(int(a), int(b)) for a, b in _FONTSIZE_RE.findall(tex)]
    if not fontsizes:
        issues.append(("HIGH", "fontsize", "no \\fontsize directives found — poster may rely on default 10pt body"))
    else:
        title_max = max(a for a, _ in fontsizes)
        if title_max < _TITLE_FONT_MIN:
            issues.append((
                "HIGH", "fontsize",
                f"largest font is {title_max}pt (<{_TITLE_FONT_MIN}pt); a 3m-readable title is needed",
            ))
        # Sizes 18-30pt cover both body and references. Anything below the
        # absolute floor (_REF_FONT_MIN = 20pt) is unreadable at 1m.
        small = [a for a, _ in fontsizes if 18 <= a <= 30]
        if small and min(small) < _REF_FONT_MIN:
            issues.append((
                "HIGH", "fontsize",
                f"smallest font {min(small)}pt < {_REF_FONT_MIN}pt absolute floor; "
                f"unreadable at 1m",
            ))
        # Body paragraphs (line height >=30pt suggests body, not refs) must be ≥24pt.
        body_sizes = [a for a, b in fontsizes if 18 <= a <= 30 and b >= 30]
        if body_sizes and min(body_sizes) < _BODY_FONT_MIN:
            issues.append((
                "MEDIUM", "fontsize",
                f"body-paragraph font {min(body_sizes)}pt < {_BODY_FONT_MIN}pt; "
                f"references can stay smaller, but body text should be 24pt+",
            ))

    # 2) Figure dominance per contentbox
    contentbox_count = len(_CONTENTBOX_RE.findall(tex))
    fig_count = len(_INCLUDEGRAPHICS_RE.findall(tex))
    if contentbox_count == 0:
        issues.append(("MEDIUM", "structure", "no contentbox detected — poster lacks chunked structure"))
    elif fig_count < contentbox_count:
        issues.append((
            "MEDIUM", "figure_dominance",
            f"only {fig_count} figures for {contentbox_count} contentbox (expect ≥1 per box)",
        ))

    # 3) Column structure
    multicols = _MULTICOLS_RE.findall(tex)
    if not multicols:
        issues.append(("MEDIUM", "columns", "no \\begin{multicols} — single-column poster has poor reading flow"))
    elif len(multicols) > 1:
        issues.append((
            "MEDIUM", "columns",
            f"{len(multicols)} multicols regions — fragments the Z-pattern reading flow",
        ))

    # 4) Section count
    section_count = len(_SECTIONBOX_RE.findall(tex))
    if section_count < 2:
        issues.append((
            "MEDIUM", "section_count",
            f"only {section_count} sectionbox — typical poster has 3-5 (intro / body 1-3 / summary)",
        ))
    elif section_count > 6:
        issues.append((
            "LOW", "section_count",
            f"{section_count} sectionbox — risks overcrowding for A1 panel",
        ))

    # 5) Reference distribution — count refs inside enumerate blocks
    ref_blocks = _REFLIST_RE.findall(tex)
    total_refs = sum(len(_CITE_RE.findall(b)) for b in ref_blocks)
    if total_refs == 0:
        issues.append(("LOW", "references", "no [N] references detected — academic poster usually cites prior work"))
    elif len(ref_blocks) == 1 and total_refs > 3:
        issues.append((
            "LOW", "references",
            f"all {total_refs} references gathered in 1 block — prefer per-contentbox distribution",
        ))

    # Render result
    lines = [
        f"poster_lint: {src_path}",
        f"  contentbox={contentbox_count}, sectionbox={section_count}, figures={fig_count}, refs={total_refs}",
        "",
    ]
    if not issues:
        lines.append("PASS — no issues detected.")
    else:
        # Sort HIGH > MEDIUM > LOW
        order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
        issues.sort(key=lambda x: order[x[0]])
        lines.append(f"FAIL — {len(issues)} issue(s):")
        for sev, axis, msg in issues:
            lines.append(f"  [{sev}] {axis}: {msg}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 3. Compile — platex + dvipdfmx (Japanese), with PDF lock release
# ---------------------------------------------------------------------------
_PDF_VIEWER_PROCESSES = (
    "Acrobat", "AcroRd32", "SumatraPDF", "FoxitReader", "FoxitPDFReader",
    "msedge", "chrome", "firefox",
)


def _release_pdf_lock(pdf_path: pathlib.Path) -> list[str]:
    """Kill processes that hold ``pdf_path`` open. Windows only.

    Mirrors the texcompile skill's lock-release logic but scoped to a
    single PDF basename match against MainWindowTitle.
    """
    killed: list[str] = []
    if not pdf_path.exists() or os.name != "nt":
        return killed
    # Try to open exclusively; if it succeeds the file is not locked.
    try:
        with open(pdf_path, "r+b"):
            return killed
    except PermissionError:
        pass
    except OSError:
        return killed

    basename = pdf_path.name
    ps_script = (
        "Get-Process | Where-Object { $_.MainWindowTitle -ne '' } | "
        f"Where-Object {{ $_.MainWindowTitle -match '{re.escape(basename)}' "
        "-or $_.ProcessName -in @("
        + ",".join(f"'{p}'" for p in _PDF_VIEWER_PROCESSES)
        + ") } | Select-Object -ExpandProperty Id"
    )
    try:
        out = subprocess.run(
            ["powershell", "-NoProfile", "-Command", ps_script],
            capture_output=True, text=True, timeout=15,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return [f"(lock probe failed: {e})"]
    for pid in (out.stdout or "").split():
        pid = pid.strip()
        if not pid.isdigit():
            continue
        try:
            subprocess.run(
                ["powershell", "-NoProfile", "-Command", f"Stop-Process -Id {pid} -Force"],
                capture_output=True, text=True, timeout=10,
            )
            killed.append(pid)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    return killed


def poster_compile(tex_path: str, engine: str = "platex") -> str:
    """Compile a poster .tex to PDF.

    Two engines are supported:

    - ``platex``: ``platex`` (twice for refs) → ``dvipdfmx``. This is the
      default for Japanese posters built on ``bxjsarticle`` like
      Kelvin.tex.
    - ``pdflatex``: single-engine pipeline (``pdflatex`` twice) for
      English / Unicode-engine posters.

    The output PDF lock (held by Acrobat / Edge / SumatraPDF / etc.) is
    released before the build, mirroring the ``texcompile`` skill.

    Args:
        tex_path: Path to the ``.tex`` source.
        engine: ``"platex"`` (default) or ``"pdflatex"``.

    Returns:
        Multi-line status string with command output tail + final PDF
        size + page count (when extractable from the log).
    """
    src = pathlib.Path(tex_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    if src.suffix.lower() != ".tex":
        return f"Error: not a .tex file: {src}"

    workdir = src.parent
    stem = src.stem
    pdf_path = workdir / f"{stem}.pdf"
    killed = _release_pdf_lock(pdf_path)

    # Clean stale aux files so failed prior runs don't poison the build.
    for ext in (".aux", ".out", ".dvi"):
        p = workdir / f"{stem}{ext}"
        if p.exists():
            try:
                p.unlink()
            except OSError:
                pass

    if engine == "platex":
        if not shutil.which("platex") or not shutil.which("dvipdfmx"):
            return ("Error: platex/dvipdfmx not on PATH. Install TeXLive and "
                    "ensure C:\\texlive\\<year>\\bin\\windows is in PATH, "
                    "or pass engine=\"pdflatex\".")
        commands = [
            ["platex", "-interaction=nonstopmode", src.name],
            ["platex", "-interaction=nonstopmode", src.name],
            ["dvipdfmx", stem],
        ]
    elif engine == "pdflatex":
        if not shutil.which("pdflatex"):
            return "Error: pdflatex not on PATH."
        commands = [
            ["pdflatex", "-interaction=nonstopmode", src.name],
            ["pdflatex", "-interaction=nonstopmode", src.name],
        ]
    else:
        return f"Error: unknown engine {engine!r}; expected 'platex' or 'pdflatex'."

    log_tail: list[str] = []
    for cmd in commands:
        try:
            proc = subprocess.run(
                cmd, cwd=workdir, capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=300,
            )
        except subprocess.TimeoutExpired:
            return f"Error: {cmd[0]} timed out after 300s."
        if proc.returncode != 0:
            tail = (proc.stdout or "").splitlines()[-25:]
            return (f"Error: {cmd[0]} returned {proc.returncode}.\n"
                    "--- log tail ---\n" + "\n".join(tail))
        # Keep the last few lines from the final invocation for the report.
        log_tail = (proc.stdout or "").splitlines()[-5:]

    if not pdf_path.exists():
        return "Error: pipeline completed but no PDF was produced."

    size_kb = pdf_path.stat().st_size / 1024
    # Try to extract page count from the log (best-effort).
    pages = "?"
    log_path = workdir / f"{stem}.log"
    if log_path.exists():
        try:
            log = log_path.read_text(encoding="utf-8", errors="replace")
            m = re.search(r"Output written on .*?\((\d+) pages?", log)
            if m:
                pages = m.group(1)
        except OSError:
            pass

    lines = [
        f"poster_compile: {pdf_path}  ({size_kb:.1f} KB, {pages} page(s))",
        f"  engine: {engine}",
    ]
    if killed:
        lines.append(f"  released PDF lock from PIDs: {','.join(killed)}")
    if log_tail:
        lines.append("  log tail:")
        lines.extend(f"    {ln}" for ln in log_tail)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 4. Figures audit
# ---------------------------------------------------------------------------
def poster_figures_audit(tex_path: str,
                         figures_dir: str | None = None) -> str:
    """Check that every ``\\includegraphics{path}`` resolves to a file.

    Also flags unused files in ``figures_dir`` so cleanup is obvious.

    Args:
        tex_path: Path to the poster ``.tex`` source.
        figures_dir: Optional figures directory. If omitted, defaults to
                     ``<tex_path>/../figures``. Resolution rules follow
                     the LaTeX ``\\graphicspath`` default: relative paths
                     in ``\\includegraphics`` are looked up against the
                     ``.tex`` directory, and a stem-only argument is
                     resolved against ``figures_dir`` with common
                     extensions (.pdf, .png, .jpg, .eps).

    Returns:
        Multi-line audit report — missing references first, then unused
        files, then a summary line.
    """
    src = pathlib.Path(tex_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"

    tex = src.read_text(encoding="utf-8")
    refs = _INCLUDEGRAPHICS_RE.findall(tex)

    figdir = pathlib.Path(figures_dir) if figures_dir else (src.parent / "figures")
    if not figdir.is_absolute():
        figdir = pathlib.Path.cwd() / figdir

    missing: list[str] = []
    found: set[pathlib.Path] = set()
    extensions = (".pdf", ".png", ".jpg", ".jpeg", ".eps")

    for ref in refs:
        candidates = []
        rp = pathlib.Path(ref)
        # 1) relative to tex dir
        candidates.append(src.parent / rp)
        # 2) relative to figures_dir
        candidates.append(figdir / rp.name)
        # 3) stem-only lookup against figures_dir with common extensions
        if rp.suffix == "":
            for ext in extensions:
                candidates.append(src.parent / rp.with_suffix(ext))
                candidates.append(figdir / (rp.name + ext))
        hit = next((c for c in candidates if c.exists()), None)
        if hit is None:
            missing.append(ref)
        else:
            found.add(hit.resolve())

    unused: list[str] = []
    if figdir.exists():
        for f in figdir.iterdir():
            if f.is_file() and f.suffix.lower() in extensions:
                if f.resolve() not in found:
                    unused.append(f.name)

    lines = [
        f"poster_figures_audit: {src}",
        f"  figures dir: {figdir} ({'exists' if figdir.exists() else 'MISSING'})",
        f"  \\includegraphics references: {len(refs)}",
    ]
    if missing:
        lines.append(f"  MISSING ({len(missing)}):")
        lines.extend(f"    - {m}" for m in missing)
    if unused:
        lines.append(f"  UNUSED in figures dir ({len(unused)}):")
        lines.extend(f"    - {u}" for u in sorted(unused))
    if not missing and not unused:
        lines.append("  PASS — all references resolved, no orphan files.")
    return "\n".join(lines)
