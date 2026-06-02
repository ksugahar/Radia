"""radia_mcp.figure._lab_api -- misuse-proof slide / paper figure API.

Born from a real incident (CEFC 2026 oral slides): figures were authored
at 150 mm, embedded with ``\\includegraphics[width=\\linewidth]`` into an
8 cm column (a 0.5x downscale -> on-page font 6 pt, not 10), carried
in-figure titles, and would have silently fallen back to DejaVu Serif on a
machine without Times New Roman.  This module makes the COMPLIANT path the
easy path and makes every violation FAIL LOUD (No-Fallback policy):

    from radia_mcp.figure import lab_figure, save_lab_figure
    fig, ax = lab_figure(embed_width_cm=8.0)     # authored AT 8 cm, 10 pt, TNR verified
    ax.plot(...)
    info = save_lab_figure(fig, "out", embed_width_cm=8.0)
    print(info["latex"])   # \\includegraphics[width=8.00cm]{out}  (embed at 100%)

Gates run at save (each RAISES on violation):
  * no in-figure title (the LaTeX caption / beamer frametitle carries it),
  * Times New Roman ACTUALLY used (no silent DejaVu fallback) -- pre
    (findfont, no fallback) AND post (scan the saved PDF bytes),
  * TrueType font embedding (pdf.fonttype == 42),
  * no Japanese text / no CJK font in the PDF.

The embed width is stashed on the figure so save_lab_figure can emit the
exact ``\\includegraphics[width=<cm>]`` snippet -- embed at 100% and the
on-page font is exactly 10 pt at that width.
"""
from __future__ import annotations

import os
import re

from ._paper_figure import (
    _get_mpl,
    _check_no_in_figure_title,
    _check_no_japanese_text,
    _check_no_japanese_font_in_pdf,
    _check_colors_are_cvd_safe,
)


# ------------------------------------------------------------------
# Times New Roman strictness (no silent DejaVu fallback)
# ------------------------------------------------------------------
def _assert_times_new_roman(strict: bool = True) -> None:
    """Raise if Times New Roman is REQUESTED in rcParams but matplotlib would
    silently fall back to another serif (DejaVu).  No-op when TNR is not the
    requested serif (e.g. a deliberately sans figure)."""
    plt, _ = _get_mpl()
    serifs = [str(s).lower() for s in plt.rcParams.get("font.serif", [])]
    fam = plt.rcParams.get("font.family", [])
    fam = [fam] if isinstance(fam, str) else list(fam)
    requests_tnr = ("serif" in [str(f).lower() for f in fam]
                    and "times new roman" in serifs)
    if not requests_tnr:
        return
    import matplotlib.font_manager as fm
    try:
        path = fm.findfont("Times New Roman", fallback_to_default=False)
    except Exception:
        path = None
    if not path or "times" not in os.path.basename(str(path)).lower():
        msg = (
            "Times New Roman is requested (font.serif[0]) but matplotlib "
            "cannot resolve it -- it would SILENTLY render DejaVu Serif. "
            "Install Times New Roman (it ships with Windows) and rebuild the "
            "cache:\n"
            "    import matplotlib.font_manager as fm\n"
            "    fm._load_fontmanager(try_read_cache=False)\n"
            f"(findfont resolved to: {path!r})"
        )
        if strict:
            raise RuntimeError(msg)


def _check_tnr_in_pdf(pdf_path: str) -> list:
    """Belt-and-suspenders post-render check: if the saved PDF embeds DejaVu
    Serif (when a serif/Times figure was intended), the font silently fell
    back.  Returns a list of violation strings (empty = clean)."""
    try:
        with open(pdf_path, "rb") as f:
            data = f.read()
    except OSError:
        return []
    bad = []
    if b"DejaVuSerif" in data or b"DejaVu Serif" in data:
        bad.append(
            f"{os.path.basename(pdf_path)} embeds DejaVu Serif "
            "-- Times New Roman silently fell back."
        )
    return bad


# ------------------------------------------------------------------
# Misuse-proof figure creation + save
# ------------------------------------------------------------------
def lab_figure(embed_width_cm, aspect: float = 0.62, nrows: int = 1, ncols: int = 1,
               use_times_roman: bool = True, legend_small=None):
    """Create ``(fig, axes)`` authored AT the on-page embed width with verified
    Times New Roman + 10 pt.

    Embed the result with ``\\includegraphics[width=<embed_width_cm>cm]`` at
    100% and the on-page body font is exactly 10 pt at that width -- the
    lab rule.  Do NOT scale with ``width=\\linewidth`` unless \\linewidth
    equals ``embed_width_cm``.

    Args:
        embed_width_cm: width the figure will OCCUPY on the page / slide.
        aspect: height / width.   nrows, ncols: subplot grid.
        use_times_roman: enforce TNR (raises if unavailable).
        legend_small: shrink legend to 8 pt (default True for ncols >= 2).

    Raises RuntimeError (No-Fallback) if TNR is requested but unavailable.
    """
    from .tools import apply_lab_style
    plt, _ = _get_mpl()
    if legend_small is None:
        legend_small = ncols >= 2
    target = ("digest_double_column_side_by_side" if legend_small
              else "digest_single_column")
    w_in, h_in = apply_lab_style(target=target, embed_width_cm=float(embed_width_cm),
                                 aspect=aspect, use_times_roman=use_times_roman)
    if use_times_roman:
        _assert_times_new_roman()
    fig, axes = plt.subplots(nrows, ncols, figsize=(w_in, h_in))
    # Stash the embed width so save_lab_figure can emit the exact LaTeX width.
    try:
        fig._lab_embed_width_cm = float(embed_width_cm)
    except Exception:
        pass
    return fig, axes


def save_lab_figure(fig, path_no_ext, embed_width_cm=None, *,
                    allow_in_figure_title: bool = False, check_cvd: bool = False,
                    save_pdf: bool = True, save_png: bool = True, dpi: int = 300,
                    tighten=True) -> dict:
    """Save a figure with FAIL-LOUD lab gates; return a dict whose ``latex``
    key is the exact ``\\includegraphics`` snippet (embed at 100% -> 10 pt)
    and whose ``axes_fraction`` reports how much of the figure the axes fill.

    Gates (each raises on violation): no in-figure title (unless
    ``allow_in_figure_title``), Times New Roman actually used (pre + post),
    TrueType embedding, no Japanese, and optionally CVD-safe colors.

    ``tighten`` (True / a float target, default 0.80): push the axes box to
    the limit WITHIN the fixed figure size (``auto_tighten`` -- overhang-safe,
    so labels never clip).  This is the supported "bbox ぎりぎり" path: the
    figure SIZE is preserved (embed width unchanged -> on-page font stays
    10 pt).  We deliberately do NOT use ``bbox_inches='tight'``, which would
    change the width and break that guarantee.
    """
    plt, _ = _get_mpl()
    ew = embed_width_cm if embed_width_cm is not None else getattr(fig, "_lab_embed_width_cm", None)
    base = os.path.basename(path_no_ext)

    # ---- pre-save gates ----
    if not allow_in_figure_title:
        v = _check_no_in_figure_title(fig)
        if v:
            raise ValueError(
                "save_lab_figure: in-figure title detected -- the LaTeX "
                "\\caption{} / beamer frametitle carries the title, not the "
                "figure.  Remove:\n  " + "\n  ".join(v))
    _assert_times_new_roman()
    jv = _check_no_japanese_text(fig)
    if jv:
        raise ValueError(
            "save_lab_figure: Japanese / CJK text in the figure (papers and "
            "slides use English labels):\n  " + "\n  ".join(str(x) for x in jv))
    if check_cvd:
        cv = _check_colors_are_cvd_safe(fig)
        if cv:
            raise ValueError(
                "save_lab_figure: non-colorblind-safe colors:\n  "
                + "\n  ".join(str(x) for x in cv))

    # ---- maximize the axes within the FIXED figure size (push the bbox to
    #     the limit, overhang-safe).  Deliberately NOT bbox_inches='tight',
    #     which would change the width and break the 10 pt @ width guarantee. ----
    eff = None
    if tighten:
        from ._paper_figure import auto_tighten, measure_figure_efficiency
        fig.tight_layout(pad=0.3)
        tgt = tighten if isinstance(tighten, float) else 0.80
        try:
            eff = auto_tighten(fig, target_axes_fraction=tgt)["final_efficiency"]
        except Exception:
            eff = measure_figure_efficiency(fig)["axes_area_fraction"]
    plt.rcParams["pdf.fonttype"] = 42
    wrote = []
    if save_pdf:
        p = path_no_ext + ".pdf"; fig.savefig(p); wrote.append(p)
    if save_png:
        p = path_no_ext + ".png"; fig.savefig(p, dpi=dpi); wrote.append(p)

    # ---- post-save gates (final bytes) ----
    if save_pdf:
        pdf = path_no_ext + ".pdf"
        tnr_bad = _check_tnr_in_pdf(pdf)
        if tnr_bad:
            raise RuntimeError("save_lab_figure: " + "; ".join(tnr_bad))
        jf = _check_no_japanese_font_in_pdf(pdf)
        if jf:
            raise ValueError("save_lab_figure: CJK font embedded in PDF:\n  "
                             + "\n  ".join(jf))

    if ew:
        latex = r"\includegraphics[width=%.2fcm]{%s}" % (ew, base)
    else:
        latex = (r"\includegraphics[width=\linewidth]{%s}  %% WARN: author width "
                 r"unknown -- pass embed_width_cm to guarantee 10 pt @ width" % base)
    plt.close(fig)
    return {"wrote": wrote, "embed_width_cm": ew, "latex": latex,
            "axes_fraction": eff, "gates": "passed"}


def legend_no_overlap(ax, *, frameon=False, outside_threshold: float = 0.03,
                      axes_shrink: float = 0.74, **legend_kw):
    """Place the legend by the standard rule -- and never on top of the data.

    Escalation (the established practice):
      1. ``loc='best'`` -- matplotlib's own minimise-overlap placement
         (considers lines, patches AND collections).
      2. If a line still passes under the legend (> ``outside_threshold`` of
         its sampled points), try the lab's exhaustive best-of-six in-axes
         search (``find_best_legend_loc``).
      3. If even that overlaps, move the legend OUTSIDE to the right and
         shrink the axes so the legend fits within the SAME figure width
         (no figure resize -> on-page font unchanged).

    Returns the placement used: ``'inside:best'`` / ``'inside:<loc>'`` /
    ``'outside-right'`` / ``'none'`` (no labelled artists).  For a legend
    that simply will not fit, prefer direct labelling
    (``label_curve_endpoints``) over an outside legend.
    """
    from .tools import check_legend_overlap, find_best_legend_loc
    handles, labels = ax.get_legend_handles_labels()
    if not labels:
        return "none"

    def _worst():
        try:
            return max((o["fraction"] for o in check_legend_overlap(ax)), default=0.0)
        except Exception:
            return 0.0

    ax.legend(loc="best", frameon=frameon, **legend_kw)
    if _worst() <= outside_threshold:
        return "inside:best"

    try:
        best, summary = find_best_legend_loc(ax)
        worst_best = max((o["fraction"] for o in summary[best]["overlaps"]), default=0.0)
    except Exception:
        best, worst_best = "best", 1.0
    if worst_best <= outside_threshold:
        ax.legend(loc=best, frameon=frameon, **legend_kw)
        return f"inside:{best}"

    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width * axes_shrink, pos.height])
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=frameon,
              borderaxespad=0.0, **legend_kw)
    return "outside-right"


# ------------------------------------------------------------------
# Figure-compliance audit of a .tex
# ------------------------------------------------------------------
_INCLUDE_RE = re.compile(r"\\includegraphics(?:\[([^\]]*)\])?\{([^}]+)\}")


def audit_tex_figures(tex_path: str) -> dict:
    """Lint every ``\\includegraphics`` in a .tex for embeds that cannot
    guarantee on-page 10 pt @ 8 cm: HEIGHT-constrained embeds (font size
    uncontrolled), ``\\linewidth`` embeds (depends on the column width /
    authoring), figures NOT found, and figure PDFs that embed DejaVu instead
    of Times New Roman.  Returns a dict {tex, n_figures, n_flagged, figures}.
    """
    with open(tex_path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    tex_dir = os.path.dirname(os.path.abspath(tex_path))
    sub = ""
    gm = re.search(r"\\graphicspath\{\{([^}]*)\}\}", text)
    if gm:
        sub = gm.group(1)

    rows = []
    for mm in _INCLUDE_RE.finditer(text):
        opts = mm.group(1) or ""
        fname = mm.group(2).strip()
        risks = []
        if re.search(r"height\s*=", opts) and not re.search(r"width\s*=", opts):
            risks.append("HEIGHT-constrained embed -> on-page font size is "
                         "uncontrolled; author at a width and use width=<cm>.")
        if "linewidth" in opts.lower():
            risks.append("width=\\linewidth -> on-page font depends on the column "
                         "width; verify the figure was authored AT that exact "
                         "width (else it is silently up/down-scaled).")
        has_fixed = bool(re.search(r"width\s*=\s*[\d.]+\s*(cm|mm|in|pt|bp)", opts))

        found = None
        for cand in (fname, fname + ".pdf", fname + ".png"):
            for d in [tex_dir] + ([os.path.join(tex_dir, sub)] if sub else []):
                p = os.path.join(d, cand)
                if os.path.isfile(p):
                    found = p
                    break
            if found:
                break
        if not found:
            risks.append("FILE NOT FOUND next to the .tex.")
        elif found.lower().endswith(".pdf"):
            try:
                with open(found, "rb") as fh:
                    if b"DejaVu" in fh.read():
                        risks.append("figure PDF embeds DejaVu (NOT Times New Roman).")
            except OSError:
                pass
        rows.append({"figure": fname, "options": opts,
                     "fixed_cm_width": has_fixed, "risks": risks})

    return {"tex": os.path.abspath(tex_path), "n_figures": len(rows),
            "n_flagged": sum(1 for r in rows if r["risks"]), "figures": rows}
