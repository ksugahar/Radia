"""Tier 3 — doc_convert_classify.

Distinguish *figure* pptx (embedded into paper/poster as graphics)
from *presentation* pptx (multi-slide deck shown live). The two have
opposite font-size rules and warrant different lints.

**Lab convention (single rule)**:

    - pptx whose parent directory is ``figures/`` → figure
    - any other pptx → presentation

That's it — no filename / slide-count guessing. Keep figure pptx in
``figures/`` subdirectories and the rest elsewhere, and classification
is deterministic.
"""
from __future__ import annotations

import pathlib


def _load_pptx():
    try:
        from pptx import Presentation  # type: ignore
        return Presentation
    except ImportError:
        return None


_FIGURE_PARENT_DIRS = {"figures", "figs", "figure"}


def _classify(path: pathlib.Path, slide_count: int) -> str:
    """Return ``"figure"`` iff the parent directory is ``figures/``."""
    parent = path.parent.name.lower()
    if parent in _FIGURE_PARENT_DIRS:
        return "figure"
    return "presentation"


def doc_convert_classify(pptx_path: str) -> str:
    """Classify a pptx as 'figure', 'presentation', or 'ambiguous'.

    Args:
        pptx_path: Source pptx.

    Returns:
        Multi-line report with kind + slide_count + reasons.
    """
    Pres = _load_pptx()
    if Pres is None:
        return "Error: python-pptx not installed."
    p = pathlib.Path(pptx_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    try:
        prs = Pres(p)
        n = len(prs.slides)
    except Exception as e:
        return f"Error: cannot open pptx: {type(e).__name__}: {e}"
    kind = _classify(p, n)
    parent = p.parent.name
    lines = [f"doc_convert_classify: {p}",
             f"  slide_count: {n}",
             f"  parent_dir: {parent!r}",
             f"  rule: parent in {sorted(_FIGURE_PARENT_DIRS)} → figure, else → presentation",
             f"  -> kind: {kind}"]
    if kind == "figure":
        lines.append("  → use doc_convert_figure_readability_check(... target_width_cm=8)")
    else:
        lines.append("  → use doc_convert_slide_lint (lab rule: ≥20pt everywhere)")
    return "\n".join(lines)


def classify_impl(pptx_path: str) -> str:
    """Internal helper used by health_report dispatch. Returns just the kind."""
    Pres = _load_pptx()
    if Pres is None:
        return "ambiguous"
    p = pathlib.Path(pptx_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return "ambiguous"
    try:
        prs = Pres(p)
        n = len(prs.slides)
    except Exception:
        return "ambiguous"
    return _classify(p, n)
