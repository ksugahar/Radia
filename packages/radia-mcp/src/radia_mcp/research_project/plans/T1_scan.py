"""Tier 1 — research_project_scan.

Walk a directory and classify artifacts:
    - paper.tex / *.tex (paper drafts)
    - poster.tex / *_poster.tex (poster drafts)
    - slides.pptx / *.pptx (presentation)
    - grant proposal indicators (申請書, application, proposal)
    - .bib files (citations)
    - PDFs (final outputs)
    - .cst (CST projects), .mph (COMSOL), .gid (GiD) — engineering source
"""
from __future__ import annotations

import pathlib


_HEURISTICS = [
    ("paper",        lambda n, t: ("paper" in n or "manuscript" in n or "tmag" in n) and t == ".tex"),
    ("poster",       lambda n, t: "poster" in n and t == ".tex"),
    ("grant",        lambda n, t: ("grant" in n or "kaken" in n or "申請" in n or "proposal" in n) and t == ".tex"),
    ("pptx",         lambda n, t: t in (".pptx", ".ppt", ".pptm")),  # split downstream by doc_convert_classify
    ("bib",          lambda n, t: t == ".bib"),
    ("output_pdf",   lambda n, t: t == ".pdf"),
    ("excel",        lambda n, t: t in (".xlsx", ".xlsm")),
    ("cst_project",  lambda n, t: t == ".cst"),
    ("comsol",       lambda n, t: t == ".mph"),
    ("matlab",       lambda n, t: t in (".m", ".mat")),
    ("readme",       lambda n, t: n in ("readme.md", "readme.txt", "claude.md")),
]


def _split_pptx_kinds(paths):
    """Sub-classify each .pptx as 'figure_pptx' or 'presentation_pptx'.

    Rule: parent directory is ``figures`` / ``figs`` / ``figure`` →
    figure_pptx; everything else → presentation_pptx.
    """
    from ...doc_convert.plans.T12_classify import classify_impl
    figs = []
    pres = []
    for p in paths:
        if classify_impl(str(p)) == "figure":
            figs.append(p)
        else:
            pres.append(p)
    return figs, pres


def research_project_scan(directory: str, recursive: bool = True) -> str:
    """Classify every file in ``directory`` by artifact kind.

    Args:
        directory: Project root.
        recursive: If True (default), descend.

    Returns:
        Per-kind file count + sample paths.
    """
    root = pathlib.Path(directory)
    if not root.is_absolute():
        root = pathlib.Path.cwd() / root
    if not root.is_dir():
        return f"Error: not a directory: {root}"

    bucket: dict[str, list[pathlib.Path]] = {}
    iterator = root.rglob("*") if recursive else root.glob("*")
    for p in iterator:
        if not p.is_file():
            continue
        name = p.name.lower()
        ext = p.suffix.lower()
        kind = next((k for k, fn in _HEURISTICS if fn(name, ext)), None)
        if kind:
            bucket.setdefault(kind, []).append(p)

    # Sub-classify pptx by parent directory ('figures' → figure, else → presentation).
    if "pptx" in bucket:
        figs, pres = _split_pptx_kinds(bucket.pop("pptx"))
        if figs:
            bucket["figure_pptx"] = figs
        if pres:
            bucket["presentation_pptx"] = pres

    lines = [f"research_project_scan: {root}",
             f"  recursive: {recursive}",
             f"  kinds found: {len(bucket)}"]
    for k in ("grant", "paper", "poster",
              "presentation_pptx", "figure_pptx",
              "bib", "output_pdf", "excel", "cst_project", "comsol",
              "matlab", "readme"):
        if k not in bucket:
            continue
        files = bucket[k]
        lines.append(f"  [{k}] {len(files)} file(s):")
        for f_ in files[:5]:
            lines.append(f"    - {f_.relative_to(root)}")
        if len(files) > 5:
            lines.append(f"    ... ({len(files) - 5} more)")
    return "\n".join(lines)
