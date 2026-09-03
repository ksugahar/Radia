"""Image-based PDF layout verification for radia_mcp.paper_writing.

Added 2026-05-26 (radia-mcp v0.88.0).  Companion to the existing
``paper_writing_check_pdf_*`` tools that operate on the PDF binary;
this module RENDERS pages to PNG so that a vision-capable agent
(Claude vision, GPT-4V, human reviewer) can visually inspect the
final TYPESET layout.

Why this matters (lab experience):

  * Overfull \\hbox warnings and PDF binary inspection catch obvious
    text overflow, but they do NOT catch SEMANTIC layout failures:
      - figure pushed to wrong page (\\ref appears 3 pages before the
        figure itself);
      - half-empty page caused by a stubborn [!h] specifier forcing
        a float to its own page;
      - figure spans both columns when it should be single-column;
      - subfigure panels squashed because total widths exceed
        \\textwidth.

  * The cheapest reliable detection is to RENDER the PDF page to a
    thumbnail and LOOK AT IT (or to algorithmically count white
    pixels per page).  This module does both.

Web-research basis (2026-05-26):
  * pymupdf get_pixmap(dpi=150) renders any PDF page to a PIL-
    compatible bitmap (Artifex PyMuPDF docs).
  * Skillnad (hgustafsson/skillnad on GitHub) overlays color-coded
    rectangles on PDF visual diffs from SyncTeX position info --
    inspired the thumbnail-strip + bbox-overlay approach below.
  * LaTeX overfull \\hbox warning catches text overflow but not
    float-placement disasters (Overleaf docs); image-based check
    is the missing piece.

Implementation:
  * pymupdf for PDF -> PNG rendering (already a radia-mcp[rag] dep
    via literature_index; for paper_writing it lives behind a lazy
    import so users without the optional dep see a clear ImportError).
  * algorithmic checks: whitespace fraction per page, multi-column
    overflow, float-far-from-reference.
  * Vision-pass recipe: dump N pages as PNGs into a folder, then the
    agent / Claude vision opens each PNG and reports anomalies.
"""

from __future__ import annotations

import json
import os
import pathlib
from typing import Optional


def _require_pymupdf():
    """Import pymupdf lazily with a clear error message."""
    try:
        import pymupdf  # type: ignore
        return pymupdf
    except ImportError as e:
        raise ImportError(
            "paper_writing image-based layout check requires pymupdf: "
            "pip install pymupdf  (or pip install radia-mcp[rag])."
        ) from e


def _require_pil():
    """Import PIL lazily."""
    try:
        from PIL import Image  # type: ignore
        return Image
    except ImportError as e:
        raise ImportError(
            "paper_writing layout-thumbnail strip requires Pillow: "
            "pip install Pillow"
        ) from e


# ============================================================
# render_pages_to_png
# ============================================================


def paper_writing_render_pages_to_png(
    pdf_path: str,
    out_dir: Optional[str] = None,
    dpi: int = 150,
    page_range: str = "",
    prefix: str = "page",
) -> dict:
    """Render PDF pages to PNG files for visual inspection.

    The fastest way to verify a typeset PDF layout is to LOOK AT THE
    PAGES.  This tool produces one PNG per page so that a
    vision-capable agent (Claude vision, GPT-4V) or a human can scan
    them for figure-placement anomalies, half-empty pages, column
    overflow, or other defects that text-based PDF inspection misses.

    Args:
        pdf_path: absolute path to the typeset PDF.
        out_dir: directory to write PNGs into.  Default: PDF parent
            dir + "_pages".  Created if missing.
        dpi: rendering resolution.  Default 150 = readable but small
            (~1100x1600 px at A4); use 300 for paper-zoom inspection.
        page_range: empty for all pages, or "1-5", "3", "1,3,5-7".
            Page numbers are 1-indexed (matches paper / TeX
            convention; converted to 0-indexed pymupdf internally).
        prefix: filename prefix.  Default "page" -> page_001.png,
            page_002.png, etc.

    Returns:
        dict with keys:
          ``"out_dir"``: where PNGs were written
          ``"pages_rendered"``: list of {"page": N, "path": ..., "size_bytes": ...}
          ``"total_pages"``: N (in source PDF)
          ``"hint"``: how to ask Claude vision to inspect them
    """
    pymupdf = _require_pymupdf()
    pdf = pathlib.Path(pdf_path)
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}
    doc = pymupdf.open(str(pdf))
    total = doc.page_count

    # Parse page_range
    if not page_range or page_range.strip() == "":
        page_nums = list(range(1, total + 1))
    else:
        page_nums = []
        for chunk in page_range.split(","):
            chunk = chunk.strip()
            if "-" in chunk:
                a, b = chunk.split("-", 1)
                page_nums.extend(range(int(a), int(b) + 1))
            else:
                page_nums.append(int(chunk))
        page_nums = sorted(set(p for p in page_nums if 1 <= p <= total))

    if not page_nums:
        doc.close()
        return {"error": f"page_range {page_range!r} matched no pages "
                          f"(PDF has {total} pages)"}

    if out_dir is None:
        out_dir = str(pdf.parent / (pdf.stem + "_pages"))
    os.makedirs(out_dir, exist_ok=True)

    rendered = []
    width = max(3, len(str(total)))
    for pn in page_nums:
        page = doc[pn - 1]
        pix = page.get_pixmap(dpi=dpi)
        out_path = os.path.join(out_dir, f"{prefix}_{pn:0{width}d}.png")
        pix.save(out_path)
        rendered.append({
            "page": pn,
            "path": out_path,
            "size_bytes": os.path.getsize(out_path),
        })
    doc.close()

    return {
        "out_dir": out_dir,
        "pages_rendered": rendered,
        "total_pages": total,
        "n_rendered": len(rendered),
        "hint": (
            "To verify layout: open each PNG and scan for half-empty pages, "
            "figures pushed to wrong pages, column-overflow, or caption-far-"
            "from-figure.  When using Claude in this MCP, feed the PNG paths "
            "back as inline images for vision inspection.  See "
            "paper_writing_layout_visual_recipe() for the full prompt template."
        ),
    }


# ============================================================
# detect_page_whitespace_anomalies
# ============================================================


def paper_writing_detect_page_whitespace_anomalies(
    pdf_path: str,
    whitespace_threshold: float = 0.95,
    dpi: int = 100,
) -> dict:
    """Flag pages that are mostly whitespace (sign of bad float placement).

    Algorithmic check (no agent needed): a page that is >95% white
    is usually a sign that a stubborn [!h] or [H] specifier forced
    a float to its own page when LaTeX could not satisfy the
    constraint with running text.  Classic lab pattern:

        \\begin{figure}[!h]            % bad: !h is rigid
            \\includegraphics{...}
        \\end{figure}
        ... 200 words of text ...
        % LaTeX gives up, dumps figure on its own page

    The fix is usually to use [htbp] or to drop the bang and let
    LaTeX choose.  This tool tells you WHICH pages are wasteful.

    Args:
        pdf_path: typeset PDF.
        whitespace_threshold: 0.0 to 1.0 white-fraction. Default 0.95.
            Normal two-column pages are often 88--90% white at pixel level,
            so the former 0.75 default reported virtually every page.
        dpi: render resolution.  Default 100 = fast (~700x1000 px);
            higher DPI doesn't change the whitespace ratio meaningfully.

    Returns:
        dict with:
          "flagged_pages": list of {page, whitespace_fraction}
          "all_pages": list of {page, whitespace_fraction}
          "total_pages": N
          "advice": text hint for fixing
    """
    pymupdf = _require_pymupdf()
    pdf = pathlib.Path(pdf_path)
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}
    doc = pymupdf.open(str(pdf))

    # We avoid PIL dep by inspecting the raw pixmap samples buffer.
    # pixmap.samples is RGB bytes (3 bytes per pixel).  "White-ish"
    # = all three channels above 240.
    all_pages = []
    flagged = []
    threshold = max(0.0, min(1.0, whitespace_threshold))
    for i, page in enumerate(doc):
        pix = page.get_pixmap(dpi=dpi)
        n = pix.width * pix.height
        samp = pix.samples
        # Stride: 3 (RGB) or 4 (RGBA).  Default get_pixmap is RGB w/o alpha.
        stride = pix.n  # 3 or 4
        white = 0
        for k in range(0, len(samp), stride):
            if samp[k] >= 240 and samp[k+1] >= 240 and samp[k+2] >= 240:
                white += 1
        frac = white / n if n else 1.0
        rec = {"page": i + 1, "whitespace_fraction": round(frac, 3)}
        all_pages.append(rec)
        if frac > threshold:
            flagged.append(rec)
    doc.close()

    advice = (
        f"Pages with whitespace > {threshold:.0%} usually come from rigid "
        f"float placement (e.g. [!h] or [H] forcing a figure to its own "
        f"page).  Loosen to [htbp] OR move the figure earlier in the source "
        f"OR group multiple figures via \\begin{{figure}}\\subfloat[...]"
        f"\\subfloat[...]\\end{{figure}}.  Also check that captions don't "
        f"\\newpage themselves -- multi-line captions can spill onto a "
        f"following blank page."
    )

    return {
        "total_pages": len(all_pages),
        "flagged_count": len(flagged),
        "flagged_pages": flagged,
        "all_pages": all_pages,
        "threshold": threshold,
        "advice": advice,
    }


# ============================================================
# make_layout_thumbnail_strip
# ============================================================


def paper_writing_layout_thumbnail_strip(
    pdf_path: str,
    out_path: Optional[str] = None,
    dpi: int = 80,
    cols: int = 4,
    gap_px: int = 8,
    bg_rgb: tuple = (240, 240, 240),
) -> dict:
    """Render every page as a thumbnail tile into ONE composite PNG.

    Useful for a single-glance visual scan: feed the strip image to
    Claude vision and ask "which pages look wasteful or have bad
    figure placement?"  Much faster than per-page review.

    Args:
        pdf_path: typeset PDF.
        out_path: composite PNG output.  Default: {pdf_path}_strip.png.
        dpi: per-thumbnail render resolution.  Default 80 = tiny
            (~565x800 px each at A4); enough to spot empty pages,
            bad figure placement, page-spanning floats.
        cols: thumbnails per row.  4 fits a typical paper on a screen.
        gap_px: pixels between thumbnails.
        bg_rgb: background colour as (R, G, B).

    Returns:
        dict with:
          "out_path", "total_pages", "thumbnail_w_h": (w, h)
          "strip_w_h": (W, H) of the composite
          "hint": how to use with vision agent
    """
    pymupdf = _require_pymupdf()
    Image = _require_pil()
    pdf = pathlib.Path(pdf_path)
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}
    doc = pymupdf.open(str(pdf))
    total = doc.page_count
    if total == 0:
        doc.close()
        return {"error": "PDF has 0 pages"}

    # Render every page to a PIL Image
    tiles = []
    for page in doc:
        pix = page.get_pixmap(dpi=dpi)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        tiles.append(img)
    doc.close()

    tw, th = tiles[0].size
    n_rows = (total + cols - 1) // cols
    W = cols * tw + (cols + 1) * gap_px
    H = n_rows * th + (n_rows + 1) * gap_px
    strip = Image.new("RGB", (W, H), tuple(bg_rgb))
    for i, t in enumerate(tiles):
        r = i // cols
        c = i % cols
        x = gap_px + c * (tw + gap_px)
        y = gap_px + r * (th + gap_px)
        strip.paste(t, (x, y))

    if out_path is None:
        out_path = str(pdf.with_name(pdf.stem + "_strip.png"))
    strip.save(out_path, optimize=True)

    return {
        "out_path": out_path,
        "total_pages": total,
        "thumbnail_w_h": (tw, th),
        "strip_w_h": (W, H),
        "size_bytes": os.path.getsize(out_path),
        "hint": (
            "Open the strip in Claude vision and ask: 'For each numbered "
            "thumbnail, list any page that looks abnormal (mostly white, "
            "figure pushed to wrong page, caption stranded across page "
            "break, single-column figure spanning both columns).'"
        ),
    }


# ============================================================
# check_floats_far_from_reference
# ============================================================


def paper_writing_check_floats_far_from_reference(
    tex_path: str,
    pdf_path: str,
    max_pages_apart: int = 1,
    aux_path: str = "",
) -> dict:
    """Detect figures whose \\ref{} appears far from the actual float.

    Reader's rule of thumb: the first \\ref{fig:foo} should sit on
    THE SAME or ADJACENT page as the figure body.  When a [!h] forces
    the float far away (or [t] dumps it 3 pages later), readers must
    flip back and forth.

    This tool reads the .tex to find label-ref pairs, then renders
    the PDF text per page (pymupdf get_text) to locate which page
    contains the \\ref{label} occurrence vs which page contains the
    figure environment containing \\label{label}.  Pages with
    |float_page - first_ref_page| > max_pages_apart are flagged.

    Args:
        tex_path: source .tex (single-file paper).  Multi-file:
            concatenate first, OR call once per file (the tool
            scans label/ref textually).
        pdf_path: typeset PDF of the same document.
        max_pages_apart: tolerance in pages.  Default 1 (figure must
            be within 1 page of its first reference).

    Returns:
        dict with:
          "flagged": list of {"label", "ref_pages", "float_page", "delta_pages"}
          "labels_found": total label count
          "advice": text hint
    """
    pymupdf = _require_pymupdf()
    tex = pathlib.Path(tex_path)
    pdf = pathlib.Path(pdf_path)
    if not tex.exists():
        return {"error": f"tex_path not found: {tex_path}"}
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}

    # Quick TeX scan for figure-block label and its corresponding refs.
    import re as _re
    src = tex.read_text(encoding="utf-8", errors="replace")

    # Find all \begin{figure}...\end{figure} blocks and capture their
    # \label{...} (if any).
    fig_re = _re.compile(
        r"\\begin\{figure[*]?\}.*?\\end\{figure[*]?\}",
        _re.DOTALL,
    )
    label_re = _re.compile(r"\\label\{([^}]+)\}")
    fig_labels = []
    for m in fig_re.finditer(src):
        block = m.group(0)
        lm = label_re.search(block)
        if lm:
            fig_labels.append(lm.group(1))

    aux = pathlib.Path(aux_path) if aux_path else tex.with_suffix(".aux")
    if not aux.is_file():
        return {
            "applicable": False,
            "status": "not_applicable",
            "labels_found": len(fig_labels),
            "flagged_count": 0,
            "flagged": [],
            "reason": f"compiled aux file not found: {aux}",
            "advice": "Compile the TeX document, preserve its .aux file, and rerun.",
        }

    aux_text = aux.read_text(encoding="utf-8", errors="strict")
    aux_labels: dict[str, tuple[str, int]] = {}
    for match in _re.finditer(
        r"\\newlabel\{([^{}]+)\}\{\{([^{}]+)\}\{([^{}]+)\}", aux_text
    ):
        label, number, page_token = match.groups()
        if label in fig_labels and page_token.isdigit():
            aux_labels[label] = (number, int(page_token))

    # Build per-page text from the PDF
    with pymupdf.open(str(pdf)) as doc:
        page_texts = [page.get_text() for page in doc]

    flagged = []
    unresolved_labels = []
    for lab in fig_labels:
        if lab not in aux_labels:
            unresolved_labels.append(lab)
            continue
        figure_number, float_page = aux_labels[lab]
        # Exclude caption-shaped lines so the caption does not count as its
        # own prose reference. The authoritative float page comes from aux.
        mention = _re.compile(
            rf"\bFig(?:ure)?s?\.?[~\s]*{_re.escape(figure_number)}\b",
            _re.IGNORECASE,
        )
        caption = _re.compile(
            rf"^\s*Fig(?:ure)?\.?[~\s]*{_re.escape(figure_number)}\s*[.:-]",
            _re.IGNORECASE,
        )
        ref_pages = []
        for page_no, page_text in enumerate(page_texts, start=1):
            if any(
                mention.search(line) and not caption.search(line)
                for line in page_text.splitlines()
            ):
                ref_pages.append(page_no)
        if not ref_pages:
            continue
        first_ref = min(ref_pages)
        delta = abs(float_page - first_ref)
        if delta > max_pages_apart:
            flagged.append({
                "label": lab,
                "figure_number": figure_number,
                "ref_pages": ref_pages,
                "first_ref_page": first_ref,
                "float_page": float_page,
                "delta_pages": delta,
            })

    advice = (
        "First \\ref{fig:X} should be within 1 page of the figure itself. "
        "When delta_pages > 1, the float was likely pushed far by a rigid "
        "[!h] or [t] specifier.  Fixes: (a) move \\begin{figure} earlier in "
        "source (LaTeX never pulls a float backward); (b) loosen to [htbp] "
        "or [tbp]; (c) add \\usepackage[section]{placeins} so floats cannot "
        "drift past section boundaries.  This heuristic is approximate; for "
        "definitive answers compile with -synctex=1 and use the Skillnad "
        "visual-diff approach.  See bem_sommerfeld_layered for an example "
        "of how the lab uses companion .vol.json metadata for similar "
        "cross-reference verification."
    )

    return {
        "applicable": True,
        "status": "warning" if flagged or unresolved_labels else "clean",
        "labels_found": len(fig_labels),
        "labels_resolved_from_aux": len(aux_labels),
        "unresolved_labels": unresolved_labels,
        "aux_path": str(aux),
        "flagged_count": len(flagged),
        "flagged": flagged,
        "advice": advice,
    }


# ============================================================
# layout_visual_recipe knowledge string
# ============================================================


def paper_writing_layout_visual_recipe() -> str:
    """Return the lab recipe for image-based PDF layout verification.

    Three-step recipe combining the algorithmic checks above with
    Claude-vision inspection:

      1. paper_writing_layout_thumbnail_strip(pdf, out='strip.png')
         -> one composite PNG for a glance over the whole document.
      2. paper_writing_detect_page_whitespace_anomalies(pdf)
         -> algorithmically flag mostly-white pages.
      3. paper_writing_render_pages_to_png(pdf, page_range=...)
         -> per-page PNGs for vision-agent / human deep inspection.
    """
    return _LAYOUT_VISUAL_RECIPE


_LAYOUT_VISUAL_RECIPE = r"""
# Lab recipe: image-based PDF layout verification (radia-mcp v0.88.0)

The text-based ``paper_writing_check_pdf_*`` tools catch text
overflow and obvious font / metadata problems, but they MISS
semantic layout failures that only become visible when you LOOK at
the rendered page.  This recipe combines algorithmic + visual checks.

## Step 1: Glance the whole document via thumbnail strip

```python
from radia_mcp.paper_writing import (
    paper_writing_layout_thumbnail_strip,
)
r = paper_writing_layout_thumbnail_strip(
    pdf_path='draft_v3.pdf',
    out_path='draft_v3_strip.png',
    dpi=80,
    cols=4,
)
# Open r['out_path'] in any image viewer (or feed to Claude vision).
```

What to look for at a glance:
  * **Mostly-white pages** -- bad float placement.
  * **Pages with one figure surrounded by white** -- [H] or [!h]
    overruled the running text.
  * **Caption stranded on its own page** -- caption is too long for
    \\begin{figure} block; consider splitting.
  * **Section header at bottom of page** with no body text following
    -- need \\nopagebreak or \\samepage around the header.

## Step 2: Algorithmic whitespace flag

```python
from radia_mcp.paper_writing import (
    paper_writing_detect_page_whitespace_anomalies,
)
r = paper_writing_detect_page_whitespace_anomalies(
    pdf_path='draft_v3.pdf',
    whitespace_threshold=0.75,
)
for p in r['flagged_pages']:
    print(f"  page {p['page']}: {p['whitespace_fraction']:.1%} white")
```

Threshold 0.75 catches typical bad-float pages without false-
positive on title pages.  Raise to 0.85 if your document has many
section-opening pages with intentional space.

## Step 3: Per-page deep inspection via Claude vision

```python
from radia_mcp.paper_writing import paper_writing_render_pages_to_png
r = paper_writing_render_pages_to_png(
    pdf_path='draft_v3.pdf',
    out_dir='draft_v3_pages',
    dpi=150,
    page_range='3-7',          # focus on flagged pages
)
# Then ask Claude (this MCP host) something like:
#   "Open draft_v3_pages/page_003.png through page_007.png and report
#    any figure that is too small to read, caption that doesn't match
#    the figure, axis labels that are clipped, or color choices that
#    fail the colorblind-safe palette."
```

## Step 4: Cross-reference proximity

```python
from radia_mcp.paper_writing import paper_writing_check_floats_far_from_reference
r = paper_writing_check_floats_far_from_reference(
    tex_path='paper.tex',
    pdf_path='draft_v3.pdf',
    max_pages_apart=1,
)
for f in r['flagged']:
    print(f"  {f['label']}: refs span pages {f['first_ref_page']}-"
          f"{f['last_ref_page']} ({f['span_pages']} pages apart)")
```

If a figure's references span more than 1 page, the figure body is
probably stranded -- check the float specifier (see
``paper_writing_tex_figure_placement`` for the canonical advice).

## When to use each step

  * Step 1 (thumbnail strip): ALWAYS, as the first check after any
    new build.
  * Step 2 (whitespace algorithm): when the strip looks suspicious
    but you want numbers to back the gut feel.
  * Step 3 (per-page PNG + vision): when step 1/2 flagged
    something -- focus on the flagged pages, do not bulk-render.
  * Step 4 (ref proximity): pre-submission sanity check; reviewers
    HATE flipping back and forth.

## Web-research basis (2026-05-26)

  * pymupdf get_pixmap(dpi=N) -- Artifex PyMuPDF docs.
  * hgustafsson/skillnad on GitHub -- color-coded visual PDF diffs
    via SyncTeX (inspired the thumbnail-strip pattern).
  * Overleaf "Understanding underfull and overfull box warnings" --
    why the text-only checks are insufficient.
"""
