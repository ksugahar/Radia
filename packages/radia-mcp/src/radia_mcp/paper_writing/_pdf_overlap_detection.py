"""Pixel-accurate PDF layout overlap/overflow detection.

Added 2026-05-26 (radia-mcp v0.92.0).  Complements the
``_pdf_layout_visual.py`` module: where layout_visual renders pages
to PNG for vision-agent / human inspection, THIS module is the
algorithmic detector that walks every text block + every image
bounding box on every page and flags:

  * **text-on-image overlap** -- caption straddling a figure,
    paragraph text running through a chart background, page number
    on top of a header logo.
  * **text-overflow off page** -- bbox extending past CropBox or
    MediaBox; common at margins for long equations, code blocks,
    or wide tables.
  * **text-on-text overlap** -- two text blocks visually
    overlapping (caption text vs body text, columns colliding in
    two-column layout, multi-line caption straddling page break).

Web-research basis (2026-05-26 GitHub + arxiv survey):

  * **PyMuPDF** ``Rect.intersects()`` / ``Rect.intersect()`` --
    primitive bbox math.
  * **VILA** (Visual Layout Groups, Shen-Lo-Wang 2021,
    arXiv:2106.00676) -- vision-based structure extraction.
  * **pdfminer.six** char_margin / line_overlap -- grouping.
  * **PaperOrchestra** (arXiv:2604.05018) -- multi-agent paper-writing
    framework that flags visual-defect issues separately from
    content issues.
  * **PaperDebugger** (arXiv:2512.02589) -- plugin in-editor
    multi-agent paper checker.

The overlap detection here is OFFLINE: requires only the typeset
PDF.  No LaTeX source or compile log needed.  Pairs well with the
existing image-based vision recipe (when overlap is detected on
page N, render page N to PNG and let the user / vision agent
verify visually).
"""

from __future__ import annotations

import os
import pathlib
from typing import Optional


def _require_pymupdf():
    try:
        import pymupdf  # type: ignore
        return pymupdf
    except ImportError as e:
        raise ImportError(
            "paper_writing PDF overlap detection requires pymupdf: "
            "pip install pymupdf  (or pip install radia-mcp[rag])."
        ) from e


def _iou(a: tuple, b: tuple) -> float:
    """Intersection-over-Union of two (x0, y0, x1, y1) rectangles.

    Returns 0.0 if no overlap or either rect has zero area.
    """
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    ix0, iy0 = max(ax0, bx0), max(ay0, by0)
    ix1, iy1 = min(ax1, bx1), min(ay1, by1)
    if ix1 <= ix0 or iy1 <= iy0:
        return 0.0
    inter = (ix1 - ix0) * (iy1 - iy0)
    a_area = max(0.0, (ax1 - ax0)) * max(0.0, (ay1 - ay0))
    b_area = max(0.0, (bx1 - bx0)) * max(0.0, (by1 - by0))
    union = a_area + b_area - inter
    if union <= 0:
        return 0.0
    return inter / union


def _area(rect: tuple) -> float:
    x0, y0, x1, y1 = rect
    return max(0.0, (x1 - x0)) * max(0.0, (y1 - y0))


def _intersection_area(a: tuple, b: tuple) -> float:
    """Area of intersection of two rects in pt^2."""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    ix0, iy0 = max(ax0, bx0), max(ay0, by0)
    ix1, iy1 = min(ax1, bx1), min(ay1, by1)
    if ix1 <= ix0 or iy1 <= iy0:
        return 0.0
    return (ix1 - ix0) * (iy1 - iy0)


# ============================================================
# paper_writing_detect_text_image_overlap
# ============================================================


def paper_writing_detect_text_image_overlap(
    pdf_path: str,
    iou_threshold: float = 0.02,
    min_intersection_area_pt2: float = 50.0,
    ignore_full_page_images: bool = True,
    full_page_image_area_ratio: float = 0.9,
) -> dict:
    """Detect text blocks overlapping with images on every page.

    Reviewer-visible defect: a caption / page number / paragraph
    runs across a figure background (often happens when an
    ``\\includegraphics[width=\\textwidth]`` extends past the
    column edge in two-column layout, or when the clip+viewport
    trap leaves the full image overlapping body text -- see
    https://en.wikibooks.org/wiki/LaTeX/Floats,_Figures_and_Captions ).

    Args:
        pdf_path: typeset PDF.
        iou_threshold: minimum IoU to flag as overlap.  Default
            0.02 = 2% of the union (catches even small caption-image
            overlaps).
        min_intersection_area_pt2: minimum intersection area in
            pt^2 to flag.  Default 50 pt^2 ~= 1mm^2 (catches even
            tiny visible overlaps).  Tighten to 200+ to suppress
            anti-aliasing-level false positives.
        ignore_full_page_images: True (default) -- if an image
            covers > 90% of the page (likely a journal-template
            background or watermark), don't flag text on top.
        full_page_image_area_ratio: threshold for "full-page image".

    Returns:
        dict with:
          "total_pages": int
          "n_overlaps": int (total flagged across all pages)
          "overlaps": list of {
              "page": int (1-indexed),
              "text_snippet": first ~40 chars of overlapping text,
              "text_bbox": (x0, y0, x1, y1) in points,
              "image_bbox": (x0, y0, x1, y1),
              "iou": float,
              "intersection_area_pt2": float,
          }
          "advice": text hint for fixing
    """
    pymupdf = _require_pymupdf()
    pdf = pathlib.Path(pdf_path)
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}
    doc = pymupdf.open(str(pdf))

    overlaps = []
    for page_idx, page in enumerate(doc):
        page_area = abs(page.rect.width * page.rect.height)
        # Collect image bboxes via get_image_info (more reliable than
        # get_images + per-name lookup).
        try:
            img_infos = page.get_image_info(xrefs=True)
        except Exception:  # noqa: BLE001
            img_infos = []
        # Filter out full-page images
        img_bboxes = []
        for info in img_infos:
            bbox = info.get("bbox")
            if not bbox or len(bbox) != 4:
                continue
            if ignore_full_page_images:
                if _area(bbox) / page_area > full_page_image_area_ratio:
                    continue
            img_bboxes.append(tuple(bbox))

        if not img_bboxes:
            continue

        # Walk text blocks
        text_dict = page.get_text("dict")
        for block in text_dict.get("blocks", []):
            if block.get("type") != 0:  # 0 = text, 1 = image
                continue
            bbox = block.get("bbox")
            if not bbox or len(bbox) != 4:
                continue
            text_bbox = tuple(bbox)
            # Extract a tiny snippet for the report
            snippet_chars = []
            for line in block.get("lines", []):
                for span in line.get("spans", []):
                    snippet_chars.append(span.get("text", ""))
            snippet = " ".join(snippet_chars).strip()
            snippet = (snippet[:40] + "...") if len(snippet) > 43 else snippet

            for ib in img_bboxes:
                iou = _iou(text_bbox, ib)
                inter_area = _intersection_area(text_bbox, ib)
                if (iou >= iou_threshold
                        and inter_area >= min_intersection_area_pt2):
                    overlaps.append({
                        "page": page_idx + 1,
                        "text_snippet": snippet,
                        "text_bbox": text_bbox,
                        "image_bbox": ib,
                        "iou": round(iou, 4),
                        "intersection_area_pt2": round(inter_area, 2),
                    })
    _n_pages = doc.page_count
    doc.close()

    advice = (
        "Text-on-image overlap = a caption / page number / paragraph "
        "is running across an image.  Likely causes: (a) \\includegraphics"
        " with clip+viewport that doesn't actually crop (the full image "
        "reappears under text), (b) figure placed with [H] forcing "
        "absolute position into a paragraph, (c) header/footer "
        "interfering with a full-width image, (d) two-column "
        "\\includegraphics width=\\textwidth in figure (NOT figure*) -- "
        "the image spans both columns but the figure environment is "
        "single-column, so text runs through the image's right half.  "
        "Inspect flagged pages via paper_writing_render_pages_to_png "
        "to confirm visually."
    )
    return {
        "total_pages": _n_pages,
        "n_overlaps": len(overlaps),
        "overlaps": overlaps,
        "iou_threshold": iou_threshold,
        "min_intersection_area_pt2": min_intersection_area_pt2,
        "advice": advice,
    }


# ============================================================
# paper_writing_detect_text_overflow_page
# ============================================================


def paper_writing_detect_text_overflow_page(
    pdf_path: str,
    use_cropbox: bool = True,
    overflow_tolerance_pt: float = 1.0,
    ignore_blocks_smaller_than_pt2: float = 25.0,
) -> dict:
    """Detect text blocks extending past the page CropBox / MediaBox.

    LaTeX's Overfull \\hbox warning catches text overflowing the
    TEXT WIDTH, but it doesn't catch:
      * tables wider than the column,
      * \\verbatim long lines spilling off the page,
      * \\includegraphics that extends past the column edge,
      * footers / watermarks past the trim box.
    This tool inspects the typeset PDF directly so the above
    cases are visible.

    Args:
        pdf_path: typeset PDF.
        use_cropbox: True (default) compares to CropBox (the visible
            area); False compares to MediaBox (paper edge).  Most
            journals use CropBox = trim area.
        overflow_tolerance_pt: ignore overflows smaller than this
            in points.  Default 1.0 pt ~ 0.35 mm (catches all
            visible overflows; raise to 3.0 for laxer check).
        ignore_blocks_smaller_than_pt2: skip tiny blocks (header
            decorations, page numbers).

    Returns:
        dict with:
          "total_pages": int
          "n_overflows": int
          "overflows": list of {
              "page": int (1-indexed),
              "side": "left" | "right" | "top" | "bottom",
              "overflow_pt": float,
              "text_snippet": str,
              "text_bbox": (x0, y0, x1, y1),
              "page_box": (x0, y0, x1, y1),
              "box_type": "normalized_cropbox" | "normalized_mediabox",
          }
          "advice": text hint
    """
    pymupdf = _require_pymupdf()
    pdf = pathlib.Path(pdf_path)
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}
    overflows = []
    with pymupdf.open(str(pdf)) as doc:
        for page_idx, page in enumerate(doc):
            # Text bboxes are in MuPDF's normalized page coordinate system.
            # Raw CropBox / MediaBox rectangles can retain a non-zero PDF
            # origin, so compare only after applying the page transform.
            if use_cropbox:
                pbox = page.rect
                box_type = "normalized_cropbox"
            else:
                pbox = page.mediabox * page.transformation_matrix
                box_type = "normalized_mediabox"
            page_box = (pbox.x0, pbox.y0, pbox.x1, pbox.y1)

            text_dict = page.get_text("dict")
            for block in text_dict.get("blocks", []):
                if block.get("type") != 0:
                    continue
                bbox = block.get("bbox")
                if not bbox or len(bbox) != 4:
                    continue
                if _area(bbox) < ignore_blocks_smaller_than_pt2:
                    continue
                x0, y0, x1, y1 = bbox

                # Snippet for report
                snippet_chars = []
                for line in block.get("lines", []):
                    for span in line.get("spans", []):
                        snippet_chars.append(span.get("text", ""))
                snippet = " ".join(snippet_chars).strip()
                snippet = (snippet[:40] + "...") if len(snippet) > 43 else snippet

                # Check each side
                sides = [
                    ("left",   page_box[0] - x0),
                    ("right",  x1 - page_box[2]),
                    ("top",    page_box[1] - y0),
                    ("bottom", y1 - page_box[3]),
                ]
                for side, overflow in sides:
                    if overflow > overflow_tolerance_pt:
                        overflows.append({
                            "page": page_idx + 1,
                            "side": side,
                            "overflow_pt": round(overflow, 2),
                            "text_snippet": snippet,
                            "text_bbox": tuple(bbox),
                            "page_box": page_box,
                            "box_type": box_type,
                        })
        _n_pages = doc.page_count

    advice = (
        "Text overflow past the page box = visible at print scale "
        "and almost always flagged by reviewers.  Common causes: "
        "(a) wide \\verbatim or code block (use \\fvset{breaklines} "
        "or listings); (b) wide table (split into two tables, or use "
        "\\resizebox{\\columnwidth}{!}{...}); (c) full-width "
        "\\includegraphics in a single-column figure (move to "
        "figure* or shrink width=\\columnwidth); (d) long inline "
        "equation -- consider displaying it.  This check complements "
        "paper_writing_check_overfull_hbox (which reads the .log)."
    )
    return {
        "total_pages": _n_pages,
        "n_overflows": len(overflows),
        "overflows": overflows,
        "use_cropbox": use_cropbox,
        "overflow_tolerance_pt": overflow_tolerance_pt,
        "advice": advice,
    }


# ============================================================
# paper_writing_detect_overlapping_text_blocks
# ============================================================


def paper_writing_detect_overlapping_text_blocks(
    pdf_path: str,
    iou_threshold: float = 0.05,
    min_intersection_area_pt2: float = 25.0,
) -> dict:
    """Detect text-on-text overlap (block-vs-block IoU on every page).

    Reviewer-visible defect: caption text running into body text,
    columns colliding in two-column layout, multi-line caption
    straddling a page break, equation overlapping with subscript
    text from another line.

    Args:
        pdf_path: typeset PDF.
        iou_threshold: minimum IoU to flag.  Default 0.05.
        min_intersection_area_pt2: minimum intersection area to flag.

    Returns:
        dict with:
          "total_pages", "n_overlaps", "overlaps", "advice"
    """
    pymupdf = _require_pymupdf()
    pdf = pathlib.Path(pdf_path)
    if not pdf.exists():
        return {"error": f"pdf_path not found: {pdf_path}"}
    doc = pymupdf.open(str(pdf))

    overlaps = []
    for page_idx, page in enumerate(doc):
        text_dict = page.get_text("dict")
        text_blocks = []
        for block in text_dict.get("blocks", []):
            if block.get("type") != 0:
                continue
            bbox = block.get("bbox")
            if not bbox or len(bbox) != 4:
                continue
            snippet_chars = []
            for line in block.get("lines", []):
                for span in line.get("spans", []):
                    snippet_chars.append(span.get("text", ""))
            snippet = " ".join(snippet_chars).strip()
            snippet = (snippet[:40] + "...") if len(snippet) > 43 else snippet
            letters = sum(char.isalpha() for char in snippet)
            if letters / max(len(snippet), 1) < 0.4:
                # Display equations are often emitted as overlapping glyph
                # fragments; they are not prose-on-prose collisions.
                continue
            text_blocks.append((tuple(bbox), snippet))

        n = len(text_blocks)
        for i in range(n):
            for j in range(i + 1, n):
                a, snip_a = text_blocks[i]
                b, snip_b = text_blocks[j]
                iou = _iou(a, b)
                inter = _intersection_area(a, b)
                if (iou >= iou_threshold
                        and inter >= min_intersection_area_pt2):
                    overlaps.append({
                        "page": page_idx + 1,
                        "block_a_snippet": snip_a,
                        "block_a_bbox": a,
                        "block_b_snippet": snip_b,
                        "block_b_bbox": b,
                        "iou": round(iou, 4),
                        "intersection_area_pt2": round(inter, 2),
                    })
    _n_pages = doc.page_count
    doc.close()

    advice = (
        "Text-on-text overlap suggests (a) caption straddling a page "
        "break (split with \\caption*{} on second part, or move the "
        "figure earlier); (b) two-column collision (a stray "
        "\\onecolumn followed by \\twocolumn can desync); (c) "
        "negative \\vspace pushing one paragraph onto another."
    )
    return {
        "total_pages": _n_pages,
        "n_overlaps": len(overlaps),
        "overlaps": overlaps,
        "iou_threshold": iou_threshold,
        "min_intersection_area_pt2": min_intersection_area_pt2,
        "advice": advice,
    }


# ============================================================
# Combined recipe knowledge
# ============================================================


_OVERLAP_RECIPE = r"""
# PDF overlap/overflow detection recipe (radia_mcp.paper_writing v0.92.0)

Three offline pymupdf-based checks complement the existing
image-render-and-look approach:

  1. paper_writing_detect_text_image_overlap(pdf_path)
        - Walks every text block + image bbox on every page.
        - Flags IoU >= 0.02 with intersection > 50 pt^2.
        - Catches: caption-on-figure, page-number-on-header-logo,
          full-image-reappearing-under-text (clip+viewport trap).
        - Ignores full-page images (journal templates) by default.

  2. paper_writing_detect_text_overflow_page(pdf_path)
        - Walks every text block, compares to CropBox (or MediaBox).
        - Flags overflow > 1.0 pt past any side.
        - Catches: wide \\verbatim, wide table, full-width
          \\includegraphics in single-column figure, long inline
          equation.
        - Complements paper_writing_check_overfull_hbox (which
          reads the .log) -- this catches OVERFLOW that LaTeX
          itself didn't warn about (typically image-related).

  3. paper_writing_detect_overlapping_text_blocks(pdf_path)
        - Walks every pair of text blocks on every page.
        - Flags IoU >= 0.05 with intersection > 25 pt^2.
        - Catches: caption straddling page break, two-column
          collision, negative \\vspace pushing paragraphs together.

## Decision tree

  Reviewer complaint about layout
      -> render thumbnail strip (paper_writing_layout_thumbnail_strip)
      -> visually identify flagged page
      -> drill in: which of the 3 overlap types?
            caption / page no / paragraph runs through image
                -> detect_text_image_overlap
            text / equation past page edge
                -> detect_text_overflow_page
            two paragraphs / caption blocks visually merged
                -> detect_overlapping_text_blocks
      -> fix root cause per advice

## Web research basis (2026-05-26)

  * arXiv:2106.00676 (VILA, Shen-Lo-Wang 2021) -- visual layout
    groups for structure extraction.
  * arXiv:2604.05018 (PaperOrchestra) -- multi-agent paper-writing.
  * arXiv:2512.02589 (PaperDebugger) -- plugin in-editor multi-
    agent paper checker.
  * pdfminer.six char_margin / line_overlap -- grouping algorithm.
  * pymupdf Rect.intersects() / Rect.intersect() -- primitive bbox
    math used by all three tools above.
"""


def paper_writing_pdf_overlap_recipe() -> str:
    """Return the recipe for PDF overlap/overflow detection."""
    return _OVERLAP_RECIPE
