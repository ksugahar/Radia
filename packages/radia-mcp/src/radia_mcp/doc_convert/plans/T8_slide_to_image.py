"""Tier 2 — doc_convert_slide_to_image.

Render specific slides of a pptx as PNG / JPG / PDF via PowerPoint
COM ``Slide.Export``. Useful when you need an image-only deliverable
(e.g. journal supplementary materials) rather than a full PDF.
"""
from __future__ import annotations

import pathlib
import sys


def _load_pptx_com():
    if sys.platform != "win32":
        return None
    try:
        import win32com.client  # type: ignore
        return win32com.client
    except ImportError:
        return None


# Slide.Export format codes (PowerPoint PpFixedFormatType + image filters).
# For Slide.Export: argument is filter name "PNG", "JPG", "GIF", "BMP", "TIF",
# "WMF", "EMF". For SaveAs we use the SaveAs format codes.


def doc_convert_slide_to_image(pptx_path: str,
                           slide_index: int,
                           out_path: str,
                           fmt: str = "PNG",
                           width: int = 1920,
                           height: int = 1080) -> str:
    """Export one slide as an image file via PowerPoint COM.

    Args:
        pptx_path: Source pptx.
        slide_index: 1-based slide index.
        out_path: Destination image file (extension should match fmt).
        fmt: ``PNG`` (default), ``JPG``, ``GIF``, ``BMP``, ``TIF``,
             ``WMF``, ``EMF``.
        width: Pixel width (default 1920 — Full HD).
        height: Pixel height (default 1080).
    """
    com = _load_pptx_com()
    if com is None:
        return "Error: PowerPoint COM unavailable."
    src = pathlib.Path(pptx_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    out = pathlib.Path(out_path)
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)

    app = com.Dispatch("PowerPoint.Application")
    presentation = None
    try:
        presentation = app.Presentations.Open(str(src))
        n = presentation.Slides.Count
        if not (1 <= slide_index <= n):
            return f"Error: slide_index {slide_index} out of range [1, {n}]"
        slide = presentation.Slides(slide_index)
        slide.Export(str(out), fmt.upper(), width, height)
        return f"doc_convert_slide_to_image: wrote {out} ({fmt} {width}x{height})"
    except Exception as e:
        return f"Error: Export failed: {type(e).__name__}: {e}"
    finally:
        if presentation is not None:
            try:
                presentation.Close()
            except Exception:
                pass
        try:
            app.Quit()
        except Exception:
            pass
