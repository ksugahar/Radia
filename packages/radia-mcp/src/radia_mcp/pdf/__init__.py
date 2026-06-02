"""pdf sub-module — generic PDF preflight, inspect, edit, audit.

Generalization of poster-specific PDF tools (poster_font_embed_check
etc.) into a domain-neutral toolkit. Used at journal submission time,
print preparation, and supplementary-material packaging.

Tools (prefix ``pdf_``):

    Tier 1 — Inspect / convert:
        pdf_inspect(pdf_path, verbose=False)
        pdf_font_embed_check(pdf_path)
        pdf_extract_text(pdf_path, mode="plain"|"layout"|"dict", page_range=None)

    Tier 2 — Surgery:
        pdf_merge(pdf_paths, out_path)
        pdf_split(pdf_path, out_dir=None)
        pdf_extract_pages(pdf_path, page_range, out_path)
        pdf_crop_whitespace(pdf_path, out_path=None, margin_pt=5)
        pdf_compress(pdf_path, out_path=None, image_dpi=150, image_quality=75)
        pdf_watermark_add(pdf_path, text="DRAFT", ...)

    Tier 3 — Metadata / outline:
        pdf_get_metadata(pdf_path)
        pdf_set_metadata(pdf_path, title=, author=, subject=, keywords=, out_path=None)
        pdf_list_bookmarks(pdf_path)
        pdf_set_bookmarks(pdf_path, outline, out_path=None)

    Tier 4 — Intelligence:
        pdf_health_report(pdf_path, max_size_mb=10.0)
"""
from . import tools as _tools


def register(mcp) -> int:
    """Register pdf_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("pdf_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
