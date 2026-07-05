"""doc_convert sub-module — document format conversion via Office / Adobe COM.

Originally ``ppt2pdf`` (PowerPoint → PDF only). Renamed to ``doc_convert``
2026-05-19 and extended with PDF → JPG via Adobe Acrobat COM.

Seed: public-safe curated corpus
H行列@菅原/原稿/Figures/ppt2pdf.py (16 lines, 2026).

Original behaviour: glob ``*.ppt*`` in the current directory and SaveAs(32)
each one. Promoted version adds:

    - explicit file / directory targeting (no implicit CWD scan)
    - PDF lock release (kill Acrobat / Edge / SumatraPDF holding the
      output file)
    - single-slide extraction helper (write one slide as its own .pptx)
    - bulk extraction by index range

Tools (prefix ``doc_convert_``):
    - doc_convert_pptx_to_pdf(pptx_path, output_dir=None)
        Convert a single .pptx to .pdf via PowerPoint COM.
    - doc_convert_pptx_to_pdf_dir(directory, recursive=False)
        Convert every .ppt/.pptx in a directory.
    - doc_convert_pdf_to_jpg(pdf_path, output_dir=None)
        Convert a single .pdf to .jpg via Adobe Acrobat COM
        (one JPG per page for multi-page PDFs).
    - doc_convert_pdf_to_jpg_dir(directory, recursive=False)
        Convert every .pdf in a directory.
    - doc_convert_extract_slide(src_pptx, slide_index, out_pptx)
        Save slide ``slide_index`` (1-based) of ``src_pptx`` as a
        single-slide pptx at ``out_pptx``.
    - doc_convert_extract_slides(src_pptx, mapping)
        Bulk variant: ``mapping = {slide_index: out_pptx, ...}``.
"""
from . import tools as _tools


def register(mcp) -> int:
    """Register doc_convert_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("doc_convert_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
