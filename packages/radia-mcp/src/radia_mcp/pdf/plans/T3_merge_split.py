"""Tier 1 — pdf_merge / pdf_split / pdf_extract_pages.

Mechanical PDF surgery: combine multiple PDFs, split a PDF into
single-page files, extract a range. All via pymupdf.
"""
from __future__ import annotations

import pathlib


def _load_pymupdf():
    try:
        import fitz  # type: ignore
        return fitz
    except ImportError:
        return None


def pdf_merge(pdf_paths: list[str], out_path: str) -> str:
    """Merge multiple PDFs into one.

    Args:
        pdf_paths: Ordered list of PDFs to merge.
        out_path: Destination path. Parent dir is created.

    Returns:
        Status string with merged size + page count.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed. Install with `pip install -e .[pdf]`."
    out = pathlib.Path(out_path)
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    dst = fitz.open()
    total_pages = 0
    for src_path in pdf_paths:
        sp = pathlib.Path(src_path)
        if not sp.is_absolute():
            sp = pathlib.Path.cwd() / sp
        if not sp.exists():
            dst.close()
            return f"Error: file not found: {sp}"
        with fitz.open(sp) as src:
            dst.insert_pdf(src)
            total_pages += len(src)
    dst.save(out, deflate=True)
    dst.close()
    size_kb = out.stat().st_size / 1024
    return f"pdf_merge: wrote {out} ({total_pages} page(s), {size_kb:.1f} KB)"


def pdf_split(pdf_path: str, out_dir: str | None = None) -> str:
    """Split a PDF into one file per page.

    Output filenames follow ``<stem>_p001.pdf``, ``<stem>_p002.pdf``...

    Args:
        pdf_path: Source PDF.
        out_dir: Output directory. Defaults to the source's directory.

    Returns:
        Status string with the count of files written.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    src_path = pathlib.Path(pdf_path)
    if not src_path.is_absolute():
        src_path = pathlib.Path.cwd() / src_path
    if not src_path.exists():
        return f"Error: file not found: {src_path}"
    out_d = pathlib.Path(out_dir) if out_dir else src_path.parent
    if not out_d.is_absolute():
        out_d = pathlib.Path.cwd() / out_d
    out_d.mkdir(parents=True, exist_ok=True)
    src = fitz.open(src_path)
    written = []
    for i in range(len(src)):
        page_doc = fitz.open()
        page_doc.insert_pdf(src, from_page=i, to_page=i)
        out_file = out_d / f"{src_path.stem}_p{i+1:03d}.pdf"
        page_doc.save(out_file)
        page_doc.close()
        written.append(out_file.name)
    src.close()
    return f"pdf_split: {src_path} → {len(written)} file(s) in {out_d}"


def pdf_extract_pages(pdf_path: str, page_range: str, out_path: str) -> str:
    """Extract a page range to a new PDF.

    Args:
        pdf_path: Source PDF.
        page_range: 1-based range, ``"3-7"`` or ``"1,3,5"`` mix.
        out_path: Destination PDF.

    Returns:
        Status string with page count + file size.
    """
    fitz = _load_pymupdf()
    if fitz is None:
        return "Error: pymupdf not installed."
    src_path = pathlib.Path(pdf_path)
    if not src_path.is_absolute():
        src_path = pathlib.Path.cwd() / src_path
    if not src_path.exists():
        return f"Error: file not found: {src_path}"
    pages = set()
    for chunk in page_range.replace(" ", "").split(","):
        if "-" in chunk:
            a, b = chunk.split("-", 1)
            try:
                pages.update(range(int(a), int(b) + 1))
            except ValueError:
                return f"Error: bad range chunk {chunk!r}"
        elif chunk:
            try:
                pages.add(int(chunk))
            except ValueError:
                return f"Error: bad page number {chunk!r}"
    out = pathlib.Path(out_path)
    if not out.is_absolute():
        out = pathlib.Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    src = fitz.open(src_path)
    dst = fitz.open()
    n_total = len(src)
    for p in sorted(pages):
        if 1 <= p <= n_total:
            dst.insert_pdf(src, from_page=p - 1, to_page=p - 1)
    dst.save(out, deflate=True)
    dst.close()
    src.close()
    return f"pdf_extract_pages: {src_path} ({page_range}) → {out} ({len(pages)} page(s))"
