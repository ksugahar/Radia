"""Tier 1 — doc_convert_pdf_to_jpg.

Convert PDF to JPG via Adobe Acrobat COM automation. Mirrors the pattern
established by the seed script ``public-safe curated corpus``:

    avDoc = win32com.client.Dispatch("AcroExch.AVDoc")
    avDoc.Open(pdf_path, "")
    pdDoc = avDoc.GetPDDoc()
    jsObject = pdDoc.GetJSObject()
    jsObject.SaveAs(out_path, "com.adobe.acrobat.jpeg")

Multi-page PDFs: Acrobat splits the output into one JPG per page,
naming them ``<stem>_ページ_1.jpg``, ``<stem>_ページ_2.jpg``, ... (the
suffix depends on Acrobat's UI language). The function returns the
files Acrobat actually produced rather than guessing the names.

Companion ``doc_convert_pdf_to_jpg_dir`` converts every PDF in a
directory.

Requirements:
    - Windows
    - Adobe Acrobat (Pro or Standard — the JS SaveAs API is not exposed
      by Acrobat Reader)
    - pywin32

This is the only Acrobat-COM tool in mcp-server-document. PowerPoint
tools in this same sub-module use ``PowerPoint.Application`` COM.
"""
from __future__ import annotations

import glob
import os
import pathlib
import sys
import time


_ACROBAT_JPEG_FORMAT_ID = "com.adobe.acrobat.jpeg"


def _load_adobe_com():
    """Import the Adobe COM dispatch helper or return None.

    Returns the ``win32com.client`` module on success. Acrobat must be
    installed; the Dispatch call below is what actually fails when the
    AcroExch.AVDoc CLSID is missing (e.g. Reader-only install).
    """
    if sys.platform != "win32":
        return None
    try:
        import win32com.client  # type: ignore
        import winerror  # type: ignore
        from win32com.client.dynamic import ERRORS_BAD_CONTEXT  # type: ignore
        if winerror.E_NOTIMPL not in ERRORS_BAD_CONTEXT:
            ERRORS_BAD_CONTEXT.append(winerror.E_NOTIMPL)
        return win32com.client
    except ImportError:
        return None


def _convert_one(com, pdf_path: pathlib.Path,
                 out_dir: pathlib.Path) -> tuple[bool, str, list[pathlib.Path]]:
    """Open ``pdf_path`` in Acrobat AVDoc, SaveAs JPEG into ``out_dir``.

    Returns ``(ok, message, produced_files)``. Acrobat owns the JPG
    filenames; we discover them by snapshotting ``out_dir`` before and
    after the SaveAs call.
    """
    target = out_dir / f"{pdf_path.stem}.jpg"
    before = set(out_dir.glob(f"{pdf_path.stem}*.jpg"))

    av_doc = com.Dispatch("AcroExch.AVDoc")
    try:
        if not av_doc.Open(str(pdf_path), ""):
            return False, f"AVDoc.Open failed for {pdf_path}", []
        pd_doc = av_doc.GetPDDoc()
        js = pd_doc.GetJSObject()
        n_pages = pd_doc.GetNumPages()
        try:
            js.SaveAs(str(target), _ACROBAT_JPEG_FORMAT_ID)
        except Exception as e:  # noqa: BLE001
            return False, f"SaveAs failed: {type(e).__name__}: {e}", []
        finally:
            try:
                pd_doc.Close()
            except Exception:
                pass
    finally:
        try:
            av_doc.Close(1)
        except Exception:
            pass

    # Give Acrobat a beat to flush, then look for what it wrote.
    time.sleep(0.2)
    after = set(out_dir.glob(f"{pdf_path.stem}*.jpg"))
    produced = sorted(after - before)
    if not produced:
        # SaveAs may overwrite an existing single-file output without
        # changing the set; fall back to the canonical target.
        if target.exists():
            produced = [target]
        else:
            return False, ("SaveAs reported success but no JPG appeared "
                           f"matching {pdf_path.stem}*.jpg in {out_dir}"), []
    pages_msg = f"{n_pages} page(s) -> {len(produced)} JPG(s)"
    return True, pages_msg, produced


def doc_convert_pdf_to_jpg(pdf_path: str,
                           output_dir: str | None = None) -> str:
    """Convert a single PDF to JPG via Adobe Acrobat COM.

    Args:
        pdf_path: Path to the source .pdf.
        output_dir: Optional output directory. Defaults to the same
                    directory as the input.

    Returns:
        Status string listing every JPG Acrobat produced (one per page
        for multi-page PDFs).
    """
    com = _load_adobe_com()
    if com is None:
        return ("Error: Adobe Acrobat COM unavailable. doc_convert_pdf_to_jpg "
                "requires Windows + pywin32 + an installed Acrobat (Pro/Standard).")
    src = pathlib.Path(pdf_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    if src.suffix.lower() != ".pdf":
        return f"Error: not a PDF file: {src}"

    out_dir = pathlib.Path(output_dir) if output_dir else src.parent
    if not out_dir.is_absolute():
        out_dir = pathlib.Path.cwd() / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    ok, msg, produced = _convert_one(com, src, out_dir)
    if not ok:
        return f"doc_convert_pdf_to_jpg: {src} -> FAIL — {msg}"
    lines = [f"doc_convert_pdf_to_jpg: {src} -> {msg}"]
    for p in produced:
        size = p.stat().st_size if p.exists() else -1
        lines.append(f"  - {p}  ({size} bytes)")
    return "\n".join(lines)


def doc_convert_pdf_to_jpg_dir(directory: str,
                               recursive: bool = False,
                               output_dir: str | None = None) -> str:
    """Convert every .pdf in ``directory`` to JPG via Adobe Acrobat COM.

    Args:
        directory: Folder to scan.
        recursive: If True, descend into subdirectories.
        output_dir: Optional unified output directory. Defaults to each
                    PDF's own directory.

    Returns:
        Per-file status + summary count.
    """
    com = _load_adobe_com()
    if com is None:
        return ("Error: Adobe Acrobat COM unavailable. doc_convert_pdf_to_jpg_dir "
                "requires Windows + pywin32 + an installed Acrobat (Pro/Standard).")
    root = pathlib.Path(directory)
    if not root.is_absolute():
        root = pathlib.Path.cwd() / root
    if not root.is_dir():
        return f"Error: not a directory: {root}"

    pattern = "**/*.pdf" if recursive else "*.pdf"
    # Case-insensitive on Windows but explicit pattern keeps macOS / Linux
    # behaviour predictable when the same .vol pipeline runs outside Cubit.
    targets = sorted(p for p in root.glob(pattern) if p.is_file())
    if not targets:
        return f"doc_convert_pdf_to_jpg_dir: no PDFs in {root}"

    out_root = pathlib.Path(output_dir) if output_dir else None
    if out_root and not out_root.is_absolute():
        out_root = pathlib.Path.cwd() / out_root
    if out_root:
        out_root.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    ok_count = 0
    total_jpgs = 0
    for src in targets:
        dst_dir = out_root or src.parent
        ok, msg, produced = _convert_one(com, src, dst_dir)
        if ok:
            ok_count += 1
            total_jpgs += len(produced)
            lines.append(f"  [OK] {src.name}: {msg}")
            for p in produced:
                lines.append(f"      - {p.name}")
        else:
            lines.append(f"  [FAIL] {src.name}: {msg}")

    head = (f"doc_convert_pdf_to_jpg_dir: {ok_count}/{len(targets)} converted "
            f"({total_jpgs} JPGs) in {root}")
    return "\n".join([head, *lines])
