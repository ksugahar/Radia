"""Tier 1 — business-card OCR + rename workflow.

Two MCP tools intended to be used together by a calling assistant
(Claude) to rename scanned-business-card JPGs to ``<company>_<name>.jpg``.

Workflow:

    1. ``doc_convert_business_card_ocr_dir("D:/")`` -- runs EasyOCR
       (ja+en) on every .jpg in the directory and returns a mapping
       ``{jpg_path: ocr_text}``.

    2. The caller (Claude) inspects the OCR text per card and decides
       what the (company, name) pair is. Layout heuristics, multilingual
       transliteration, and corporate-suffix logic (株式会社 / Co., Ltd. /
       LLC / etc.) are left to the caller because OCR text alone is too
       ambiguous for a pure regex pass.

    3. ``doc_convert_rename_business_card(jpg_path, company, name)``
       performs the actual filesystem rename with Windows-safe
       sanitization and collision-safe suffixing.

This split (extract -> human-judgement -> rename) matches the
``existing ocr_extract_text + Claude parse`` decision made on
2026-05-19.
"""
from __future__ import annotations

import pathlib
import re

# Windows-forbidden filename characters (NTFS / Win32).
_FORBIDDEN_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')

# Trailing dots and spaces are stripped by Win32 -- enforce explicitly.
_TRAIL_STRIP = re.compile(r"[\.\s]+$")


def _sanitize_segment(s: str, segment_name: str) -> tuple[str, list[str]]:
    """Make ``s`` safe to use as part of a Windows filename.

    Returns ``(cleaned, warnings)``. ``warnings`` lists transformations
    applied so the caller can surface them.
    """
    warnings: list[str] = []
    if s is None:
        return "", [f"{segment_name} was None"]
    cleaned = str(s).strip()
    if not cleaned:
        return "", [f"{segment_name} was empty after strip"]

    # Remove forbidden chars.
    sub = _FORBIDDEN_CHARS.sub("", cleaned)
    if sub != cleaned:
        warnings.append(
            f"{segment_name}: removed forbidden chars ({cleaned!r} -> {sub!r})"
        )
    cleaned = sub

    # Collapse whitespace runs.
    collapsed = re.sub(r"\s+", " ", cleaned)
    if collapsed != cleaned:
        warnings.append(f"{segment_name}: collapsed whitespace")
    cleaned = collapsed

    # Reject the underscore that we ourselves use as the separator,
    # so the resulting filename is unambiguous to split back later.
    if "_" in cleaned:
        replaced = cleaned.replace("_", "-")
        warnings.append(
            f"{segment_name}: replaced '_' with '-' to keep the "
            f"company_name separator unambiguous ({cleaned!r} -> {replaced!r})"
        )
        cleaned = replaced

    # Strip trailing dots / spaces.
    cleaned = _TRAIL_STRIP.sub("", cleaned)

    if not cleaned:
        warnings.append(f"{segment_name}: nothing left after sanitization")
    return cleaned, warnings


def _unique_path(path: pathlib.Path) -> pathlib.Path:
    """Return ``path`` if free, else ``<stem>_2.jpg``, ``_3.jpg`` ..."""
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    for i in range(2, 1000):
        candidate = parent / f"{stem}_{i}{suffix}"
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"could not find a free filename near {path}")


_IMAGE_EXTS = (".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp")
_PDF_EXTS = (".pdf",)


def doc_convert_business_card_ocr_dir(
    directory: str,
    recursive: bool = False,
    gpu: bool = False,
    include_pdfs: bool = True,
    pdf_render_dpi: int = 300,
) -> dict:
    """Run EasyOCR on every business-card file in ``directory``.

    Accepted formats: raster images (.jpg/.jpeg/.png/.tif/.tiff/.bmp)
    and .pdf (rendered to numpy arrays via pymupdf). For multi-page
    PDFs, each page is OCR'd separately and reported as its own item
    (``page`` key indicates the 1-based page index).

    Intended to be called once per batch so the assistant can read all
    cards in a single thought-pass and call
    ``doc_convert_rename_business_card`` per file.

    Args:
        directory: folder containing scanned business-card files.
        recursive: if True, descend into subdirectories.
        gpu: enable EasyOCR GPU mode.
        include_pdfs: if True (default), also OCR .pdf files. Disable to
            scan images only.
        pdf_render_dpi: DPI at which to render PDF pages before OCR.

    Returns:
        dict with::

            {
              "directory": str,
              "count":     int,                # items (1 per image, 1 per PDF page)
              "items": [
                {"file": <path>, "page": int or None,
                 "lines": [...], "text": "<joined>"},
                ...
              ]
            }

        On a hard error per-file the entry has ``{"file": ..., "error": ...}``.
        For backwards compatibility, image entries also carry ``"jpg":
        <path>`` (same as ``file``) so callers written against the
        original JPG-only version continue to work.
    """
    root = pathlib.Path(directory)
    if not root.is_absolute():
        root = pathlib.Path.cwd() / root
    if not root.is_dir():
        return {"directory": str(root), "error": f"not a directory: {root}"}

    accepted = _IMAGE_EXTS + (_PDF_EXTS if include_pdfs else ())
    glob_pattern = "**/*" if recursive else "*"
    candidates = sorted(
        p for p in root.glob(glob_pattern)
        if p.is_file() and p.suffix.lower() in accepted
    )

    if not candidates:
        return {"directory": str(root), "count": 0, "items": []}

    # Lazy import to avoid pulling EasyOCR / pymupdf on module load.
    from ...ocr import core as _ocr_core
    import numpy as np
    from PIL import Image

    reader = _ocr_core._get_reader(gpu=gpu)
    items: list[dict] = []
    for path in candidates:
        suf = path.suffix.lower()
        if suf in _IMAGE_EXTS:
            try:
                # EasyOCR delegates path opening to OpenCV imread, which
                # fails silently on Windows for non-ASCII paths (e.g.
                # names containing 'ページ'). Open via PIL and pass a
                # numpy array to sidestep this.
                with Image.open(str(path)) as im:
                    arr = np.array(im.convert("RGB"))
                raw = reader.readtext(arr, detail=1, paragraph=False)
                lines = [item[1] for item in raw]
                items.append({
                    "file": str(path),
                    "jpg": str(path),  # legacy compat
                    "page": None,
                    "lines": lines,
                    "text": "\n".join(lines),
                })
            except Exception as e:  # noqa: BLE001
                items.append({"file": str(path), "jpg": str(path),
                              "error": f"{type(e).__name__}: {e}"})
        else:  # pdf
            try:
                import pymupdf  # type: ignore
            except ImportError:
                items.append({"file": str(path),
                              "error": "pymupdf not installed; "
                                       "cannot OCR PDF business cards"})
                continue
            try:
                with pymupdf.open(str(path)) as doc:
                    for page_idx, page in enumerate(doc, start=1):
                        pix = page.get_pixmap(dpi=pdf_render_dpi)
                        # pix.samples is bytes; build a PIL image then ndarray.
                        mode = "RGB" if pix.alpha == 0 else "RGBA"
                        img = Image.frombytes(mode, (pix.width, pix.height),
                                              pix.samples)
                        arr = np.array(img.convert("RGB"))
                        raw = reader.readtext(arr, detail=1, paragraph=False)
                        lines = [item[1] for item in raw]
                        items.append({
                            "file": str(path),
                            "page": page_idx,
                            "lines": lines,
                            "text": "\n".join(lines),
                        })
            except Exception as e:  # noqa: BLE001
                items.append({"file": str(path),
                              "error": f"{type(e).__name__}: {e}"})
    return {
        "directory": str(root),
        "count": len(items),
        "items": items,
    }


def doc_convert_rename_business_card(
    jpg_path: str,
    company: str,
    name: str,
    dry_run: bool = False,
) -> dict:
    """Rename a scanned-business-card file to ``<company>_<name><ext>``.

    Accepted source extensions: .jpg / .jpeg / .png / .tif / .tiff /
    .bmp / .pdf.  The source extension is preserved (so a .pdf input
    becomes ``<company>_<name>.pdf``).

    Sanitization (Windows-safe):
        - strip leading/trailing whitespace
        - remove ``< > : " / \\ | ? *`` and control chars
        - collapse internal whitespace runs to single spaces
        - replace ``_`` inside either segment with ``-`` so the
          ``company_name`` separator stays unambiguous
        - strip trailing dots / spaces

    If the target filename already exists, a numeric suffix is appended
    (``<company>_<name>_2.jpg``, ``_3.jpg``, ...).

    Args:
        jpg_path: source JPG path (absolute or relative).
        company: company name from the business card (Claude-parsed).
        name: person's name from the business card (Claude-parsed).
        dry_run: if True, return the resolved target path without
            renaming the file.

    Returns:
        dict with::

            {
              "src":      <input path>,
              "dst":      <target path>,
              "renamed":  True / False,
              "warnings": ["sanitization note", ...]
            }

        On error returns ``{"src": ..., "error": "..."}``.
    """
    src = pathlib.Path(jpg_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return {"src": str(src), "error": f"file not found: {src}"}
    # Accept any common business-card file format. The rename operation
    # itself is format-agnostic; the source extension is preserved.
    _ACCEPTED = (".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp", ".pdf")
    if src.suffix.lower() not in _ACCEPTED:
        return {"src": str(src),
                "error": f"unsupported extension {src.suffix!r}; "
                         f"accepted: {_ACCEPTED}"}

    clean_company, w1 = _sanitize_segment(company, "company")
    clean_name, w2 = _sanitize_segment(name, "name")
    warnings = w1 + w2
    if not clean_company:
        return {"src": str(src),
                "error": "company name is empty after sanitization",
                "warnings": warnings}
    if not clean_name:
        return {"src": str(src),
                "error": "person name is empty after sanitization",
                "warnings": warnings}

    target_name = f"{clean_company}_{clean_name}{src.suffix.lower()}"
    target = src.parent / target_name
    target = _unique_path(target)

    if dry_run:
        return {
            "src": str(src),
            "dst": str(target),
            "renamed": False,
            "dry_run": True,
            "warnings": warnings,
        }

    if target == src:
        return {
            "src": str(src),
            "dst": str(target),
            "renamed": False,
            "warnings": warnings + ["target == source; no rename needed"],
        }

    src.rename(target)
    return {
        "src": str(src),
        "dst": str(target),
        "renamed": True,
        "warnings": warnings,
    }
