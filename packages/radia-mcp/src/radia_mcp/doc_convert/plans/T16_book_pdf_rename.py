"""Tier 1 - book-cover title extraction + PDF rename workflow.

Two MCP tools that together rename a scanned-book PDF to a filename
derived from its visual cover (first page).  The extract step shells
out to the local ``claude`` CLI with ``--allowedTools Read`` so the
title judgement runs in a normal Claude Code session against the
rendered cover image; the rename step is a pure filesystem op.

Workflow (peer to ``doc_convert_business_card_*`` in T15):

    1. ``doc_convert_extract_pdf_title_by_cover(pdf_path)`` -- render
       page 1 of the PDF to a temp PNG, ask Claude to read it and emit
       a single-line title.  Returns the extracted title (or ``None``
       on failure).

    2. ``doc_convert_rename_pdf_by_title(pdf_path, title,
       output_dir=None)`` -- Windows-safe sanitize the title and rename
       the source PDF in place.  Collision-safe (``" (2)"`` suffix).

Originally these helpers lived inline in
``O:\\BookScanner\\kindle_scanner.py``.  They are promoted to MCP tools
so any post-OCR workflow (Kindle scan, paper scan re-OCR, ad-hoc
title-from-cover) can reuse them without copy-paste.
"""
from __future__ import annotations

import pathlib
import re
import shutil
import subprocess

_FORBIDDEN_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_TRAIL_STRIP = re.compile(r"[\.\s]+$")
_TITLE_QUOTE_STRIP = '"\'「」『』 　'  # ASCII + JP brackets + ideographic space
_MAX_TITLE_LEN = 100  # filename safety; some FSes cap at 255 bytes UTF-8


def _sanitize_title(title: str) -> tuple[str, list[str]]:
    """Make ``title`` safe to use as a Windows filename stem.

    Returns ``(cleaned, warnings)``.  Mirrors the policy from
    ``T15_business_card_rename._sanitize_segment`` minus the
    underscore-separator escape (book titles are a single segment so
    underscores are fine).
    """
    warnings: list[str] = []
    if title is None:
        return "", ["title was None"]
    cleaned = str(title).strip(_TITLE_QUOTE_STRIP).strip()
    if not cleaned:
        return "", ["title was empty after strip"]
    sub = _FORBIDDEN_CHARS.sub("", cleaned)
    if sub != cleaned:
        warnings.append(
            f"removed forbidden chars ({cleaned!r} -> {sub!r})"
        )
    cleaned = sub
    collapsed = re.sub(r"\s+", " ", cleaned)
    if collapsed != cleaned:
        warnings.append("collapsed whitespace")
    cleaned = _TRAIL_STRIP.sub("", collapsed)
    if len(cleaned) > _MAX_TITLE_LEN:
        # Try to break on a space near the cap so we don't slice mid-word.
        head = cleaned[:_MAX_TITLE_LEN]
        if " " in head:
            head = head.rsplit(" ", 1)[0]
        warnings.append(f"truncated from {len(cleaned)} to {len(head)} chars")
        cleaned = head
    if not cleaned:
        warnings.append("nothing left after sanitization")
    return cleaned, warnings


def _unique_path(path: pathlib.Path) -> pathlib.Path:
    """Return ``path`` if free, else ``<stem> (2).<suffix>``, ``(3)``..."""
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    for i in range(2, 1000):
        candidate = parent / f"{stem} ({i}){suffix}"
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"could not find a free filename near {path}")


def _render_first_page(pdf_path: pathlib.Path, out_png: pathlib.Path,
                       dpi: int) -> bool:
    """Render page 1 of ``pdf_path`` to ``out_png`` using pymupdf.

    Returns True on success, False on failure.  Kept as a separate
    helper so the title-extract tool stays focused on the Claude
    interaction.
    """
    try:
        import pymupdf  # type: ignore
    except ImportError:
        return False
    try:
        with pymupdf.open(str(pdf_path)) as doc:
            if doc.page_count == 0:
                return False
            page = doc[0]
            mat = pymupdf.Matrix(dpi / 72.0, dpi / 72.0)
            pix = page.get_pixmap(matrix=mat)
            pix.save(str(out_png))
        return True
    except Exception:  # noqa: BLE001
        return False


def _ask_claude_for_title(image_path: pathlib.Path,
                          timeout_s: int) -> tuple[str | None, str, str | None]:
    """Run ``claude -p '<prompt>' --allowedTools Read`` against an image.

    Returns ``(title, raw_stdout, error)`` where exactly one of ``title``
    or ``error`` is non-None on a successful subprocess return.
    """
    claude_cmd = shutil.which("claude")
    if not claude_cmd:
        return None, "", "claude CLI not found on PATH"

    abs_image = str(image_path.resolve())
    prompt = (
        f"Read ツールで次の画像ファイルを読み込んでください: {abs_image}\n"
        "これは本の表紙または最初のページです。本のタイトルだけを1行で"
        "出力してください。\n"
        "説明文・前置き・引用符・括弧は一切不要。"
        "タイトル文字列のみ出力すること。"
    )

    try:
        result = subprocess.run(
            [claude_cmd, "-p", prompt, "--allowedTools", "Read"],
            capture_output=True, text=True,
            timeout=timeout_s, encoding="utf-8",
        )
    except subprocess.TimeoutExpired:
        return None, "", f"claude subprocess timed out after {timeout_s}s"
    except Exception as e:  # noqa: BLE001
        return None, "", f"{type(e).__name__}: {e}"

    raw = result.stdout or ""
    if result.returncode != 0:
        err = (f"claude exit {result.returncode}: "
               f"{(result.stderr or '').strip()[:500]}")
        return None, raw, err

    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
    if not lines:
        return None, raw, "claude stdout was empty"
    candidate = lines[-1].strip(_TITLE_QUOTE_STRIP)
    if not candidate or len(candidate) <= 1:
        return None, raw, "title candidate too short or empty"
    return candidate, raw, None


_IMAGE_EXTS_FOR_TITLE = (".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff")


def doc_convert_extract_title_from_image(
    image_path: str,
    timeout_s: int = 180,
) -> dict:
    """Ask Claude Code to read a cover image and return the book title.

    Spawns ``claude -p '<prompt>' --allowedTools Read`` against the
    given image.  Intended for scanner workflows where the cover is
    already a raster (the first scanned PNG) — no PDF intermediate
    needed.

    Args:
        image_path: path to a raster image (.png/.jpg/.jpeg/.bmp/.tif).
        timeout_s: kill the ``claude`` subprocess after this many seconds.

    Returns:
        dict with::

            {
              "image":   <input path>,
              "title":   "<extracted title>"     # or None on failure
              "raw":     "<full Claude stdout>"  # for inspection
              "error":   "<reason>"              # only on failure
            }

    Companion to ``doc_convert_extract_pdf_title_by_cover`` (which
    renders page 1 of a PDF to a PNG and delegates here).
    """
    src = pathlib.Path(image_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return {"image": str(src), "title": None,
                "error": f"file not found: {src}"}
    if src.suffix.lower() not in _IMAGE_EXTS_FOR_TITLE:
        return {"image": str(src), "title": None,
                "error": f"unsupported image extension: {src.suffix!r}"}

    title, raw, err = _ask_claude_for_title(src, timeout_s)
    out: dict = {"image": str(src), "title": title, "raw": raw}
    if err:
        out["error"] = err
    return out


def doc_convert_extract_pdf_title_by_cover(
    pdf_path: str,
    render_dpi: int = 144,
    timeout_s: int = 180,
    keep_image: bool = False,
) -> dict:
    """Render page 1 of a PDF and ask Claude Code to read its title.

    Thin wrapper: renders the cover to a temp PNG via pymupdf, then
    delegates to ``doc_convert_extract_title_from_image``.  Designed
    for re-OCR pipelines that already have a merged PDF; scanners that
    have the cover as a raw PNG should call
    ``doc_convert_extract_title_from_image`` directly.

    Args:
        pdf_path: source PDF (absolute or relative).
        render_dpi: resolution for the temp cover image (144 is enough
            for title legibility while staying small).
        timeout_s: kill the ``claude`` subprocess if it has not
            produced stdout within this many seconds.
        keep_image: if True, do not delete the temp PNG after the
            extraction (useful for debugging).  Default deletes it.

    Returns:
        dict with::

            {
              "pdf":     <input path>,
              "title":   "<extracted title>"     # or None on failure
              "raw":     "<full Claude stdout>"  # for inspection
              "image":   "<temp PNG path>"       # only if keep_image=True
              "error":   "<reason>"              # only on failure
            }
    """
    src = pathlib.Path(pdf_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return {"pdf": str(src), "title": None,
                "error": f"file not found: {src}"}
    if src.suffix.lower() != ".pdf":
        return {"pdf": str(src), "title": None,
                "error": f"not a .pdf: {src.suffix!r}"}

    img_path = src.with_suffix(".cover.png")
    if not _render_first_page(src, img_path, render_dpi):
        return {"pdf": str(src), "title": None,
                "error": "failed to render first page (pymupdf)"}

    title, raw, err = _ask_claude_for_title(img_path, timeout_s)
    out: dict = {"pdf": str(src), "title": title, "raw": raw}
    if err:
        out["error"] = err
    if keep_image:
        out["image"] = str(img_path)
    else:
        try:
            img_path.unlink()
        except OSError:
            pass
    return out


def doc_convert_rename_pdf_by_title(
    pdf_path: str,
    title: str,
    output_dir: str | None = None,
    dry_run: bool = False,
) -> dict:
    """Rename a PDF to ``<sanitized-title>.pdf`` (collision-safe).

    Sanitization (Windows-safe):
        - strip outer ASCII quotes + 「」 / 『』 / 全角空白 / spaces
        - remove ``< > : " / \\ | ? *`` and control chars
        - collapse internal whitespace runs to single spaces
        - strip trailing dots / spaces
        - cap length at 100 chars on a word boundary if possible

    If the target filename already exists, a parenthesised numeric
    suffix is appended (``<title> (2).pdf``, ``(3).pdf``, ...).

    Args:
        pdf_path: source PDF path.
        title: title string (typically from
            ``doc_convert_extract_pdf_title_by_cover``).
        output_dir: optional destination directory.  Default keeps the
            file in the source directory.
        dry_run: if True, return the resolved target path without
            actually renaming.

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
    src = pathlib.Path(pdf_path)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return {"src": str(src), "error": f"file not found: {src}"}
    if src.suffix.lower() != ".pdf":
        return {"src": str(src),
                "error": f"not a .pdf: {src.suffix!r}"}

    clean, warnings = _sanitize_title(title)
    if not clean:
        return {"src": str(src),
                "error": "title empty after sanitization",
                "warnings": warnings}

    dest_dir = pathlib.Path(output_dir) if output_dir else src.parent
    if not dest_dir.is_absolute():
        dest_dir = pathlib.Path.cwd() / dest_dir
    dest_dir.mkdir(parents=True, exist_ok=True)

    target = dest_dir / f"{clean}.pdf"
    target = _unique_path(target)

    if dry_run:
        return {"src": str(src), "dst": str(target),
                "renamed": False, "dry_run": True,
                "warnings": warnings}

    if target == src:
        return {"src": str(src), "dst": str(target),
                "renamed": False,
                "warnings": warnings + ["target == source; no rename needed"]}

    src.rename(target)
    return {"src": str(src), "dst": str(target),
            "renamed": True, "warnings": warnings}
