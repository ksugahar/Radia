"""pdf2ppt_pdfgear -- Drive PDFgear PDF->PPTX conversion via UIAutomation.

Reconstructed 2026-06-14 from SKILL.md after the original was lost to a
parallel-session race condition. The functional spec is unchanged:

  1. Write conv_docs.json with ConvertType=PDFToPPT
  2. Launch:   pdfconverter.exe "app1" "PDFToPPT" "<conv_docs.json>"
  3. Find the window by pdfconverter.exe PID (NOT by title -- the window
     title carries the document name, which is unreliable).
  4. Ensure the ``上級モード`` (Advanced Mode) CheckBox is ON via the UIA
     TogglePattern. This is what makes the output EDITABLE (text boxes +
     vector shapes, not one rasterised picture per slide).
  5. Find the ``変換`` button and fire it via the UIA InvokePattern.
     Invoke works even when the control reports visible=False, and no
     mouse hijack is needed.
  6. Poll the PDFgear output dir for ``<stem> conv.pptx`` (note the
     literal space before ``conv``).
  7. If --output was given, move the result there.

Caveats inherited from SKILL.md:

  - 上級モード OFF -> one PICTURE per slide (non-editable image-only).
  - The conversion POSTs the PDF to apiw.pdfgear.com -- do NOT use for
    confidential or unpublished documents.
  - If PDFgear's free cloud tier needs an interactive sign-in, the
    unattended click stalls and we time out with a clear message.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# PDFgear paths and constants
# ---------------------------------------------------------------------------
PDFCONVERTER_EXE = r"C:\Program Files\PDFgear\pdfconverter.exe"

# PDFgear's "出力先" setting. The default install drops converted files
# into %USERPROFILE%\OneDrive\PDFgear; override with --outdir if changed.
DEFAULT_OUTDIR = Path(os.environ.get("USERPROFILE", "")) / "OneDrive" / "PDFgear"

# PDFgear appends " conv" before the extension (literal space).
OUTPUT_SUFFIX = " conv.pptx"

# 240 s covers a typical multi-page deck through the cloud round-trip.
DEFAULT_TIMEOUT = 240.0

# .NET DateTime.Ticks = 100-nanosecond units since 0001-01-01 UTC.
_NET_EPOCH = datetime(1, 1, 1, tzinfo=timezone.utc)
_TICKS_PER_SECOND = 10_000_000


def _net_ticks_now() -> int:
    delta = datetime.now(timezone.utc) - _NET_EPOCH
    return int(delta.total_seconds() * _TICKS_PER_SECOND)


def _write_conv_docs(pdf_path: Path, json_dir: Path) -> Path:
    """Write conv_docs.json in the schema PDFgear's standalone converter wants."""
    out = json_dir / "conv_docs.json"
    data = {
        "ConvertType": "PDFToPPT",
        "FilesPath": [
            {
                "FilePath": str(pdf_path),
                "OpenDate": str(_net_ticks_now()),
                "Password": None,
            }
        ],
    }
    out.write_text(json.dumps(data), encoding="utf-8")
    return out


def _connect_top_window(pid: int, timeout: float = 30.0):
    """Wait for the top-level pdfconverter window owned by ``pid``.

    Match by PID -- not by window title -- per SKILL.md note: the title
    carries the document name and changes during loading.
    """
    from pywinauto import Application

    end = time.monotonic() + timeout
    last_err: Exception | None = None
    while time.monotonic() < end:
        try:
            app = Application(backend="uia").connect(process=pid)
            top = app.top_window()
            top.wait("visible exists", timeout=2.0)
            return app, top
        except Exception as e:  # window not up yet
            last_err = e
            time.sleep(0.5)
    raise RuntimeError(
        "Could not connect to pdfconverter.exe "
        f"(PID {pid}) within {timeout:.0f}s: {last_err}"
    )


def _set_advanced_mode(top, on: bool) -> bool:
    """Toggle the ``上級モード`` (Advanced Mode) CheckBox via TogglePattern.

    Returns True if a matching CheckBox was found and set, False otherwise.
    Looks for either the Japanese or English label.
    """
    for ck in top.descendants(control_type="CheckBox"):
        try:
            name = (ck.window_text() or "").strip()
        except Exception:
            continue
        if not name:
            continue
        nl = name.lower()
        if "上級" in name or "advanced" in nl:
            try:
                tp = ck.iface_toggle  # UIA TogglePattern
                state = tp.CurrentToggleState  # 0=Off, 1=On, 2=Indeterminate
                if (state == 1) != bool(on):
                    tp.Toggle()
                return True
            except Exception:
                continue
    return False


def _click_convert(top, retries: int = 20, sleep: float = 0.5) -> None:
    """Find the ``変換`` button and fire it via InvokePattern.

    InvokePattern works even when the control reports visible=False, so
    we do NOT need to bring it into view or move the mouse.
    """
    last_err: Exception | None = None
    for _ in range(retries):
        for btn in top.descendants(control_type="Button"):
            try:
                name = (btn.window_text() or "").strip()
            except Exception:
                continue
            if name in ("変換", "Convert"):
                try:
                    btn.iface_invoke.Invoke()
                    return
                except Exception as e:
                    last_err = e
                    continue
        time.sleep(sleep)
    raise RuntimeError(
        "Could not find or invoke the 変換 / Convert button. "
        f"Last error: {last_err}"
    )


def _poll_output(outdir: Path, expected_name: str, timeout: float) -> Path:
    """Wait for ``outdir/expected_name`` to appear and finish being written."""
    target = outdir / expected_name
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        if target.exists():
            try:
                size0 = target.stat().st_size
                if size0 > 0:
                    time.sleep(1.0)  # let the write finalise
                    if target.stat().st_size == size0:
                        return target
            except OSError:
                pass
        time.sleep(1.0)
    raise TimeoutError(
        f"Output {target} did not appear within {timeout:.0f}s.\n"
        "Likely causes:\n"
        "  - PDFgear cloud sign-in required (open the GUI once to sign in)\n"
        "  - Network down or apiw.pdfgear.com unreachable\n"
        "  - PDFgear UI changed (上級モード checkbox or 変換 button not found)\n"
        "  - Output dir setting differs from the default; pass --outdir"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Drive PDFgear's standalone converter for PDF -> PPTX via "
            "UIAutomation. Output is EDITABLE pptx (text boxes + vector "
            "shapes) when --no-advanced is NOT set."
        ),
    )
    ap.add_argument("pdf", help="input PDF path")
    ap.add_argument(
        "-o", "--output",
        help="destination .pptx path; default = leave at PDFgear output dir",
    )
    ap.add_argument(
        "--outdir",
        help=f"PDFgear output dir override (default: {DEFAULT_OUTDIR})",
    )
    ap.add_argument(
        "--timeout", type=float, default=DEFAULT_TIMEOUT,
        help=f"seconds to wait for converted .pptx (default {DEFAULT_TIMEOUT:g})",
    )
    ap.add_argument(
        "--no-advanced", action="store_true",
        help=(
            "do NOT enable 上級モード -- produces one PICTURE per slide "
            "(image-only, non-editable). Default is to enable it."
        ),
    )
    args = ap.parse_args(argv)

    pdf_path = Path(args.pdf).resolve()
    if not pdf_path.is_file():
        print(f"ERROR: PDF not found: {pdf_path}", file=sys.stderr)
        return 2
    if not Path(PDFCONVERTER_EXE).is_file():
        print(f"ERROR: PDFgear not installed at {PDFCONVERTER_EXE}", file=sys.stderr)
        return 3

    outdir = Path(args.outdir).resolve() if args.outdir else DEFAULT_OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)
    expected_name = pdf_path.stem + OUTPUT_SUFFIX

    # Remove any stale prior output so the poll only catches the new one.
    stale = outdir / expected_name
    if stale.exists():
        try:
            stale.unlink()
        except OSError as e:
            print(f"WARN: could not delete stale {stale}: {e}", file=sys.stderr)

    proc: subprocess.Popen | None = None
    with tempfile.TemporaryDirectory(prefix="pdf2ppt_") as tmp_dir:
        json_path = _write_conv_docs(pdf_path, Path(tmp_dir))

        # Launch the standalone converter pre-loaded in PDFToPPT mode.
        # The "app1" positional and the PDFToPPT mode keyword are the
        # canonical sequence PDFgear's editor uses internally.
        proc = subprocess.Popen(
            [PDFCONVERTER_EXE, "app1", "PDFToPPT", str(json_path)],
            cwd=str(Path(PDFCONVERTER_EXE).parent),
        )
        try:
            _app, top = _connect_top_window(proc.pid, timeout=30.0)

            # 上級モード ON -> editable pptx. Default unless --no-advanced.
            if not args.no_advanced:
                _set_advanced_mode(top, on=True)

            # Fire the 変換 button via InvokePattern.
            _click_convert(top)

            # Poll for the produced file.
            produced = _poll_output(outdir, expected_name, timeout=args.timeout)
        finally:
            # Close the converter window. PDFgear also pops an Explorer
            # window at the output dir (its "変換後エクスプローラーで表示"
            # setting); we leave that alone -- it is the user's setting,
            # not ours to override.
            if proc is not None:
                try:
                    proc.terminate()
                except OSError:
                    pass

    # If the user asked for a specific destination, move the result.
    if args.output:
        dst = Path(args.output).resolve()
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            dst.unlink()
        shutil.move(str(produced), str(dst))
        produced = dst

    print(str(produced))
    return 0


if __name__ == "__main__":
    sys.exit(main())
