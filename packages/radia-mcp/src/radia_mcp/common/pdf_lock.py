"""Release a PDF held open by a dedicated viewer before it is rewritten (Windows)."""

from __future__ import annotations

import os
import pathlib
import subprocess

# Dedicated viewers only. A browser shows a PDF as one tab of a shared process,
# so stopping it would close every unrelated window; a browser-held PDF stays
# locked and the caller's write fails loudly instead.
PDF_VIEWER_PROCESSES = (
    "Acrobat", "AcroRd32", "SumatraPDF", "FoxitReader", "FoxitPDFReader",
)

# The file name reaches PowerShell through the environment, never through the
# script text, so quotes or metacharacters in it cannot alter the command.
_PROBE = (
    "$name = $env:RADIA_MCP_PDF_LOCK_NAME; "
    "Get-Process -Name " + ",".join(PDF_VIEWER_PROCESSES) + " -ErrorAction SilentlyContinue | "
    "Where-Object { $_.MainWindowTitle.IndexOf($name, "
    "[System.StringComparison]::OrdinalIgnoreCase) -ge 0 } | "
    "Select-Object -ExpandProperty Id"
)


def _is_locked(pdf_path: pathlib.Path) -> bool:
    try:
        with open(pdf_path, "r+b"):
            return False
    except PermissionError:
        return True
    except OSError:
        return False


def release_pdf_lock(pdf_path: pathlib.Path) -> list[str]:
    """Stop dedicated PDF viewers whose window shows ``pdf_path``; return their PIDs."""
    if os.name != "nt" or not pdf_path.exists() or not _is_locked(pdf_path):
        return []
    env = dict(os.environ, RADIA_MCP_PDF_LOCK_NAME=pdf_path.name)
    try:
        probe = subprocess.run(
            ["powershell", "-NoProfile", "-NonInteractive", "-Command", _PROBE],
            capture_output=True, text=True, timeout=15, env=env,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return []
    stopped: list[str] = []
    for pid in (probe.stdout or "").split():
        if not pid.isdigit():
            continue
        try:
            result = subprocess.run(
                ["taskkill", "/PID", pid, "/F"],
                capture_output=True, text=True, timeout=10,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        if result.returncode == 0:
            stopped.append(pid)
    return stopped
