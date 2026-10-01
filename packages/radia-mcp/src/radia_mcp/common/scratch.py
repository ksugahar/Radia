"""Disposable working-directory root for generated MCP artifacts."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


def scratch_root() -> Path:
    """Return ``RADIA_MCP_TEMP``, else ``C:\\temp`` on Windows, else the OS temp dir."""
    configured = os.environ.get("RADIA_MCP_TEMP", "").strip()
    if configured:
        return Path(configured)
    return Path(r"C:\temp") if os.name == "nt" else Path(tempfile.gettempdir())
