"""Shared result metadata for readable FEM validation JSON."""

from __future__ import annotations

import datetime as dt
import importlib.metadata as metadata
import platform
import sys
from pathlib import Path


def _pkg_version(name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def add_result_metadata(summary: dict, script_file: str) -> dict:
    """Return ``summary`` with timestamp and runtime version metadata."""
    out = dict(summary)
    out.setdefault("schema", f"radia.validation.fem_readable.{Path(script_file).stem}.v1")
    out["generated_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(
        timespec="seconds"
    ).replace("+00:00", "Z")
    out["versions"] = {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "radia_mcp_version": _pkg_version("radia-mcp"),
        "radia_version": _pkg_version("radia"),
    }
    return out
