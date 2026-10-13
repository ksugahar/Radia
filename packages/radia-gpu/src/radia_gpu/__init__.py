"""Locate and verify the independently installed MATLAB GPU payload."""
from __future__ import annotations

import hashlib
import json
from importlib.resources import files
from pathlib import Path

__version__ = "0.1.3"


def matlab_path() -> Path:
    """Return installed resources only; never borrow an editable solver tree."""
    return Path(str(files(__package__).joinpath("matlab")))


def verify() -> dict:
    """Verify the installed payload without importing MATLAB, CUDA or Radia."""
    root = Path(str(files(__package__)))
    manifest_bytes = (root / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest.get("schema") != "radia.gpu.distribution.v1" or manifest.get("version") != __version__:
        raise RuntimeError("Unsupported or mismatched radia-gpu manifest")
    for relative, expected in manifest["files"].items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file():
            raise RuntimeError(f"Missing or invalid payload path: {relative}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Payload hash mismatch: {relative}")
    return dict(passed=True, runtime_checked=False, version=__version__, matlab_path=str(matlab_path()),
                source_commit=manifest["source_commit"], capabilities=manifest["capabilities"],
                manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                runtime_requirements=manifest["runtime_requirements"])
