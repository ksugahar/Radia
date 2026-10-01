"""Every third-party module the package imports is declared somewhere in its metadata."""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
SOURCE = PACKAGE / "src" / "radia_mcp"

# Import name -> requirement that provides it (distribution, optionally with an extra).
PROVIDED_BY = {
    "anyio": "mcp", "pydantic": "mcp", "certifi": "mcp",   # hard dependencies of mcp
    "OCP": "build123d",
    "cycler": "matplotlib",
    "PIL": "Pillow",
    "cv2": "opencv-python-headless",
    "fitz": "pymupdf",
    "pptx": "python-pptx",
    "win32com": "pywin32", "winerror": "pywin32",
    "google": "google-cloud-vision",
    "netgen": "radia", "ngsolve": "radia",                # radia depends on NGSolve
    "torch": "radia[urn]",
    "youtube_transcript_api": "youtube-transcript-api",
    "edge_tts": "edge-tts",
    "cubit_mesh_export": "cubit-mesh-export",
}


def _requirements() -> set[str]:
    """Normalized ``name`` and ``name[extra]`` of every declared requirement."""
    text = (PACKAGE / "pyproject.toml").read_text(encoding="utf-8")
    declared = set()
    for body in re.findall(r"^[\w-]+\s*=\s*\[(.*?)^\]", text, re.MULTILINE | re.DOTALL):
        for requirement in re.findall(r'"([^"]+)"', body):
            match = re.match(r"([A-Za-z0-9_.-]+)(\[[^\]]+\])?", requirement)
            if not match:
                continue
            name = match.group(1).lower().replace("_", "-")
            declared.add(name)
            for extra in (match.group(2) or "").strip("[]").split(","):
                if extra.strip():
                    declared.add(f"{name}[{extra.strip()}]")
    return declared


def _third_party_imports() -> dict[str, str]:
    stdlib = set(sys.stdlib_module_names)
    found: dict[str, str] = {}
    for path in SOURCE.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            else:
                continue
            for name in names:
                top = name.split(".")[0]
                if top not in stdlib and top not in {"radia_mcp", "__future__"}:
                    found.setdefault(top, path.relative_to(SOURCE).as_posix())
    return found


def test_every_imported_distribution_is_declared():
    declared = _requirements()
    missing = {}
    for module, where in sorted(_third_party_imports().items()):
        requirement = PROVIDED_BY.get(module, module).lower().replace("_", "-")
        if requirement not in declared:
            missing[module] = f"{requirement} (first import: {where})"
    assert missing == {}
