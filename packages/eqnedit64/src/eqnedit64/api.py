"""Checked wrappers around the bundled standalone Eqnedit64 executable."""
from __future__ import annotations

from importlib.resources import files
from pathlib import Path
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ET


_COPY_TARGETS = {
    "office": "office",
    "powerpoint": "office",
    "google-slides": "slides",
    "png": "clipboard-png",
}
_RENDER_SUFFIXES = {".png", ".emf"}


def export_powerpoint_math(tex: str, output: str | Path,
                           timeout_s: float = 30.0) -> Path:
    """Export checked native PowerPoint OOXML math, with no clipboard/images.

    Unsupported TeX raises. The destination is replaced only after successful
    CLI execution and namespace validation; an existing file survives errors.
    """
    if not 0 < float(timeout_s) <= 120:
        raise ValueError("timeout_s must be in (0, 120]")
    if not isinstance(tex, str) or not tex.strip():
        raise ValueError("tex must be a non-empty string")
    destination = Path(output).expanduser().resolve()
    if destination.suffix.lower() != ".pptxmath":
        raise ValueError("output extension must be .pptxmath")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = _tex_file(tex)
    try:
        with tempfile.TemporaryDirectory(prefix="eqnedit64-math-", dir=destination.parent) as scratch:
            candidate = Path(scratch) / "equation.pptxmath"
            result = _invoke([str(source), str(candidate)], timeout_s)
            if result.returncode:
                raise RuntimeError(
                    f"Eqnedit64 math export exited {result.returncode}: "
                    f"{(result.stderr or result.stdout).strip()}"
                )
            if not candidate.is_file() or not 0 < candidate.stat().st_size <= 10_000_000:
                raise RuntimeError("Eqnedit64 did not create a bounded math output")
            root = ET.fromstring(candidate.read_bytes())
            math = "http://schemas.openxmlformats.org/officeDocument/2006/math"
            drawing = "http://schemas.openxmlformats.org/drawingml/2006/main"
            if root.tag != f"{{{math}}}oMathPara" or root.find(f"{{{math}}}oMath") is None:
                raise RuntimeError("Invalid PowerPoint math root")
            if any(not node.tag.startswith((f"{{{math}}}", f"{{{drawing}}}"))
                   or node.tag == f"{{{drawing}}}t" for node in root.iter()):
                raise RuntimeError("Invalid PowerPoint math dialect")
            os.replace(candidate, destination)
    finally:
        source.unlink(missing_ok=True)
    return destination


def backend_path() -> Path:
    """Return the bundled standalone executable path."""
    path = Path(str(files("eqnedit64").joinpath("Eqnedit64.exe")))
    if not path.is_file():
        raise FileNotFoundError(f"bundled Eqnedit64.exe is missing: {path}")
    return path


def web_asset(name: str = "equation-editor.js") -> Path:
    """Return a checked path to a bundled browser-editor asset."""
    if name not in {"equation-editor.js", "equation-editor.fragment.html"}:
        raise ValueError("unknown Eqnedit64 Web asset")
    path = Path(str(files("eqnedit64").joinpath("web", name)))
    if not path.is_file():
        raise FileNotFoundError(f"bundled Eqnedit64 Web asset is missing: {path}")
    return path


def _invoke(arguments: list[str], timeout_s: float) -> subprocess.CompletedProcess[str]:
    creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    return subprocess.run(
        [str(backend_path()), *arguments],
        capture_output=True,
        check=False,
        creationflags=creationflags,
        text=True,
        timeout=timeout_s,
    )


def _tex_file(tex: str) -> Path:
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tex", encoding="utf-8", newline="\n", delete=False
    )
    try:
        handle.write(tex)
        return Path(handle.name)
    finally:
        handle.close()


def copy_equation(tex: str, target: str = "office", timeout_s: float = 15.0) -> None:
    """Publish TeX to the Windows clipboard using a named target contract."""
    normalized = target.strip().lower().replace("_", "-")
    if normalized not in _COPY_TARGETS:
        raise ValueError("target must be office, powerpoint, google-slides, or png")
    source = _tex_file(tex)
    try:
        result = _invoke([str(source), _COPY_TARGETS[normalized]], timeout_s)
    finally:
        source.unlink(missing_ok=True)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"Eqnedit64 clipboard command exited {result.returncode}: {detail}"
        )


def render_equation(tex: str, output: str | Path, timeout_s: float = 30.0) -> Path:
    """Render TeX to a checked PNG or EMF output file."""
    destination = Path(output).expanduser().resolve()
    suffix = destination.suffix.lower()
    if suffix not in _RENDER_SUFFIXES:
        raise ValueError("output extension must be .png or .emf")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = _tex_file(tex)
    try:
        result = _invoke([str(source), str(destination)], timeout_s)
    finally:
        source.unlink(missing_ok=True)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"Eqnedit64 render command exited {result.returncode}: {detail}"
        )
    if not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError("Eqnedit64 did not create a non-empty output file")
    return destination
