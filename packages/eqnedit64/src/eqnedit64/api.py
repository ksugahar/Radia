"""Checked wrappers around the bundled standalone Eqnedit64 executable."""
from __future__ import annotations

from importlib.resources import files
from pathlib import Path
import subprocess
import tempfile


_COPY_SWITCHES = {
    "office": "--copy-tex-file",
    "powerpoint": "--copy-tex-file",
    "google-slides": "--copy-google-slides-file",
    "png": "--copy-png-file",
}
_RENDER_SWITCHES = {
    ".png": "--render-png-file",
    ".emf": "--render-emf-file",
}


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
    if normalized not in _COPY_SWITCHES:
        raise ValueError("target must be office, powerpoint, google-slides, or png")
    source = _tex_file(tex)
    try:
        result = _invoke([_COPY_SWITCHES[normalized], str(source)], timeout_s)
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
    if suffix not in _RENDER_SWITCHES:
        raise ValueError("output extension must be .png or .emf")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = _tex_file(tex)
    try:
        result = _invoke(
            [_RENDER_SWITCHES[suffix], str(source), str(destination)], timeout_s
        )
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
