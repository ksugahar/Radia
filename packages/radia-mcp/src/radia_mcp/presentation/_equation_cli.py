"""Eqnedit64 command-line bridge and equation-product policy for presentations.

The native and Web/JavaScript editions are both maintained under Radia's
``tools/eqnedit64`` tree.  The laboratory homepage is the Web publication
surface; the native executable remains the Windows rendering and clipboard
backend for ``radia_mcp.presentation``.
"""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable


_COPY_TARGETS = {
    "office": ("--copy-tex-file", [
        "HTML Format", "CF_UNICODETEXT", "LaTeX",
        "CF_ENHMETAFILE", "CF_DIBV5",
    ]),
    "powerpoint": ("--copy-tex-file", [
        "HTML Format", "CF_UNICODETEXT", "LaTeX",
        "CF_ENHMETAFILE", "CF_DIBV5",
    ]),
    "google-slides": ("--copy-google-slides-file", ["HTML Format", "PNG"]),
    "png": ("--copy-png-file", ["PNG", "CF_DIBV5"]),
}


def presentation_equation_policy() -> dict:
    """Return the canonical source, publication, and format responsibilities."""
    return {
        "schema": "radia-mcp.presentation-equation-policy.v2",
        "source_of_truth": {
            "native": "tools/eqnedit64/src",
            "web": "tools/eqnedit64/web/equation-editor.js",
            "web_mount": "tools/eqnedit64/web/equation-editor.fragment.html",
            "parity": "tools/eqnedit64/docs/PRODUCT_PARITY.md",
        },
        "publication": {
            "web": "Sugahara Laboratory homepage",
            "native": "ksugahar/Radia GitHub Releases",
        },
        "web_publication_contract": {
            "mode": "build-time import from Radia checkout",
            "repository": "https://github.com/ksugahar/Radia",
            "checkout_environment_variable": "RADIA_REPOSITORY",
            "integrity": "SHA-256 equality after copy",
            "homepage_source_copy": False,
            "release_gate": "site_builder/tools/run_eqnedit64_release_qa.ps1",
            "test_scope": "Eqnedit64 page, browser contract, and hidden PowerPoint render only",
        },
        "presentation_paths": {
            "batch_pptx": "radia.equation.markdown_to_pptx (native OMML)",
            "interactive_clipboard": "presentation_copy_equation (Eqnedit64)",
            "image_file": "presentation_render_equation (Eqnedit64)",
        },
        "powerpoint_normal_paste": {
            "command": 'Application.CommandBars.ExecuteMso("Paste")',
            "clipboard_carrier": "HTML Format with inline MathML",
            "office_math": "editable inline m:oMath",
            "alignment": "left",
            "font_points": 18,
            "font_scope": "standard blank PowerPoint presentation",
            "priority": "left alignment over forced 24 pt",
            "native_web_parity": True,
            "forbidden_acceptance_probe": "slide.Shapes.Paste()",
            "reason": (
                "Shapes.Paste preserves 24 pt through a different object-model "
                "path; ordinary Ctrl+V applies PowerPoint destination formatting"
            ),
        },
        "canonical_input": "TeX",
        "retired_formats": ["MTEF", ".eqn"],
        "homepage_is_source_of_truth": False,
    }


def _candidate_executables(explicit: str | None) -> Iterable[Path]:
    if explicit:
        yield Path(explicit).expanduser()
        return
    configured = os.environ.get("EQNEDIT64_EXE")
    if configured:
        yield Path(configured).expanduser()
    found = shutil.which("Eqnedit64.exe") or shutil.which("Eqnedit64")
    if found:
        yield Path(found)
    # Editable Radia checkout: find tools/eqnedit64 without depending on a
    # drive letter or on the laboratory archive path.
    for parent in Path(__file__).resolve().parents:
        component = parent / "tools" / "eqnedit64"
        yield component / "dist" / "Eqnedit64.exe"
        yield component / "build" / "Eqnedit64.exe"


def _resolve_executable(explicit: str | None) -> Path:
    checked: list[str] = []
    for candidate in _candidate_executables(explicit):
        resolved = candidate.resolve()
        text = str(resolved)
        if text in checked:
            continue
        checked.append(text)
        if resolved.is_file():
            return resolved
    raise FileNotFoundError(
        "Eqnedit64.exe was not found. Pass executable=..., set "
        "EQNEDIT64_EXE, place it on PATH, or build tools/eqnedit64."
    )


def _invoke(command: list[str], timeout_s: float) -> subprocess.CompletedProcess[str]:
    creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
    return subprocess.run(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout_s,
        check=False,
        creationflags=creationflags,
    )


def _tex_file(tex: str) -> Path:
    if not isinstance(tex, str) or not tex.strip():
        raise ValueError("tex must be a non-empty string")
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="\n", suffix=".tex",
        prefix="radia-mcp-eqnedit64-", delete=False,
    )
    try:
        handle.write(tex)
    finally:
        handle.close()
    return Path(handle.name)


def _failure(operation: str, exc: Exception) -> dict:
    return {
        "ok": False,
        "backend": "Eqnedit64",
        "operation": operation,
        "error": f"{type(exc).__name__}: {exc}",
    }


def presentation_equation_backend(executable: str | None = None) -> dict:
    """Locate the standalone Eqnedit64 backend used by presentation tools."""
    try:
        path = _resolve_executable(executable)
        return {"ok": True, "backend": "Eqnedit64", "executable": str(path)}
    except (OSError, ValueError) as exc:
        return _failure("locate", exc)


def presentation_copy_equation(
    tex: str,
    target: str = "office",
    executable: str | None = None,
    timeout_s: float = 30.0,
) -> dict:
    """Copy TeX using Eqnedit64's native clipboard contract.

    ``target="office"`` (or ``"powerpoint"``) publishes editable, left-aligned
    18 pt Office Math through normal PowerPoint paste, plus TeX, EMF, and opaque
    DIBV5 fallbacks. ``"google-slides"``
    publishes the 300 dpi/24 pt PNG+HTML contract. ``"png"`` publishes only
    PNG and DIBV5 image formats. The user's foreground window is never
    activated.
    """
    operation = "copy"
    input_path: Path | None = None
    try:
        key = str(target).strip().lower().replace("_", "-")
        if key not in _COPY_TARGETS:
            raise ValueError("target must be office, powerpoint, google-slides, or png")
        if not 0 < float(timeout_s) <= 120:
            raise ValueError("timeout_s must be in (0, 120]")
        app = _resolve_executable(executable)
        input_path = _tex_file(tex)
        switch, formats = _COPY_TARGETS[key]
        completed = _invoke([str(app), switch, str(input_path)], float(timeout_s))
        if completed.returncode != 0:
            raise RuntimeError(
                f"Eqnedit64 exited {completed.returncode}: {completed.stderr.strip()}"
            )
        return {
            "ok": True,
            "backend": "Eqnedit64",
            "operation": operation,
            "target": key,
            "formats": formats,
            "executable": str(app),
        }
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        return _failure(operation, exc)
    finally:
        if input_path is not None:
            input_path.unlink(missing_ok=True)


def presentation_render_equation(
    tex: str,
    output_path: str,
    image_format: str = "png",
    executable: str | None = None,
    timeout_s: float = 30.0,
) -> dict:
    """Render TeX to a PNG or EMF file through the Eqnedit64 CLI."""
    operation = "render"
    input_path: Path | None = None
    try:
        fmt = str(image_format).strip().lower()
        if fmt not in {"png", "emf"}:
            raise ValueError("image_format must be png or emf")
        if not 0 < float(timeout_s) <= 120:
            raise ValueError("timeout_s must be in (0, 120]")
        app = _resolve_executable(executable)
        output = Path(output_path).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        input_path = _tex_file(tex)
        switch = "--render-png-file" if fmt == "png" else "--render-emf-file"
        completed = _invoke(
            [str(app), switch, str(input_path), str(output)], float(timeout_s)
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"Eqnedit64 exited {completed.returncode}: {completed.stderr.strip()}"
            )
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError("Eqnedit64 reported success without a non-empty output file")
        result = {
            "ok": True,
            "backend": "Eqnedit64",
            "operation": operation,
            "format": fmt,
            "output_path": str(output),
            "bytes": output.stat().st_size,
            "executable": str(app),
        }
        if fmt == "png":
            fields = completed.stdout.strip().split()
            if len(fields) == 2 and all(field.isdigit() for field in fields):
                result["pixel_size"] = [int(fields[0]), int(fields[1])]
                result["dpi"] = 300
                result["font_points"] = 24
        return result
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        return _failure(operation, exc)
    finally:
        if input_path is not None:
            input_path.unlink(missing_ok=True)
