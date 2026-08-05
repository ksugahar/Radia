"""Environment capability detection (detect_matlab_toolboxes twin).

One call answers "what can gmsh do on THIS machine": package presence,
gmsh version, build features (OCC, FLTK, MED, ...), and -- crucially --
whether an FLTK graphics context can actually be created (probed for
real in a subprocess), which decides upfront if gmsh_render /
gmsh_export_animation will work.
"""

from __future__ import annotations

import sys
from importlib import util as _importlib_util
from typing import Any

from ._gmsh_subprocess import gmsh_available, run_gmsh_json_subprocess

_DETECT_SCRIPT = r"""
import json
import sys

out_path = sys.argv[1]
result = {"ok": False, "ran": False}
try:
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        info = {
            "gmsh_version": gmsh.option.getString("General.Version"),
            "build_options": gmsh.option.getString("General.BuildOptions"),
        }
        graphics = {"available": False}
        try:
            gmsh.fltk.initialize()
            gmsh.fltk.finalize()
            graphics["available"] = True
        except Exception as exc:
            graphics["error"] = f"{type(exc).__name__}: {exc}"
        result.update({"ok": True, "ran": True, "info": info,
                       "graphics": graphics})
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""

_BUILD_FEATURE_TOKENS = (
    ("opencascade", "OpenCASCADE"),
    ("fltk", "Fltk"),
    ("med", "Med"),
    ("netgen", "Netgen"),
    ("cairo", "Cairo"),
    ("png", "Png"),
    ("jpeg", "Jpeg"),
    ("mpeg", "Mpeg"),
)


def detect_capabilities() -> dict[str, Any]:
    """Report gmsh availability, version, build features, graphics."""
    caps: dict[str, Any] = {
        "ok": True,
        "gmsh_package_installed": gmsh_available(),
        "pillow_installed": _importlib_util.find_spec("PIL") is not None,
        "python_version": sys.version.split()[0],
        "lanes": {
            "one_shot_gating": ["gmsh_inspect_msh", "gmsh_validate_msh",
                                "gmsh_validate_geo", "gmsh_field_stats",
                                "gmsh_diff_msh", "gmsh_audit_msh_directory",
                                "gmsh_verify"],
            "rendering": ["gmsh_render", "gmsh_export_animation"],
            "session": ["gmsh_exec", "gmsh_run_file",
                        "gmsh_session_status", "gmsh_session_shutdown"],
        },
    }
    if caps["gmsh_package_installed"]:
        probe = run_gmsh_json_subprocess(_DETECT_SCRIPT, [],
                                         timeout_s=180.0,
                                         prefix="radia_mcp_gmsh_detect_")
        if probe.get("ran"):
            build = probe["info"].get("build_options", "")
            caps["gmsh_version"] = probe["info"].get("gmsh_version")
            caps["graphics"] = probe["graphics"]
            caps["rendering_available"] = bool(
                probe["graphics"].get("available"))
            caps["build_features"] = {
                key: token in build for key, token in _BUILD_FEATURE_TOKENS}
        else:
            caps["gmsh_probe_error"] = probe.get("error")
            caps["rendering_available"] = False
    else:
        caps["rendering_available"] = False

    # Session state without side effects (import here to avoid a cycle
    # at module import time).
    from .session import GmshSession
    session = GmshSession.peek()
    caps["session_running"] = bool(session is not None and session.alive())
    return caps
