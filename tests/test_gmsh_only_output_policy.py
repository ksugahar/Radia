from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SCANNED_ROOTS = (
    ROOT / "src" / "core",
    ROOT / "src" / "lib",
    ROOT / "src" / "radia",
    ROOT / "src" / "ext",
    ROOT / "matlab",
    ROOT / "validation_test",
    ROOT / "packages" / "radia-mcp" / "src" / "radia_mcp" / "ih",
    ROOT / "packages" / "radia-mcp" / "src" / "radia_mcp" / "mor",
)

CUBIT_VTK_BOUNDARY = (
    "src/radia/cubit_toolbar_smoke.py",
    "src/radia/install_panels.py",
    "src/radia/panels/calc_common.py",
    "src/radia/panels/calc_mesh_eval.py",
    "src/radia/panels/cubit_toolbar/",
    "src/radia/panels/cubit_toolbar_probe.py",
    "src/radia/panels/radia_export_menu.py",
    "src/radia/panels/register_toolbar.py",
    "validation_test/cubit/",
    "validation_test/panels/test_radia_export_menu.py",
)

SOURCE_PATHS = tuple(
    str(root.relative_to(ROOT).as_posix()) + "/**"
    for root in SCANNED_ROOTS
) + ("docs/**/*.py",)
FORBIDDEN_OUTPUT = r"VTKOutput|[\"'][^\"']*\.(vtk|vtu|vts)[\"']"
IMPLEMENTATION_SUFFIXES = {
    ".bat", ".c", ".cc", ".cpp", ".h", ".hpp", ".m", ".ps1", ".py",
}


def _is_cubit_vtk_boundary(relative: str) -> bool:
    return any(
        relative == allowed or relative.startswith(allowed)
        for allowed in CUBIT_VTK_BOUNDARY
    )


def test_non_cubit_code_has_no_vtk_output_path() -> None:
    completed = subprocess.run(
        [
            "git",
            "grep",
            "-n",
            "-I",
            "-i",
            "-E",
            FORBIDDEN_OUTPUT,
            "--",
            *SOURCE_PATHS,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert completed.returncode in (0, 1), completed.stderr

    violations: list[str] = []
    for match in completed.stdout.splitlines():
        relative = match.split(":", 1)[0].replace("\\", "/")
        if Path(relative).suffix.lower() not in IMPLEMENTATION_SUFFIXES:
            continue
        if not _is_cubit_vtk_boundary(relative):
            violations.append(match)

    assert not violations, (
        "Only cubit-mesh-export may produce VTK. Use GmshPostExport and checked "
        ".msh v4.1 output elsewhere:\n" + "\n".join(violations)
    )


def test_docs_do_not_track_paraview_state_files() -> None:
    state_files = sorted(
        path.relative_to(ROOT).as_posix() for path in (ROOT / "docs").rglob("*.pvsm")
    )
    assert not state_files, (
        "Docs use saved notebook output and GMSH post-processing artifacts; "
        "ParaView state files are stale machine-specific UI state:\n"
        + "\n".join(state_files)
    )
