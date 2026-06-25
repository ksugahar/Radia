"""Render the cuboid mesh as a PNG using a standalone GMSH executable.

Mesh files live in C:\\temp (ASCII-safe for GMSH on Windows).  This script
does not import the pip ``gmsh`` runtime; it writes a small ``.geo`` display
script and asks standalone ``gmsh.exe`` to print the PNG.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

SAFE = Path(r"C:\temp\cuboid_521_h020.msh")
SAFE_GEO = Path(r"C:\temp\cuboid_521_h020_render.geo")
SAFE_PNG = Path(r"C:\temp\cuboid_521_h020_mesh.png")
PNG_OUT_DIR = Path(__file__).parent.parent  # CLN root, where the .tex lives
PNG_OUT = PNG_OUT_DIR / "cuboid_521_h020_mesh.png"


def _gmsh_executable() -> str:
    """Return a standalone gmsh.exe path, rejecting pip's gmsh.bat wrapper."""
    candidates = [
        Path(r"C:\gmsh.exe"),
        Path(r"C:\tools\gmsh.exe"),
        Path(r"C:\Program Files\Gmsh\gmsh.exe"),
        Path(r"C:\Program Files (x86)\Gmsh\gmsh.exe"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    path_hit = shutil.which("gmsh.exe")
    if path_hit and Path(path_hit).suffix.lower() == ".exe":
        return path_hit

    raise FileNotFoundError(
        "standalone gmsh.exe was not found. Install standalone GMSH or set it "
        "on PATH as gmsh.exe; the pip gmsh.bat wrapper is intentionally not used."
    )


def _geo_path(path: Path) -> str:
    return path.resolve().as_posix()


def _write_render_geo() -> None:
    SAFE_GEO.write_text(
        "\n".join(
            [
                f'Merge "{_geo_path(SAFE)}";',
                "General.Trackball = 0;",
                "General.RotationX = -75;",
                "General.RotationY = -10;",
                "General.RotationZ = -45;",
                "Mesh.SurfaceFaces = 1;",
                "Mesh.SurfaceEdges = 1;",
                "Mesh.VolumeEdges = 0;",
                "Mesh.ColorCarousel = 2;",
                "General.GraphicsWidth = 1000;",
                "General.GraphicsHeight = 700;",
                f'Print "{_geo_path(SAFE_PNG)}";',
                "Exit;",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    if not SAFE.exists():
        raise FileNotFoundError(
            f"mesh not found: {SAFE}. Run save_mesh_and_view.py --no-view first."
        )

    gmsh_exe = _gmsh_executable()
    _write_render_geo()
    subprocess.run([gmsh_exe, str(SAFE_GEO), "-nopopup"], check=True)

    if not SAFE_PNG.exists():
        raise RuntimeError(f"gmsh did not write expected PNG: {SAFE_PNG}")

    shutil.copy(SAFE_PNG, PNG_OUT)
    print(f"Wrote {PNG_OUT}")


if __name__ == "__main__":
    main()
