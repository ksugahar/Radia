"""Generate the cuboid mesh used in tree-cotree validation, save as .vol and .msh,
then optionally launch GMSH GUI to view it.

Usage:
    python save_mesh_and_view.py            # save + view in GMSH
    python save_mesh_and_view.py --no-view  # save only
"""
import sys
import subprocess
from pathlib import Path

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from ngsolve import TaskManager

AX, AY, AZ = 5e-3, 2e-3, 1e-3
H_COND = 0.20e-3

OUT_DIR = Path(r"C:\temp")  # ASCII-safe for GMSH on Windows (cp932 path issues)
VOL_PATH = OUT_DIR / "cuboid_521_h020.vol"
MSH_PATH = OUT_DIR / "cuboid_521_h020.msh"


def main(view: bool):
    print(f"Building {AX*1000}x{AY*1000}x{AZ*1000} mm Cu cuboid (h={H_COND*1000} mm)")
    box = Box(Pnt(0, 0, 0), Pnt(AX, AY, AZ))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    geo = OCCGeometry(box)

    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=H_COND))
        print(f"  ne = {mesh.ne}, nv = {mesh.nv}, nedge = {mesh.nedge}, nface = {mesh.nface}")

        mesh.ngmesh.Save(str(VOL_PATH))
        print(f"  Wrote .vol: {VOL_PATH}")

        # Convert to gmsh format via NGSolve (ngmesh.Export)
        try:
            mesh.ngmesh.Export(str(MSH_PATH), "Gmsh2 Format")
            print(f"  Wrote .msh: {MSH_PATH}")
        except Exception as e:
            print(f"  .msh export failed: {e} — using gmsh.merge fallback in viewer")

        if view:
            print("\nLaunching GMSH GUI...")
            target = MSH_PATH if MSH_PATH.exists() else VOL_PATH
            subprocess.Popen(["gmsh", str(target)])


if __name__ == "__main__":
    main(view="--no-view" not in sys.argv)
