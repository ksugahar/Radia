"""Render the cuboid mesh as a PNG using the GMSH Python API headless.

Mesh files live in C:\\temp (ASCII-safe — GMSH on Windows can't open paths
containing Japanese chars). PNG copied back to ngsolve_validation/ for tex use.
"""
import shutil
import gmsh
from pathlib import Path

SAFE = Path(r"C:\temp\cuboid_521_h020.msh")
SAFE_PNG = Path(r"C:\temp\cuboid_521_h020_mesh.png")
PNG_OUT_DIR = Path(__file__).parent.parent  # CLN root, where the .tex lives
PNG_OUT = PNG_OUT_DIR / "cuboid_521_h020_mesh.png"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.merge(str(SAFE))

gmsh.option.setNumber("General.Trackball", 0)
gmsh.option.setNumber("General.RotationX", -75)
gmsh.option.setNumber("General.RotationY", -10)
gmsh.option.setNumber("General.RotationZ", -45)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
gmsh.option.setNumber("Mesh.VolumeEdges", 0)
gmsh.option.setNumber("Mesh.ColorCarousel", 2)
gmsh.option.setNumber("General.GraphicsWidth", 1000)
gmsh.option.setNumber("General.GraphicsHeight", 700)
gmsh.option.setColor("General.Background", 255, 255, 255)
gmsh.option.setColor("General.Foreground", 0, 0, 0)
gmsh.option.setColor("General.Text", 0, 0, 0)

gmsh.fltk.initialize()
gmsh.write(str(SAFE_PNG))
shutil.copy(SAFE_PNG, PNG_OUT)
print(f"Wrote {PNG_OUT}")
gmsh.finalize()
