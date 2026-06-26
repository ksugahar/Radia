"""Generate the hollow-cylinder conductor tube fixture for the Volume-3D
stream-function coil designer golden test.

Mirrors the validated geometry of
``examples/vim/foliated_clebsch_solenoid.py`` (R0=0.05, R1=0.10, ZL=0.10,
maxh=0.03) and saves a FIXED .vol so the golden is reproducible across machines
(Netgen meshing is not byte-identical across versions/platforms; freezing the
.vol locks the mesh).

Run once to (re)create the fixture:
    python validation_test/panels/fixtures/generate_streamfunction_volume_tube.py
"""
import os

from netgen.csg import CSGeometry, Cylinder, Plane, Pnt, Vec
from ngsolve import Mesh, TaskManager

R0, R1, ZL, MAXH = 0.05, 0.10, 0.10, 0.03
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "streamfunction_volume_tube.vol")


def main():
    geo = CSGeometry()
    outer = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R1)
    inner = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R0)
    top = Plane(Pnt(0, 0, ZL), Vec(0, 0, 1))
    bot = Plane(Pnt(0, 0, -ZL), Vec(0, 0, -1))
    tube = (outer - inner) * top * bot
    geo.Add(tube.mat("cond"))
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=MAXH))
    mesh.ngmesh.Save(OUT)
    print(f"wrote {OUT}: ne={mesh.ne}, nv={mesh.nv}")


if __name__ == "__main__":
    main()
