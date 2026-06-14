"""Generate cube_5mm.vol for testing."""
from netgen.occ import Box, OCCGeometry, Pnt
from ngsolve import Mesh

L = 5e-3
box = Box(Pnt(0, 0, 0), Pnt(L, L, L))
box.mat("Cu")
for f in box.faces:
    f.name = "outer"
geo = OCCGeometry(box)
ng_mesh = geo.GenerateMesh(maxh=L / 6)
mesh = Mesh(ng_mesh)
print(f"ne = {mesh.ne}, nv = {mesh.nv}")
ng_mesh.Save("cube_5mm.vol")
print("Saved cube_5mm.vol")
