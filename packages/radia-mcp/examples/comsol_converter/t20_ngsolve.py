"""NGSolve TEAM-20 lifting force, LINEAR steel mur=1000, circular
(cylindrical-annulus) coil with an azimuthal current density. Cross-check:
Fz and |B|@pole vs a stored independent reference
(Fz=-7.76 N, |B|@pole=0.681 T)."""
import os, sys, time
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src")))
import numpy as np
from ngsolve import *
from netgen.occ import Box, Cylinder, Pnt, Dir, Glue, OCCGeometry
from netgen.meshing import MeshingParameters
from radia_mcp.radia_ngsolve.force import MU0
from radia_mcp.radia_ngsolve.solve import solve_magnetostatic_Aform

NU0 = 1.0 / MU0
mm = 1e-3
mur = 1000.0
cxc, cyc = 63.5*mm, 12.5*mm
J0 = 575109.0
e = 1e-4

def _yoke():
    return (Box(Pnt(0,0,0), Pnt(127*mm,25*mm,150*mm))
            - Box(Pnt(25*mm,-e,25*mm), Pnt(102*mm,25*mm+e,125*mm))
            - Box(Pnt(50*mm,-e,125*mm), Pnt(77*mm,25*mm+e,150*mm)))
def _pole():
    return Box(Pnt(51*mm,7.5*mm,26.5*mm), Pnt(76*mm,17.5*mm,125*mm))
def _coil():
    co = Cylinder(Pnt(cxc,cyc,26.7*mm), Dir(0,0,1), r=37.5*mm, h=96.6*mm)
    ci = Cylinder(Pnt(cxc,cyc,26.6*mm), Dir(0,0,1), r=19.5*mm, h=96.8*mm)
    return co - ci

t0 = time.time()
yoke = _yoke(); yoke.mat('yoke')
pole = _pole(); pole.mat('pole')
for f in pole.faces: f.name = 'pole_surface'
coil = _coil(); coil.mat('coil')
big = Box(Pnt(-80*mm,-105*mm,-80*mm), Pnt(207*mm,130*mm,230*mm))
for f in big.faces: f.name = 'outer'
air = big - _yoke() - _pole() - _coil(); air.mat('air')
yoke.faces.maxh = 6*mm; pole.faces.maxh = 4*mm; coil.faces.maxh = 5*mm; air.maxh = 45*mm
mesh = Mesh(OCCGeometry(Glue([yoke,pole,coil,air])).GenerateMesh(
    MeshingParameters(maxh=45*mm, grading=0.3))); mesh.Curve(2)
print(f"mesh: {mesh.ne} elems  mats={mesh.GetMaterials()}  ({time.time()-t0:.1f}s)")
assert "pole_surface" in set(mesh.GetBoundaries())

nu = mesh.MaterialCF({'yoke':NU0/mur, 'pole':NU0/mur}, default=NU0)
rrel = sqrt((x-cxc)**2 + (y-cyc)**2 + 1e-12)
Jcoil = CoefficientFunction((-(y-cyc)/rrel*J0, (x-cxc)/rrel*J0, 0.0))
J = mesh.MaterialCF({'coil': Jcoil}, default=CoefficientFunction((0.0,0.0,0.0)))

t1 = time.time()
gfu = solve_magnetostatic_Aform(mesh, nu, source=J, order=2, dirichlet='outer', reg=1e-6)
B = curl(gfu)
print(f"solve {time.time()-t1:.1f}s")

bp = sqrt(InnerProduct(B,B))(mesh(cxc,cyc,75*mm))
fl = L2(mesh, order=1, definedon=mesh.Materials('air'))
gB = [GridFunction(fl) for _ in range(3)]
for k in range(3): gB[k].Set(B[k])
Ba = CoefficientFunction(tuple(BoundaryFromVolumeCF(gB[k]) for k in range(3)))
n = specialcf.normal(3)
Tz = (1.0/MU0)*(Ba[2]*InnerProduct(Ba,n) - 0.5*InnerProduct(Ba,Ba)*n[2])
Fz = -Integrate(Tz*ds(definedon=mesh.Boundaries('pole_surface')), mesh)
print(f"NGSolve TEAM20 (linear mur=1000, circular coil): Fz={Fz:.4f} N  |B|@pole={bp:.4f} T")
print(f"stored reference: Fz=-7.76 N  |B|@pole=0.681 T")
