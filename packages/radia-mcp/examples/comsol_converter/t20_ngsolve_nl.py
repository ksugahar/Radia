"""NGSolve NONLINEAR TEAM-20 (B-H steel) with the circular coil + azimuthal
current -- the relay partner for the COMSOL nonlinear solve. Same geometry as
t20_ngsolve.py but solve_magnetostatic_newton with the 38-pt B-H curve."""
import sys, time
sys.path.insert(0, r"W:\00_CAE\Radia\01_GitHub\packages\radia-mcp\src")
import numpy as np
from ngsolve import *
from netgen.occ import Box, Cylinder, Pnt, Dir, Glue, OCCGeometry
from netgen.meshing import MeshingParameters
from radia_mcp.radia_ngsolve.force import MU0
from radia_mcp.radia_ngsolve.solve import solve_magnetostatic_newton

NU0 = 1.0 / MU0
mm = 1e-3
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
yoke.faces.maxh = 7*mm; pole.faces.maxh = 4*mm; coil.faces.maxh = 6*mm; air.maxh = 50*mm
mesh = Mesh(OCCGeometry(Glue([yoke,pole,coil,air])).GenerateMesh(
    MeshingParameters(maxh=50*mm, grading=0.3))); mesh.Curve(2)
print(f"mesh: {mesh.ne} elems  ({time.time()-t0:.1f}s)", flush=True)
assert "pole_surface" in set(mesh.GetBoundaries())

# nonlinear steel B-H -> co-energy density phi(B)=int H dB
Bd = [0,0.01,0.025,0.05,0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.0,
      1.1,1.2,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,
      2.0,2.05,2.1,2.15,2.2,2.25,2.3]
Hd = [0,27,58,100,153,185,205,233,255,285,320,355,405,470,555,673,836,1065,
      1220,1420,1720,2130,2670,3480,4500,5950,7650,10100,13000,15900,21100,
      26300,32900,42700,61700,84300,110000,135000]
Ms = 2.16
for bb in (2.4,2.6,2.8,3.0,3.5,4.0): Bd.append(bb); Hd.append((bb-Ms)/MU0)
energy_density = BSpline(2, [0]+list(Bd), list(Hd)).Integrate()

# azimuthal coil current
rrel = sqrt((x-cxc)**2 + (y-cyc)**2 + 1e-12)
jx, jy = -(y-cyc)/rrel*J0, (x-cxc)/rrel*J0
Jv = mesh.MaterialCF({"coil": CoefficientFunction((jx,jy,0.0))},
                     default=CoefficientFunction((0.0,0.0,0.0)))

t1 = time.time()
gfu = solve_magnetostatic_newton(mesh, Jv, energy_density, "yoke|pole",
                                 load_steps=3, max_iter=25, tol=1e-7)
B = curl(gfu)
print(f"newton solve {time.time()-t1:.1f}s", flush=True)

bp = sqrt(InnerProduct(B,B))(mesh(cxc,cyc,75*mm))
fl = L2(mesh, order=1, definedon=mesh.Materials('air'))
gB = [GridFunction(fl) for _ in range(3)]
for k in range(3): gB[k].Set(B[k])
Ba = CoefficientFunction(tuple(BoundaryFromVolumeCF(gB[k]) for k in range(3)))
n = specialcf.normal(3)
Tz = (1.0/MU0)*(Ba[2]*InnerProduct(Ba,n) - 0.5*InnerProduct(Ba,Ba)*n[2])
Fz = -Integrate(Tz*ds(definedon=mesh.Boundaries('pole_surface')), mesh)
print(f"NGSolve TEAM20 (NONLINEAR B-H, circular coil): Fz={Fz:.4f} N  |B|@pole={bp:.4f} T", flush=True)
