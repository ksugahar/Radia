"""COMSOL-class #8 -- HALBACH CYLINDER (accelerator / undulator dipole magnet).

An ideal k=1 (dipole) Halbach permanent-magnet cylinder: remanence Br with the easy
axis rotating as 2*phi over the annulus Ri<r<Ro produces a UNIFORM transverse field in
the bore (and ZERO field outside). Exact (2D / infinite length):

    B_bore = Br * ln(Ro / Ri)

Method: 2D planar A_z magnetostatics (solve_planar_magnetostatic with magnet_cf) and
the rotating magnetization CF from halbach_dipole_magnetization. Validates |B| at
several bore points (uniformity) against the closed form. radia's home turf.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, sqrt
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           halbach_dipole_magnetization,
                                           halbach_bore_field, NU0)

Br, Ri, Ro, Rbox = 1.2, 0.02, 0.04, 0.12      # NdFeB-ish; bore 20 mm, magnet OD 40 mm

geo = SplineGeometry()
geo.AddCircle((0, 0), Rbox, leftdomain=3, rightdomain=0, bc="outer")
geo.AddCircle((0, 0), Ro, leftdomain=2, rightdomain=3)
geo.AddCircle((0, 0), Ri, leftdomain=1, rightdomain=2)
geo.SetMaterial(1, "bore"); geo.SetMaterial(2, "magnet"); geo.SetMaterial(3, "air")
mesh = Mesh(geo.GenerateMesh(maxh=0.004))
print(f"mesh: {mesh.ne} elems")

M_cf = halbach_dipole_magnetization(mesh, Br, region="magnet")
gfu = solve_planar_magnetostatic(mesh, NU0, magnet_cf=M_cf, order=3, dirichlet="outer")
B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))           # (Bx, By)
Bmag = sqrt(B[0] * B[0] + B[1] * B[1])

B_ref = halbach_bore_field(Br, Ri, Ro)
print(f"B_exact = Br ln(Ro/Ri) = {B_ref:.4f} T")
print("bore point (mm)    |B|[T]    err%")
worst = 0.0
for px, py in ((0, 0), (0.01, 0), (0, 0.01), (0.007, 0.007), (-0.008, 0.004)):
    val = float(Bmag(mesh(px, py)))
    err = abs(val - B_ref) / B_ref * 100; worst = max(worst, err)
    print(f"  ({px*1e3:6.1f},{py*1e3:6.1f})   {val:.4f}   {err:5.1f}")
print("OK" if worst < 5 else f"OFF (worst {worst:.1f}%)")
