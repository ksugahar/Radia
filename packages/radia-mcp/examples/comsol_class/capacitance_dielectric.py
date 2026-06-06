"""COMSOL-class #7 -- DIELECTRIC capacitance (AC/DC "Computing Capacitance" with
dielectric materials). Spherical capacitor whose gap is filled with TWO concentric
dielectric shells (eps_r1 for a<r<c, eps_r2 for c<r<b) -- series dielectrics.

Exact (spherical symmetry, no fringing):
    C = 4 pi eps0 / [ (1/eps_r1)(1/a - 1/c) + (1/eps_r2)(1/c - 1/b) ]

Validates the piecewise-permittivity (mesh.MaterialCF) path of
radia_ngsolve.electrostatic3d -- single uniform dielectric is the c->b limit
(C = 4 pi eps0 eps_r a b/(b-a)).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh
from netgen.csg import CSGeometry, Sphere, Pnt
from radia_mcp.radia_ngsolve.electrostatic3d import (
    solve_electrostatic_3d, capacitance_from_energy, EPS0)

V0 = 1.0
print("a,c,b[mm]  er1 er2   C_fem[pF]  C_exact[pF]  err%")
worst = 0.0
for a, c, b, er1, er2 in ((0.04, 0.06, 0.10, 2.0, 5.0),
                          (0.04, 0.06, 0.10, 1.0, 1.0),     # -> uniform vacuum check
                          (0.03, 0.05, 0.09, 4.0, 2.0)):
    sph_a = Sphere(Pnt(0, 0, 0), a).bc("inner")
    sph_c = Sphere(Pnt(0, 0, 0), c)
    sph_b = Sphere(Pnt(0, 0, 0), b).bc("outer")
    geo = CSGeometry()
    geo.Add((sph_c - sph_a).mat("diel1").maxh(0.008))
    geo.Add((sph_b - sph_c).mat("diel2").maxh(0.010))
    mesh = Mesh(geo.GenerateMesh(maxh=0.012, grading=0.3)); mesh.Curve(3)

    eps_cf = mesh.MaterialCF({"diel1": EPS0 * er1, "diel2": EPS0 * er2})
    gfV = solve_electrostatic_3d(mesh, eps_cf, {"inner": V0, "outer": 0.0}, order=2)
    C_fem = capacitance_from_energy(gfV, eps_cf, V0, mesh)
    C_ref = 4 * math.pi * EPS0 / ((1/er1)*(1/a - 1/c) + (1/er2)*(1/c - 1/b))
    err = abs(C_fem - C_ref) / C_ref * 100; worst = max(worst, err)
    print(f"{a*1e3:.0f},{c*1e3:.0f},{b*1e3:.0f}    {er1:.0f}   {er2:.0f}    "
          f"{C_fem*1e12:7.3f}    {C_ref*1e12:7.3f}    {err:5.1f}")
print("OK" if worst < 3 else f"OFF (worst {worst:.1f}%)")
