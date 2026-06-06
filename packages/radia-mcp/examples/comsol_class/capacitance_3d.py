"""COMSOL-class #6 -- 3D CAPACITANCE (AC/DC blog "Computing Capacitance").

Spherical capacitor: inner conductor radius a at V0, outer conductor radius b
grounded, vacuum gap. Exact:  C = 4 pi eps0 a b / (b - a).

Method: 3D electrostatics on H1 with the conductors as Dirichlet boundaries, then
capacitance from the field ENERGY  C = 2W/V0^2 (radia_ngsolve.electrostatic3d) --
the same energy route COMSOL uses for the capacitance matrix. Sweeps the gap ratio
and cross-checks every point against the closed form.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh
from netgen.csg import CSGeometry, Sphere, Pnt
from radia_mcp.radia_ngsolve.electrostatic3d import (
    solve_electrostatic_3d, capacitance_from_energy, spherical_capacitor_C, EPS0)

V0 = 1.0
print("a[mm] b[mm]  C_fem[pF]  C_exact[pF]  err%")
worst = 0.0
for a, b in ((0.05, 0.10), (0.05, 0.075), (0.04, 0.12)):
    sph_a = Sphere(Pnt(0, 0, 0), a).bc("inner")
    sph_b = Sphere(Pnt(0, 0, 0), b).bc("outer")
    geo = CSGeometry(); geo.Add((sph_b - sph_a).mat("gap").maxh(0.01))
    mesh = Mesh(geo.GenerateMesh(maxh=0.012, grading=0.3)); mesh.Curve(3)

    gfV = solve_electrostatic_3d(mesh, EPS0, {"inner": V0, "outer": 0.0}, order=2)
    C_fem = capacitance_from_energy(gfV, EPS0, V0, mesh)
    C_ref = spherical_capacitor_C(a, b)
    err = abs(C_fem - C_ref) / C_ref * 100; worst = max(worst, err)
    print(f"{a*1e3:4.0f} {b*1e3:4.0f}   {C_fem*1e12:7.3f}    {C_ref*1e12:7.3f}    {err:5.1f}")
print("OK" if worst < 3 else f"OFF (worst {worst:.1f}%)")
