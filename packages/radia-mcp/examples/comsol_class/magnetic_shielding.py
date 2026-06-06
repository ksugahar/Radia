"""COMSOL-class #5 -- MAGNETIC SHIELDING (AC/DC blog "How to Model Magnetic Shielding").

A mu-metal SPHERICAL SHELL (inner radius a, outer b, relative permeability mu_r) in a
uniform applied field B0. The shell reduces the field inside the cavity by the
shielding factor S = B0 / B_cavity. Exact analytic reference (uniform cavity field):

    B_in/B0 = 9 mu_r / [ (2 mu_r+1)(mu_r+2) - 2 (a/b)^3 (mu_r-1)^2 ]

Method: scattered-field A-form magnetostatics driven by the permeability jump (no
coil) -- radia_ngsolve.solve.solve_scattered_uniform_field. The cavity field is read
as a VOLUME AVERAGE (robust vs point eval). Reproduces the COMSOL parametric
shielding-factor-vs-mu_r study and cross-checks each point against the closed form.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, curl, dx, Integrate
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt
from radia_mcp.radia_ngsolve.solve import (solve_scattered_uniform_field,
                                           shell_shielding_factor)

a, b, L, B0 = 0.04, 0.06, 0.30, 1.0          # cavity r<a, shell a<r<b, box half-len L

sph_a = Sphere(Pnt(0, 0, 0), a); sph_b = Sphere(Pnt(0, 0, 0), b)
box = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
geo = CSGeometry()
geo.Add(sph_a.mat("cavity").maxh(0.012))
geo.Add((sph_b - sph_a).mat("shell").maxh(0.006))
geo.Add((box - sph_b).mat("air"))
mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.3)); mesh.Curve(2)
vol_cav = Integrate(CoefficientFunction(1.0) * dx("cavity"), mesh)
print(f"mesh: {mesh.ne} elems")

print("mu_r    S_fem    S_exact    err%")
worst = 0.0
for mu_r in (50.0, 100.0, 200.0):
    gfA = solve_scattered_uniform_field(mesh, {"shell": mu_r}, (0.0, 0.0, B0))
    Bz_in = Integrate((curl(gfA)[2] + B0) * dx("cavity"), mesh) / vol_cav
    S_fem = B0 / Bz_in
    Sx = shell_shielding_factor(mu_r, a, b, "sphere")
    err = abs(S_fem - Sx) / Sx * 100; worst = max(worst, err)
    print(f"{mu_r:5.0f}   {S_fem:7.2f}   {Sx:7.2f}   {err:5.1f}")
print("OK" if worst < 8 else f"OFF (worst {worst:.1f}%)")
