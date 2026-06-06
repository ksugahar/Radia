"""COMSOL-class #11 -- CYLINDER magnetic shielding (2D transverse field).

The 2D / infinite-length twin of #5 (sphere): a mu-metal cylindrical SHELL (inner radius
a, outer b) in a uniform TRANSVERSE field B0. Exact shielding factor S = B0/B_cavity:

    S = [ (mu_r+1)^2 - (a/b)^2 (mu_r-1)^2 ] / (4 mu_r)

Method: 2D planar A_z, scattered field. The permeability jump is a MAGNETIZATION source
M = nu_pert * B0 on the shell (nu_pert = nu0(1-1/mu_r)) -- reusing the magnet_cf hook of
solve_planar_magnetostatic. Total B = B0 + curl(A_z); cavity field read as a VOLUME
AVERAGE. shell_shielding_factor(.., "cylinder") is the closed form.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, dx, Integrate
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           shell_shielding_factor, NU0)

B0, a, b, Rbox = 1.0, 0.04, 0.06, 0.60     # box at 10*b: 2D dipole decays ~1/r, needs room

geo = SplineGeometry()
geo.AddCircle((0, 0), Rbox, leftdomain=3, rightdomain=0, bc="outer")
geo.AddCircle((0, 0), b, leftdomain=2, rightdomain=3)
geo.AddCircle((0, 0), a, leftdomain=1, rightdomain=2)
geo.SetMaterial(1, "bore"); geo.SetMaterial(2, "shell"); geo.SetMaterial(3, "air")
mesh = Mesh(geo.GenerateMesh(maxh=0.01))
vol = Integrate(CoefficientFunction(1.0) * dx("bore"), mesh)
print(f"mesh: {mesh.ne} elems")

print("mu_r    S_fem    S_exact    err%")
worst = 0.0
for mu_r in (50.0, 100.0, 200.0):
    nu_cf = mesh.MaterialCF({"shell": NU0 / mu_r}, default=NU0)
    ind = mesh.MaterialCF({"shell": 1.0}, default=0.0)
    magnet_cf = ind * (NU0 * (1.0 - 1.0 / mu_r)) * CoefficientFunction((B0, 0.0))
    gfu = solve_planar_magnetostatic(mesh, nu_cf, magnet_cf=magnet_cf, order=3)
    Bx_in = (Integrate((grad(gfu)[1] + B0) * dx("bore"), mesh)) / vol   # B_x = dA/dy + B0
    S_fem = B0 / Bx_in
    Sx = shell_shielding_factor(mu_r, a, b, "cylinder")
    err = abs(S_fem - Sx) / Sx * 100; worst = max(worst, err)
    print(f"{mu_r:5.0f}   {S_fem:7.2f}   {Sx:7.2f}   {err:5.1f}")
print("OK" if worst < 5 else f"OFF (worst {worst:.1f}%)")
