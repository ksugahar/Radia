"""COMSOL-class #17 -- HELMHOLTZ coil pair (uniform-field generation + uniformity).

Two coaxial N-turn loops (radius a) at z = +/- a/2 carrying the SAME current produce a
maximally uniform field at the centre (d^2B/dz^2 = 0 -- the defining Helmholtz
condition). The on-axis field:

    B_z(z) = (mu0 N I a^2/2)[((z-a/2)^2+a^2)^-3/2 + ((z+a/2)^2+a^2)^-3/2],
    centre  B0 = (4/5)^(3/2) mu0 N I / a = 0.7155 mu0 N I / a.

Validated outputs: the on-axis field vs the closed form (smooth -> clean L2 extraction),
the centre value vs (4/5)^(3/2) mu0 NI/a, and the FIELD UNIFORMITY (max |dB/B0| over the
central |z| <= a/4 working region) vs the analytic flatness -- the figure of merit for a
uniform-field magnet (NMR shim, sensor calibration, beam optics).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import math
from ngsolve import (Mesh, CoefficientFunction, TaskManager, L2, GridFunction,
                     grad, x as r_cf)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic, coil_pair_axial_field, NU0, MU0

a, N, I, w = 0.10, 100.0, 10.0, 0.004      # loop radius, turns, current, loop side [m]
s = a                                       # Helmholtz: separation = radius
Rbox, Zbox, nearR, nearZ = 0.80, 0.80, 0.20, 0.40
J_loop = N * I / (w * w)                     # azimuthal current density in each loop

# --- axisymmetric (r,z): two loop squares at (a, +/- a/2) + graded near/far air
def loop(zc):
    f = MoveTo(a - w / 2, zc - w / 2).Rectangle(w, w).Face(); f.maxh = 0.0015
    return f
top, bot = loop(+s / 2), loop(-s / 2)
top.faces.name = "ltop"; bot.faces.name = "lbot"
near = MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face() - top - bot
near.faces.name = "air"; near.maxh = 0.004
near.edges.Min(X).name = "axis"; near.edges.Min(X).maxh = 0.0006
far = MoveTo(0, -Zbox).Rectangle(Rbox, 2 * Zbox).Face() \
    - MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face()
far.faces.name = "air"; far.maxh = 0.06
far.edges.Min(X).name = "axis"
far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
mesh = Mesh(OCCGeometry(Glue([far, near, top, bot]), dim=2).GenerateMesh(maxh=0.06))
print(f"mesh: {mesh.ne} elems, materials={mesh.GetMaterials()}")

with TaskManager():
    Jr = mesh.MaterialCF({"ltop": J_loop, "lbot": J_loop}, default=0.0)
    gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0), Jr=Jr, order=2)
bzgf = GridFunction(L2(mesh, order=2)); bzgf.Set(grad(gfu)[0] + gfu / r_cf)
radii = [0.001, 0.0015, 0.002, 0.0025, 0.003]
# average over near-axis radii AND a small +/- z pair: avoids both the L2 inter-element
# jump (radii) and the exact-symmetry-plane sampling artifact at z=0 (the dz pair).
def bz0(z):
    return sum(bzgf(mesh(r, zz)) for r in radii for zz in (z - 2e-4, z + 2e-4)) / (2 * len(radii))

z_samples = [0.0, a / 8, a / 4, a / 2, 0.7 * a]
print(f"\n z/a     B_z_fem[T]   B_z_exact[T]   err%")
worst = 0.0
for z in z_samples:
    bf, be = bz0(z), coil_pair_axial_field(N, I, a, s, z)
    err = abs(bf - be) / abs(be) * 100
    worst = max(worst, err)
    print(f"{z/a:5.3f}   {bf:9.6f}   {be:9.6f}   {err:5.2f}")

B0 = bz0(0.0); B0_exact = (4.0 / 5.0) ** 1.5 * MU0 * N * I / a
# field uniformity over the central working region |z| <= a/4
zc = [k * a / 40 for k in range(-10, 11)]
unif_fem = max(abs(bz0(z) / B0 - 1.0) for z in zc) * 100
unif_ana = max(abs(b / coil_pair_axial_field(N, I, a, s, 0.0) - 1.0)
               for b in coil_pair_axial_field(N, I, a, s, zc)) * 100
print(f"\ncentre B0 = {B0:.6f} T  exact (4/5)^1.5 mu0 NI/a = {B0_exact:.6f}  "
      f"({abs(B0-B0_exact)/B0_exact*100:.2f}%)")
print(f"uniformity over |z|<=a/4:  FEM {unif_fem:.3f}%   analytic {unif_ana:.3f}%  "
      f"(Helmholtz flatness; FEM ~ analytic + numerical scatter)")
# the Helmholtz property = the field stays flat to <~1 % over the working volume;
# the analytic ideal is 0.42 %, the FEM adds sub-% numerical scatter on top.
ok = worst < 1.5 and abs(B0 - B0_exact) / B0_exact < 0.015 and unif_fem < 1.5
print("OK" if ok else f"OFF (field {worst:.2f}%, B0 {abs(B0-B0_exact)/B0_exact*100:.2f}%, unif {unif_fem:.2f}%)")
