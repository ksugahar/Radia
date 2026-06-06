"""COMSOL-class #15 -- finite air-core SOLENOID (on-axis field + self-inductance).

The air-cored-coil design workhorse. A solenoid (radius a, length L, n turns/m,
current I) modelled as an axisymmetric azimuthal current sheet K = n I (winding
current density J_phi = n I / t over a thin annulus). Two validated outputs:

  * ON-AXIS B_z(z) vs the exact closed form (smooth on axis -- no pole corner, so
    the L2-projection / near-axis radial-mean extraction recovers it cleanly):
        B_z(z) = (mu0 n I/2)[(L/2-z)/sqrt((L/2-z)^2+a^2) + (L/2+z)/sqrt((L/2+z)^2+a^2)]
  * SELF-INDUCTANCE via the field ENERGY (robust, no axis sampling):
        L_ind = 2 W / I^2,  W = 2 pi int nu/2 |B|^2 r dr dz   (inductance_axi)
    The ratio L_ind / L_long (L_long = mu0 n^2 pi a^2 L) is the Nagaoka coefficient
    k_N(2a/L) < 1; checked against Wheeler's continuous formula k_N ~ 1/(1+0.9 a/L).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import math
from ngsolve import (Mesh, CoefficientFunction, TaskManager, L2, GridFunction,
                     grad, x as r_cf)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_axi_magnetostatic,
    solenoid_axial_field, solenoid_inductance_long, NU0)
from radia_mcp.radia_ngsolve.force import inductance_axi

a, L, t = 0.05, 0.30, 0.001              # coil radius, length, winding thickness [m]
n, I = 1000.0, 2.0                        # turns/m, current [A]
Rbox, Zbox = 1.20, 1.20                   # big box: capture the external flux ENERGY (-> L)
J_phi = n * I / t                         # azimuthal winding current density [A/m^2]

# --- axisymmetric (r,z) geometry: winding annulus + graded near/far air -------
coil = MoveTo(a - t / 2, -L / 2).Rectangle(t, L).Face()
coil.faces.name = "coil"; coil.maxh = 0.0008
nearR, nearZ = 0.12, 0.50
near = MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face() - coil
near.faces.name = "air"; near.maxh = 0.004
near.edges.Min(X).name = "axis"; near.edges.Min(X).maxh = 0.001
far = MoveTo(0, -Zbox).Rectangle(Rbox, 2 * Zbox).Face() \
    - MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face()
far.faces.name = "air"; far.maxh = 0.06
far.edges.Min(X).name = "axis"
far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
mesh = Mesh(OCCGeometry(Glue([far, near, coil]), dim=2).GenerateMesh(maxh=0.05))
print(f"mesh: {mesh.ne} elems, materials={mesh.GetMaterials()}")

# --- solve: azimuthal coil current J_phi ; air-core (nu = NU0 everywhere) ------
with TaskManager():
    Jr = mesh.MaterialCF({"coil": J_phi}, default=0.0)
    gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0), Jr=Jr, order=2)
B = CoefficientFunction((grad(gfu)[0] + gfu / r_cf, -grad(gfu)[1]))   # (B_z, B_r)

# --- on-axis B_z: L2-project then near-axis radial mean (the #14 recipe) -------
bzgf = GridFunction(L2(mesh, order=2)); bzgf.Set(grad(gfu)[0] + gfu / r_cf)
radii = [0.001, 0.0015, 0.002, 0.0025, 0.003]
def bz0(z): return sum(bzgf(mesh(r, z)) for r in radii) / len(radii)

z_samples = [0.0, L / 4, 0.45 * L, 0.5 * L, 0.6 * L, 0.8 * L, L]
print(f"\n z/L     B_z_fem[T]   B_z_exact[T]   err%")
worst = 0.0
for z in z_samples:
    bf, be = bz0(z), solenoid_axial_field(n, I, a, L, z)
    err = abs(bf - be) / abs(be) * 100 if be else 0.0
    worst = max(worst, err)
    print(f"{z/L:5.2f}   {bf:9.6f}   {be:9.6f}   {err:5.2f}")

# --- self-inductance via field energy (robust; per-turn current I) ------------
L_fem = inductance_axi(B, mesh, NU0, I)
L_long = solenoid_inductance_long(n, a, L)
kN_fem = L_fem / L_long
kN_wheeler = 1.0 / (1.0 + 0.9 * a / L)            # Wheeler continuous (~1 % accurate)
print(f"\nL_fem      = {L_fem*1e3:.4f} mH")
print(f"L_long     = {L_long*1e3:.4f} mH  (mu0 n^2 pi a^2 L, Nagaoka->1 limit)")
print(f"k_Nagaoka  = L_fem/L_long = {kN_fem:.4f}   Wheeler 1/(1+0.9a/L) = {kN_wheeler:.4f}")
# interior (z <= 0.6 L) is the design-relevant bore field; the far external field
# (z >= 0.8 L) is tiny and mesh-sensitive (as in #14), reported but not gated.
interior_worst = max(abs(bz0(z) - solenoid_axial_field(n, I, a, L, z))
                     / abs(solenoid_axial_field(n, I, a, L, z)) * 100
                     for z in (0.0, L / 4, 0.45 * L, 0.5 * L, 0.6 * L))
b0_err = abs(bz0(0.0) - solenoid_axial_field(n, I, a, L, 0.0)) / solenoid_axial_field(n, I, a, L, 0.0) * 100
kN_err = abs(kN_fem - kN_wheeler) / kN_wheeler * 100
print(f"\ncentre-field err {b0_err:.2f}% ; bore (z<=0.6L) worst {interior_worst:.2f}% ; "
      f"full-profile worst {worst:.2f}% (far field tiny) ; inductance vs Wheeler {kN_err:.2f}%")
print("OK" if (interior_worst < 3.0 and kN_err < 3.0) else f"OFF (bore {interior_worst:.2f}%, L {kN_err:.2f}%)")
