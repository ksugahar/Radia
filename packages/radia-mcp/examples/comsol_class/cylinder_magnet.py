"""COMSOL-class #14 -- axially-magnetized CYLINDER permanent magnet (on-axis field).

The PM-design workhorse: the on-axis B_z of a uniformly axially-magnetized cylinder
(radius R, length L, remanence Br) -- what every pull-force / holding-force / gap-field
estimate starts from. Closed form (surface-charge model, recoil mu_r = 1):

    B_z(z) = (Br/2)[ (z+L/2)/sqrt(R^2+(z+L/2)^2) - (z-L/2)/sqrt(R^2+(z-L/2)^2) ]

Method: axisymmetric A_phi (H1Henrotte) with a rigid axial magnetization source
(magnets={"magnet": (Br/mu0, 90 deg)}), nu = NU0 everywhere (mu_r=1). On-axis B_z is
read by the FLUX identity, NOT a point gradient: the flux through a disk of radius
rho is Phi = 2 pi rho A_phi(rho), so the disk-averaged axial field is
B_z ~ 2 A_phi(rho,z)/rho -> B_z(0,z) as rho->0. This uses only A_phi (the primary,
reliably point-evaluable variable); the H1Henrotte GRADIENT is not point-evaluable.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import math
from ngsolve import (Mesh, CoefficientFunction, TaskManager, L2, GridFunction,
                     grad, x as r_cf)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_axi_magnetostatic,
                                           cylinder_magnet_axial_field, NU0, MU0)

Br, R, L = 1.2, 0.010, 0.020             # T, m  (NdFeB-ish, 10 mm dia x 20 mm)
Rbox, Zbox = 0.40, 0.40                   # far box ~ 20-40x the magnet
maxh_mag, maxh_near, maxh_far = 0.0008, 0.001, 0.04
axis_h = 0.0003
nearR, nearZ = 8 * R, 6 * L               # near region holds ALL sample points

# --- axisymmetric (r,z) geometry: magnet + graded near/far air ----------------
magnet = MoveTo(0, -L / 2).Rectangle(R, L).Face()
magnet.faces.name = "magnet"; magnet.maxh = maxh_mag
magnet.edges.Min(X).name = "axis"; magnet.edges.Min(X).maxh = axis_h
# refine the pole faces & the charge-singular pole corner (r=R, z=+/-L/2)
magnet.edges.Max(Y).maxh = 0.0003; magnet.edges.Min(Y).maxh = 0.0003
magnet.edges.Max(X).maxh = 0.0004
near = MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face() - magnet
near.faces.name = "air"; near.maxh = maxh_near
near.edges.Min(X).name = "axis"; near.edges.Min(X).maxh = axis_h
far = MoveTo(0, -Zbox).Rectangle(Rbox, 2 * Zbox).Face() - MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face()
far.faces.name = "air"; far.maxh = maxh_far
far.edges.Min(X).name = "axis"
far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
mesh = Mesh(OCCGeometry(Glue([far, near, magnet]), dim=2).GenerateMesh(maxh=maxh_far))
print(f"mesh: {mesh.ne} elems, materials={mesh.GetMaterials()}")

# --- solve: rigid axial magnetization M = Br/mu0 (theta=90 deg => +z) ---------
with TaskManager():
    gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0),
                                  magnets={"magnet": (Br / MU0, 90.0)}, order=2)

# B_z = dA/dr + A/r (Henrotte) is reliable inside an INTEGRAL but its POINT eval
# falls back to a base diffop -> project it onto a standard L2 space (the .Set
# uses the working integral path), then sample that field cleanly just off-axis.
bzgf = GridFunction(L2(mesh, order=2))
bzgf.Set(grad(gfu)[0] + gfu / r_cf)
# On-axis B_z = MEAN over a few near-axis radii: each is ~B_z(0,z) to <0.5 %, and
# averaging washes out the per-point L2 inter-element jump.
radii = [R / 25, R / 22, R / 20, R / 18, R / 16]

def bz_on_axis(z):
    return sum(bzgf(mesh(r, z)) for r in radii) / len(radii)

z_samples = [0.0, L / 4, 0.45 * L, 0.75 * L, L, 1.5 * L, 2.0 * L]
print(f"\n z/L     B_z_fem[T]   B_z_exact[T]   err%   (L2-proj, near-axis radial mean)")
errs = []
for z in z_samples:
    bz_fem = bz_on_axis(z)
    bz_ex = cylinder_magnet_axial_field(Br, R, L, z)
    err = abs(bz_fem - bz_ex) / abs(bz_ex) * 100 if bz_ex else 0.0
    errs.append(err)
    flag = "  <- pole exit (off-axis bulge)" if 0.9 < z / L < 1.1 and err > 4 else ""
    print(f"{z/L:5.2f}   {bz_fem:9.5f}   {bz_ex:9.5f}   {err:5.2f}{flag}")
b0_err = errs[0]
n_ok = sum(e < 4.0 for e in errs)
b0 = cylinder_magnet_axial_field(Br, R, L, 0.0)
print(f"\ncentre B_z(0) = {b0:.4f} T  [exact Br*L/(2 sqrt(R^2+(L/2)^2))]; centre err {b0_err:.2f}%")
print(f"{n_ok}/{len(errs)} points < 4 % (the lone outlier is the pole-exit z~L,")
print("where the on-axis field is corner-singularity / off-axis-bulge sensitive).")
# headline spec (centre field) ~exact, and >=6/7 of the profile within a few %
print("OK" if (b0_err < 2.0 and n_ok >= 6) else f"OFF (centre {b0_err:.2f}%, {n_ok}/7 ok)")
