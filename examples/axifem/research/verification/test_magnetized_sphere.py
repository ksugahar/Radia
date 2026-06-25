"""Permanent-magnet SOURCE term for the axisymmetric A_phi formulation,
validated against the analytical uniformly-magnetized linear sphere.

This is the FEMM `prob3big.cpp` magnetization line-integral (Hc * MagDir)
ported to the H1Henrotte axisymmetric FESpace as an NGSolve LinearForm.

Physics
-------
A sphere (radius a) of linear material (relative permeability mu_r,
coercivity Hc, axial magnetization) in vacuum has a UNIFORM internal
flux density and an external dipole field:

    B_in = 2 * mu0 * mu_r * Hc / (mu_r + 2)            (uniform, +z, inside)
    B_z(equator, r) = -(1/2) * B_in * (a / r)^3        (dipole, outside)

For mu_r = 2, Hc = 3e5 A/m  ->  B_in = 0.376991 T.

Magnet source (Galerkin weak form  RHS += int nu * B_rem . (curl v) dV)
----------------------------------------------------------------------
With nu*B_rem = Hc*(cos th, sin th) and th measured from the r-axis
(th = 90 deg => axial magnetization), in the calc_axisym_volumetric
A_phi convention (B_z = grad(u)[0] + u/r, B_r = -grad(u)[1]):

    f += [ Hc*sin(th)*(r*grad(v)[0] + v) - Hc*cos(th)*r*grad(v)[1] ] dx_magnet

Stiffness (static, sigma = 0):
    a += nu*(1/r)*(r*grad(u)[0]+u)*(r*grad(v)[0]+v)*dx + nu*r*grad(u)[1]*grad(v)[1]*dx

Reference: D. Meeker, FEMM 4.2 `fkn/prob3big.cpp` StaticAxisymmetric,
magnetization edge loop (the `H_c * (cos t * dr + sin t * dz)` contribution).
"""
import math
import numpy as np
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction, grad,
                     dx, x as r_cf, TaskManager, CoefficientFunction, Integrate)
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia.axifem import H1Henrotte

MU0 = 4e-7 * math.pi
A_SPHERE = 1.0
MU_R = 2.0
HC = 3.0e5
R_FAR = 50.0
B_IN_EXACT = 2.0 * MU0 * MU_R * HC / (MU_R + 2.0)   # 0.376991 T


def build_mesh(maxh_mag=0.05, maxh_air=0.5):
    """Half-disk magnet (r in [0, a]) embedded in an air box; r=0 -> 'axis'."""
    disk = WorkPlane().Circle(0, 0, A_SPHERE).Face()
    half = MoveTo(0, -R_FAR).Rectangle(R_FAR, 2 * R_FAR).Face()
    magnet = disk * half
    magnet.faces.name = "magnet"
    magnet.maxh = maxh_mag
    magnet.edges.Min(X).name = "axis"
    air = MoveTo(0, -R_FAR).Rectangle(R_FAR, 2 * R_FAR).Face() - magnet
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "outer"
    air.edges.Max(Y).name = "outer"
    air.edges.Min(Y).name = "outer"
    shape = Glue([air, magnet])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_air))


def add_magnet_source(f, v, mesh, Hc, theta_deg, region="magnet"):
    """FEMM prob3big.cpp magnetization source as an NGSolve LinearForm term."""
    th = math.radians(theta_deg)
    reg = mesh.MaterialCF({region: 1.0}, default=0.0)
    f += reg * Hc * math.sin(th) * (r_cf * grad(v)[0] + v) * dx
    f += -reg * Hc * math.cos(th) * (r_cf * grad(v)[1]) * dx


def solve_magnet(order=2, theta_deg=90.0):
    mesh = build_mesh()
    nu_cf = CoefficientFunction(
        1.0 / (MU0 * mesh.MaterialCF({"magnet": MU_R}, default=1.0)))
    fes = H1Henrotte(mesh, order=order, dirichlet="axis|outer")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * (1.0 / r_cf) * (r_cf * grad(u)[0] + u) * (r_cf * grad(v)[0] + v) * dx
    a += nu_cf * r_cf * grad(u)[1] * grad(v)[1] * dx
    f = LinearForm(fes)
    add_magnet_source(f, v, mesh, HC, theta_deg)
    gfu = GridFunction(fes)
    with TaskManager():
        a.Assemble(); f.Assemble()
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return mesh, gfu


def main():
    mesh, gfu = solve_magnet(order=2, theta_deg=90.0)
    print(f"Mesh: {mesh.nv} vertices, {mesh.ne} elements; ndof check done")
    B_z = grad(gfu)[0] + gfu / r_cf
    B_r = -grad(gfu)[1]

    # ---- internal field: volume-averaged B_z over the magnet (robust) ----
    vol = Integrate(r_cf, mesh, definedon=mesh.Materials("magnet"))
    bz_avg = Integrate(r_cf * B_z, mesh, definedon=mesh.Materials("magnet")) / vol
    err_in = (bz_avg - B_IN_EXACT) / B_IN_EXACT
    print(f"  <B_z>_magnet = {bz_avg:.6f} T  (exact {B_IN_EXACT:.6f})  err = {100*err_in:+.3f}%")

    # interior B_r must be ~0 (purely axial inside)
    br_max = max(abs(B_r(mesh(rr, zz)))
                 for rr in (0.3, 0.5, 0.7) for zz in (-0.3, 0.0, 0.3))
    print(f"  max |B_r| interior off-axis = {br_max:.3e} T (expect ~0)")

    # ---- external dipole at the equator (off-axis, reliable recovery) ----
    print("  external equatorial B_z vs dipole -(1/2)B_in(a/r)^3:")
    ext_ok = True
    for rr in (3.0, 5.0):
        bz = B_z(mesh(rr, 0.0))
        ref = -0.5 * B_IN_EXACT * (A_SPHERE / rr) ** 3
        rel = abs(bz - ref) / abs(ref)
        print(f"    r={rr:.1f}: B_z={bz:+.6e}  ref={ref:+.6e}  rel={100*rel:.2f}%")
        ext_ok = ext_ok and rel < 0.03

    assert abs(err_in) < 5e-3, f"internal B_z off by {100*err_in:.3f}% (>0.5%)"
    assert br_max < 1e-3, f"interior B_r should be ~0, got {br_max:.3e}"
    assert ext_ok, "external dipole field off by >3%"
    print("\n[OK] magnetized-sphere source validated on H1Henrotte "
          f"(<B_z> within {100*abs(err_in):.3f}% of analytical).")


if __name__ == "__main__":
    main()
