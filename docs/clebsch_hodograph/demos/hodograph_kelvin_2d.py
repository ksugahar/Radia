"""Hodograph + Kelvin in 2-D Cartesian: an EXACT open boundary, no air box.

RESEARCH example (track A, rung 1 of "Kelvin in the hodograph").  The dipole
FEM rung truncates the exterior with a Dirichlet AIR BOX; this replaces it with
the **2-D Cartesian Kelvin transformation**, which is especially clean:

  In 2-D the Dirichlet energy INT |grad u|^2 is CONFORMALLY INVARIANT, and the
  Kelvin inversion z -> R^2 / z-bar is conformal -- so the Kelvin exterior is
  WEIGHT-FREE (mu' = mu0, no (R/rho')^2 factor, unlike the axisymmetric case in
  hodograph_kelvin_axisym.py).  The open boundary at infinity maps to a single
  interior point (GND), and the field is exact.

This is the classical complex-potential pole-design picture made into an FEM
realisation with an exact open boundary.  In a current-free 2-D region the flux
function A_z and the scalar potential V are CAUCHY-RIEMANN CONJUGATES:
``W = A_z + i mu0 V`` is analytic, so the orthogonal nets {A_z = const} (flux
lines) and {V = const} (equipotentials) ARE the hodograph grid.  The Kelvin
inversion is conformal, so W extends analytically across the open boundary --
"Kelvin in the hodograph" is exactly analytic continuation of W.

Geometry (netgen.occ 2-D): an inner disk (magnetisable cylinder + air, radius
R_K) glued to an offset Kelvin exterior disk (radius R_K), the two boundary
circles identified periodically (kelvin_int <-> kelvin_ext), GND at the Kelvin
centre = the image of infinity.  Excitation: a uniform applied field H0 y-hat
(reduced scalar potential; the Kelvin exterior carries the 2-D reduced-potential
background -F_s, via radia.kelvin_material.make_reduced_potential_background_cf
with dim=2).

Verified (mu_r = 100, cylinder a = 0.3, Kelvin R = 1.0, order 3, maxh 0.04):
  interior B = 2 mu_r/(mu_r+1) B0,  **field_error ~2e-8** (machine precision)
      -- vs ~3e-3 for a truncated air box at r/a = 6 (the Kelvin win);
  Bx interior ~1e-17 (no transverse field by symmetry);
  hodograph consistency B(from A_z) vs B(from V) ~1e-4 (air, off the cylinder).

run:  python hodograph_kelvin_2d.py
"""
import math
import os

from numpy import pi
from ngsolve import (Mesh, H1, Periodic, GridFunction, grad, InnerProduct, dx,
                     CF, x, y, IfPos, BilinearForm, LinearForm, TaskManager,
                     Integrate)
from netgen.occ import (Circle, Vertex, Pnt, Glue, OCCGeometry,
                        IdentificationType)
from radia.kelvin_material import make_reduced_potential_background_cf

MU0 = 4 * pi * 1e-7


def _kelvin_geometry(a, R_K, offset):
    """Inner disk (cylinder + air) + offset Kelvin disk, periodic circles."""
    iron = Circle((0, 0), a).Face()
    iron.faces.name = "magnetic"
    for e in iron.edges:
        e.name = "iron_bnd"
    innerdisk = Circle((0, 0), R_K).Face()
    for e in innerdisk.edges:                    # name BEFORE the boolean (a full
        e.name = "kelvin_int"                    # circle's edge-centroid is the origin)
    air_inner = innerdisk - iron
    air_inner.faces.name = "air_inner"
    outerdisk = Circle((offset, 0), R_K).Face()
    outerdisk.faces.name = "air_outer"
    for e in outerdisk.edges:
        e.name = "kelvin_ext"
    gnd = Vertex(Pnt(offset, 0, 0))
    gnd.name = "GND"
    shape = Glue([air_inner, iron, outerdisk, gnd])
    for ie in [e for e in shape.edges if e.name == "kelvin_int"]:
        for ee in [e for e in shape.edges if e.name == "kelvin_ext"]:
            ie.Identify(ee, "kperiodic", IdentificationType.PERIODIC)
    return OCCGeometry(shape, dim=2)


def _truncated_geometry(a, R_trunc):
    """Plain air box (Dirichlet) at radius R_trunc -- the thing Kelvin replaces."""
    iron = Circle((0, 0), a).Face()
    iron.faces.name = "magnetic"
    for e in iron.edges:
        e.name = "iron_bnd"
    disk = Circle((0, 0), R_trunc).Face()
    for e in disk.edges:
        e.name = "outer"
    air = disk - iron
    air.faces.name = "air_inner"
    return OCCGeometry(Glue([air, iron]), dim=2)


def _interior_B(mesh, B):
    m = mesh.MaterialCF({"magnetic": 1.0}, default=0.0)
    area = Integrate(m, mesh)
    return (Integrate(m * B[0], mesh) / area, Integrate(m * B[1], mesh) / area)


def solve(mu_r=100.0, order=3, maxh=0.04, a=0.3, R_K=1.0, offset=3.0, H0=1.0,
          x_eval=0.45, plot=False):
    """Reduced-Omega + 2-D Kelvin solve, the hodograph net, and the air-box
    contrast.  Returns the field error (exact open boundary) + hodograph
    consistency."""
    with TaskManager():
        mesh = Mesh(_kelvin_geometry(a, R_K, offset).GenerateMesh(maxh=maxh))
        mesh.Curve(order)
        Mu = mesh.MaterialCF({"magnetic": mu_r * MU0}, default=MU0)
        # uniform applied H0 y-hat; the Kelvin exterior carries the 2-D
        # reduced-potential background -F_s (sign flip for periodic matching).
        Hs3 = make_reduced_potential_background_cf(
            mesh, lambda xc, yc, zc: CF((0.0, H0, 0.0)),
            R_K=R_K, offset=(offset, 0.0, 0.0), kelvin_mats=("air_outer",), dim=2)
        Hs = CF((Hs3[0], Hs3[1]))

        fes = Periodic(H1(mesh, order=order, dirichlet_bbnd="GND"))
        u, v = fes.TnT()
        aO = BilinearForm(Mu * InnerProduct(grad(u), grad(v)) * dx)
        aO.Assemble()
        fO = LinearForm(Mu * InnerProduct(Hs, grad(v)) * dx)
        fO.Assemble()
        gfO = GridFunction(fes, name="Omega")
        gfO.vec.data = aO.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * fO.vec
        H = Hs - grad(gfO)
        B = Mu * H

        Bx_in, By_in = _interior_B(mesh, B)
        B0 = MU0 * H0
        By_analytic = 2.0 * mu_r / (mu_r + 1.0) * B0
        field_error = abs(By_in - By_analytic) / abs(By_analytic)

        # ---- hodograph net: flux A_z (conjugate of V) on the physical air ----
        # B = (dA_z/dy, -dA_z/dx) -> grad(A_z) = (-B_y, B_x).  Recover A_z by a
        # Poisson projection on air_inner.  A_z is fixed only up to a CONSTANT
        # (gauge) -- pin it with a tiny L2 term (the RHS is orthogonal to the
        # constant, so the gradient is unaffected); do NOT Dirichlet a whole
        # boundary (A_z = the flux, varies there).
        fesA = H1(mesh, order=order, definedon="air_inner")
        uA, vA = fesA.TnT()
        aA = BilinearForm((InnerProduct(grad(uA), grad(vA))
                          + 1e-8 * uA * vA) * dx("air_inner"), symmetric=True)
        rotB = CF((-B[1], B[0]))
        fA = LinearForm(InnerProduct(rotB, grad(vA)) * dx("air_inner"))
        aA.Assemble()
        fA.Assemble()
        gA = GridFunction(fesA, name="A_z")
        gA.vec.data = aA.mat.Inverse(fesA.FreeDofs(),
                                     inverse="sparsecholesky") * fA.vec
        B_from_A = CF((grad(gA)[1], -grad(gA)[0]))     # flux-conjugate field
        # scalar V: total scalar potential Phi = Omega - H0 y (H = -grad Phi)
        Phi = gfO - H0 * y

        airmask = mesh.MaterialCF({"air_inner": 1.0}, default=0.0)
        offax = airmask * IfPos(x * x + y * y - x_eval * x_eval, 1.0, 0.0)
        den = Integrate(offax * InnerProduct(B, B), mesh)
        consistency = math.sqrt(Integrate(
            offax * InnerProduct(B_from_A - B, B_from_A - B), mesh) / den)

        # ---- the air-box contrast (no Kelvin): truncated Dirichlet disk ----
        mesh2 = Mesh(_truncated_geometry(a, 6.0 * a).GenerateMesh(maxh=maxh))
        mesh2.Curve(order)
        Mu2 = mesh2.MaterialCF({"magnetic": mu_r * MU0}, default=MU0)
        fes2 = H1(mesh2, order=order, dirichlet="outer")
        u2, v2 = fes2.TnT()
        a2 = BilinearForm(Mu2 * InnerProduct(grad(u2), grad(v2)) * dx)
        a2.Assemble()
        f2 = LinearForm(Mu2 * InnerProduct(CF((0.0, H0)), grad(v2)) * dx)
        f2.Assemble()
        gf2 = GridFunction(fes2)
        gf2.vec.data = a2.mat.Inverse(fes2.FreeDofs(),
                                      inverse="sparsecholesky") * f2.vec
        B2 = Mu2 * (CF((0.0, H0)) - grad(gf2))
        _, By2 = _interior_B(mesh2, B2)
        airbox_error = abs(By2 - By_analytic) / abs(By_analytic)

        if plot:
            _plot_net(mesh, gA, Phi, a, R_K)

    return {
        "mu_r": mu_r, "order": order, "ne": int(mesh.ne),
        "By_in": float(By_in), "Bx_in": float(Bx_in),
        "By_analytic": float(By_analytic),
        "field_error": float(field_error),          # ~2e-8 (Kelvin, exact)
        "airbox_error": float(airbox_error),        # ~3e-3 (truncated r/a=6)
        "consistency": float(consistency),          # hodograph B(A_z) vs B(V)
    }


def _plot_net(mesh, gA, Phi, a, R_K):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = 220
    g = np.linspace(-R_K * 0.98, R_K * 0.98, n)
    XX, YY = np.meshgrid(g, g)
    AZ = np.full(XX.shape, np.nan)
    PH = np.full(XX.shape, np.nan)
    for i in range(n):
        for j in range(n):
            if math.hypot(g[j], g[i]) < R_K * 0.98:
                mip = mesh(g[j], g[i])
                PH[i, j] = Phi(mip)
                if math.hypot(g[j], g[i]) > a:        # A_z defined on air_inner
                    AZ[i, j] = gA(mip)
    fig, ax = plt.subplots(figsize=(5.4, 5.4), dpi=150)
    ax.contour(XX, YY, AZ, levels=18, colors="C0", linewidths=0.7)
    ax.contour(XX, YY, PH, levels=18, colors="C3", linewidths=0.7)
    th = np.linspace(0, 2 * math.pi, 120)
    ax.fill(a * np.cos(th), a * np.sin(th), color="lightblue", alpha=0.6)
    ax.plot(R_K * np.cos(th), R_K * np.sin(th), "g--", lw=1.2)
    ax.set_aspect("equal")
    ax.set_xlim(-R_K, R_K)
    ax.set_ylim(-R_K, R_K)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("2-D hodograph net on the EXACT (Kelvin) field\n"
                 "flux lines $\\{A_z\\}$ (blue) $\\perp$ equipotentials "
                 "$\\{V\\}$ (red)")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Hodograph + Kelvin in 2-D Cartesian (exact open boundary, no air box)\n")
    r = solve(plot=True)
    print(f"  mu_r={r['mu_r']:.0f}  ne={r['ne']}  order={r['order']}")
    print(f"  interior B = {r['By_in']:.6e}  (analytic 2mu_r/(mu_r+1)B0 = "
          f"{r['By_analytic']:.6e}),  Bx = {r['Bx_in']:.1e}")
    print(f"  Kelvin field_error  = {r['field_error']:.2e}   (EXACT open boundary)")
    print(f"  air-box field_error = {r['airbox_error']:.2e}   (truncated r/a=6 "
          f"Dirichlet -- the thing Kelvin replaces)")
    print(f"  -> Kelvin is ~{r['airbox_error']/r['field_error']:.0e}x more "
          f"accurate; the 2-D Kelvin is conformal => WEIGHT-FREE + exact.")
    print(f"  hodograph consistency B(from A_z) vs B(from V) = "
          f"{r['consistency']:.2e}")
    print("\n  => 'Kelvin in the hodograph' = analytic continuation of "
          "W = A_z + i mu0 V")
    print("     across the conformal inversion -- the open boundary is a point.")


if __name__ == "__main__":
    main()
