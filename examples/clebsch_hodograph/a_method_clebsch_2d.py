"""A-method Clebsch hodograph on the 2-D cylinder (vector potential primary).

RESEARCH example (NOT a panel mode).  The verified forward "Clebsch
hodograph" panel mode (radia-electromagnet) solves the reduced SCALAR
potential Omega as primary; this script is its DUAL -- it solves the
VECTOR potential as primary, the way the user asked for: "express A with
Clebsch potentials and formulate in the A-method".

Clebsch representation
----------------------
  B = grad(phi) x grad(psi),   A = phi grad(psi).
In 2-D with the ignorable z, take phi = A_z (the flux function) and
psi = z, so  A = A_z z_hat  -- the standard 2-D vector potential, and
B = curl(A_z z_hat) = (d_y A_z, -d_x A_z).

A-method solve (PRIMARY = A_z)
------------------------------
  div( (1/mu) grad A_z ) = 0,   A_z = B0 y (+ dipole) on the far boundary.
Then the CONJUGATE scalar potential V (H = grad V, current-free air) is
recovered from the SAME field.  (A_z, mu0 V) are Cauchy-Riemann conjugates
-- W = A_z + i mu0 V is analytic -- which is the classical complex-
potential / conformal pole-design picture (Rogowski, Halbach), here as an
FEM realisation.  The orthogonal nets {A_z = const} (flux lines) and
{V = const} (equipotentials) ARE the hodograph grid.

Hodograph self-consistency
--------------------------
  consistency = || B(from A_z) - B(from V) || / || B ||   (air)
tends to 0 as the FEM field becomes exact -- the A-method analogue of the
panel mode's a-posteriori field-quality metric.

Verified (mu0 = B0 = a = 1, mu_r = 1000, order 3, maxh 0.08):
  interior B = 2 mu_r/(mu_r+1) B0 = 1.998  (err ~3e-5)
  A_z / V accuracy ~5e-5,  consistency ~4e-4.

Exact 2-D cylinder (conjugate pair, OPPOSITE dipole signs):
  A_z(air) = B0 y (1 + beta2 a^2/r^2),   V(air) = B0 x (1 - beta2 a^2/r^2),
  beta2 = (mu_r - 1)/(mu_r + 1).
"""
import math
import os

from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx,
                     x, y, CF, BilinearForm, LinearForm, TaskManager, Integrate)
from netgen.geom2d import CSG2d, Circle


def solve(mu_r=1000.0, order=3, maxh=0.08, a_cyl=1.0, r_out=6.0, B0=1.0,
          plot=False):
    """Solve the A-method Clebsch problem and return the verification dict.

    Returns a dict with interior-field accuracy, A_z / V accuracy vs the
    exact 2-D cylinder, and the hodograph bidirectional consistency.
    """
    beta2 = (mu_r - 1.0) / (mu_r + 1.0)
    r2 = x * x + y * y
    Az_exact = B0 * y * (1.0 + beta2 * a_cyl**2 / r2)    # vector potential (air)
    V_exact = B0 * x * (1.0 - beta2 * a_cyl**2 / r2)     # scalar potential (air)

    geo = CSG2d()
    iron = Circle(center=(0, 0), radius=a_cyl, bc="iface").Mat("iron")
    air = (Circle(center=(0, 0), radius=r_out, bc="outer").Mat("air") - iron)
    geo.Add(iron)
    geo.Add(air)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        mu = mesh.MaterialCF({"iron": mu_r, "air": 1.0}, default=1.0)
        ironmask = mesh.MaterialCF({"iron": 1.0, "air": 0.0}, default=0.0)
        airmask = mesh.MaterialCF({"iron": 0.0, "air": 1.0}, default=0.0)

        # --- A-METHOD: vector potential A_z is the PRIMARY unknown ---
        fesA = H1(mesh, order=order, dirichlet="outer")
        u, v = fesA.TnT()
        aA = BilinearForm((1.0 / mu) * InnerProduct(grad(u), grad(v)) * dx,
                          symmetric=True)
        aA.Assemble()
        fA = LinearForm(fesA)
        fA.Assemble()
        gA = GridFunction(fesA, name="A_z")
        gA.Set(CF(Az_exact), definedon=mesh.Boundaries("outer"))
        rA = gA.vec.CreateVector()
        rA.data = fA.vec - aA.mat * gA.vec
        gA.vec.data += aA.mat.Inverse(fesA.FreeDofs(),
                                      inverse="sparsecholesky") * rA
        dA = grad(gA)
        B = CF((dA[1], -dA[0]))                  # B = curl(A_z z_hat)

        vol = Integrate(ironmask, mesh)
        Bx_in = Integrate(ironmask * B[0], mesh) / vol
        Bin_analytic = 2.0 * mu_r / (mu_r + 1.0) * B0
        err_Bin = abs(Bx_in - Bin_analytic) / abs(Bin_analytic)
        errAz = math.sqrt(Integrate(airmask * (gA - Az_exact)**2, mesh)
                          / Integrate(airmask * Az_exact**2, mesh))

        # --- CONJUGATE scalar V recovered from the SAME field ---
        fesV = H1(mesh, order=order, definedon=mesh.Materials("air"),
                  dirichlet="outer")
        uV, vV = fesV.TnT()
        aV = BilinearForm(InnerProduct(grad(uV), grad(vV)) * dx, symmetric=True)
        aV.Assemble()
        fV = LinearForm(InnerProduct(B, grad(vV)) * dx)
        fV.Assemble()
        gV = GridFunction(fesV, name="V")
        gV.Set(CF(V_exact), definedon=mesh.Boundaries("outer"))
        rV = gV.vec.CreateVector()
        rV.data = fV.vec - aV.mat * gV.vec
        gV.vec.data += aV.mat.Inverse(fesV.FreeDofs(),
                                      inverse="sparsecholesky") * rV
        B_from_V = grad(gV)

        den = Integrate(airmask * InnerProduct(B, B), mesh)
        consistency = math.sqrt(
            Integrate(airmask * InnerProduct(B_from_V - B, B_from_V - B), mesh)
            / den)
        errV = math.sqrt(Integrate(airmask * (gV - V_exact)**2, mesh)
                         / Integrate(airmask * V_exact**2, mesh))

        if plot:
            _plot_hodograph_net(mesh, gA, gV, a_cyl, r_out)

    return {
        "mu_r": mu_r, "order": order, "ne": mesh.ne,
        "ndof": int(fesA.ndof + fesV.ndof),
        "B_in": float(Bx_in), "B_in_analytic": float(Bin_analytic),
        "field_error": float(err_Bin),
        "Az_error": float(errAz), "V_error": float(errV),
        "consistency": float(consistency),
    }


def _plot_hodograph_net(mesh, gA, gV, a_cyl, r_out):
    """Save the orthogonal hodograph net: flux lines {A_z=const} and
    equipotentials {V=const} next to this script."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = 220
    xs = np.linspace(-3.0, 3.0, n)
    ys = np.linspace(-3.0, 3.0, n)
    XX, YY = np.meshgrid(xs, ys)
    AZ = np.full(XX.shape, np.nan)
    VV = np.full(XX.shape, np.nan)
    for i in range(n):
        for j in range(n):
            rr = math.hypot(xs[j], ys[i])
            if a_cyl + 1e-3 < rr < r_out - 1e-3:
                mip = mesh(xs[j], ys[i])
                AZ[i, j] = gA(mip)
                VV[i, j] = gV(mip)
    fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=150)
    ax.contour(XX, YY, AZ, levels=25, colors="C0", linewidths=0.7)
    ax.contour(XX, YY, VV, levels=25, colors="C3", linewidths=0.7)
    th = np.linspace(0, 2 * math.pi, 200)
    ax.fill(a_cyl * np.cos(th), a_cyl * np.sin(th), color="0.85", zorder=3)
    ax.set_aspect("equal")
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Hodograph net: flux $A_z$ (blue) $\\perp$ scalar $V$ (red)")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    res = solve(plot=True)
    print("A-method Clebsch hodograph (2-D cylinder)")
    print(f"  mu_r={res['mu_r']:.0f}  ne={res['ne']}  ndof={res['ndof']}")
    print(f"  interior B = {res['B_in']:.5f}  analytic = "
          f"{res['B_in_analytic']:.5f}  field_error = {res['field_error']:.2e}")
    print(f"  A_z accuracy (air) = {res['Az_error']:.2e}")
    print(f"  V   accuracy (air) = {res['V_error']:.2e}")
    print(f"  hodograph consistency B(A_z) vs B(V) = {res['consistency']:.2e}")


if __name__ == "__main__":
    main()
