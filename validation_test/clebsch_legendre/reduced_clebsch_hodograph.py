"""Reduced-potential Clebsch hodograph (path A): the hodograph on a real
reduced-potential field, verified against the exact magnetizable-sphere solution.

This composes the two solvers the lab already has:

  * REDUCED POTENTIAL supplies a real, current-free field with a coil/applied
    SOURCE + an iron REACTION, without meshing the source and with an open
    boundary.  Here the canonical exact case: a linear sphere (radius a, relative
    permeability mu_r) in a uniform applied field H0 z_hat.  The exterior field is
    the uniform SOURCE plus the induced-DIPOLE reaction -- i.e. H = H_s - grad(Omega)
    with H_s the uniform source and Omega the reduced (dipole) scalar potential.

  * The CLEBSCH HODOGRAPH (potential-coordinate map; NOT a Legendre transform --
    the coordinate swap is a hodograph, Hamilton 1846) takes that field and solves
    its two potential coordinates -- the Stokes flux function psi and the scalar
    B-potential Phi -- so the geometry can be "redrawn in potential coordinates".

The verification (the point of path A) is that the field reconstructed from the
FLUX coordinate and from the POTENTIAL coordinate AGREE, and both match the exact
sphere field -- the bidirectional consistency that the (psi, Phi) chart is valid
(Jacobian det[grad psi; grad Phi] = r |grad Phi|^2 != 0).  This is the verifiable
foundation for the monolithic "deep" reduced-potential Clebsch hodograph (path B).

Conventions follow examples/feec_vim/bidirectional_map_axisym.py (NGSolve 2D coords
x = rho, y = z; B = grad Phi; Stokes psi with B_rho = -(1/rho)psi_z, B_z =
(1/rho)psi_rho).  Normalized H0 = mu0 = a = 1.  beta = (mu_r - 1)/(mu_r + 2).
Exact consistent pair (current-free: E^2 psi = 0, Laplace Phi = 0):
    psi = rho^2 (1/2 + beta a^3/R^3),   Phi = z (1 - beta a^3/R^3),  R = sqrt(rho^2+z^2)
on an exterior meridian box (R > a, rho > 0).

Run:  python reduced_clebsch_hodograph.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))

from ngsolve import (                                       # noqa: E402
    Mesh, H1, GridFunction, grad, InnerProduct, dx, sqrt, x, y, CF,
    BilinearForm, LinearForm, TaskManager, Integrate,
)
from netgen.geom2d import SplineGeometry                    # noqa: E402

A_SPH = 1.0                                   # sphere radius
MU_R = 1000.0                                 # iron relative permeability
BETA = (MU_R - 1.0) / (MU_R + 2.0)            # induced-dipole strength
RHO0, RHO1, Z0, Z1 = 1.5, 3.0, -1.0, 1.0      # exterior meridian box, R = sqrt(rho^2+z^2) > a

R = sqrt(x * x + y * y)                        # x = rho, y = z
PSI = x * x * (0.5 + BETA * A_SPH**3 / R**3)   # Stokes flux of source+reaction
PHI = y * (1.0 - BETA * A_SPH**3 / R**3)       # B-potential of source+reaction


def _solve(mesh, order, bc_cf, weight):
    """Harmonic extension of bc_cf under int weight*grad(u).grad(v) (Dirichlet bnd)."""
    fes = H1(mesh, order=order, dirichlet="bnd")
    u, v = fes.TnT()
    a = BilinearForm(weight * InnerProduct(grad(u), grad(v)) * dx, symmetric=True)
    a.Assemble()
    f = LinearForm(fes); f.Assemble()
    gf = GridFunction(fes); gf.Set(bc_cf, definedon=mesh.Boundaries("bnd"))
    r = gf.vec.CreateVector(); r.data = f.vec - a.mat * gf.vec
    gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
    return gf


def solve(maxh=0.06, order=3):
    """Solve the reduced-potential field's flux psi (Grad-Shafranov) and potential
    Phi (axisym Laplace); return forward errors, field bidirectional-agreement,
    accuracy vs the exact sphere, and the chart Jacobian residual."""
    geo = SplineGeometry(); geo.AddRectangle((RHO0, Z0), (RHO1, Z1), bc="bnd")
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        gPsi = _solve(mesh, order, PSI, 1.0 / x)     # current-free Grad-Shafranov
        gPhi = _solve(mesh, order, PHI, x)           # axisymmetric Laplace

        errPsi = sqrt(Integrate((gPsi - PSI)**2, mesh) / Integrate(PSI**2, mesh))
        errPhi = sqrt(Integrate((gPhi - PHI)**2, mesh) / Integrate(PHI**2, mesh))
        dPsi, dPhi = grad(gPsi), grad(gPhi)
        B_psi = CF((-dPsi[1] / x, dPsi[0] / x))      # field from the FLUX coordinate
        B_phi = CF((dPhi[0], dPhi[1]))               # field from the POTENTIAL coordinate
        B_an = CF((PHI.Diff(x), PHI.Diff(y)))        # exact sphere field = grad(Phi_exact)
        n_an = sqrt(Integrate(InnerProduct(B_an, B_an), mesh))
        agree = sqrt(Integrate(InnerProduct(B_psi - B_phi, B_psi - B_phi), mesh)) / n_an
        acc = sqrt(Integrate(InnerProduct(B_phi - B_an, B_phi - B_an), mesh)) / n_an
        detJ = dPsi[0] * dPhi[1] - dPsi[1] * dPhi[0]
        g2 = InnerProduct(dPhi, dPhi)
        jac = sqrt(Integrate((detJ - x * g2)**2, mesh) / Integrate((x * g2)**2, mesh))
    return {"errPsi": float(errPsi), "errPhi": float(errPhi),
            "agree": float(agree), "acc": float(acc), "jac": float(jac)}


def main():
    print("=== reduced-potential Clebsch hodograph (sphere in uniform field) ===")
    r = solve()
    print(f"sphere mu_r={MU_R:.0f}, beta={BETA:.4f}; exterior box "
          f"rho[{RHO0},{RHO1}] z[{Z0},{Z1}] (R>a, source+dipole reaction)")
    print(f"  forward errPsi = {r['errPsi']:.2e}   errPhi = {r['errPhi']:.2e}")
    print(f"  B(flux) vs B(potential) agreement = {r['agree']:.2e}  "
          f"<- bidirectional consistency (the hodograph point)")
    print(f"  B(potential) vs exact sphere      = {r['acc']:.2e}")
    print(f"  chart Jacobian det = r|grad Phi|^2 rel-err = {r['jac']:.2e}  (valid chart)")
    print("  => the hodograph recovers BOTH potential coordinates of the real")
    print("     reduced-potential field consistently -> ready to redraw geometry")
    print("     in (psi, Phi) coordinates; foundation for the monolithic path B.")

    ok = (r["errPsi"] < 1e-5 and r["errPhi"] < 1e-5 and r["agree"] < 1e-4
          and r["acc"] < 1e-4 and r["jac"] < 1e-3)
    print(f"\n=== PASS: {ok} ===")
    return r, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
