"""Reduced-potential Clebsch hodograph, path B increment 1: the reduced field is
SOLVED in NGSolve (not analytic), then the hodograph is applied -- both potential
coordinates, with the DUAL boundary conditions made explicit.

2-D cylinder (radius a, mu_r) in a uniform applied field H0 x_hat. The reduced
scalar potential phi (H = H_s - grad phi, H_s = H0 x_hat) is solved by FEM:

    int mu grad(phi).grad(v) dx = int mu H_s.grad(v) dx,     mu = mu_r (r<a), 1 (r>a)

(the iron-air interface is NATURAL -- B_n and H_t continuous). The air field
B = H_s - grad(phi) is curl-free, so it has BOTH hodograph coordinates:

    scalar  Phi = (H_s.x) - phi = x - phi      (B = grad Phi)      directly from the solve
    flux    psi   with   grad(psi) = R90 B      (B_x = psi_y, B_y = -psi_x)

DUAL BOUNDARY CONDITIONS (the key point): Phi and psi are conjugate, so their
Dirichlet/Neumann roles SWAP.  The flux psi is solved with Dirichlet ONLY on the
clean far boundary (to fix the additive constant); the iron interface is NATURAL
(Neumann) -- it is ~an equipotential where the flux ENTERS and psi VARIES (where
Phi would be Dirichlet).  Imposing Dirichlet psi on the interface over-constrains
it and breaks the solve (an instructive bug: 34% vs 8e-4 once fixed).

Verified vs the exact 2-D solution (H0 = mu0 = a = 1, beta = (mu_r-1)/(mu_r+1)):
  reduced phi -> exact; air field -> exact 2-D dipole; and B(flux psi) == B(potential Phi).

Run:  python reduced_clebsch_hodograph_fem.py
"""
import math

from netgen.geom2d import Circle, CSG2d
from ngsolve import (
    CF,
    H1,
    BilinearForm,
    GridFunction,
    InnerProduct,
    Integrate,
    LinearForm,
    Mesh,
    TaskManager,
    dx,
    grad,
    sqrt,
    x,
    y,
)

A_CYL, MU_R, R_OUT = 1.0, 1000.0, 6.0
BETA2 = (MU_R - 1.0) / (MU_R + 1.0)
r2 = x * x + y * y
PHI_in = BETA2 * x                                  # reduced potential, iron
PHI_out = BETA2 * A_CYL**2 * x / r2                 # reduced potential, air (2-D dipole)
PSI_air = y * (1.0 + BETA2 * A_CYL**2 / r2)         # exact air flux (Dirichlet anchor)


def solve(maxh=0.10, order=3):
    geo = CSG2d()
    iron = Circle(center=(0, 0), radius=A_CYL, bc="iface").Mat("iron")
    geo.Add(iron)
    geo.Add(Circle(center=(0, 0), radius=R_OUT, bc="outer").Mat("air") - iron)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        mu = mesh.MaterialCF({"iron": MU_R, "air": 1.0}, default=1.0)
        Hs = CF((1.0, 0.0))

        # --- reduced scalar potential solve (interface NATURAL) ---
        fes = H1(mesh, order=order, dirichlet="outer")
        u, v = fes.TnT()
        a = BilinearForm(mu * InnerProduct(grad(u), grad(v)) * dx, symmetric=True)
        f = LinearForm(mu * InnerProduct(Hs, grad(v)) * dx)
        a.Assemble(); f.Assemble()
        gf = GridFunction(fes); gf.Set(CF(PHI_out), definedon=mesh.Boundaries("outer"))
        rr = gf.vec.CreateVector(); rr.data = f.vec - a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * rr

        phi_ex = mesh.MaterialCF({"iron": PHI_in, "air": PHI_out}, default=PHI_out)
        err_phi = sqrt(Integrate((gf - phi_ex)**2, mesh) / Integrate(phi_ex**2, mesh))

        airmask = mesh.MaterialCF({"iron": 0.0, "air": 1.0}, default=0.0)
        B_phi = Hs - grad(gf)                            # field from the potential (= grad Phi)
        H_ex = Hs - CF((PHI_out.Diff(x), PHI_out.Diff(y)))
        den = Integrate(airmask * InnerProduct(H_ex, H_ex), mesh)
        err_field = math.sqrt(Integrate(airmask * InnerProduct(B_phi - H_ex, B_phi - H_ex),
                                        mesh) / den)

        # --- hodograph flux coordinate (DUAL BCs: Dirichlet far only, interface Neumann) ---
        fes2 = H1(mesh, order=order, definedon=mesh.Materials("air"), dirichlet="outer")
        u2, v2 = fes2.TnT()
        a2 = BilinearForm(InnerProduct(grad(u2), grad(v2)) * dx, symmetric=True)
        f2 = LinearForm(InnerProduct(CF((-B_phi[1], B_phi[0])), grad(v2)) * dx)  # R90 B
        a2.Assemble(); f2.Assemble()
        gpsi = GridFunction(fes2); gpsi.Set(PSI_air, definedon=mesh.Boundaries("outer"))
        r2v = gpsi.vec.CreateVector(); r2v.data = f2.vec - a2.mat * gpsi.vec
        gpsi.vec.data += a2.mat.Inverse(fes2.FreeDofs(), inverse="sparsecholesky") * r2v
        dpsi = grad(gpsi)
        B_psi = CF((dpsi[1], -dpsi[0]))
        agree = math.sqrt(Integrate(airmask * InnerProduct(B_psi - B_phi, B_psi - B_phi),
                                    mesh) / den)

    return {"err_phi": float(err_phi), "err_field": err_field, "agree": agree}


def main():
    print("=== reduced-potential Clebsch hodograph, FEM (path B inc.1): 2-D cylinder ===")
    r = solve()
    print(f"cylinder mu_r={MU_R:.0f}, beta={BETA2:.4f}, R_out={R_OUT}")
    print(f"  reduced-Omega FEM solve:  err_phi = {r['err_phi']:.2e}")
    print(f"  air field B = H_s - grad(phi) vs exact 2-D dipole = {r['err_field']:.2e}")
    print(f"  hodograph: B(flux psi) vs B(potential Phi) = {r['agree']:.2e}  "
          f"(bidirectional; dual BCs: psi Dirichlet far only, interface Neumann)")
    print("  => the reduced field is SOLVED in NGSolve and admits both hodograph")
    print("     coordinates consistently -> ready to redraw geometry in (psi, Phi).")
    ok = (r["err_phi"] < 1e-2 and r["err_field"] < 1e-2 and r["agree"] < 1e-2)
    print(f"\n=== PASS: {ok} ===")
    return r, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
