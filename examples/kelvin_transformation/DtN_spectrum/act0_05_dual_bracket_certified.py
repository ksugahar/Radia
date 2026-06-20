# -*- coding: utf-8 -*-
# DEMO 2 (verified, verdict=confirms): complementary (Prager-Synge) dual bracketing.
# -div grad u = f on (0,1)^2, u=0; u_ex=sin(pi x)sin(pi y), E_ex=pi^2/4.
# Primal H1 energy (lower) and equilibrated RT-flux complementary energy (upper)
# bracket the consistent-data energy E(f_h): E_primal <= E(f_h) <= E_comp on every mesh.
# CAVEAT: strict bound needs EXACT equilibration (div sigma + f = 0); mixed RT/L2 gives
# div sigma = -Pi_h f, so bracket the projected-data energy E(f_h) -> E_ex (data oscillation).
import numpy as np
from ngsolve import *
from netgen.geom2d import unit_square
pi = np.pi
E_EX = pi ** 2 / 4.0


def run(h, order=1):
    mesh = Mesh(unit_square.GenerateMesh(maxh=h))
    f = 2 * pi * pi * sin(pi * x) * sin(pi * y)
    Q = L2(mesh, order=order); fh = GridFunction(Q); fh.Set(f)   # consistent data Pi_h f
    Vp = H1(mesh, order=order, dirichlet=".*"); up, vp = Vp.TnT()
    ap = BilinearForm(grad(up) * grad(vp) * dx); ap.Assemble()
    lp = LinearForm(fh * vp * dx); lp.Assemble()
    uh = GridFunction(Vp)
    uh.vec.data = ap.mat.Inverse(Vp.FreeDofs(), inverse="sparsecholesky") * lp.vec
    E_primal = 0.5 * Integrate(grad(uh) * grad(uh), mesh)
    fes = HDiv(mesh, order=order, RT=True) * L2(mesh, order=order)
    (sig, u), (tau, w) = fes.TnT()
    a = BilinearForm((sig * tau + div(sig) * w + div(tau) * u) * dx); a.Assemble()
    L = LinearForm(-fh * w * dx); L.Assemble()
    gf = GridFunction(fes); gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * L.vec
    sig_h = gf.components[0]
    equil_res = sqrt(Integrate((div(sig_h) + fh) ** 2, mesh))
    E_comp = 0.5 * Integrate(sig_h * sig_h, mesh)
    V6 = H1(mesh, order=6, dirichlet=".*"); u6, v6 = V6.TnT()
    a6 = BilinearForm(grad(u6) * grad(v6) * dx); a6.Assemble()
    l6 = LinearForm(fh * v6 * dx); l6.Assemble()
    g6 = GridFunction(V6); g6.vec.data = a6.mat.Inverse(V6.FreeDofs(), inverse="sparsecholesky") * l6.vec
    E_ref_fh = 0.5 * Integrate(grad(g6) * grad(g6), mesh)
    return dict(h=h, ndof=Vp.ndof, E_primal=E_primal, E_comp=E_comp,
                equil_res=equil_res, E_ref_fh=E_ref_fh)


if __name__ == "__main__":
    print(f"E_ex = pi^2/4 = {E_EX:.10f}")
    print(f"{'h':>5} {'E_primal':>11} {'E(f_h)':>11} {'E_comp':>11} {'gap':>10} {'equil':>9}")
    rows = []
    for h in [0.4, 0.2, 0.1, 0.05]:
        r = run(h); g = r["E_comp"] - r["E_primal"]; rows.append((r, g))
        print(f"{h:>5.3f} {r['E_primal']:>11.7f} {r['E_ref_fh']:>11.7f} {r['E_comp']:>11.7f} {g:>10.4e} {r['equil_res']:>9.1e}")
    print("bracket E_primal<=E(f_h)<=E_comp every mesh:",
          all(r["E_primal"] <= r["E_ref_fh"] + 1e-9 <= r["E_comp"] + 1e-9 for r, _ in rows))
