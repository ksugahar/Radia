"""Vector-T convex volume inverse design + ACA+TSVD  (Stage A / #2).

The MOST GENERAL volume stream function: J = curl T with the vector potential
T in H(curl).  Any divergence-free volume current is J = curl T, so -- unlike
the two-scalar Clebsch J = grad(lambda) x grad(mu) -- there is NO helicity
restriction (verified sympy 2026-06-11: T = lambda*grad(mu) forces
T.curl T = 0 = zero helicity; the full vector T has no such constraint).

CONVEXITY (the key result): J = curl T is LINEAR in T, so B = BiotSavart(curl T)
is linear in T and the inverse design

    A T = B_target    (M axis targets, N H(curl) edge DOFs)

is a LINEAR (convex) least-norm problem the EXISTING radia.stream_function
(ACA+)+TSVD engine solves UNCHANGED.  The GAUGE null space (T -> T + grad(chi),
curl(grad chi) = 0 -> same J -> same field, verified) lies in ker(A); TSVD
TRUNCATES it automatically (the ACA+ rank k_aca collapses to <= M << N, and the
min-norm T is orthogonal to gradients) -- NO explicit tree-cotree gauge needed.

Response row (target p_i, Bz):  [curl(v) x K_i]_z,  K_i = (p_i - x)/|p_i-x|^3.

Run:  python validation_test/feec/vim_legacy/vector_t_inverse.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, HCurl, H1, GridFunction, CF, grad, curl, InnerProduct, dx, x, y, z,
    BilinearForm, LinearForm, TaskManager, Integrate, Cross,
)
from netgen.csg import CSGeometry, Cylinder, Plane, Pnt, Vec

from radia.stream_function import aca_tsvd, pseudo_inverse_solve

MU0 = 4e-7 * math.pi

R0, R1 = 0.05, 0.10
ZL = 0.10
B0 = 1.0e-3


def build_mesh(maxh=0.03):
    geo = CSGeometry()
    outer = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R1)
    inner = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R0)
    top = Plane(Pnt(0, 0, ZL), Vec(0, 0, 1))
    bot = Plane(Pnt(0, 0, -ZL), Vec(0, 0, -1))
    geo.Add(((outer - inner) * top * bot).mat("cond"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def assemble_response(fes, targets):
    """Dense M x N response A[i,e] = Bz at target i from basis current curl(N_e)."""
    u, v = fes.TnT()
    X = CF((x, y, z))
    M, N = len(targets), fes.ndof
    A = np.zeros((M, N))
    for i, (px, py, pz) in enumerate(targets):
        P = CF((px, py, pz))
        d = P - X
        K = d / (InnerProduct(d, d) ** 1.5)
        lf = LinearForm(fes)
        lf += (MU0 / (4 * math.pi)) * Cross(curl(v), K)[2] * dx
        lf.Assemble()
        A[i, :] = lf.vec.FV().NumPy()
    return A


def bz_from_T(fes, T_vec, targets, mesh):
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = T_vec
    J = curl(gf)
    X = CF((x, y, z))
    out = []
    for (px, py, pz) in targets:
        P = CF((px, py, pz))
        d = P - X
        K = d / (InnerProduct(d, d) ** 1.5)
        out.append(float((MU0 / (4 * math.pi)) *
                         Integrate(Cross(J, K)[2] * dx, mesh)))
    return np.array(out), J


def weak_div_rel(mesh, J, scale):
    """||div J||_{H^-1} / (|J| * scale) -- current conservation (J=curl T)."""
    fesq = H1(mesh, order=2, dirichlet=".*")
    q, p = fesq.TnT()
    Kf = BilinearForm(fesq, symmetric=True)
    Kf += InnerProduct(grad(q), grad(p)) * dx
    bf = LinearForm(fesq)
    bf += InnerProduct(J, grad(p)) * dx
    Kf.Assemble()
    bf.Assemble()
    w = bf.vec.CreateVector()
    w.data = Kf.mat.Inverse(fesq.FreeDofs(), inverse="sparsecholesky") * bf.vec
    hm1 = math.sqrt(abs(InnerProduct(bf.vec, w)))
    jL2 = math.sqrt(Integrate(InnerProduct(J, J) * dx, mesh))
    return hm1 / (jL2 * scale + 1e-30)


def main():
    print("=== vector-T convex volume inverse + ACA+TSVD ===")
    with TaskManager():
        mesh = build_mesh()
        fes = HCurl(mesh, order=1)
        ndof = fes.ndof
        zt = np.linspace(-0.6 * ZL, 0.6 * ZL, 9)
        targets = [(0.0, 0.0, float(zz)) for zz in zt]
        M = len(targets)

        A = assemble_response(fes, targets)
        Btar = np.full(M, B0)

        res = aca_tsvd(M, ndof, lambda i, j: A[i, j],
                       modes=M, kmax=M, aca_eps=1e-10, method=3)
        Tco = pseudo_inverse_solve(res, Btar)

        fit_res = np.linalg.norm(A @ Tco - Btar) / np.linalg.norm(Btar)
        Bz_check, J = bz_from_T(fes, Tco, targets, mesh)
        check_res = np.linalg.norm(Bz_check - Btar) / np.linalg.norm(Btar)
        div_rel = weak_div_rel(mesh, J, 2 * R1)

        # GAUGE: a pure-gradient T = grad(chi) produces ZERO current/field.
        # Response applied to the H(curl)-embedded gradient of an H1 bump ~ 0.
        h1 = H1(mesh, order=1)
        gchi = GridFunction(h1)
        gchi.Set(x * y + z * z)            # arbitrary smooth chi
        gg = GridFunction(fes)
        gg.Set(grad(gchi))                 # grad(chi) embedded in H(curl)
        gauge_field = np.abs(A @ gg.vec.FV().NumPy()).max() / B0   # rel to B0

    print(f"   mesh ne={mesh.ne}, T ndof={ndof}, M={M}, k_aca={res.k_aca} (<= M: gauge truncated)")
    print(f"   ACA+TSVD fit residual      = {fit_res:.2e}")
    print(f"   direct Bz(curl T) residual = {check_res:.2e}")
    print(f"   div J / |J| (rel)          = {div_rel:.2e}   (J = curl T, div-free)")
    print(f"   gauge field |A grad(chi)|  = {gauge_field:.2e}   (grad-only T -> ~0 field)")
    print(f"   mean Bz_check              = {Bz_check.mean():.4e} T  (target {B0:.1e})")

    ok = (fit_res < 1e-3 and check_res < 5e-2 and div_rel < 1e-3 and
          gauge_field < 5e-2 and res.k_aca <= M)
    print(f"\n=== PASS: {ok} ===")
    return {"fit_res": fit_res, "check_res": check_res, "div_rel": div_rel,
            "gauge_field": gauge_field, "k_aca": int(res.k_aca),
            "ndof": ndof, "M": M}, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
