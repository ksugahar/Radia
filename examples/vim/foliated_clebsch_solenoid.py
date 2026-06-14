"""Foliated current-Clebsch volume stream function + ACA+TSVD  (Stage A / Phase 3).

A 3D VOLUME stream function on a solid conductor: fix a foliation mu (here the
nested cylinders mu = r = sqrt(x^2+y^2)) and put a single scalar lambda on the
volume.  The current is

    J = grad(lambda) x grad(mu)        (div J = 0 identically, verified sympy)

With mu = r this is the azimuthal current  J = (d lambda/dz) phi_hat, so each
mu = const cylinder carries a surface stream function lambda -- a STACK of
surface stream functions = a true 3D bulk winding.  Wires = the level-set
intersections {lambda = const} INT {mu = const} = lambda-contours on each
cylinder.

Because mu is FIXED, J is LINEAR in lambda, so the inverse design

    A lambda = B_target   (M axis targets, N volume H1 DOFs)

is a LINEAR least-norm problem solved by the EXISTING kernel-agnostic
radia.stream_function (ACA+)+TSVD engine UNCHANGED -- only the per-DOF basis
current is the volume grad(phi_j) x grad(mu) instead of the surface n x grad.

Biot-Savart of the basis current (target p_i, Bz component): since grad(mu)_z=0,
    [(grad v x grad mu) x K_i]_z = -(dv/dz)(grad mu . K_i),   K_i=(p_i-x)/|p_i-x|^3
so the response row is a single LinearForm in (dv/dz).

Run:  python examples/vim/foliated_clebsch_solenoid.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, H1, GridFunction, CF, grad, InnerProduct, dx, sqrt,
    x, y, z, BilinearForm, LinearForm, TaskManager, Integrate, Cross,
)
from netgen.csg import CSGeometry, Cylinder, Plane, Pnt, Vec

from stream_function import aca_tsvd, pseudo_inverse_solve

MU0 = 4e-7 * math.pi

# conductor tube  (a thick solenoid winding region)
R0, R1 = 0.05, 0.10     # inner / outer radius [m]
ZL = 0.10               # half-length [m]
B0 = 1.0e-3             # target axial field [T]


def build_mesh(maxh=0.03):
    geo = CSGeometry()
    outer = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R1)
    inner = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R0)
    top = Plane(Pnt(0, 0, ZL), Vec(0, 0, 1))
    bot = Plane(Pnt(0, 0, -ZL), Vec(0, 0, -1))
    tube = (outer - inner) * top * bot
    geo.Add(tube.mat("cond"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def grad_mu_cf():
    r = sqrt(x * x + y * y)
    return CF((x / r, y / r, 0.0))


def assemble_response(mesh, fes, targets):
    """Dense M x N Biot-Savart response A[i,j] = Bz at target i from basis j."""
    gmu = grad_mu_cf()
    u, v = fes.TnT()
    X = CF((x, y, z))
    M, N = len(targets), fes.ndof
    A = np.zeros((M, N))
    for i, (px, py, pz) in enumerate(targets):
        P = CF((px, py, pz))
        d = P - X
        r3 = (InnerProduct(d, d)) ** 1.5
        K = d / r3
        # [(grad v x grad mu) x K]_z = -(dv/dz)(grad mu . K)   since grad(mu)_z=0
        lf = LinearForm(fes)
        lf += -(MU0 / (4 * math.pi)) * grad(v)[2] * InnerProduct(gmu, K) * dx
        lf.Assemble()
        A[i, :] = lf.vec.FV().NumPy()
    return A


def bz_from_lambda(mesh, fes, lam_vec, targets):
    """Recompute Bz at targets directly from a solved lambda (cross-check)."""
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = lam_vec
    gmu = grad_mu_cf()
    J = Cross(grad(gf), gmu)                      # volume current density
    X = CF((x, y, z))
    out = []
    for (px, py, pz) in targets:
        P = CF((px, py, pz))
        d = P - X
        r3 = (InnerProduct(d, d)) ** 1.5
        K = d / r3
        Bz = (MU0 / (4 * math.pi)) * Integrate(Cross(J, K)[2] * dx, mesh)
        out.append(float(Bz))
    return np.array(out), J


def main():
    print("=== foliated current-Clebsch volume SF + ACA+TSVD (solenoid) ===")
    with TaskManager():
        mesh = build_mesh()
        fes = H1(mesh, order=2)            # lambda free everywhere (gauge -> TSVD)
        ndof = fes.ndof
        # axial design targets: uniform Bz = B0 along the bore
        zt = np.linspace(-0.6 * ZL, 0.6 * ZL, 9)
        targets = [(0.0, 0.0, float(zz)) for zz in zt]
        M = len(targets)

        A = assemble_response(mesh, fes, targets)
        Btar = np.full(M, B0)

        # ACA+TSVD least-norm solve (wrap dense A as the entry callback)
        res = aca_tsvd(M, ndof, lambda i, j: A[i, j],
                       modes=M, kmax=M, aca_eps=1e-10, method=3)
        lam = pseudo_inverse_solve(res, Btar)

        fit = A @ lam
        fit_res = np.linalg.norm(fit - Btar) / np.linalg.norm(Btar)

        Bz_check, J = bz_from_lambda(mesh, fes, lam, targets)
        check_res = np.linalg.norm(Bz_check - Btar) / np.linalg.norm(Btar)

        # div J = div(grad(lam) x grad(mu)) = 0 analytically (sympy-verified).
        # Numerical confirmation via the WEAK divergence: for q in H1_0,
        # int J.grad(q) = -int (div J) q = 0.  ||div J||_{H^-1} relative to
        # |J| * length scale (avoids the L2-projection artifact of HDiv.Set).
        fesq = H1(mesh, order=2, dirichlet=".*")
        q, p = fesq.TnT()
        Kf = BilinearForm(fesq, symmetric=True)
        Kf += InnerProduct(grad(q), grad(p)) * dx
        bf = LinearForm(fesq)
        bf += InnerProduct(J, grad(p)) * dx          # p = test function
        Kf.Assemble()
        bf.Assemble()
        w = bf.vec.CreateVector()
        w.data = Kf.mat.Inverse(fesq.FreeDofs(), inverse="sparsecholesky") * bf.vec
        hminus1 = math.sqrt(abs(InnerProduct(bf.vec, w)))
        jL2 = math.sqrt(Integrate(InnerProduct(J, J) * dx, mesh))
        div_rel = hminus1 / (jL2 * (2 * R1) + 1e-30)

    print(f"   mesh ne={mesh.ne}, lambda ndof={ndof}, M targets={M}, k_aca={res.k_aca}")
    print(f"   ACA+TSVD fit residual      = {fit_res:.2e}")
    print(f"   direct Bz(lambda) residual = {check_res:.2e}")
    print(f"   div J / |J| (rel)          = {div_rel:.2e}   (current conservation)")
    print(f"   mean Bz_check              = {Bz_check.mean():.4e} T  (target {B0:.1e})")

    ok = (fit_res < 1e-3 and check_res < 5e-2 and div_rel < 1e-3)
    print(f"\n=== PASS: {ok} ===")
    return {"fit_res": fit_res, "check_res": check_res, "div_rel": div_rel,
            "k_aca": int(res.k_aca), "ndof": ndof, "ne": mesh.ne}, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
