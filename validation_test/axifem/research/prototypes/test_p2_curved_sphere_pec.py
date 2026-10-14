"""Curved P2 sphere PEC benchmark.

Cu sphere half-disk with A_phi = 0 Dirichlet on the curved sphere surface
(equivalent to PEC at r = R_sphere). The leading eigenvalue tau_1 of
K v = lambda M v gives the slowest-decaying interior diffusion mode.
For axisymmetric sphere this matches Stoll's analytical Bessel formula
tau_1 = mu sigma a^2 / pi^2 = 738.48 us at R = 10 mm, sigma = 5.8e7.

Compares straight-edge baseline (Curve(0)) vs Curve(2) vs Curve(3) on
the SAME mesh -- isolates the geometry-fidelity improvement from
mesh.Curve(p) propagating through the curved-P2 FE.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from math import pi
import time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

from netgen.occ import OCCGeometry, WorkPlane
from ngsolve import Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

R_SPHERE = 10.0e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
STOLL_TAU_1 = MU0 * SIGMA_CU * R_SPHERE ** 2 / (pi ** 2)


def build_sphere_only_geometry():
    """Half-disk of conducting sphere in (r, z) plane (axisym), no air box.
    OCC-built so mesh.Curve(p) projects mid-edge nodes onto the true circle.
    A_phi = 0 on the curved arc (PEC); axis at r = 0."""
    wp = WorkPlane().MoveTo(0.0, -R_SPHERE).Direction(1.0, 0.0)
    wp = wp.Arc(R_SPHERE, 180.0).LineTo(0.0, -R_SPHERE)
    face = wp.Face()
    face.name = "conductor"
    # name the curved arc edge and the axis line for Dirichlet
    # OCC edges in 2D: distinguish by point distance — arc is the longest,
    # the axis is the straight segment at r=0.
    for e in face.edges:
        # Sample midpoint of the edge in (r, z) — arc midpoint is at (R, 0),
        # axis midpoint is at (0, 0).
        # Use a simple heuristic: distance of edge centre from origin.
        try:
            c = e.center
            r_mid = c[0]
        except Exception:
            r_mid = 0.0
        if r_mid > R_SPHERE * 0.4:
            e.name = "sphereBND"
        else:
            e.name = "axis"
    return OCCGeometry(face, dim=2)


def to_csr(mat, n):
    rs, cs, vs = mat.COO()
    return sp.coo_matrix(
        (np.array(vs), (np.array(rs), np.array(cs))),
        shape=(n, n),
    ).toarray()


def run_one(maxh, curve_order):
    ngsglobals.msg_level = 0
    geo = build_sphere_only_geometry()
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    mesh.Curve(curve_order if curve_order > 0 else 1)

    fes = H1Henrotte(mesh, order=2, dirichlet="sphereBND")
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()

    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m_bf.Assemble()

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)
    K = to_csr(a.mat, fes.ndof)[np.ix_(free, free)]
    M = to_csr(m_bf.mat, fes.ndof)[np.ix_(free, free)]
    K = 0.5 * (K + K.T)
    M = 0.5 * (M + M.T)

    diag_M = np.diag(M); diag_K = np.diag(K)
    keep = (diag_M > 1e-30 * (np.max(diag_M) + 1e-30)) | \
           (diag_K > 1e-30 * (np.max(diag_K) + 1e-30))
    K = K[np.ix_(keep, keep)]
    M = M[np.ix_(keep, keep)]
    keep_M = np.diag(M) > 1e-30 * np.max(np.diag(M))
    K = K[np.ix_(keep_M, keep_M)]
    M = M[np.ix_(keep_M, keep_M)]

    eigvals = sla.eigvalsh(K, M, subset_by_index=[0, min(3, K.shape[0] - 1)])
    pos = eigvals[eigvals > 0]
    tau = 1.0 / pos[0] if len(pos) > 0 else np.nan
    gap = (tau - STOLL_TAU_1) / STOLL_TAU_1 * 100
    return mesh.ne, fes.ndof, tau, gap


def main():
    print("=" * 70)
    print(f"Curved-P2 sphere PEC benchmark vs Stoll {STOLL_TAU_1*1e6:.4f} us")
    print("=" * 70)

    levels = [5e-3, 3e-3, 2e-3, 1.5e-3, 1e-3]
    for maxh in levels:
        print(f"\nmaxh = {maxh:.1e}")
        print(f"  {'curve':>8s}  {'ne':>5s}  {'ndof':>6s}  {'tau [us]':>10s}  {'gap':>9s}")
        for curve in [0, 2, 3, 5]:
            try:
                ne, ndof, tau, gap = run_one(maxh, curve)
                print(f"  Curve({curve})  {ne:>5d}  {ndof:>6d}  {tau*1e6:>10.4f}  {gap:>+9.4f}%")
            except Exception as e:
                print(f"  Curve({curve})  FAIL: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
