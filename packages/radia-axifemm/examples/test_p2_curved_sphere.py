"""Curved-P2 sphere benchmark: OCC sphere half (axisym half-disk in (r,z))
with mesh.Curve(p) so the conductor surface follows the analytical sphere
to higher order.

This tests whether the AxiHenrotteFE_P2_Triangle picks up the curved
mid-edge node positions correctly. Compares straight-edge baseline
(test_p2_cpp_sphere.py) vs Curve(2) and Curve(3) on the same mesh.

Cu sphere R = 10 mm, air box R_air = 50 mm. Stoll tau_1 = 738.48 us.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from math import pi
import time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

from netgen.geom2d import SplineGeometry
from netgen.meshing import meshsize
from ngsolve import Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

R_SPHERE = 10.0e-3
R_AIR = 50.0e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
STOLL_TAU_1 = MU0 * SIGMA_CU * R_SPHERE ** 2 / (pi ** 2)


def build_axisym_sphere_geometry():
    """2D half-disk in (x, y) = (r, z) for axisymmetric sphere + air, via
    SplineGeometry. Two domains: conductor (sphere half-disk) and air
    (annulus between R_SPHERE and R_AIR), separated by the sphere boundary.
    """
    geo = SplineGeometry()
    # 6 points: axis ends (top/bottom) of inner+outer + sphere/air apex
    p_inner_bot = geo.AppendPoint(0.0, -R_SPHERE)
    p_inner_apex = geo.AppendPoint(R_SPHERE, 0.0)
    p_inner_top = geo.AppendPoint(0.0, R_SPHERE)
    p_outer_bot = geo.AppendPoint(0.0, -R_AIR)
    p_outer_apex = geo.AppendPoint(R_AIR, 0.0)
    p_outer_top = geo.AppendPoint(0.0, R_AIR)

    # Sphere boundary arcs (left = conductor, right = air). Quadratic spline
    # arcs use 3 control points: (start, control, end). Approximate circle
    # with control point chosen so the spline interpolates the circle.
    # SplineGeometry "spline3" = quadratic Bezier; for a quarter circle the
    # ideal control offset is r * (sqrt(2) - 1) but quadratic only approximates.
    # Use multiple line segments instead for exact axis-aligned bound.
    # Actually SplineGeometry supports "spline3" type with 3 pts that traces
    # the conic through them. Let's just use that.

    # Inner arc bottom-half: from p_inner_bot through (R, -0?) skip -- use
    # quarter arcs: bot to apex, apex to top.
    geo.Append(["spline3", p_inner_bot, p_inner_apex, p_inner_top],
               bc="sphereBND", leftdomain=1, rightdomain=2)
    # Outer arc: bot to apex to top (going right then up-right).
    geo.Append(["spline3", p_outer_bot, p_outer_apex, p_outer_top],
               bc="outer", leftdomain=2, rightdomain=0)
    # Inner axis bottom (conductor): from p_outer_bot to p_inner_bot (axis r=0)
    geo.Append(["line", p_outer_bot, p_inner_bot],
               bc="axis", leftdomain=2, rightdomain=0)
    geo.Append(["line", p_inner_bot, p_inner_top],
               bc="axis_in", leftdomain=1, rightdomain=0)
    geo.Append(["line", p_inner_top, p_outer_top],
               bc="axis", leftdomain=2, rightdomain=0)

    geo.SetMaterial(1, "conductor")
    geo.SetMaterial(2, "air")
    return geo


def to_csr(mat, n):
    rs, cs, vs = mat.COO()
    return sp.coo_matrix(
        (np.array(vs), (np.array(rs), np.array(cs))),
        shape=(n, n),
    ).toarray()


def run_one(maxh, curve_order):
    ngsglobals.msg_level = 0
    geo = build_axisym_sphere_geometry()
    mesh = Mesh(geo.GenerateMesh(maxh=maxh, quad_dominated=False))
    if curve_order > 0:
        mesh.Curve(curve_order)

    fes = H1Henrotte(mesh, order=2, dirichlet="outer")
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
    cond_mask = diag_M > 1e-30 * (np.max(diag_M) + 1e-30)
    nonzero_K = diag_K > 1e-30 * (np.max(diag_K) + 1e-30)
    air_mask = (~cond_mask) & nonzero_K
    cond_idx = np.where(cond_mask)[0]
    air_idx = np.where(air_mask)[0]

    if len(air_idx) > 0:
        K_cc = K[np.ix_(cond_idx, cond_idx)]
        K_ca = K[np.ix_(cond_idx, air_idx)]
        K_aa = K[np.ix_(air_idx, air_idx)]
        X = sla.solve(K_aa, K_ca.T, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K[np.ix_(cond_idx, cond_idx)]
    M_cc = M[np.ix_(cond_idx, cond_idx)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)

    eigvals = sla.eigvalsh(S, M_cc, subset_by_index=[0, min(3, S.shape[0] - 1)])
    pos = eigvals[eigvals > 0]
    tau = 1.0 / pos[0] if len(pos) > 0 else np.nan
    gap = (tau - STOLL_TAU_1) / STOLL_TAU_1 * 100
    return mesh.ne, fes.ndof, tau, gap


def main():
    print("=" * 70)
    print(f"Curved-P2 OCC sphere benchmark vs Stoll {STOLL_TAU_1*1e6:.4f} us")
    print(f"  Cu sphere R={R_SPHERE*1e3} mm, air box R_air={R_AIR*1e3} mm")
    print("=" * 70)

    levels = [5e-3, 3e-3, 2e-3, 1.5e-3]
    for maxh in levels:
        print(f"\nmaxh = {maxh:.1e}")
        print(f"  {'curve':>8s}  {'ne':>6s}  {'ndof':>6s}  {'tau [us]':>10s}  {'gap':>8s}")
        for curve in [0, 2, 3]:
            try:
                ne, ndof, tau, gap = run_one(maxh, curve)
                print(f"  Curve({curve})  {ne:>6d}  {ndof:>6d}  {tau*1e6:>10.4f}  {gap:>+8.3f}%")
            except Exception as e:
                print(f"  Curve({curve})  FAIL: {e}")


if __name__ == "__main__":
    main()
