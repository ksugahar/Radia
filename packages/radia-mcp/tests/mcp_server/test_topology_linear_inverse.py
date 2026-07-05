# -*- coding: utf-8 -*-
"""TSVD / Tikhonov regularized linear inversion (field-synthesis solver), verified against
the SVD pseudo-inverse and the filter-factor / L-curve identities on a controlled
ill-conditioned, rank-deficient system with a known SVD.

Grounds the lab's truncated-SVD magnet/coil field-synthesis method (public-safe curated corpus
30_Optimization\\2019_10_16_TSVD) as the linear-inverse core of topology_opt_applications
topic `field_synthesis`.  Exercises radia_mcp.topology_optimization.linear_inverse.
"""
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.topology_optimization.linear_inverse import (
    tsvd_solve, tikhonov_solve, filter_factors, lcurve, lcurve_corner)


def _build():
    """Controlled A = U diag(s) V^T : 8x5, rank 3, singular values [5, 1, 0.05]."""
    rng = np.random.default_rng(0)
    U, _ = np.linalg.qr(rng.standard_normal((8, 3)))      # 8x3 orthonormal columns
    Vt, _ = np.linalg.qr(rng.standard_normal((5, 3)))     # 5x3 -> V^T is 3x5
    s = np.array([5.0, 1.0, 0.05])
    A = U @ np.diag(s) @ Vt.T
    x_true = Vt @ rng.standard_normal(3)                  # truth in row space (recoverable)
    b = A @ x_true + 1e-6 * rng.standard_normal(8)        # mild noise
    return A, b, s, x_true


def test_lcurve_corner_selects_interior_regularization():
    A, b, _, _ = _build()
    lams = np.array([1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0])
    corner = lcurve_corner(A, b, lams)

    assert 0 < corner["index"] < len(lams) - 1
    assert corner["lambda"] == lams[corner["index"]]
    assert corner["curvature"][0] == 0.0
    assert corner["curvature"][-1] == 0.0
    assert corner["curvature"][corner["index"]] > 0.0
    assert corner["residual_norm"] == corner["points"][corner["index"]][0]
    assert corner["solution_norm"] == corner["points"][corner["index"]][1]


def main():
    A, b, s_true, x_true = _build()
    pinv = np.linalg.pinv(A) @ b

    # 1) TSVD at k=rank == minimum-norm least-squares (pseudo-inverse)
    x3 = tsvd_solve(A, b, 3)
    e_pinv = np.linalg.norm(x3 - pinv) / np.linalg.norm(pinv)
    print(f"  TSVD(k=3) vs pinv         rel = {e_pinv:.2e}")
    assert e_pinv < 1e-10, e_pinv

    # 2) Tikhonov lam->0 -> pinv
    xt = tikhonov_solve(A, b, 1e-7)
    e_tik = np.linalg.norm(xt - pinv) / np.linalg.norm(pinv)
    print(f"  Tikhonov(lam=1e-7) vs pinv rel = {e_tik:.2e}")
    assert e_tik < 1e-4, e_tik

    # 3) filter factors phi_n = s^2/(s^2+lam^2) (and the recovered singular values)
    s_meas = np.linalg.svd(A, compute_uv=False)
    assert np.allclose(np.sort(s_meas)[::-1][:3], s_true, rtol=1e-9), s_meas
    lam = 0.2
    phi = filter_factors(s_true, lam)
    assert np.allclose(phi, s_true**2 / (s_true**2 + lam**2)), phi
    print(f"  filter factors @lam=0.2   = {np.round(phi,4)}")

    # 4) TSVD monotonicity: residual non-increasing, solution-norm non-decreasing in k
    res = [np.linalg.norm(A @ tsvd_solve(A, b, k) - b) for k in (1, 2, 3)]
    nrm = [np.linalg.norm(tsvd_solve(A, b, k)) for k in (1, 2, 3)]
    print(f"  TSVD residual(k=1,2,3) = {np.round(res,6)}")
    print(f"  TSVD solnorm(k=1,2,3)  = {np.round(nrm,4)}")
    assert res[0] >= res[1] >= res[2] - 1e-12, res
    assert nrm[0] <= nrm[1] + 1e-12 <= nrm[2] + 1e-12, nrm

    # 5) Tikhonov L-curve: residual increases, solution norm decreases as lam grows
    pts = lcurve(A, b, [1e-3, 1e-1, 1.0])
    r = [p[0] for p in pts]; xn = [p[1] for p in pts]
    assert r[0] <= r[1] <= r[2], r            # bigger lam -> worse fit
    assert xn[0] >= xn[1] >= xn[2], xn        # bigger lam -> smaller solution
    print(f"  L-curve residuals = {np.round(r,5)}  solnorms = {np.round(xn,4)}")

    corner = lcurve_corner(A, b, [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0])
    assert 0 < corner["index"] < 4, corner
    print(f"  L-curve corner lambda = {corner['lambda']:.1e}  "
          f"curvature = {corner['curvature'][corner['index']]:.3g}")

    print("\n[OK] TSVD/Tikhonov field-synthesis inversion: TSVD(k=rank)==pinv (1e-10), "
          "Tikhonov(lam->0)==pinv, filter factors, TSVD + L-curve monotonicity verified.")


if __name__ == "__main__":
    main()
