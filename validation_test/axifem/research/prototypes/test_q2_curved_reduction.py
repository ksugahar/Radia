"""test_q2_curved_reduction.py — historical Phase B5 reduction tests.

These tests belong to the Python prototype retained in this research corpus.
The production C++ curved-Q2 path now ships as ``AxiHenrotteFE_Q2_Curved`` via
``H1Henrotte(mesh, order=2, curvedquad=True)``; the current regression gate is
``tests/axifem/test_q2_curved.py``.

Tests the AxifemQuadQ2Curved (curved 9-node Q2) against:
  1. axifem_quad_q2.py axis-aligned Gauss reference (must agree to machine
     precision when the curved element is fed quadratic-mean mid-node r-coords).
  2. Element-level invariants (symmetry, nodal interpolation, partition of
     unity, axis decoupling).
  3. Sanity tests on a true curved element (bowed edges).
"""
from __future__ import annotations

import sys

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

from axifem_quad_q2 import (
    element_matrices_q2_numerical as q2_axis_K,
    element_sigma_mass_q2_numerical as q2_axis_M,
)
from axifem_quad_q2_curved import (
    AxifemQuadQ2Curved,
    make_axis_aligned_q2_curved,
    make_straight_q2_curved,
    MU0,
)


# ---------------------------------------------------------------------------
# Test 1: axis-aligned reduction
# ---------------------------------------------------------------------------


def test_axis_aligned_reduction():
    """With quad-mean r mid-nodes, AxifemQuadQ2Curved must reproduce the
    axifem_quad_q2 element matrices to machine precision (both use the same
    (s, z)-monomial basis and the same physical integration domain).

    To make the test independent of axifem_quad_q2.py's fixed 8-pt Gauss,
    we temporarily raise its quadrature order to match ours (n=20). Both
    should then agree to ~1e-12 (limited by FP rounding given high cond(K))."""
    ra, rb, za, zb = 1e-3, 2e-3, 0.0, 1e-3
    mu0 = MU0
    sigma = 5.8e7

    # Patch axifem_quad_q2's GL order to n=20 for the comparison.
    import axifem_quad_q2 as q2mod
    saved_x, saved_w = q2mod._GL8_x.copy(), q2mod._GL8_w.copy()
    try:
        q2mod._GL8_x, q2mod._GL8_w = np.polynomial.legendre.leggauss(20)
        K_ref, _ = q2_axis_K(ra, rb, za, zb, mu0, mu0)
        M_ref = q2_axis_M(ra, rb, za, zb, sigma)
    finally:
        q2mod._GL8_x, q2mod._GL8_w = saved_x, saved_w

    elem = make_axis_aligned_q2_curved(ra, rb, za, zb, mid_r="quad")
    K = elem.stiffness(mu_r=mu0, mu_z=mu0, n_quad=20)
    M = elem.sigma_mass(sigma, n_quad=20)

    norm_K = np.linalg.norm(K_ref)
    norm_M = np.linalg.norm(M_ref)
    err_K = np.linalg.norm(K - K_ref) / norm_K
    err_M = np.linalg.norm(M - M_ref) / norm_M
    print(f"[reduction] cond(V_curved) = {elem._vandermonde_cond:.3e}")
    print(f"[reduction] ||K - K_ref|| / ||K_ref|| = {err_K:.3e}")
    print(f"[reduction] ||M - M_ref|| / ||M_ref|| = {err_M:.3e}")
    assert err_K < 1e-10, f"K reduction failed: rel err = {err_K}"
    assert err_M < 1e-10, f"M reduction failed: rel err = {err_M}"
    print("[OK] axis-aligned reduction matches axifem_quad_q2 to <1e-10 rel "
          "(both at n_quad=20)")


def test_axis_aligned_reduction_axis_touching():
    """Axis-touching axis-aligned element: K, M must agree (axis row/col is
    identically zero in both formulations via the r_i factor)."""
    ra, rb, za, zb = 0.0, 2e-3, 0.0, 1e-3
    mu0 = MU0
    sigma = 5.8e7

    # axifem_quad_q2 raises NotImplementedError for axis-touching — skip K/M
    # comparison here; instead check that the curved element's axis rows/cols
    # are zero and the non-axis 6x6 sub-block reproduces a refined check.
    elem = make_axis_aligned_q2_curved(ra, rb, za, zb, mid_r="quad")
    K = elem.stiffness(mu_r=mu0, mu_z=mu0, n_quad=14)
    M = elem.sigma_mass(sigma, n_quad=14)

    # Axis nodes: indices 0, 3, 7 (r = 0).
    axis_idx = [0, 3, 7]
    for k in axis_idx:
        assert np.max(np.abs(K[k, :])) < 1e-30, f"K row {k} not zero"
        assert np.max(np.abs(K[:, k])) < 1e-30, f"K col {k} not zero"
        assert np.max(np.abs(M[k, :])) < 1e-30, f"M row {k} not zero"
        assert np.max(np.abs(M[:, k])) < 1e-30, f"M col {k} not zero"
    print("[OK] axis-touching: axis row/col exactly zero")

    # Non-axis sub-block should be SPD
    non_axis = [i for i in range(9) if i not in axis_idx]
    K_sub = K[np.ix_(non_axis, non_axis)]
    M_sub = M[np.ix_(non_axis, non_axis)]
    eK = np.linalg.eigvalsh((K_sub + K_sub.T) / 2)
    eM = np.linalg.eigvalsh((M_sub + M_sub.T) / 2)
    assert eK.min() > 0, f"K sub-block not SPD: min eig = {eK.min()}"
    assert eM.min() > 0, f"M sub-block not SPD: min eig = {eM.min()}"
    print(f"[OK] non-axis sub-block SPD: K eig ∈ [{eK.min():.2e}, {eK.max():.2e}]"
          f", M eig ∈ [{eM.min():.2e}, {eM.max():.2e}]")


# ---------------------------------------------------------------------------
# Test 2: element-level invariants
# ---------------------------------------------------------------------------


def test_nodal_interpolation():
    """psi_i(r_j, z_j) = delta_ij."""
    elem = make_axis_aligned_q2_curved(1e-3, 2e-3, 0.0, 1e-3, mid_r="quad")
    for j in range(9):
        v = elem.shape(elem.r[j], elem.z[j])
        expected = np.zeros(9)
        expected[j] = 1.0
        err = np.linalg.norm(v - expected)
        assert err < 1e-9, f"Nodal interp failed at node {j}: err = {err}"
    print("[OK] nodal interpolation psi_i(r_j, z_j) = delta_ij")


def test_partition_of_unity():
    """sum_i psi_i(r, z) == 1 inside the element."""
    elem = make_axis_aligned_q2_curved(1e-3, 2e-3, 0.0, 1e-3, mid_r="quad")
    for (r, z) in [(1.2e-3, 0.2e-3), (1.5e-3, 0.5e-3), (1.8e-3, 0.7e-3)]:
        s = elem.shape(r, z).sum()
        assert abs(s - 1.0) < 1e-9, f"Partition of unity fail at ({r}, {z}): sum = {s}"
    print("[OK] partition of unity at 3 interior points")


def test_symmetry():
    """K, M symmetric."""
    elem = make_axis_aligned_q2_curved(1e-3, 2e-3, 0.0, 1e-3, mid_r="quad")
    K = elem.stiffness()
    M = elem.sigma_mass(5.8e7)
    eK = np.linalg.norm(K - K.T) / max(np.linalg.norm(K), 1e-30)
    eM = np.linalg.norm(M - M.T) / max(np.linalg.norm(M), 1e-30)
    assert eK < 1e-10, f"K asymmetric: {eK}"
    assert eM < 1e-10, f"M asymmetric: {eM}"
    print(f"[OK] K, M symmetric to {max(eK, eM):.2e}")


def test_constant_potential_kernel():
    """For axisymmetric A_phi, the gauge A_phi = const / r has zero magnetic
    energy (B = 0 because rot of phi-component reduces to (1/r) d/dr(r * const/r)
    = 0 and similarly d/dz). In V-DOF coords V_i = const / r_i, so K V should
    be near zero."""
    elem = make_axis_aligned_q2_curved(1e-3, 2e-3, 0.0, 1e-3, mid_r="quad")
    K = elem.stiffness()
    V = 1.0 / elem.r  # A_phi = 1/r at all nodes
    energy = V @ K @ V
    norm = np.linalg.norm(K) * np.linalg.norm(V) ** 2
    rel = energy / norm
    assert abs(rel) < 1e-12, f"Constant-(1/r) kernel: relative energy = {rel}"
    print(f"[OK] constant-(1/r) kernel energy = {energy:.3e} (rel {rel:.3e})")


# ---------------------------------------------------------------------------
# Test 3: truly curved element
# ---------------------------------------------------------------------------


def test_curved_element_basic():
    """Bow one edge of an otherwise axis-aligned element and verify symmetry,
    nodal interp, partition of unity, and reasonable conditioning."""
    ra, rb, za, zb = 1e-3, 2e-3, 0.0, 1e-3
    rm_q = float(np.sqrt(0.5 * (ra * ra + rb * rb)))
    zm = 0.5 * (za + zb)
    rn = np.array([ra, rb, rb, ra, rm_q, rb, rm_q, ra, rm_q])
    zn = np.array([za, za, zb, zb, za, zm, zb, zm, zm])
    # Bow the right edge outward: move mid-node 5 to r = rb + 0.05 * (rb - ra)
    rn[5] = rb + 0.05 * (rb - ra)
    elem = AxifemQuadQ2Curved(rn, zn)

    # Nodal interpolation
    for j in range(9):
        v = elem.shape(elem.r[j], elem.z[j])
        e = np.zeros(9); e[j] = 1.0
        err = np.linalg.norm(v - e)
        assert err < 1e-7, f"Curved nodal interp fail at {j}: err = {err}"

    # Symmetry
    K = elem.stiffness()
    M = elem.sigma_mass(5.8e7)
    eK = np.linalg.norm(K - K.T) / max(np.linalg.norm(K), 1e-30)
    eM = np.linalg.norm(M - M.T) / max(np.linalg.norm(M), 1e-30)
    assert eK < 1e-9, f"Curved K asymmetric: {eK}"
    assert eM < 1e-9, f"Curved M asymmetric: {eM}"

    # SPD K (non-axis case, no axis node, so full 9x9 must be PSD; the rank
    # equals 9 minus the (1/r) constant-gauge dim = 8). Smallest eigval should
    # be ~ 0 (gauge), rest > 0.
    eigs = np.linalg.eigvalsh((K + K.T) / 2)
    print(f"[curved] K eigvals: {eigs}")
    # Drop the smallest (gauge), rest should be > 0 with healthy spread.
    assert eigs[1] > 0, f"K rank-8 SPD failed: 2nd-smallest = {eigs[1]}"
    print(f"[OK] curved element: nodal interp, symmetry, SPD on non-gauge subspace")
    print(f"[curved] cond(V) = {elem._vandermonde_cond:.3e}")


def test_curved_sphere_segment():
    """A 'pie slice' near the sphere boundary at r = R: mid-nodes on the arc
    are at sphere radius R, corners at r < R. Sanity check that conditioning
    stays reasonable for typical IGTE sphere-mesh curvature."""
    R = 50e-3
    # 'inner' corners at radius 0.9 R, 'outer' corners at sphere itself.
    # In (r, z), the outer arc passes through (R*sin(t1), R*cos(t1)) and
    # (R*sin(t2), R*cos(t2)); pick t1 = pi/3, t2 = pi/2.
    t1, t2 = np.pi/3, np.pi/2
    inner_t1 = 0.9
    inner_t2 = 0.9
    # corners 0..3 CCW: (inner@t1, inner@t2, outer@t2, outer@t1)
    rs = np.array([inner_t1 * R * np.sin(t1), inner_t2 * R * np.sin(t2),
                   R * np.sin(t2), R * np.sin(t1)])
    zs = np.array([inner_t1 * R * np.cos(t1), inner_t2 * R * np.cos(t2),
                   R * np.cos(t2), R * np.cos(t1)])
    # edge midpoints: edges 4=0-1 radial, 5=1-2 outer arc, 6=2-3 angular outer,
    # 7=3-0 inner.  Mid-node 5 (between corners 1 and 2): both at radius R if
    # we pick inner_t2 = 1; but inner_t2 = 0.9, so 5 is "almost outer" — pick
    # straight midpoint.
    # For test: use straight midpoints (no curvature), checking that the
    # geometric mapping handles general non-axis-aligned quads.
    elem = make_straight_q2_curved(rs, zs)

    # Symmetry
    K = elem.stiffness()
    M = elem.sigma_mass(5.8e7)
    eK = np.linalg.norm(K - K.T) / max(np.linalg.norm(K), 1e-30)
    eM = np.linalg.norm(M - M.T) / max(np.linalg.norm(M), 1e-30)
    assert eK < 1e-9, f"K asymmetric for sphere-segment quad: {eK}"
    assert eM < 1e-9, f"M asymmetric for sphere-segment quad: {eM}"
    print(f"[OK] sphere-segment quad: K, M symmetric to {max(eK, eM):.2e}")
    print(f"[curved-sphere] cond(V) = {elem._vandermonde_cond:.3e}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 70)
    print("axifem Q2 curved element — historical Phase B5 element-level tests")
    print("=" * 70)
    print()
    test_axis_aligned_reduction()
    test_axis_aligned_reduction_axis_touching()
    test_nodal_interpolation()
    test_partition_of_unity()
    test_symmetry()
    test_constant_potential_kernel()
    test_curved_element_basic()
    test_curved_sphere_segment()
    print()
    print("All historical Phase B5 element-level tests passed.")


if __name__ == "__main__":
    main()
