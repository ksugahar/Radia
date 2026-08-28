"""Tests for radia.stream_function -- generic (ACA+)+TSVD least-norm solver.

The C++ core (src/core/rad_stream_function.cpp) is kernel-agnostic: it does
(ACA+)+TSVD of an M x N matrix whose entries A(i,j) come from a caller-supplied
callback.  ACA+ itself is HACApK's cHACApK_acaplus (single source of truth);
only the recompression to a truncated SVD (the standard QR-of-a-low-rank-product;
peer review JIAM-2026-36) lives in Radia.  There is ONE recompression method,
``method="aca_qr_tsvd"`` (default); ``method="dense"`` is the direct-SVD
reference.  NO backward compatibility -- the legacy manuscript Method 2/3 (and
their integer selectors) were removed 2026-07-05.

Two kernels are exercised:
  1. a coil Biot-Savart kernel (mirrored rectangular current loops) implemented
     in numpy HERE, to cross-check against the validated Fortran reference
     coil_solver.f90 (the recompression gives the same SVD by either route);
  2. Radia's OWN field computation (radia.Fld over permanent-magnet objects via
     radia_field_kernel), to prove the same machinery serves magnetic materials,
     not just coils.
"""
import os
import sys

import numpy as np
import pytest

from radia.stream_function import (
    aca_tsvd, pseudo_inverse_solve, radia_field_kernel,
    RegularizedTSVD, StreamTSVD, pseudo_inverse_solve_regularized,
    abe_nearest_field_distance_scales, abe_reduce_node_potential_scales,
    solve_abe_current_potential,
    solve_abe_bounded_current_potential,
)

# f2py reference (LAB-only network drive).
REF_DIR = r"W:\04_卒論論文関係\2025年度\046_伊藤海人\2026_01_06_f2py_matlab比較\f2py"


# --------------------------------------------------------------------------
# Coil Biot-Savart kernel (mirrored rectangular loops) -- numpy reference that
# replicates coil_solver.f90's biot_savart_entry exactly, used to drive the
# f90 cross-check.  This lives in the TEST only; the production module embeds
# no field kernel.
# --------------------------------------------------------------------------
def _seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1
    R21 = float(d @ d)
    o1 = O - P1
    o2 = O - P2
    OR1 = np.sqrt(o1 @ o1)
    OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d)
    O221 = float(o2 @ d)
    Rc12 = O121 / OR1
    Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f  # Hz


def _loop_Hz(obs, center, offsets):
    obs = np.asarray(obs, float)
    C = np.asarray(center, float) + np.asarray(offsets, float)  # (4,3) upper corners
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, C[c], C[n])
    Cl = C.copy()
    Cl[:, 2] = -Cl[:, 2]  # mirror through z
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, Cl[c], Cl[n])
    return hz


def _square_loop(a):
    return np.array([[-a / 2, -a / 2, 0.0],
                     [ a / 2, -a / 2, 0.0],
                     [ a / 2,  a / 2, 0.0],
                     [-a / 2,  a / 2, 0.0]])


def _coil_geometry(obs_z, n_side=10, m_side=5, z0=0.005, a=0.02,
                   span=0.1, obs_span=0.08):
    xs = np.linspace(-span, span, n_side)
    ys = np.linspace(-span, span, n_side)
    cx, cy = np.meshgrid(xs, ys, indexing="ij")
    centers = np.column_stack([cx.ravel(), cy.ravel(), np.full(cx.size, z0)])
    n = centers.shape[0]
    offsets = np.broadcast_to(_square_loop(a), (n, 4, 3)).copy()

    oxs = np.linspace(-obs_span, obs_span, m_side)
    oys = np.linspace(-obs_span, obs_span, m_side)
    ox, oy = np.meshgrid(oxs, oys, indexing="ij")
    obs = np.column_stack([ox.ravel(), oy.ravel(), np.full(ox.size, obs_z)])
    return obs, centers, offsets


def _coil_dense_A(obs, centers, offsets):
    M, N = obs.shape[0], centers.shape[0]
    A = np.empty((M, N))
    for i in range(M):
        for j in range(N):
            A[i, j] = _loop_Hz(obs[i], centers[j], offsets[j])
    return A


def _recon_err(A, res, k=None):
    if k is None:
        k = res.k_aca
    Alr = res.U[:, :k] @ np.diag(res.S[:k]) @ res.V[:, :k].T
    return np.linalg.norm(A - Alr) / np.linalg.norm(A)


# --------------------------------------------------------------------------
# Reconstruction vs true dense A (coil kernel)
# --------------------------------------------------------------------------
def test_reconstruction_near_field_full_rank():
    """Near field is effectively full rank (k_aca = min(M,N)); the
    (ACA+)+TSVD factorization reconstructs A to machine precision."""
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-8)
    assert res.k_aca == min(M, N)
    assert _recon_err(A, res) < 1.0e-10
    assert np.all(res.S >= 0.0)
    assert np.all(np.diff(res.S[:res.k_aca]) <= 1.0e-12)


def test_low_rank_far_field_compresses():
    """Far field is smooth -> ACA+ compresses (k_aca << min(M,N))."""
    obs, centers, offsets = _coil_geometry(obs_z=0.5)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-8)
    assert res.k_aca < min(M, N)        # compression happened
    assert res.k_aca >= 1
    assert _recon_err(A, res) < 1.0e-4


def test_truncation_monotonic():
    """Reconstruction error decreases (weakly) as more TSVD modes are kept."""
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-10)
    errs = [_recon_err(A, res, k=k) for k in (5, 10, 15, res.k_aca)]
    for a, b in zip(errs, errs[1:]):
        assert b <= a + 1.0e-12


def test_qr_and_dense_methods_agree_with_svd():
    """The two supported methods (peer review JIAM-2026-36) -- and ONLY two:
      - ``method="dense"`` IS the direct TSVD of A (materialise + numpy.linalg.svd);
      - ``method="aca_qr_tsvd"`` (default, ACA+QR+TSVD) matches it to the
        ACA tolerance, with ORTHONORMAL U / V (what RegularizedTSVD requires).
    This makes the reviewer's point concrete: QR is validated against the CORRECT
    dense baseline, not a broken 'naive' one.  NO backward compatibility: the
    legacy ints 2/3, the terse "qr"/"aca", and any unknown value all RAISE."""
    obs, centers, offsets = _coil_geometry(obs_z=0.1)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    s_dense = np.linalg.svd(A, compute_uv=False)

    # method="dense" == the exact direct TSVD of A
    rd = aca_tsvd(M, N, entry, modes=min(M, N), method="dense")
    assert rd.method == "dense"
    nd = rd.modes
    assert np.linalg.norm(rd.S[:nd] - s_dense[:nd]) / np.linalg.norm(s_dense) < 1.0e-12
    assert np.linalg.norm(A - (rd.U * rd.S) @ rd.V.T) / np.linalg.norm(A) < 1.0e-12

    # method="aca_qr_tsvd" (default) matches the dense SVD to the ACA tolerance
    rq = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N), aca_eps=1.0e-10)
    assert rq.method == "aca_qr_tsvd"
    nq = rq.modes
    assert np.linalg.norm(rq.S[:nq] - s_dense[:nq]) / np.linalg.norm(s_dense) < 1.0e-8
    # orthonormal factors -- required by RegularizedTSVD.from_stiffness (W = V^T V = I)
    assert np.linalg.norm(rq.U.T @ rq.U - np.eye(nq)) < 1.0e-10, "U not orthonormal"
    assert np.linalg.norm(rq.V.T @ rq.V - np.eye(nq)) < 1.0e-10, "V not orthonormal"

    # NO backward compatibility (MCP-delivered software): the legacy ints 2/3,
    # the terse "qr"/"aca", and any unknown value all RAISE.
    for bad in (2, 3, "qr", "aca", "bogus"):
        with pytest.raises(ValueError):
            aca_tsvd(M, N, entry, modes=min(M, N), method=bad)


# --------------------------------------------------------------------------
# Pseudo-inverse solve
# --------------------------------------------------------------------------
def test_pseudo_inverse_solve_recovers_range():
    """For a target B = A @ phi0 in range(A), the least-norm solve reproduces
    B (near field is full row rank, so A phi == B exactly)."""
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-10)
    rng = np.random.default_rng(0)
    phi0 = rng.standard_normal(N)
    B = A @ phi0
    phi = pseudo_inverse_solve(res, B, k_mode=res.k_aca)
    assert phi.shape == (N,)
    assert np.linalg.norm(A @ phi - B) / np.linalg.norm(B) < 1.0e-8
    assert np.linalg.norm(phi) <= np.linalg.norm(phi0) + 1.0e-9


def test_pseudo_inverse_solve_validates_B_length():
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    res = aca_tsvd(M, N, lambda i, j: A[i, j], modes=min(M, N),
                   kmax=min(M, N), aca_eps=1.0e-8)
    with pytest.raises(ValueError):
        pseudo_inverse_solve(res, np.zeros(M + 1))


def test_aca_tsvd_validates_args():
    with pytest.raises(ValueError):
        aca_tsvd(0, 5, lambda i, j: 0.0, modes=1, kmax=1)
    with pytest.raises(TypeError):
        aca_tsvd(5, 5, 123, modes=1, kmax=1)  # not callable


def test_aca_tsvd_rethrows_entry_callback_failure():
    def failing_entry(i, j):
        raise RuntimeError(f"entry failure at {i},{j}")

    with pytest.raises(RuntimeError, match="entry failure"):
        aca_tsvd(4, 4, failing_entry, modes=2, kmax=2, aca_eps=1.0e-8)


# --------------------------------------------------------------------------
# Improved DUCAS / Abe node-current-potential selection
# --------------------------------------------------------------------------
def test_abe_distance_weights_are_nearest_distance_squared_not_area():
    nodes = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]])
    fields = np.array([[0.0, 1.0], [4.0, 0.0]])
    got = abe_nearest_field_distance_scales(nodes, fields, normalize=False)
    assert np.allclose(got, [1.0, 4.0, 4.0])


def test_abe_distance_weights_reject_touching_source_and_field_meshes():
    with pytest.raises(ValueError, match="coincides"):
        abe_nearest_field_distance_scales([[0.0, 0.0]], [[0.0, 0.0]])
    got = abe_nearest_field_distance_scales(
        [[0.0, 0.0]], [[0.0, 0.0]], distance_floor=0.25,
        normalize=False)
    assert got.tolist() == pytest.approx([0.25 ** 2])


def test_abe_node_scale_reduction_matches_weighted_component_average():
    """Abe eqs. 24--25: tied full nodes share their averaged delta scale."""
    sparse = pytest.importorskip("scipy.sparse")
    R = sparse.csr_matrix(np.array([
        [1.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.0, 1.0],
    ]))
    got = abe_reduce_node_potential_scales(R, [1.0, 3.0, 2.0, 4.0])
    assert np.allclose(got, [2.0 / 3.0, 1.0])


def test_abe_weighted_initial_potential_matches_dense_formula():
    """The public solve is exactly Abe eqs. 9--16, including T0 and W_B."""
    A = np.array([
        [2.0, -1.0, 0.5, 0.0],
        [0.0, 1.5, -0.5, 1.0],
        [1.0, 0.0, 0.0, -2.0],
    ])
    R = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    delta = np.array([0.4, 1.0, 0.7])
    wb = np.array([2.0, 0.5, 1.5])
    initial = np.array([0.2, -0.1, -0.1, 0.3])
    external = np.array([0.05, -0.02, 0.01])
    target = np.array([0.8, -0.4, 0.6])
    weighted = wb[:, None] * (A @ R) * delta[None, :]
    rhs = wb * (target - external - A @ initial)
    q = np.linalg.pinv(weighted, rcond=1.0e-14) @ rhs
    expected = initial + R @ (delta * q)
    result = solve_abe_current_potential(
        A, target, reduction=R, field_weights=wb,
        independent_potential_scales=delta, initial_potential=initial,
        external_field=external, method="dense")
    assert np.linalg.norm(result.potential - expected) < 1.0e-11
    assert np.linalg.norm(result.reconstructed_field -
                          (external + A @ expected)) < 1.0e-12
    assert result.converged


def test_abe_mode_strength_selection_is_not_a_contiguous_prefix():
    """Strong mode 3 survives while weak mode 2 is omitted (Abe eq. 19)."""
    A = np.diag([9.0, 5.0, 2.0])
    target = np.array([1.0, 1.0e-3, 1.0])
    result = solve_abe_current_potential(
        A, target, minimum_mode_strength=0.1, method="dense")
    assert result.selected_modes.tolist() == [0, 2]
    assert np.allclose(result.reconstructed_field, [1.0, 0.0, 1.0],
                       atol=1.0e-13)
    assert result.diagnostics.rejection_reason[1] == \
        "mode_strength_below_threshold"


def test_abe_initial_potential_retains_unselected_high_order_content():
    """T0's unselected content remains while selected modes correct B."""
    A = np.diag([9.0, 5.0, 2.0])
    initial = np.array([0.0, 0.4, 0.0])
    target = A @ initial + np.array([9.0, 0.0, 4.0])
    result = solve_abe_current_potential(
        A, target, initial_potential=initial,
        allowed_modes=[0, 2], method="dense")
    assert result.selected_modes.tolist() == [0, 2]
    assert np.allclose(result.potential, [1.0, 0.4, 2.0], atol=1.0e-13)


def test_abe_residual_target_stops_before_unneeded_modes():
    """Residual quality, not a fixed k, determines the accumulated modes."""
    A = np.diag([9.0, 5.0, 2.0])
    target = np.array([1.0, 1.0e-2, 1.0e-3])
    result = solve_abe_current_potential(
        A, target, residual_rms=1.0e-2, method="dense")
    assert result.selected_modes.tolist() == [0]
    assert result.converged
    assert result.stop_reason == "residual_target_met"
    assert result.diagnostics.rejection_reason[1] == \
        "not_needed_after_residual_target"


def test_abe_potential_limit_fails_loudly_after_field_fit():
    """A field fit outside the engineering current limit is not 'converged'."""
    A = np.diag([9.0, 5.0, 2.0])
    result = solve_abe_current_potential(
        A, [18.0, 0.0, 0.0], maximum_abs_potential=1.0,
        method="dense")
    assert not result.converged
    assert result.stop_reason == "maximum_abs_potential_exceeded"
    assert result.peak_abs_potential == pytest.approx(2.0)


def test_abe_precomputed_factor_reuses_identical_weighted_operator():
    A = np.diag([9.0, 5.0, 2.0])
    first = solve_abe_current_potential(A, [1.0, 0.2, -0.3], method="dense")
    second = solve_abe_current_potential(
        A, [-0.4, 0.7, 0.1], precomputed_factor=first.factor,
        method="dense")
    direct = solve_abe_current_potential(A, [-0.4, 0.7, 0.1], method="dense")
    assert second.factor is first.factor
    assert np.linalg.norm(second.potential - direct.potential) < 1.0e-13


def test_abe_rejects_malformed_precomputed_factor_and_fractional_modes():
    A = np.diag([3.0, 2.0])
    malformed = StreamTSVD(
        U=np.eye(2), S=np.array([]), V=np.eye(2), k_aca=0, method="dense")
    with pytest.raises(ValueError, match="at least one mode"):
        solve_abe_current_potential(
            A, [1.0, 0.0], precomputed_factor=malformed)
    with pytest.raises(ValueError, match="integer mode indices"):
        solve_abe_current_potential(
            A, [1.0, 0.0], allowed_modes=[0.5], method="dense")


def test_abe_bounded_iteration_requires_a_physical_acceptance_target():
    with pytest.raises(ValueError, match="acceptance criterion"):
        solve_abe_bounded_current_potential(
            np.eye(2), [1.0, 0.0], lower_potential=0.0,
            method="dense")


def test_abe_positive_only_iteration_redistributes_clipped_solution():
    """Abe's positive-iron loop repeatedly solves the post-clip error."""
    A = np.array([[1.0, -1.0]])
    result = solve_abe_bounded_current_potential(
        A, [1.0], lower_potential=0.0, residual_rms=1.0e-9,
        max_iterations=64, method="dense")
    assert result.converged
    assert result.stop_reason == "bounded_residual_target_met"
    assert result.iterations > 1
    assert np.all(result.solution.potential >= 0.0)
    assert result.solution.potential[0] == pytest.approx(1.0, abs=2.0e-9)
    assert result.solution.potential[1] == pytest.approx(0.0, abs=1.0e-14)
    assert np.any(result.clipped_dof_history > 0)


def test_abe_bounded_iteration_reports_unattainable_capacity():
    result = solve_abe_bounded_current_potential(
        np.array([[1.0]]), [2.0], lower_potential=0.0,
        upper_potential=1.0, residual_rms=1.0e-10,
        max_iterations=10, method="dense")
    assert not result.converged
    assert result.stop_reason == "bounded_stagnation"
    assert result.solution.potential.tolist() == [1.0]
    assert result.solution.residual_rms == pytest.approx(1.0)


def test_abe_bounded_iteration_rejects_post_clip_constraint_violation():
    with pytest.raises(ValueError, match="already reduced response"):
        solve_abe_bounded_current_potential(
            np.eye(2), [1.0, 0.0], reduction=np.eye(2),
            lower_potential=0.0, method="dense")


# --------------------------------------------------------------------------
# RegularizedTSVD -- regularisation-folded pseudo-inverse
#   phi = (S^-1 V) . W^-1 . diag(1/Sigma) . U^T B,  W = V^T S^-1 V (k x k)
# Solves  min phi^T S phi  s.t.  A phi = B  with the factorisation cached.
# --------------------------------------------------------------------------
def _random_dense_problem(M=12, N=90, seed=0):
    """A random over-parameterised system + a target in range(A)."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((M, N)) / np.sqrt(N)
    B = rng.standard_normal(M)
    res = aca_tsvd(M, N, lambda i, j: float(A[i, j]),
                   modes=M, kmax=min(M, N), aca_eps=1.0e-13)
    return A, B, res


def _spd(N, seed):
    rng = np.random.default_rng(seed)
    G = rng.standard_normal((N, N))
    return G @ G.T + 0.1 * np.eye(N)


def test_regularized_identity_matches_standard_pseudo_inverse():
    """S = I reduces to the plain min-Euclidean-norm pseudo-inverse:
    W = V^T V = I_k (V columns orthonormal from TSVD), so the formula
    collapses to phi = V Sigma^-1 U^T B == pseudo_inverse_solve."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    phi_std = pseudo_inverse_solve(res, B)
    phi_reg = pseudo_inverse_solve_regularized(res, B, np.eye(N))
    rel = np.linalg.norm(phi_reg - phi_std) / np.linalg.norm(phi_std)
    assert rel < 1.0e-10


def test_regularized_matches_direct_lagrangian_spd():
    """For arbitrary SPD S the cached form matches the direct Lagrangian
    solution  phi = S^-1 A^T (A S^-1 A^T)^-1 B  to ACA+ tolerance."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=1)
    phi_reg = pseudo_inverse_solve_regularized(res, B, S)
    # Direct reference
    Sinv_AT = np.linalg.solve(S, A.T)
    small = A @ Sinv_AT
    lam = np.linalg.solve(small, B)
    phi_ref = Sinv_AT @ lam
    rel = np.linalg.norm(phi_reg - phi_ref) / np.linalg.norm(phi_ref)
    assert rel < 1.0e-7


def test_regularized_satisfies_constraint():
    """A phi = B holds (the constraint is exact for full-row-rank A)."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=2)
    phi = pseudo_inverse_solve_regularized(res, B, S)
    assert np.linalg.norm(A @ phi - B) / np.linalg.norm(B) < 1.0e-7


def test_regularized_minimises_S_seminorm():
    """The regularised solution has the SMALLEST phi^T S phi among all
    phi satisfying A phi = B.  Check against the standard min-Euclidean
    solution (which minimises phi^T phi, not phi^T S phi): for a
    non-identity S the regularised phi^T S phi must be <= the
    Euclidean one's phi^T S phi."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=3)
    phi_reg = pseudo_inverse_solve_regularized(res, B, S)
    phi_l2 = pseudo_inverse_solve(res, B)
    sn_reg = float(phi_reg @ S @ phi_reg)
    sn_l2 = float(phi_l2 @ S @ phi_l2)
    assert sn_reg <= sn_l2 + 1.0e-9
    # ... and the random-perturbation-in-nullspace test: adding any
    # nullspace component to phi_reg can only increase phi^T S phi.
    rng = np.random.default_rng(9)
    # build a nullspace vector: remove r's component in row(A) = range(A^T).
    # coeff solves min ||A^T coeff - r||; A (r - A^T coeff) = 0 since
    # coeff = (A A^T)^-1 A r is the least-squares solution.
    r = rng.standard_normal(N)
    coeff, *_ = np.linalg.lstsq(A.T, r, rcond=None)
    null_vec = r - A.T @ coeff
    assert np.linalg.norm(A @ null_vec) < 1.0e-8   # genuinely in nullspace
    sn_perturbed = float((phi_reg + 0.3 * null_vec) @ S @ (phi_reg + 0.3 * null_vec))
    assert sn_perturbed >= sn_reg - 1.0e-9


def test_regularized_cache_reuse_across_B():
    """RegularizedTSVD.from_stiffness factorises ONCE; .solve(B) reused
    across many right-hand sides gives the same answer as the one-shot
    helper (this is the Path-A iteration pattern)."""
    A, _, res = _random_dense_problem()
    N, M = A.shape[1], A.shape[0]
    S = _spd(N, seed=4)
    reg = RegularizedTSVD.from_stiffness(res, S)
    rng = np.random.default_rng(5)
    for _ in range(5):
        B = rng.standard_normal(M)
        phi_cached = reg.solve(B)
        phi_oneshot = pseudo_inverse_solve_regularized(res, B, S)
        assert np.linalg.norm(phi_cached - phi_oneshot) < 1.0e-12


def test_regularized_k_mode_truncation():
    """Truncating to k_mode < modes re-inverts the k_mode x k_mode block
    of W (W_inv depends on k).  The truncated solve must still satisfy
    the projected constraint and reduce to the right rank."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=6)
    reg = RegularizedTSVD.from_stiffness(res, S)
    k = max(1, res.modes // 2)
    phi_k = reg.solve(B, k_mode=k)
    assert phi_k.shape == (N,)
    # The truncated solve uses only the top-k singular directions; the
    # residual on the retained directions (U[:, :k]^T (A phi - B)) ~ 0.
    Uk = res.U[:, :k]
    proj_res = Uk.T @ (A @ phi_k - B)
    assert np.linalg.norm(proj_res) / (np.linalg.norm(B) + 1e-30) < 1.0e-6


def test_regularized_sparse_S_matches_dense():
    """The scipy-sparse S path (splu factor-once) matches the dense path
    (np.linalg.solve) to ~1e-9 on a well-conditioned S (same S, two
    representations).  The PRODUCTION folded seminorm is well-conditioned
    (cond(W) ~ 1e2..1e7); the dangerous-ill-conditioning case is locked
    separately in test_regularized_rejects_singular_and_illconditioned_S."""
    sparse = pytest.importorskip("scipy.sparse")
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=7)
    phi_dense = pseudo_inverse_solve_regularized(res, B, S)
    phi_sparse = pseudo_inverse_solve_regularized(res, B, sparse.csr_matrix(S))
    assert np.linalg.norm(phi_dense - phi_sparse) / np.linalg.norm(phi_dense) < 1.0e-9


def test_regularized_validates_B_length():
    A, _, res = _random_dense_problem()
    N = A.shape[1]
    reg = RegularizedTSVD.from_stiffness(res, np.eye(N))
    with pytest.raises(ValueError):
        reg.solve(np.zeros(res.M + 1))


def test_regularized_validates_S_shape():
    """from_stiffness rejects an S whose size does not match V's rows."""
    A, _, res = _random_dense_problem()
    N = A.shape[1]
    with pytest.raises(ValueError):
        RegularizedTSVD.from_stiffness(res, np.eye(N + 3))


def test_regularized_rejects_singular_and_illconditioned_S():
    """No-Fallback: from_stiffness must RAISE (not silently degrade) on an S
    it cannot use, and the happy path must not regress.

    (1) Singular S -- the dense path raises numpy LinAlgError and the sparse
        splu path raises RuntimeError; BOTH are normalised to ONE ValueError
        mentioning 'singular'.
    (2) A catastrophically ill-conditioned core W = V^T S^-1 V (cond > 1e12,
        where the exact-fit solve would silently lose most digits) raises a
        ValueError mentioning 'ill-conditioned'.  (Production seminorms sit at
        cond(W) ~ 1e2..1e7 and pass.)
    (3) A well-conditioned S past the guards still gives an EXACT data fit to
        machine precision via the factored W-solve (the explicit inverse is
        gone)."""
    sparse = pytest.importorskip("scipy.sparse")
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    V = res.V                                         # (N, k)

    # (1) singular S: dense LinAlgError + sparse RuntimeError -> one ValueError
    S_sing = np.zeros((N, N))
    with pytest.raises(ValueError, match="singular"):
        RegularizedTSVD.from_stiffness(res, S_sing)
    with pytest.raises(ValueError, match="singular"):
        RegularizedTSVD.from_stiffness(res, sparse.csr_matrix(S_sing))

    # (2) SPD S acting as diag(d), d in [1, 1e-13], on span(V) (so
    #     W = V^T S^-1 V has cond ~ 1e13 > 1e12) and identity off V.
    d = np.logspace(0.0, -13.0, V.shape[1])
    P = V @ V.T
    S_ill = V @ np.diag(d) @ V.T + (np.eye(N) - P)
    S_ill = 0.5 * (S_ill + S_ill.T)
    with pytest.raises(ValueError, match="ill-conditioned"):
        RegularizedTSVD.from_stiffness(res, S_ill)

    # (3) well-conditioned S: exact data fit (alpha=0) to machine precision.
    reg = RegularizedTSVD.from_stiffness(res, _spd(N, seed=11))
    phi = reg.solve(B)
    assert np.linalg.norm(A @ phi - B) / np.linalg.norm(B) < 1.0e-9


# --------------------------------------------------------------------------
# RegularizedTSVD -- folded Tikhonov (alpha > 0, soft data fit)
#   phi(alpha) = (S^-1 V) . (alpha I + Sigma^2 W)^-1 . Sigma . U^T B
#              == (A^T A + alpha S)^-1 A^T B
# Same cached factorisation as alpha=0; only an alpha I added to the core.
# --------------------------------------------------------------------------
def test_tikhonov_matches_dense_spd():
    """solve(B, alpha>0) == dense (A^T A + alpha S)^-1 A^T B for SPD S,
    across a range of alpha -- to ACA+ tolerance."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=11)
    reg = RegularizedTSVD.from_stiffness(res, S)
    AtA, AtB = A.T @ A, A.T @ B
    for alpha in (1.0e-3, 1.0e-1, 1.0, 1.0e1):
        phi_fold = reg.solve(B, alpha=alpha)
        phi_dense = np.linalg.solve(AtA + alpha * S, AtB)
        rel = np.linalg.norm(phi_fold - phi_dense) / np.linalg.norm(phi_dense)
        assert rel < 1.0e-6, f"alpha={alpha}: rel={rel:.2e}"


def test_tikhonov_identity_reduces_to_filter_factors():
    """For S = I the core collapses to the classic Tikhonov filter
    factors:  phi = V . diag(sigma/(sigma^2+alpha)) . U^T B  (exact, since
    W = V^T V = I and S^-1 V = V -- no S-solve, machine precision)."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    reg = RegularizedTSVD.from_stiffness(res, np.eye(N))
    U, Sig, V = res.U, res.S, res.V
    UtB = U.T @ B
    for alpha in (1.0e-3, 1.0e-1, 1.0, 1.0e1):
        phi_fold = reg.solve(B, alpha=alpha)
        phi_filter = V @ ((Sig / (Sig ** 2 + alpha)) * UtB)
        rel = np.linalg.norm(phi_fold - phi_filter) / np.linalg.norm(phi_filter)
        assert rel < 1.0e-10, f"alpha={alpha}: rel={rel:.2e}"


def test_tikhonov_alpha_zero_recovers_exact_fit():
    """alpha=0 returns the exact-fit min-seminorm solution (== solve(B));
    and alpha << sigma_min^2 converges to it."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=12)
    reg = RegularizedTSVD.from_stiffness(res, S)
    phi0 = reg.solve(B)
    assert np.linalg.norm(reg.solve(B, alpha=0.0) - phi0) < 1.0e-12
    alpha_tiny = (float(res.S.min()) ** 2) * 1.0e-6
    rel = np.linalg.norm(reg.solve(B, alpha=alpha_tiny) - phi0) \
        / np.linalg.norm(phi0)
    assert rel < 1.0e-4, f"alpha_tiny rel={rel:.2e}"


def test_tikhonov_lcurve_monotone():
    """Increasing alpha trades data misfit UP for seminorm DOWN -- the
    L-curve / Pareto monotonicity that makes the alpha-sweep a front."""
    A, B, res = _random_dense_problem()
    N = A.shape[1]
    S = _spd(N, seed=13)
    reg = RegularizedTSVD.from_stiffness(res, S)
    alphas = [1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1]
    misfit = [float(np.linalg.norm(A @ reg.solve(B, alpha=a) - B))
              for a in alphas]
    seminorm = [float(reg.solve(B, alpha=a) @ S @ reg.solve(B, alpha=a))
                for a in alphas]
    for i in range(len(alphas) - 1):
        assert misfit[i] <= misfit[i + 1] + 1.0e-9
        assert seminorm[i] >= seminorm[i + 1] - 1.0e-9


# --------------------------------------------------------------------------
# Generic path through Radia's OWN field computation (magnetic materials)
# --------------------------------------------------------------------------
def test_radia_field_kernel_magnets():
    """The same (ACA+)+TSVD machinery works with Radia's fixed-magnet field over
    permanent-magnet objects (not coils), via radia_field_kernel + radia.Fld.
    Proves the solver is kernel-agnostic and reuses Radia's existing kernels."""
    import radia as rad
    rad.UtiDelAll()
    try:
        # N small permanent-magnet cubes on a plane, M observation points above.
        xs = np.linspace(-0.04, 0.04, 6)
        ys = np.linspace(-0.04, 0.04, 6)
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        centers = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
        N = centers.shape[0]
        sources = [rad.magnet_box(c.tolist(), [0.01, 0.01, 0.01], [0.0, 0.0, 954930.0])
                   for c in centers]

        oxs = np.linspace(-0.03, 0.03, 4)
        oys = np.linspace(-0.03, 0.03, 4)
        ox, oy = np.meshgrid(oxs, oys, indexing="ij")
        obs = np.column_stack([ox.ravel(), oy.ravel(), np.full(ox.size, 0.06)])
        M = obs.shape[0]

        entry = radia_field_kernel(obs, sources, component=2, field="b")
        A = np.array([[entry(i, j) for j in range(N)] for i in range(M)])

        res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                       aca_eps=1.0e-8)
        assert 1 <= res.k_aca <= min(M, N)
        assert _recon_err(A, res) < 1.0e-5
    finally:
        rad.UtiDelAll()


# --------------------------------------------------------------------------
# f90 cross-check (LAB only) -- runs in a fresh subprocess to dodge conftest's
# DLL-search pollution (the f2py module bundles its own Intel/MKL DLLs and
# cannot be imported alongside conftest's preloaded ngsolve/cubit).
# Verified bit-exact: k_aca identical, ||S_f90 - S_radia||/||S_f90|| ~ 1e-15.
# --------------------------------------------------------------------------
_F90_SUBPROCESS = r"""
import os, sys
import numpy as np
import radia  # noqa: F401  (loads radia + ngsolve MKL first, as in production)
from radia.stream_function import aca_tsvd

ref = sys.argv[1]
try:
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(ref)
    if ref not in sys.path:
        sys.path.insert(0, ref)
    import coil_solver
except Exception as exc:  # noqa: BLE001
    print("SKIP", type(exc).__name__)
    sys.exit(0)

def seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1; R21 = float(d @ d)
    o1 = O - P1; o2 = O - P2
    OR1 = np.sqrt(o1 @ o1); OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d); O221 = float(o2 @ d)
    Rc12 = O121 / OR1; Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f

def loop_Hz(obs, center, offsets):
    C = center + offsets
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4; hz += seg_Hz(obs, C[c], C[n])
    Cl = C.copy(); Cl[:, 2] = -Cl[:, 2]
    for c in range(4):
        n = (c + 1) % 4; hz += seg_Hz(obs, Cl[c], Cl[n])
    return hz

n, m, z0, a, span, ospan, obs_z = 10, 5, 0.005, 0.02, 0.1, 0.08, 0.1
xs = np.linspace(-span, span, n); ys = np.linspace(-span, span, n)
cx, cy = np.meshgrid(xs, ys, indexing="ij")
centers = np.column_stack([cx.ravel(), cy.ravel(), np.full(cx.size, z0)])
N = centers.shape[0]
loop = np.array([[-a/2,-a/2,0.],[a/2,-a/2,0.],[a/2,a/2,0.],[-a/2,a/2,0.]])
off = np.broadcast_to(loop, (N, 4, 3)).copy()
oxs = np.linspace(-ospan, ospan, m); oys = np.linspace(-ospan, ospan, m)
ox, oy = np.meshgrid(oxs, oys, indexing="ij")
obs = np.column_stack([ox.ravel(), oy.ravel(), np.full(ox.size, obs_z)])
M = obs.shape[0]; kmax = min(M, N); modes = kmax; eps = 1.0e-8

A = np.empty((M, N))
for i in range(M):
    for j in range(N):
        A[i, j] = loop_Hz(obs[i], centers[j], off[j])
entry = lambda i, j: A[i, j]

fails = []
# radia has ONE recompression (aca_qr_tsvd) -> the SAME SVD as BOTH f90 methods
# (all recompress the same ACA product), so compute it once.
r = aca_tsvd(M, N, entry, modes=modes, kmax=kmax, aca_eps=eps)
for f90name in ("method_aca_tsvd_2", "method_aca_tsvd_1"):
    fn = getattr(coil_solver, f90name)
    Uf, Sf, Vf, kf, _t1, _t2, _pk = fn(
        np.ascontiguousarray(obs[:, 0]), np.ascontiguousarray(obs[:, 1]),
        np.ascontiguousarray(obs[:, 2]), np.ascontiguousarray(centers[:, 0]),
        np.ascontiguousarray(centers[:, 1]), np.ascontiguousarray(centers[:, 2]),
        np.ascontiguousarray(off[:, :, 0]), np.ascontiguousarray(off[:, :, 1]),
        np.ascontiguousarray(off[:, :, 2]), modes, kmax, eps)
    ns = min(int(kf), r.k_aca)
    rel = np.linalg.norm(Sf[:ns] - r.S[:ns]) / np.linalg.norm(Sf[:ns])
    if int(kf) != r.k_aca or rel > 1.0e-9:
        fails.append(f"{f90name}: kf={kf} k_aca={r.k_aca} relS={rel:.2e}")

if fails:
    print("FAIL", "; ".join(fails)); sys.exit(1)
print("PASS"); sys.exit(0)
"""


@pytest.mark.skipif(not os.path.isdir(REF_DIR),
                    reason="f2py coil_solver reference (W: drive) not available")
def test_matches_f90_reference():
    """radia C++ port matches the f2py coil_solver.f90 reference exactly
    (same k_aca, same singular values) for BOTH method pairs, using a coil
    Biot-Savart kernel.  Runs in a fresh subprocess to dodge conftest's
    DLL-search pollution."""
    import subprocess
    proc = subprocess.run([sys.executable, "-c", _F90_SUBPROCESS, REF_DIR],
                          capture_output=True, text=True, timeout=300)
    out = (proc.stdout + "\n" + proc.stderr).strip()
    if proc.stdout.startswith("SKIP"):
        pytest.skip(f"f2py coil_solver could not be imported: {proc.stdout.strip()}")
    assert "PASS" in proc.stdout, f"f90 cross-check failed (rc={proc.returncode}):\n{out}"
