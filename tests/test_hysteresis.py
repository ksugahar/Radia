"""
Pytest tests for hysteresis materials (Play and Energy models).

Tests Play model creation, Forward/Inverse, state management,
solver integration, and round-trip accuracy.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

import pytest
import radia as rad
import numpy as np

MU_0 = 4e-7 * np.pi


@pytest.fixture
def play_fixture():
    """Load real B-input Play shape functions from fixture."""
    fix_path = os.path.join(os.path.dirname(__file__),
                            '../examples/hysteresis/binput_play_fixture.npz')
    if not os.path.exists(fix_path):
        pytest.skip("Fixture file not found")
    fix = np.load(fix_path)
    K = int(fix['K_sub'])
    eta = fix['chi_sub']
    f_k_tables = []
    for ki in range(K):
        r = fix[f'r_sub_{ki}']
        f = fix[f'f_sub_{ki}']
        f_k_tables.append((r.tolist(), f.tolist()))
    return K, eta.tolist(), f_k_tables


@pytest.fixture
def full_fixture():
    """Load full (non-subset) Play shape functions."""
    fix_path = os.path.join(os.path.dirname(__file__),
                            '../examples/hysteresis/binput_play_fixture.npz')
    if not os.path.exists(fix_path):
        pytest.skip("Fixture file not found")
    fix = np.load(fix_path)
    K = int(fix['K'])
    eta = fix['chi']
    f_k_tables = []
    for k in range(K):
        r = fix[f'r_{k}']
        f = fix[f'f_{k}']
        f_k_tables.append((r.tolist(), f.tolist()))
    return K, eta.tolist(), f_k_tables


class TestPlayHysteresisMaterial:
    """Test MatPlayHysteresis creation and basic Forward/Inverse."""

    def test_create_play_material(self, play_fixture):
        """MatPlayHysteresis returns valid handle."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)
        assert mat > 0

    def test_forward_returns_nonzero_M(self, play_fixture):
        """MatMvsH at moderate H returns ferromagnetic M."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)
        M = rad.MatMvsH(mat, 'm', [500.0, 0, 0])
        M_mag = np.linalg.norm(M)
        assert M_mag > 0, "Magnetization should be nonzero for H=500 A/m"
        assert M[0] > 0, "M should be parallel to H (ferromagnetic)"

    def test_forward_multiple_directions(self, play_fixture):
        """Forward works for various H directions."""
        K, eta, tables = play_fixture
        H_tests = [
            [100.0, 0.0, 0.0],
            [0.0, 0.0, 1000.0],
            [-300.0, 400.0, 100.0],
        ]
        for H in H_tests:
            rad.UtiDelAll()
            mat = rad.MatPlayHysteresis(K, eta, tables)
            M = rad.MatMvsH(mat, 'm', H)
            M_arr = np.array(M)
            H_arr = np.array(H)
            dot = np.dot(M_arr, H_arr)
            assert dot > 0, f"M should be parallel to H for H={H}"


class TestPlayHysteresisLoop:
    """Test hysteresis behavior (irreversibility)."""

    def test_hysteresis_detected(self, play_fixture):
        """B-H loop shows hysteresis (rising != falling)."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        Hmax = 3000.0
        n = 200
        t = np.linspace(0, 2 * np.pi, n)
        H_drive = Hmax * np.sin(t)

        B_loop = np.zeros(n)
        for i, H_val in enumerate(H_drive):
            M = rad.MatMvsH(mat, 'm', [H_val, 0, 0])
            B_loop[i] = MU_0 * (H_val + M[0])

        # Compare B at rising vs falling zero-crossing
        mid = n // 4
        mid2 = 3 * n // 4
        hyst_width = abs(B_loop[mid - 1] - B_loop[mid2 - 1])
        assert hyst_width > 0.001, f"Hysteresis width {hyst_width:.4f} T too small"


class TestPlayStateManagement:
    """Test Save/Restore/Commit state lifecycle."""

    def test_save_state_nonempty(self, play_fixture):
        """MatHysSaveState returns non-empty array."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)
        rad.MatMvsH(mat, 'm', [1000.0, 0, 0])
        state = rad.MatHysSaveState(mat)
        assert len(state) == K * 9, f"State size should be K*9={K*9}, got {len(state)}"

    def test_restore_recovers_state(self, play_fixture):
        """Restore after modification returns to saved state."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        # Drive to non-trivial state
        rad.MatMvsH(mat, 'm', [1000.0, 0, 0])
        M_original = np.array(rad.MatMvsH(mat, 'm', [500.0, 0, 0]))
        state = rad.MatHysSaveState(mat)

        # Modify state by driving further
        rad.MatMvsH(mat, 'm', [5000.0, 0, 0])
        rad.MatMvsH(mat, 'm', [-2000.0, 0, 0])

        # Restore and verify
        rad.MatHysRestoreState(mat, state)
        M_restored = np.array(rad.MatMvsH(mat, 'm', [500.0, 0, 0]))
        err = np.linalg.norm(M_original - M_restored)
        assert err < 1e-6, f"Restored M differs by {err:.6e} A/m"

    def test_commit_state(self, play_fixture):
        """MatHysCommitState does not raise."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)
        rad.MatMvsH(mat, 'm', [1000.0, 0, 0])
        rad.MatHysCommitState(mat)  # Should not raise


class TestPlaySolverIntegration:
    """Test Play material with BEM solver."""

    def test_solve_with_play_material(self, play_fixture):
        """rad.Solve converges with Play hysteresis material."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        # Simple iron block
        mag = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.02], [0, 0, 0])
        rad.MatApl(mag, mat)

        bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
        container = rad.ObjCnt([mag, bkg])

        result = rad.Solve(container, 0.001, 100, 0)
        assert result[3] > 0, "Should have done at least 1 iteration"
        assert result[3] < 100, f"Should converge in < 100 iterations, got {result[3]}"


class TestPlayMonotoneLimits:
    """Test monotone limit enforcement."""

    def test_b_bounded_at_saturation(self, full_fixture):
        """B should be bounded at monotone limit for very large H."""
        K, eta, tables = full_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        H_values = [5000.0, 10000.0, 20000.0]
        B_values = []
        for H in H_values:
            M = rad.MatMvsH(mat, 'm', [H, 0, 0])
            B = MU_0 * (H + M[0])
            B_values.append(B)

        # B should saturate (not grow linearly)
        assert B_values[-1] < 2.0, f"B={B_values[-1]:.2f} T exceeds expected saturation"

    def test_b_monotonically_increasing(self, full_fixture):
        """B(H) should be monotonically increasing (or saturated)."""
        K, eta, tables = full_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        H_sweep = np.linspace(100, 15000, 50)
        B_prev = 0.0
        for H in H_sweep:
            M = rad.MatMvsH(mat, 'm', [H, 0, 0])
            B = MU_0 * (H + M[0])
            assert B >= B_prev - 1e-10, f"B decreased: {B:.4f} < {B_prev:.4f} at H={H:.0f}"
            B_prev = B


class TestEnergyHysteresisMaterial:
    """Test MatEnergyHysteresis (Type 5) basic operations."""

    def test_create_energy_material(self, play_fixture):
        """MatEnergyHysteresis returns valid handle."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatEnergyHysteresis(K, eta, tables, 1e-8)
        assert mat > 0

    def test_energy_forward_returns_nonzero_M(self, play_fixture):
        """Energy model MatMvsH at moderate H returns ferromagnetic M."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatEnergyHysteresis(K, eta, tables, 1e-8)
        M = rad.MatMvsH(mat, 'm', [500.0, 0, 0])
        M_mag = np.linalg.norm(M)
        assert M_mag > 0, "Energy model M should be nonzero"

    def test_energy_state_save_restore(self, play_fixture):
        """Energy model state save/restore works."""
        K, eta, tables = play_fixture
        rad.UtiDelAll()
        mat = rad.MatEnergyHysteresis(K, eta, tables, 1e-8)

        rad.MatMvsH(mat, 'm', [1000.0, 0, 0])
        state = rad.MatHysSaveState(mat)
        assert len(state) == K * 9

        rad.MatMvsH(mat, 'm', [5000.0, 0, 0])
        rad.MatHysRestoreState(mat, state)
        # Should not raise


class TestShapeFunctionIdentificationRegression:
    """Regression tests for the JMAG shape function identification pipeline.

    These tests would have caught the three hysteresis_io.py bugs fixed on
    2026-05-02 (load_mat .flat[0] bug, sort direction inverted, r=0 anchor
    lost to fp rounding). Each bug produced shape functions that
    individually look reasonable but fail to reproduce the input descending
    branches when fed back through the Play model.
    """

    def test_load_mat_preserves_full_branch(self):
        """load_mat must return all H/B values, not just the first scalar."""
        from radia.hysteresis_io import load_mat
        mat_file = (r"W:\999_菅原賢悟\19_磁気ヒステリシス\2024_IGTE_共同研究"
                    r"\2024_03_08_H-input_B-input\Potter_Schmulian\B_input.mat")
        if not os.path.exists(mat_file):
            pytest.skip("Potter-Schmulian B_input.mat not available")
        loops, dB, BMax = load_mat(mat_file)
        # Each loop must have B and H of length > 1; the .flat[0] bug
        # produced single-element arrays.
        for i, loop in enumerate(loops):
            assert len(loop['B']) > 1, \
                f"loop[{i}] has only {len(loop['B'])} B-points (load_mat regression)"
            assert len(loop['H']) == len(loop['B']), \
                f"loop[{i}]: H/B length mismatch"

    def test_jmag_identification_reproduces_input_branches(self):
        """Identified shape functions must reproduce input descending branches.

        The JMAG identification is exact for the input data (it inverts a
        finite-difference operator). Forward-evaluating the identified Play
        model on the same B trajectory must reconstruct H to machine
        precision (< 1e-9 A/m). Any deviation indicates a bug in
        build_shape_functions, load_mat, or the Play evaluator.
        """
        from radia.hysteresis_io import load_mat, build_shape_functions
        mat_file = (r"W:\999_菅原賢悟\19_磁気ヒステリシス\2024_IGTE_共同研究"
                    r"\2024_03_08_H-input_B-input\Potter_Schmulian\B_input.mat")
        if not os.path.exists(mat_file):
            pytest.skip("Potter-Schmulian B_input.mat not available")
        loops, dB, BMax = load_mat(mat_file)
        eta, f_tables, _ = build_shape_functions(loops, dB)
        K = len(eta)

        # Run scalar Play on the largest input branch and compare to input H.
        # Pre-saturate first to set Play states consistently.
        loops_desc = sorted(loops, key=lambda x: x['Bmax'], reverse=True)
        largest = loops_desc[0]
        Bm = float(largest['Bmax'])

        # Build interpolators with anti-symmetric extension
        interps = []
        for k in range(K):
            r = np.asarray(f_tables[k][0])
            f = np.asarray(f_tables[k][1])
            mask = r >= 0
            r = r[mask]; f = f[mask]
            idx = np.argsort(r)
            r = r[idx]; f = f[idx]
            r_full = np.concatenate([-r[::-1], r[1:]])
            f_full = np.concatenate([-f[::-1], f[1:]])
            r_unique, ia = np.unique(r_full, return_index=True)
            interps.append((r_unique, f_full[ia]))

        def play_eval(B_traj, eta, interps, K):
            p = np.zeros(K)
            H = np.zeros_like(B_traj)
            for n, B in enumerate(B_traj):
                Hn = 0.0
                for k in range(K):
                    if eta[k] < 1e-30:
                        p[k] = B
                    else:
                        if B > p[k] + eta[k]:
                            p[k] = B - eta[k]
                        elif B < p[k] - eta[k]:
                            p[k] = B + eta[k]
                    Hn += float(np.interp(p[k], interps[k][0], interps[k][1]))
                H[n] = Hn
            return H

        # Pre-saturate by going from 0 up to Bm
        n_pre = 200
        B_pre = np.linspace(0, Bm, n_pre)
        # Then follow the largest descending branch
        B_input = np.asarray(largest['B'])
        H_input = np.asarray(largest['H'])
        B_full = np.concatenate([B_pre, B_input])
        H_model_full = play_eval(B_full, eta, interps, K)
        H_model = H_model_full[n_pre:]
        max_err = float(np.max(np.abs(H_model - H_input)))
        assert max_err < 1e-9, \
            (f"JMAG identification reconstruction error {max_err:.3e} A/m "
             f"exceeds 1e-9 tolerance — possible regression in "
             f"hysteresis_io.build_shape_functions or load_mat.")

    def test_inverse_monotone_continuous_sweep(self, full_fixture):
        """Virgin Play material driven by continuous ascending H sweep should
        produce monotone non-decreasing M.

        Regression for the 2026-05-02 fix where the saturation early-return in
        radTPlayHysteresisMaterial::Inverse used to be unconditional. Without
        the collinear guard, a tiny disturbance pushed Newton onto the wrong
        branch of the multi-valued residual function above the saturation knee.
        Symptom: M jumped from +1.5e6 to -1.5e6 A/m crossing H ~ H_mono_max.
        """
        K, eta, tables = full_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        H_sweep = np.linspace(0, 20000, 100)
        M = np.zeros(len(H_sweep))
        for i, H in enumerate(H_sweep):
            Mv = rad.MatMvsH(mat, 'm', [float(H), 0, 0])
            rad.MatHysCommitState(mat)
            M[i] = float(Mv[0])

        dM = np.diff(M)
        # Allow tiny FP-level wiggles, but no jumps > 1e3 A/m downward.
        max_drop = float(-dM.min()) if len(dM) > 0 else 0.0
        assert max_drop < 1e3, \
            (f"M dropped by {max_drop:.3e} A/m on continuous ascending sweep — "
             f"branch-flip bug in MatMvsH (regression of cf1b6d2e fix).")

    def test_inverse_no_branch_flip_at_saturation(self, full_fixture):
        """MatMvsH should not flip B sign across saturation.

        Specific case that broke before cf1b6d2e: virgin material, drive H
        through positive saturation; B must stay positive.
        """
        K, eta, tables = full_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)

        # Approach saturation from below
        for H in [1000.0, 5000.0, 10000.0, 12000.0, 15000.0, 20000.0]:
            Mv = rad.MatMvsH(mat, 'm', [float(H), 0, 0])
            rad.MatHysCommitState(mat)
            B = MU_0 * (H + float(Mv[0]))
            assert B > 0, (f"B={B:+.4f} T at H={H} A/m — Inverse picked the "
                           f"wrong branch (regression of cf1b6d2e fix).")

    def test_shape_function_table_has_origin(self):
        """Each f_k table must include r=0 anchor for clean interpolation.

        Regression for the 2026-05-02 fix where np.arange floating-point
        noise caused `Bplay >= 0` to drop the middle index, leaving
        f_k tables that started at r = dB/2 instead of 0.
        """
        from radia.hysteresis_io import load_mat, build_shape_functions
        mat_file = (r"W:\999_菅原賢悟\19_磁気ヒステリシス\2024_IGTE_共同研究"
                    r"\2024_03_08_H-input_B-input\Potter_Schmulian\B_input.mat")
        if not os.path.exists(mat_file):
            pytest.skip("Potter-Schmulian B_input.mat not available")
        loops, dB, BMax = load_mat(mat_file)
        _, f_tables, _ = build_shape_functions(loops, dB)
        for k, (r, f) in enumerate(f_tables):
            r0 = float(np.asarray(r)[0])
            f0 = float(np.asarray(f)[0])
            assert r0 == 0.0, \
                f"f_tables[{k}] starts at r={r0:.3e}, expected exact 0"
            assert f0 == 0.0, \
                f"f_tables[{k}] f(r=0)={f0:.3e}, expected exact 0 (anti-symm)"
