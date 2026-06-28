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
                            '../validation_test/hysteresis/binput_play_fixture.npz')
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
                            '../validation_test/hysteresis/binput_play_fixture.npz')
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
        # State = K*9 (Play states) + 7 (warm-start cache, fixed 2026-05-02)
        assert len(state) == K * 9 + 7, \
            f"State size should be K*9+7={K*9+7}, got {len(state)}"

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
        mag = rad.magnet_box([0, 0, 0], [0.02, 0.02, 0.02], [0, 0, 0])
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
        # K*9 (per-particle J states) + 7 (warm-start cache, fixed 2026-05-02)
        assert len(state) == K * 9 + 7

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

    def test_type5_isotropy_3d(self):
        """Type 5 (Energy) must be 3D isotropic to machine precision.

        Regression for the 2026-05-02 fixes:
          - 00a13857: L_inf -> L2 norm in U_k argument
          - 2b6db1cf: auto-scale eps to (max_chi * 1e-4)^2 floor

        Synthetic isotropic shape functions f_k(r) = a_k * r driven by
        the SAME H amplitude trajectory along (a) +x and (b) (1,2,3)/||.||
        must produce B related by exact rotation. Pre-fix: 37% relative
        error (broken). Post-fix: < 1e-6 relative (machine-precision).
        """
        K = 10
        a_k = list(np.geomspace(1e3, 5e4, K))
        chi_k = list(np.geomspace(20.0, 800.0, K))
        r_max = 0.20
        r_grid = list(np.linspace(0, r_max, 200))
        tables = [(r_grid, [a * r for r in r_grid]) for a in a_k]

        n = 30
        Hs = 800.0 * np.sin(np.linspace(0, 4 * np.pi, n))

        def drive(direction):
            rad.UtiDelAll()
            mat = rad.MatEnergyHysteresis(K, chi_k, tables, 1e-8)
            B = np.zeros((n, 3))
            for i, h in enumerate(Hs):
                Hv = h * np.asarray(direction)
                Mv = rad.MatMvsH(mat, 'm', list(Hv))
                rad.MatHysCommitState(mat)
                B[i] = MU_0 * (Hv + np.array(Mv))
            return B

        B_x = drive([1.0, 0, 0])

        # Test multiple out-of-plane directions
        for axis in [(1, 1, 0), (1, 2, 3), (3, 1, 7), (1, 1, 1)]:
            e = np.array(axis, dtype=float); e /= np.linalg.norm(e)
            B_e = drive(e)

            # Build rotation matrix sending +x to e (Rodrigues' formula)
            a = np.array([1.0, 0, 0])
            v = np.cross(a, e)
            sn = np.linalg.norm(v); cs = a @ e
            if sn > 1e-30:
                nx = v / sn
                Vx = np.array([[0, -nx[2], nx[1]],
                                [nx[2], 0, -nx[0]],
                                [-nx[1], nx[0], 0]])
                R = np.eye(3) + sn * Vx + (1 - cs) * Vx @ Vx
            else:
                R = np.eye(3)

            B_x_rot = (R @ B_x.T).T
            err = float(np.max(np.linalg.norm(B_e - B_x_rot, axis=1)))
            scale = float(np.max(np.linalg.norm(B_x, axis=1)))
            rel = err / max(scale, 1e-30)
            assert rel < 1e-6, \
                (f"Type 5 isotropy failed for axis {axis}: rel error "
                 f"{rel:.3e} (was 1e-2 pre-fix, must be < 1e-6 post-fix).")

    def test_state_warmstart_invariance_play(self, full_fixture):
        """RestoreStateFromArray must be state-equivalent to having driven
        the material to that state (not just to the Play states alone).

        Regression for the 2026-05-02 fix that extended GetStateSize from
        K*9 to K*9+7 to include the warm-start cache (m_last_B, m_last_H,
        m_has_result). Without that cache, a Save → drive-elsewhere →
        Restore → probe sequence cold-started the next Inverse, which
        landed on a different Newton basin on hard branches and produced a
        2.83e6 A/m sign flip in M.
        """
        K, eta, tables = full_fixture

        # Path A: in-line drive
        rad.UtiDelAll()
        mat_A = rad.MatPlayHysteresis(K, eta, tables)
        rad.MatMvsH(mat_A, 'm', [2000.0, 0, 0]); rad.MatHysCommitState(mat_A)
        M_A = np.array(rad.MatMvsH(mat_A, 'm', [2000.0, 500.0, -300.0]))

        # Path B: drive elsewhere, restore, probe
        rad.UtiDelAll()
        mat_B = rad.MatPlayHysteresis(K, eta, tables)
        rad.MatMvsH(mat_B, 'm', [2000.0, 0, 0]); rad.MatHysCommitState(mat_B)
        state = rad.MatHysSaveState(mat_B)
        rad.MatMvsH(mat_B, 'm', [-30000.0, 0, 0]); rad.MatHysCommitState(mat_B)
        rad.MatMvsH(mat_B, 'm', [+30000.0, 30000.0, 0]); rad.MatHysCommitState(mat_B)
        rad.MatHysRestoreState(mat_B, state)
        M_B = np.array(rad.MatMvsH(mat_B, 'm', [2000.0, 500.0, -300.0]))

        err = float(np.linalg.norm(M_A - M_B))
        assert err < 1.0, (f"||M_A - M_B|| = {err:.3e} A/m — RestoreState is "
                            f"not state-equivalent to in-line drive (regression "
                            f"of the K*9+7 state-size fix).")

    def test_roundtrip_BHB_play_hysteresis(self, full_fixture):
        """B -> H -> B round-trip via the C++ Type 5 (Play) material.

        The IGTE'26 digest claims B -> H -> B accuracy of 1.4e-10 T from the
        MATLAB BInputEnergyModel reference. This test verifies the same
        property holds for the production C++ Play implementation in Radia
        (the Play forward is algebraically identical to the Energy forward
        by the rev/irrev separation, so accuracy must transfer).

        Forward B -> H:  H = nu_rev * B + MatHysIrreversible(B)
        Inverse H -> B:  M = MatMvsH(H);  B' = mu_0 * (H + M)

        State is saved before the forward pass and restored before the
        inverse, so both directions evaluate from the same play history.
        """
        K, eta, tables = full_fixture
        rad.UtiDelAll()
        mat = rad.MatPlayHysteresis(K, eta, tables)
        nu_rev = rad.MatHysGetNuRev(mat)

        # Sweep B over a range that includes both linear and saturation
        Bm = 1.5
        B_seq = np.linspace(0.1, Bm, 20).tolist() + [Bm, 0.5 * Bm, -0.5 * Bm]

        max_err = 0.0
        for B in B_seq:
            B_vec = [float(B), 0.0, 0.0]

            # Save play state before any evaluation
            state = rad.MatHysSaveState(mat)

            # Forward B -> H
            H_irr = rad.MatHysIrreversible(mat, B_vec)
            H = [nu_rev * B_vec[i] + H_irr[i] for i in range(3)]

            # Restore state so the inverse starts from the same history
            rad.MatHysRestoreState(mat, state)

            # Inverse H -> M -> B'
            M = rad.MatMvsH(mat, 'm', H)
            B_rec = [MU_0 * (H[i] + M[i]) for i in range(3)]

            err = abs(B - B_rec[0])
            max_err = max(max_err, err)

            # Commit and continue to the next B
            rad.MatHysCommitState(mat)

        # MATLAB reference achieves ~1.4e-10 T on this fixture; allow
        # a bit of headroom for the C++ Newton tolerance settings.
        assert max_err < 1e-8, \
            f"B->H->B round-trip max error {max_err:.3e} T exceeds 1e-8 T " \
            f"(MATLAB reference: 1.4e-10 T). Possible regression in " \
            f"MatMvsH/MatHysIrreversible coupling for Type 6 material."

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


class TestEnergyBasedPlayModelPython:
    """Regression tests for the Python EnergyBasedPlayModel reference.

    Covers the 4 bugs fixed in commit 6909af1c (2026-05-02) in
    src/radia/energy_play_model.py — the documented Python reference of the
    rev/irrev separation form of B-input Play.
    """

    def test_compute_nu_rev_matches_virgin_curve_scan(self, full_fixture):
        """Bug 1: nu_rev must equal max_B sum_k f_k'(p_k(B)) on the virgin curve.

        The pre-fix implementation returned `sum_k max|f_k'|`, a loose
        triangle-inequality upper bound (typically 2-3x the actual maximum).
        Symptom: Picard contraction ratio close to 1 -> slow inverse
        convergence. Verify nu_rev matches an independent virgin-curve scan.
        """
        from radia.energy_play_model import EnergyBasedPlayModel
        K, eta, f_k_tables = full_fixture
        f_k_tables = [(np.asarray(r), np.asarray(f)) for r, f in f_k_tables]
        model = EnergyBasedPlayModel(eta, f_k_tables)

        # Independent virgin-curve scan of max_B sum_k f_k'(max(B - eta_k, 0))
        eta_arr = np.asarray(eta)
        B_sat = max(float(np.max(np.abs(np.asarray(rk))))
                    for rk, _ in f_k_tables)
        df_interps = []
        for rk, fk in f_k_tables:
            r = np.asarray(rk); f = np.asarray(fk)
            mask = r >= 0
            r = r[mask]; f = f[mask]
            idx = np.argsort(r); r = r[idx]; f = f[idx]
            df_interps.append((r, np.gradient(f, r)))
        n_scan = 1000
        max_dHdB = 0.0
        for i in range(1, n_scan + 1):
            B_val = B_sat * i / n_scan
            dHdB = 0.0
            for k in range(K):
                pk = B_val - eta_arr[k]
                if pk <= 1e-30:
                    continue
                r, df = df_interps[k]
                dHdB += float(np.interp(pk, r, df))
            max_dHdB = max(max_dHdB, dHdB)

        rel_err = abs(model.nu_rev - max_dHdB) / max_dHdB
        assert rel_err < 0.05, \
            f"nu_rev = {model.nu_rev:.3e} should match virgin-curve max " \
            f"= {max_dHdB:.3e} within 5% (got {rel_err*100:.1f}%); " \
            f"loose triangle-inequality bound (old bug) is typically 2-3x larger"

    def test_play_operator_three_regimes(self, full_fixture):
        """Bug 2: _play_operator must implement the Play update, not np.clip.

        Standard B-input Play: p_new = max(B - eta_k, min(B + eta_k, p_old))
          - elastic regime |B - p_old| <= eta:  p_new = p_old
          - following up    B > p_old + eta:    p_new = B - eta
          - following down  B < p_old - eta:    p_new = B + eta

        Pre-fix used np.clip(B, p-eta, p+eta), which clamps B to that range
        and returns the wrong scalar (e.g., elastic regime returned B
        instead of p_old).
        """
        from radia.energy_play_model import EnergyBasedPlayModel
        K, eta, f_k_tables = full_fixture
        f_k_tables = [(np.asarray(r), np.asarray(f)) for r, f in f_k_tables]
        model = EnergyBasedPlayModel(eta, f_k_tables)

        # Pick a hysteron with non-trivial threshold
        k_test = next(k for k in range(K) if eta[k] > 1e-3)
        eta_k = float(eta[k_test])
        p_old = 0.5

        # Elastic regime: |B - p_old| <= eta_k -> p_new = p_old
        model._p[k_test] = p_old
        B_elastic = p_old + 0.5 * eta_k
        p_new = model._play_operator(B_elastic, k_test)
        assert abs(p_new - p_old) < 1e-12, \
            f"Elastic regime: p_new should equal p_old={p_old}, got {p_new}"

        # Following-up regime: B > p_old + eta_k -> p_new = B - eta_k
        model._p[k_test] = p_old
        B_above = p_old + 2.0 * eta_k
        p_new = model._play_operator(B_above, k_test)
        assert abs(p_new - (B_above - eta_k)) < 1e-12, \
            f"Following-up: p_new should equal B-eta_k={B_above-eta_k}, " \
            f"got {p_new}"

        # Following-down regime: B < p_old - eta_k -> p_new = B + eta_k
        model._p[k_test] = p_old
        B_below = p_old - 2.0 * eta_k
        p_new = model._play_operator(B_below, k_test)
        assert abs(p_new - (B_below + eta_k)) < 1e-12, \
            f"Following-down: p_new should equal B+eta_k={B_below+eta_k}, " \
            f"got {p_new}"

    def test_irreversible_anti_symmetric_at_negative_B(self, full_fixture):
        """Bug 3: H_irr(-Bm) must be NEGATIVE, not positive.

        Pre-fix evaluated g_k_interp(pk) with signed pk, but tables span
        only r >= 0; linear extrapolation gave wrong sign for large negative
        |pk|. Symptom: H(-Bm) came out POSITIVE (sign wrong).
        Fix: evaluate at |pk|, then multiply by sign(pk) (anti-symmetry).
        """
        from radia.energy_play_model import EnergyBasedPlayModel
        K, eta, f_k_tables = full_fixture
        f_k_tables = [(np.asarray(r), np.asarray(f)) for r, f in f_k_tables]
        model = EnergyBasedPlayModel(eta, f_k_tables)

        Bm = 1.5
        N = 200
        B_seq = np.concatenate([
            np.linspace(0.0, Bm, N),
            np.linspace(Bm, -Bm, N),
        ])
        H_seq = np.array([model.forward(B) for B in B_seq])

        # H at the final point (B = -Bm) must be negative
        assert H_seq[-1] < 0, \
            f"H(B=-Bm) = {H_seq[-1]:.3e} should be NEGATIVE, " \
            f"old bug returned positive value (sign error in g_k extrapolation)"

        # Sweep range should be roughly symmetric (|min| ~ |max|)
        H_max, H_min = float(H_seq.max()), float(H_seq.min())
        asymmetry = abs(H_max + H_min) / max(H_max, abs(H_min))
        assert asymmetry < 0.1, \
            f"H range [{H_min:.1f}, {H_max:.1f}] asymmetric: " \
            f"{asymmetry*100:.1f}% > 10% threshold (bug 3 regression)"

    def test_energy_nonzero_on_negative_half_plane(self, full_fixture):
        """Bug 4: W_irr(B) for negative B must be > 0, not zero.

        Pre-fix used `mask = rk <= pk` with rk >= 0 (g_k tables span r >= 0)
        and pk < 0; the mask was always empty, so W_irr stayed 0 on the
        entire negative half-plane. Fix: integrate over |pk|, exploiting
        evenness of G_k by anti-symmetry of g_k.
        """
        from radia.energy_play_model import EnergyBasedPlayModel
        K, eta, f_k_tables = full_fixture
        f_k_tables = [(np.asarray(r), np.asarray(f)) for r, f in f_k_tables]
        model = EnergyBasedPlayModel(eta, f_k_tables)

        Bm = 1.0
        N = 100
        B_seq = np.concatenate([
            np.linspace(0.0, Bm, N),
            np.linspace(Bm, -Bm, N),
        ])
        for B in B_seq:
            model.forward(B)  # advance Play state through the trajectory

        W_neg = model.energy(-Bm)
        W_rev = 0.5 * model.nu_rev * Bm**2
        W_irr_neg = W_neg - W_rev

        # The pre-fix bug returned W_irr = 0 EXACTLY on the entire negative
        # half-plane (empty integration mask). Post-fix: W_irr is non-zero
        # because integration is over |pk|. Sign of W_irr is not the issue
        # (G_k can be negative by design); the issue was identical-zero.
        assert abs(W_irr_neg) > 1e-3 * abs(W_rev), \
            f"W_irr(-Bm) = {W_irr_neg:.3e} should be non-trivial " \
            f"(|W_irr| >> 0); old bug had empty integration mask on " \
            f"negative B and returned exactly 0"
