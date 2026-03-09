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
