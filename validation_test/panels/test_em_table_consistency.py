"""1% fixed-T consistency gate for the 2D EM table.

For a small 5x5 table at constant sigma, asserts that
``em_table.interp_qsurf(H, T_grid_node)`` reproduces a direct
ESIMCellProblemSolver.solve(H).Z value (via q_surf = 0.5*Re(Z)*H^2)
to within 1%.  This proves the table is the same physical quantity
the existing Karl outer loop in calc_fem_kelvin converges to, just
looked up instead of solved.

Also exercises the bilinear interpolation contract (grid-node hit,
mid-cell hit, out-of-range clamp).

No Cubit, no NGSolve -- runs in ~5 seconds.
"""

from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

_panels_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia", "panels"))
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)

from em_table import load_em_table, interp_qsurf, interp_Zs  # noqa: E402


# Build one table per session so the ESIM solves are amortised across
# the assertion tests.
@pytest.fixture(scope="module")
def small_table(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("em_table")
    npz_path = out_dir / "em_table_steel_5x5.npz"
    cmd = [
        sys.executable,
        os.path.join(_panels_dir, "calc_em_table.py"),
        "--material", "steel",
        "--frequency", "50e3",
        "--H-grid-min", "1e3",
        "--H-grid-max", "1e5",
        "--n-H", "5",
        "--T-grid-min", "20",
        "--T-grid-max", "600",
        "--n-T", "5",
        "--em-table-output", str(npz_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, (
        f"calc_em_table.py failed: stdout={proc.stdout}\n"
        f"stderr={proc.stderr}")
    assert npz_path.exists(), f"Table not written to {npz_path}"
    return load_em_table(str(npz_path))


class TestTableShape:
    def test_grid_shapes(self, small_table):
        assert small_table.H_grid.shape == (5,)
        assert small_table.T_grid.shape == (5,)
        assert small_table.Zs_re.shape == (5, 5)
        assert small_table.Zs_im.shape == (5, 5)
        assert small_table.q_surf.shape == (5, 5)

    def test_grid_monotone(self, small_table):
        assert np.all(np.diff(small_table.H_grid) > 0)
        assert np.all(np.diff(small_table.T_grid) > 0)

    def test_q_surf_consistent_with_Zs(self, small_table):
        # By construction in calc_em_table: q_surf = 0.5 * Re(Z_s) * H^2
        expected = 0.5 * small_table.Zs_re * (small_table.H_grid[:, None] ** 2)
        np.testing.assert_allclose(small_table.q_surf, expected, rtol=1e-12)

    def test_meta_carries_material_and_frequency(self, small_table):
        assert small_table.meta["material"] == "steel"
        assert small_table.meta["frequency"] == 50e3
        # No --sigma-T-curve given -> sigma_T_path empty.
        assert small_table.meta["sigma_T_curve"] == ""


class TestFixedT_OnePercentGate:
    """The headline check: table lookup vs single direct ESIM solve."""

    @pytest.fixture(scope="class")
    def direct_solver(self):
        # Same construction calc_em_table uses internally.
        sys.path.insert(0, os.path.join(_panels_dir, ".."))
        from esim_cell_problem import ESIMCellProblemSolver
        from em_material import EMMaterial
        mat = EMMaterial.from_name("steel")
        return ESIMCellProblemSolver(
            bh_curve=mat.bh_curve, sigma=mat.sigma,
            frequency=50e3, n_nodes=100, num_skin_depths=10)

    @pytest.mark.parametrize("H_t", [3e3, 1e4, 3e4, 1e5])
    def test_grid_node_T20_qsurf_within_1pct(self, small_table, direct_solver, H_t):
        # H may not hit a grid node, but the constant-sigma table has
        # identical columns and bilinear interp is exact along H, so
        # the only error is the H-direction interp at the queried H.
        # We restrict the H values to the grid nodes themselves to
        # isolate the "lookup vs solve" gate from interpolation error.
        H_grid = small_table.H_grid
        # Snap H to the nearest grid node to make this an exact lookup.
        H_node = float(H_grid[int(np.argmin(np.abs(H_grid - H_t)))])
        T_node = float(small_table.T_grid[0])  # 20 degC

        q_table = float(interp_qsurf(small_table, H_node, T_node))

        Z_direct = complex(direct_solver.solve(H_node)["Z"])
        q_direct = 0.5 * Z_direct.real * H_node ** 2

        rel_err = abs(q_table - q_direct) / abs(q_direct)
        assert rel_err < 0.01, (
            f"Table q_surf={q_table:.4e}, direct={q_direct:.4e}, "
            f"rel_err={rel_err:.2%} > 1% at H={H_node:.2e} A/m, "
            f"T={T_node:.1f} C")

    @pytest.mark.parametrize("H_t", [3e3, 1e4, 3e4, 1e5])
    def test_grid_node_T20_Zs_within_1pct(self, small_table, direct_solver, H_t):
        H_grid = small_table.H_grid
        H_node = float(H_grid[int(np.argmin(np.abs(H_grid - H_t)))])
        T_node = float(small_table.T_grid[0])

        Z_table = complex(interp_Zs(small_table, H_node, T_node))
        Z_direct = complex(direct_solver.solve(H_node)["Z"])

        rel_err_re = abs(Z_table.real - Z_direct.real) / abs(Z_direct.real)
        rel_err_im = abs(Z_table.imag - Z_direct.imag) / abs(Z_direct.imag)
        assert rel_err_re < 0.01 and rel_err_im < 0.01, (
            f"Table Z_s={Z_table}, direct={Z_direct}, "
            f"rel_err Re={rel_err_re:.2%} Im={rel_err_im:.2%}")


class TestInterpolation:
    def test_grid_node_exact(self, small_table):
        # Querying a grid node must return that exact cell.
        i, j = 2, 3
        H = float(small_table.H_grid[i])
        T = float(small_table.T_grid[j])
        q = float(interp_qsurf(small_table, H, T))
        np.testing.assert_allclose(q, small_table.q_surf[i, j], rtol=1e-12)

    def test_mid_cell_bounded(self, small_table):
        # Mid-cell q_surf must lie within the 4 cell corner values
        # (monotone bilinear interp).
        i, j = 1, 1
        H = float(np.exp(0.5 * (np.log(small_table.H_grid[i])
                                + np.log(small_table.H_grid[i + 1]))))
        T = float(0.5 * (small_table.T_grid[j] + small_table.T_grid[j + 1]))
        q = float(interp_qsurf(small_table, H, T))
        corners = [small_table.q_surf[i, j], small_table.q_surf[i + 1, j],
                   small_table.q_surf[i, j + 1], small_table.q_surf[i + 1, j + 1]]
        assert min(corners) <= q <= max(corners)

    def test_out_of_range_clamps(self, small_table):
        # Below H_min and above H_max should clamp to grid edges, not raise.
        T = float(small_table.T_grid[0])
        H_min, H_max = float(small_table.H_grid[0]), float(small_table.H_grid[-1])
        q_below = float(interp_qsurf(small_table, 0.1 * H_min, T))
        q_at_min = float(interp_qsurf(small_table, H_min, T))
        np.testing.assert_allclose(q_below, q_at_min, rtol=1e-12)

        q_above = float(interp_qsurf(small_table, 10.0 * H_max, T))
        q_at_max = float(interp_qsurf(small_table, H_max, T))
        np.testing.assert_allclose(q_above, q_at_max, rtol=1e-12)

    def test_vector_input(self, small_table):
        H_arr = np.array([1e3, 1e4, 1e5])
        T_arr = np.array([20.0, 100.0, 600.0])
        q = interp_qsurf(small_table, H_arr, T_arr)
        assert q.shape == (3,)
        # Each entry must match the scalar call.
        for k in range(3):
            qk = float(interp_qsurf(small_table, float(H_arr[k]), float(T_arr[k])))
            np.testing.assert_allclose(q[k], qk, rtol=1e-12)
