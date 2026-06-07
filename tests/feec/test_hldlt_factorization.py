"""Golden test for the symmetric H-LDL^T factorization (cHACApK_harith, Phase 6).

H-LDL^T is the symmetric H-matrix factorization the HDiv-type VIM operator (symmetric A = (1/chi)
M_mass - N) will use once it is assembled as an H-matrix -- it stores only the lower triangle
(~half the storage of the non-symmetric H-LU) and uses Bunch-Kaufman 2x2 pivoting for the
indefinite case.  This locks the standalone factorization (validated in a worktree, then merged):
LDL^T solves a symmetric-indefinite block tree accurately, exercises 2x2 pivots, and stores ~half
of what H-LU does.  (Wiring it to the HDiv operator needs the HDiv-operator-as-H-matrix, a later
production step; this test guards the factorization itself.)
"""
import pytest
import radia._radia_pybind as _rp


@pytest.mark.parametrize("depth,nb", [
    (1, 8),
    (2, 16),
    (3, 24),
    (4, 8),
])
def test_hldlt_symmetric_indefinite(depth, nb):
    r = _rp._hldlt_self_test(depth, nb)
    # (1) LDL^T solves the symmetric-indefinite system accurately (vs LAPACKE_dsysv reference)
    assert 0.0 <= r["ldlt_residual"] < 1e-9, f"LDL^T residual {r['ldlt_residual']:.2e} too large"
    # (2) the indefinite Bunch-Kaufman path is genuinely exercised (2x2 pivots used)
    assert r["n_2x2_pivots"] > 0, "no 2x2 pivots -> indefinite path not exercised"
    # (3) storage: H-LU stores the lower+upper off-diagonal = exactly 2x the LDL^T lower-only
    assert r["hlu_offdiag_doubles"] == 2 * r["ldlt_lower_doubles"], (
        f"storage ratio != 2: H-LU={r['hlu_offdiag_doubles']} vs 2*LDLT={2*r['ldlt_lower_doubles']}")


def test_hldlt_more_accurate_than_nopivot_hlu_on_indefinite():
    """On a symmetric INDEFINITE matrix, the pivoted LDL^T is at least as accurate as the no-pivot
    H-LU (a robustness bonus of the symmetric factorization)."""
    r = _rp._hldlt_self_test(3, 16)
    assert r["ldlt_residual"] <= r["hlu_residual"] + 1e-12, (
        f"LDL^T {r['ldlt_residual']:.2e} worse than H-LU {r['hlu_residual']:.2e}")


# --- rk-aware H-LDL^T (low-rank / ACA-compressed off-diagonal leaves) ---------
# A real ACA-compressed H-matrix stores the admissible far-field blocks as rk
# (U V^T) leaves.  These tests lock the rk-aware factorization: the off-diagonal
# leaves are built EXACTLY rank-rk_rank so the rk arithmetic (offdiag-solve +
# the 4 trailing-update combos) is EXACT (no truncation), giving a near-machine
# residual when the rk path is correct.

@pytest.mark.parametrize("nb,rk", [
    (24, 4),
    (32, 8),
])
def test_hldlt_rk_offdiagonal_depth1(nb, rk):
    """depth-1: dense indefinite diagonal leaves + an EXACT-rank rk off-diagonal leaf."""
    r = _rp._hldlt_self_test_rk(nb, rk, 1)
    assert 0.0 <= r["ldlt_residual"] < 1e-8, (
        f"rk LDL^T residual {r['ldlt_residual']:.2e} too large (rk path bug)")


@pytest.mark.parametrize("nb,rk", [
    (24, 4),
    (32, 8),
])
def test_hldlt_rk_offdiagonal_depth2(nb, rk):
    """depth-2: rk off-diagonal leaves at TWO levels (rk*rk trailing update into an
    internal sub-tree)."""
    r = _rp._hldlt_self_test_rk(nb, rk, 2)
    assert 0.0 <= r["ldlt_residual"] < 1e-8, (
        f"rk depth-2 LDL^T residual {r['ldlt_residual']:.2e} too large (rk path bug)")


@pytest.mark.parametrize("s", [4, 5])
def test_hldlt_rk_mixed_all_four_combos(s):
    """Mixed dense/rk off-diagonal leaves: the trailing update must exercise ALL FOUR
    operand-kind combos (dense*dense, rk*dense, dense*rk, rk*rk) AND stay machine-exact."""
    r = _rp._hldlt_self_test_rk_mixed(20, 4, s)
    assert 0.0 <= r["ldlt_residual"] < 1e-8, (
        f"mixed rk LDL^T residual {r['ldlt_residual']:.2e} too large")
    # every combo path must actually fire (otherwise the test gives false coverage)
    assert r["combo_dense_dense"] > 0, "combo dense*dense never fired"
    assert r["combo_rk_dense"]    > 0, "combo rk*dense never fired"
    assert r["combo_dense_rk"]    > 0, "combo dense*rk never fired"
    assert r["combo_rk_rk"]       > 0, "combo rk*rk never fired"
