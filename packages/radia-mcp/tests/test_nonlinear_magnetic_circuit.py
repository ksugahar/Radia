r"""Nonlinear magnetic-circuit operating point (#67) -- test.

NI = H_iron(B) * l_fe + (B/mu0) * gap, solved for B. Self-consistency (residual -> 0), the
constant-mu limit reproduces the linear reluctance formula exactly, a gap de-saturates the iron,
and the saturating BH curve H = 51 B + 2.5 B^15 reaches the stored ~1.56 T operating point at
NI = 32 A-turns over a ~15.6 mm path. Tool-independent (root-find against an analytic BH curve)."""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (
    MU0,
    magnetic_circuit_bh_operating_point,
    magnetic_circuit_bh_operating_summary,
)

H_of_B = lambda B: 51.0 * B + 2.5 * B ** 15            # the iron BH curve [A/m]
dH_dB = lambda B: 51.0 + 37.5 * B ** 14                # analytic derivative [A/m/T]
L_FE = 0.0156                                          # mean iron path [m]


def test_self_consistent_and_monotone():
    prev = -1.0
    for NI in (8, 16, 32, 64, 128, 256):
        B = magnetic_circuit_bh_operating_point(NI, L_FE, H_of_B)
        # Ampere's law residual ~ 0 (closed loop, gap=0)
        assert math.isclose(H_of_B(B) * L_FE, NI, rel_tol=1e-9, abs_tol=1e-6)
        assert B > prev                                # B increases monotonically with NI
        prev = B
    # hard saturation: 8x the MMF does NOT 8x B (the BH knee)
    B1 = magnetic_circuit_bh_operating_point(32, L_FE, H_of_B)
    B8 = magnetic_circuit_bh_operating_point(256, L_FE, H_of_B)
    assert B8 < 1.4 * B1                               # far sub-linear in NI


def test_constant_mu_limit_is_exact():
    mu_r = 1000.0
    lin = lambda B: B / (MU0 * mu_r)                   # H = B/(mu0 mu_r) -> linear reluctance
    for NI, gap in ((100.0, 0.0), (100.0, 1e-3), (50.0, 2e-3)):
        B = magnetic_circuit_bh_operating_point(NI, L_FE, lin, gap=gap)
        B_exact = NI / (L_FE / (MU0 * mu_r) + gap / MU0)
        assert math.isclose(B, B_exact, rel_tol=1e-9)


def test_gap_desaturates():
    NI = 200.0
    B_closed = magnetic_circuit_bh_operating_point(NI, L_FE, H_of_B, gap=0.0)
    B_gap = [magnetic_circuit_bh_operating_point(NI, L_FE, H_of_B, gap=g)
             for g in (0.2e-3, 0.5e-3, 1.0e-3)]
    assert B_closed > B_gap[0] > B_gap[1] > B_gap[2]   # more gap -> less flux
    # even a sub-mm gap drops the iron well below the closed-loop saturation
    assert B_gap[1] < 0.6 * B_closed


def test_stored_operating_point():
    # the saturating BH curve reaches the stored ~1.56 T at NI = 32 A-turns, l_fe ~ 15.6 mm
    B = magnetic_circuit_bh_operating_point(32.0, L_FE, H_of_B)
    assert math.isclose(B, 1.56, abs_tol=0.03)
    # and H there is the saturation-regime field
    assert 1900.0 < H_of_B(B) < 2200.0


def test_operating_summary_mmf_split_and_incremental_mu():
    rows = [
        magnetic_circuit_bh_operating_summary(NI, L_FE, H_of_B, dh_db=dH_dB)
        for NI in (8.0, 16.0, 32.0, 64.0, 128.0)
    ]
    assert max(abs(row["mmf_residual_A_turn"]) for row in rows) < 1.0e-9
    assert all(row["mmf_gap_A_turn"] == 0.0 for row in rows)
    assert all(row["incremental_mu_r"] < row["apparent_mu_r"] for row in rows)
    assert all(a["B_T"] < b["B_T"] for a, b in zip(rows, rows[1:]))
    assert all(a["incremental_mu_r"] > b["incremental_mu_r"] for a, b in zip(rows, rows[1:]))


def test_operating_summary_constant_mu_limit():
    mu_r = 750.0
    lin = lambda B: B / (MU0 * mu_r)
    dlin = lambda B: 1.0 / (MU0 * mu_r)
    row = magnetic_circuit_bh_operating_summary(80.0, L_FE, lin, gap=0.8e-3, dh_db=dlin)
    B_exact = 80.0 / (L_FE / (MU0 * mu_r) + 0.8e-3 / MU0)
    assert math.isclose(row["B_T"], B_exact, rel_tol=1e-9)
    assert math.isclose(row["apparent_mu_r"], mu_r, rel_tol=1e-12)
    assert math.isclose(row["incremental_mu_r"], mu_r, rel_tol=1e-12)
    assert math.isclose(
        row["secant_permeance_T_per_A_turn"],
        row["incremental_permeance_T_per_A_turn"],
        rel_tol=1e-12,
    )


def test_operating_summary_gap_drops_but_linearizes_iron():
    rows = [
        magnetic_circuit_bh_operating_summary(96.0, L_FE, H_of_B, gap=gap, dh_db=dH_dB)
        for gap in (0.0, 0.2e-3, 0.5e-3, 1.0e-3)
    ]
    assert all(a["B_T"] > b["B_T"] for a, b in zip(rows, rows[1:]))
    assert all(a["incremental_mu_r"] < b["incremental_mu_r"] for a, b in zip(rows, rows[1:]))
    gap_rows = rows[1:]
    assert all(a["gap_mmf_fraction"] < b["gap_mmf_fraction"] for a, b in zip(gap_rows, gap_rows[1:]))
    assert rows[-1]["gap_mmf_fraction"] > 0.9


if __name__ == "__main__":
    test_self_consistent_and_monotone()
    test_constant_mu_limit_is_exact()
    test_gap_desaturates()
    test_stored_operating_point()
    test_operating_summary_mmf_split_and_incremental_mu()
    test_operating_summary_constant_mu_limit()
    test_operating_summary_gap_drops_but_linearizes_iron()
    print("[OK] nonlinear magnetic circuit: NI=H(B)l_fe+(B/mu0)gap; self-consistent, const-mu exact, "
          "gap de-saturates, 1.56 T @ NI=32.")
