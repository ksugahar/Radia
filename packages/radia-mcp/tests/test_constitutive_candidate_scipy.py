"""Optional independent SciPy oracle; keep minimal gate tests collectable."""
import math
import pytest
from radia_mcp.radia_ngsolve.field_profile_gate import build_constitutive_comparison_candidate


def test_build_comparison_candidate_matches_pchip_and_energy_identity():
    PchipInterpolator = pytest.importorskip("scipy.interpolate").PchipInterpolator

    table = [[0.0, 0.0], [100.0, 0.5], [1000.0, 1.4], [10000.0, 1.7]]
    h_values = [0.0, 50.0, 500.0, 10000.0, 20000.0]
    candidate = build_constitutive_comparison_candidate(table, h_values)
    # Radia's law: monotone PCHIP of the rising magnetization M = B/mu0 - H.
    mu_0 = 4.0e-7 * math.pi
    magnetization = PchipInterpolator(
        [row[0] for row in table], [row[1] / mu_0 - row[0] for row in table],
        extrapolate=False,
    )
    expected_inside = [
        mu_0 * (value + float(magnetization(value))) for value in h_values[:-1]
    ]
    assert candidate["identity"]["constitutive_interpolation"] == "monotone_pchip"
    assert candidate["identity"]["solver_runtime_mode"] is True
    assert candidate["B_T"][:-1] == pytest.approx(expected_inside, rel=1.0e-12)
    assert candidate["B_T"][-1] == pytest.approx(
        table[-1][1] + 4.0e-7 * math.pi * (h_values[-1] - table[-1][0])
    )
    for h_value, b_value, energy, coenergy in zip(
        candidate["H_A_per_m"],
        candidate["B_T"],
        candidate["energy_density_J_per_m3"],
        candidate["coenergy_density_J_per_m3"],
    ):
        assert energy + coenergy == pytest.approx(h_value * b_value, abs=1.0e-10)
