import copy
import json
import math

import pytest

from radia_mcp.radia_ngsolve.lossy_dielectric_power_gate import (
    lossy_dielectric_complex_power_refinement_gate,
)
from radia_mcp.radia_ngsolve.server import (
    lossy_dielectric_complex_power_refinement_gate as mcp_gate,
)


def _encoded(value: complex | float) -> dict[str, float]:
    number = complex(value)
    return {"real": number.real, "imag": number.imag, "abs": abs(number)}


def good_summary() -> dict:
    frequency = 60.0
    epsilon_0 = 8.8541878128e-12
    epsilon_r = 6.0
    sigma = 1.0e-8
    omega = 2.0 * math.pi * frequency
    loss_ratio = sigma / (omega * epsilon_0 * epsilon_r)
    rows = []
    for index, (mesh_size, elements, p) in enumerate(
        ((0.1, 6500, 8.90e-8), (0.05, 23000, 8.88e-8), (0.025, 88000, 8.864e-8))
    ):
        q = p / loss_ratio
        w = q / (2.0 * omega)
        s = complex(p * (1.0 - 2.0e-4), q * (1.0 + 5.0e-5))
        rows.append(
            {
                "mesh_size_in": mesh_size,
                "element_count": elements,
                "real_power_W": _encoded(p),
                "reactive_power_var": _encoded(complex(q, q * 1.0e-5)),
                "apparent_power_VA": _encoded(s),
                "time_average_stored_energy_J": _encoded(w),
                "voltage_drop_V": _encoded(10.0 + index * 1.0e-10),
            }
        )
    return {
        "frequency_Hz": frequency,
        "sigma_S_per_m": sigma,
        "epsilon_r": epsilon_r,
        "epsilon_0_F_per_m": epsilon_0,
        "rows": rows,
    }


def test_accepts_complex_power_and_three_level_refinement_closure():
    result = lossy_dielectric_complex_power_refinement_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["complex_power_imaginary_part_matches_reactive_power"] is True
    assert result["metrics"]["last_pair_real_power_relative_change"] < 0.002
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_rejects_treating_complex_power_real_part_as_magnitude():
    bad = copy.deepcopy(good_summary())
    for row in bad["rows"]:
        p = row["real_power_W"]["real"]
        q = row["reactive_power_var"]["real"]
        wrong = math.hypot(p, q)
        row["apparent_power_VA"] = _encoded(complex(wrong, 0.0))
    result = lossy_dielectric_complex_power_refinement_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["complex_power_real_part_matches_real_power"] is False
    assert result["checks"]["complex_power_imaginary_part_matches_reactive_power"] is False


def test_rejects_nonfinite_or_internally_inconsistent_complex_encoding():
    bad = good_summary()
    bad["rows"][0]["apparent_power_VA"]["abs"] *= 2.0
    with pytest.raises(ValueError, match="inconsistent"):
        lossy_dielectric_complex_power_refinement_gate(bad)

    nonfinite = good_summary()
    nonfinite["rows"][1]["real_power_W"]["real"] = math.nan
    with pytest.raises(ValueError, match="finite"):
        lossy_dielectric_complex_power_refinement_gate(nonfinite)
