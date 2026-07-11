import cmath
import json
import math

from radia_mcp.radia_ngsolve.acoustic_kernel_gate import helmholtz_double_layer_low_frequency_gate
from radia_mcp.radia_ngsolve.server import helmholtz_double_layer_low_frequency_gate as mcp_gate


def _summary():
    distance = 0.8
    normal_dot = 0.6
    laplace = normal_dot / (4.0 * math.pi * distance**3)
    rows = []
    for k in [0.0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]:
        z = 1j * k * distance
        stable_factor = sum(
            (1 - order) * z**order / math.factorial(order)
            for order in range(2, 9)
        )
        correction = laplace * stable_factor
        direct = laplace * cmath.exp(z) * (1.0 - z)
        series = laplace * (-z**2 / 2.0 - z**3 / 3.0)
        rows.append({
            "wavenumber_per_m": k,
            "kr_abs": abs(z),
            "laplace_real": laplace,
            "laplace_imag": 0.0,
            "correction_real": correction.real,
            "correction_imag": correction.imag,
            "split_direct_relative_error": abs(laplace + correction - direct) / abs(direct),
            "correction_series_relative_error": (
                0.0 if correction == 0.0 else abs(correction - series) / abs(correction)
            ),
            "cancellation_ratio": None if correction == 0.0 else abs(laplace) / abs(correction),
        })
    return {
        "kernel_family": "helmholtz_source_normal_double_layer",
        "time_convention": "exp(+i*k*r)",
        "distance_m": distance,
        "normal_dot_m": normal_dot,
        "rows": rows,
    }


def test_double_layer_low_frequency_gate_accepts_quadratic_regular_part():
    result = helmholtz_double_layer_low_frequency_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["correction_starts_quadratic"] is True


def test_double_layer_low_frequency_gate_rejects_normal_sign_and_stale_correction():
    summary = _summary()
    summary["normal_dot_m"] *= -1.0
    summary["rows"][2]["correction_real"] *= 0.5
    result = helmholtz_double_layer_low_frequency_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["laplace_part_matches_geometry"] is False
    assert result["checks"]["correction_starts_quadratic"] is False


def test_double_layer_low_frequency_mcp_gate_rejects_nonfinite_rows():
    summary = _summary()
    summary["rows"][1]["kr_abs"] = "nan"
    result = json.loads(mcp_gate(json.dumps(summary)))
    assert result["status"] == "invalid_input"
