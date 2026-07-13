import copy
import json
import math

from radia_mcp.radia_ngsolve.heterogeneous_current_flow_gate import (
    heterogeneous_current_flow_p1_reintegration_gate,
)
from radia_mcp.radia_ngsolve.server import (
    heterogeneous_current_flow_p1_reintegration_gate as mcp_gate,
)


POINT_ORDER = [
    "V",
    "Jx",
    "Jy",
    "Kx",
    "Ky",
    "Ex",
    "Ey",
    "eps_x_relative",
    "eps_y_relative",
    "Jdx",
    "Jdy",
    "sigma_x_S_per_m",
    "sigma_y_S_per_m",
    "Jcx",
    "Jcy",
]
ODD_INDICES = {0, 1, 2, 5, 6, 9, 10, 13, 14}


def _encoded(value: complex | float) -> dict[str, float]:
    number = complex(value)
    return {"real": number.real, "imag": number.imag, "abs": abs(number)}


def _point_values(sign: float) -> list[object]:
    values: list[object] = []
    for index in range(len(POINT_ORDER)):
        if index in ODD_INDICES:
            values.append(_encoded(sign * complex(index + 1, 0.25 * index)))
        elif index in {7, 8, 11, 12}:
            values.append(float(index + 1))
        else:
            values.append(_encoded(complex(index + 1, 0.1)))
    return values


def _case(
    name: str,
    mesh: float,
    elements: int,
    p: float,
    q: float,
    sign: float = 1.0,
) -> dict:
    omega = 2.0 * math.pi * 50.0
    energy = q / (2.0 * omega)
    return {
        "case": name,
        "mesh_size_mm": mesh,
        "high_voltage_V": sign * 10.0,
        "element_count": elements,
        "real_power_W": _encoded(p),
        "apparent_power_VA": _encoded(complex(p, q)),
        "time_average_stored_energy_J": _encoded(energy),
        "area_m2": _encoded(1.0),
        "volume_m3": _encoded(1.0),
        "hv_voltage_V": _encoded(sign * complex(10.0, 0.2)),
        "point_rows": [
            {
                "material_a": {
                    "point_mm": [0.0, 0.0],
                    "values": _point_values(sign),
                },
                "material_b": {
                    "point_mm": [1.0, 0.0],
                    "values": _point_values(sign),
                },
            }
        ],
        "maximum_local_identity_errors": {
            "total_current_split_relative_error": 1.0e-12,
            "complex_conductivity_relative_error": 1.0e-12,
            "conduction_current_relative_error": 1.0e-12,
            "displacement_current_relative_error": 1.0e-12,
        },
        "independent_anc_reintegration": {
            "element_count": elements,
            "total_area_m2": 1.0,
            "total_complex_power_VA": _encoded(complex(p, q)),
            "total_energy_J": energy,
            "two_omega_energy_var": q,
            "material_rows": {
                "material_a": {
                    "element_count": elements // 2,
                    "area_m2": 0.4,
                    "complex_power_VA": _encoded(complex(0.4 * p, 0.4 * q)),
                    "energy_J": 0.4 * energy,
                },
                "material_b": {
                    "element_count": elements - elements // 2,
                    "area_m2": 0.6,
                    "complex_power_VA": _encoded(complex(0.6 * p, 0.6 * q)),
                    "energy_J": 0.6 * energy,
                },
            },
        },
    }


def good_summary() -> dict:
    fine = _case("fine", 0.5, 400, 10.001, 20.001)
    return {
        "frequency_Hz": 50.0,
        "depth_m": 1.0,
        "postprocess_contract": {"point_value_order": POINT_ORDER},
        "anc_contract": {
            "element_order": "P1_triangle",
            "material_resolution": "element -> material row",
            "independent_power_identity": "S = 0.5 integral(conj(E) dot J) dV",
        },
        "cases": [
            _case("coarse", 2.0, 100, 9.8, 20.2),
            _case("medium", 1.0, 200, 10.0, 20.0),
            fine,
            copy.deepcopy(fine) | {"case": "fine_repeat"},
            _case("fine_negative", 0.5, 400, 10.001, 20.001, sign=-1.0),
        ],
    }


def test_accepts_p1_reintegration_replay_and_sign_covariance():
    summary = good_summary()
    result = heterogeneous_current_flow_p1_reintegration_gate(summary)
    assert result["status"] == "ok"
    assert result["checks"]["material_partition_sums_close"] is True
    assert result["checks"][
        "voltage_and_fields_are_odd_while_quadratic_outputs_are_even"
    ] is True
    assert json.loads(mcp_gate(json.dumps(summary)))["status"] == "ok"


def test_rejects_a_corrupted_independent_triangle_power_integral():
    summary = good_summary()
    summary["cases"][2]["independent_anc_reintegration"][
        "total_complex_power_VA"
    ] = _encoded(complex(15.0, 20.001))
    result = heterogeneous_current_flow_p1_reintegration_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_p1_reintegration_closes"] is False


def test_rejects_a_missing_exact_replay_case():
    summary = good_summary()
    summary["cases"].pop(3)
    try:
        heterogeneous_current_flow_p1_reintegration_gate(summary)
    except ValueError as exc:
        assert "exactly five" in str(exc)
    else:
        raise AssertionError("missing replay case was accepted")
