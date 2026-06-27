"""Validation-class coaxial capacitor electrostatic force.

Run:

    python validation_test/electrostatics/validation_coaxial_capacitor_force.py

The coaxial capacitor is the cylindrical analogue of the parallel-plate force
gate: capacitance, stored energy, Maxwell pressure, and capacitance-gradient
force all collapse to the same closed form.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.electrostatics import (  # noqa: E402
    EPS0,
    capacitance_gradient_force_summary,
    coaxial_capacitor_energy_force,
)
from result_metadata import add_result_metadata  # noqa: E402

OUT_JSON = Path(__file__).with_name("validation_coaxial_capacitor_force_summary.json")

EPS_R = 2.5
R_INNER = 0.01
R_OUTER = 0.03
LENGTH = 0.2
VOLTAGE = 120.0


def _assert_close(actual: float, expected: float, name: str, rel: float = 1.0e-12) -> float:
    scale = max(1.0, abs(expected))
    err = abs(actual - expected)
    if err > rel * scale:
        raise AssertionError(f"{name}: got {actual:.16e}, expected {expected:.16e}")
    return err


def main() -> None:
    out = coaxial_capacitor_energy_force(EPS_R, R_INNER, R_OUTER, LENGTH, VOLTAGE)
    eps = EPS0 * EPS_R
    log_ratio = math.log(R_OUTER / R_INNER)
    expected_c = 2.0 * math.pi * eps * LENGTH / log_ratio
    expected_inner_force = (
        math.pi * eps * LENGTH * VOLTAGE * VOLTAGE / (R_INNER * log_ratio * log_ratio)
    )
    expected_outer_force = (
        -math.pi * eps * LENGTH * VOLTAGE * VOLTAGE / (R_OUTER * log_ratio * log_ratio)
    )

    inner_grad = capacitance_gradient_force_summary(
        out["C"],
        out["dCdr_inner_F_per_m"],
        voltage_V=VOLTAGE,
    )
    outer_grad = capacitance_gradient_force_summary(
        out["C"],
        out["dCdr_outer_F_per_m"],
        voltage_V=VOLTAGE,
    )

    checks = {
        "capacitance_abs_error_F": _assert_close(out["C"], expected_c, "capacitance"),
        "inner_force_abs_error_N": _assert_close(
            out["inner_radius_force_N"],
            expected_inner_force,
            "inner force",
        ),
        "outer_force_abs_error_N": _assert_close(
            out["outer_radius_force_N"],
            expected_outer_force,
            "outer force",
        ),
        "inner_pressure_force_error_N": _assert_close(
            out["inner_radius_force_N"],
            out["inner_pressure_area_force_N"],
            "inner pressure force",
        ),
        "outer_pressure_force_error_N": _assert_close(
            out["outer_radius_force_N"],
            out["outer_pressure_area_force_N"],
            "outer pressure force",
        ),
        "inner_gradient_force_error_N": _assert_close(
            out["inner_radius_force_N"],
            inner_grad["fixed_voltage_force_N"],
            "inner gradient force",
        ),
        "outer_gradient_force_error_N": _assert_close(
            out["outer_radius_force_N"],
            outer_grad["fixed_voltage_force_N"],
            "outer gradient force",
        ),
    }

    summary = {
        "kind": "coaxial_capacitor_force_validation",
        "parameters": {
            "eps_r": EPS_R,
            "r_inner_m": R_INNER,
            "r_outer_m": R_OUTER,
            "length_m": LENGTH,
            "voltage_V": VOLTAGE,
        },
        "closed_form": out,
        "inner_radius_capacitance_gradient": inner_grad,
        "outer_radius_capacitance_gradient": outer_grad,
        "checks": checks,
    }
    summary = add_result_metadata(summary, __file__)
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("[Coaxial capacitor electrostatic force]")
    print(f"  C_F: {out['C']:.12e}")
    print(f"  inner_radius_force_N: {out['inner_radius_force_N']:.12e}")
    print(f"  outer_radius_force_N: {out['outer_radius_force_N']:.12e}")
    print(f"  pressure_inner_Pa: {out['pressure_inner_Pa']:.12e}")
    print(f"  wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
