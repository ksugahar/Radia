"""Validation-class electrostatic force from a capacitance gradient.

For a displacement coordinate x and capacitance C(x),

    fixed voltage: F_x = 0.5 V^2 dC/dx
    fixed charge:  F_x = 0.5 Q^2 C^-2 dC/dx

The sign is along increasing x.  If x is a gap or height, dC/dx is usually
negative, so the force is negative: attraction toward the smaller gap.

Run:

    python validation_test/electrostatics/validation_capacitance_gradient_force.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.electrostatics import (  # noqa: E402
    EPS0,
    capacitance_gradient_force_summary,
    parallel_plate_capacitor_energy_force,
    sphere_above_plane_capacitance,
)
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_capacitance_gradient_force_summary.json"

EPS_R = 2.5
AREA_M2 = 2.0e-4
GAP_M = 5.0e-4
VOLTAGE_V = 120.0
SPHERE_RADIUS_M = 0.01
SPHERE_HEIGHT_M = 2.0 * SPHERE_RADIUS_M


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-18) -> float:
    error = abs(actual - expected)
    if error > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")
    return error


def _sphere_capacitance(height_m: float) -> float:
    return sphere_above_plane_capacitance(
        SPHERE_RADIUS_M,
        height_m,
        eps_r=1.0,
        n_terms=250,
    )["C"]


def build_summary() -> dict:
    plate = parallel_plate_capacitor_energy_force(EPS_R, AREA_M2, GAP_M, VOLTAGE_V)
    dcdgap = -plate["C"] / GAP_M
    charge = plate["C"] * VOLTAGE_V
    gap_force = capacitance_gradient_force_summary(
        plate["C"],
        dcdgap,
        voltage_V=VOLTAGE_V,
        charge_C=charge,
    )
    closing_force = capacitance_gradient_force_summary(
        plate["C"],
        -dcdgap,
        voltage_V=VOLTAGE_V,
        charge_C=charge,
    )
    plate_checks = {
        "gap_force_matches_minus_pressure_force_N": _assert_close(
            gap_force["fixed_voltage_force_N"],
            -plate["force"],
        ),
        "closing_force_matches_pressure_force_N": _assert_close(
            closing_force["fixed_voltage_force_N"],
            plate["force"],
        ),
        "fixed_voltage_charge_route_abs_error_N": _assert_close(
            gap_force["fixed_voltage_force_N"],
            gap_force["fixed_charge_force_N"],
        ),
    }

    step = 1.0e-5 * SPHERE_HEIGHT_M
    c0 = _sphere_capacitance(SPHERE_HEIGHT_M)
    c_minus = _sphere_capacitance(SPHERE_HEIGHT_M - step)
    c_plus = _sphere_capacitance(SPHERE_HEIGHT_M + step)
    dcdh = (c_plus - c_minus) / (2.0 * step)
    sphere_force = capacitance_gradient_force_summary(
        c0,
        dcdh,
        voltage_V=VOLTAGE_V,
    )
    if not sphere_force["fixed_voltage_force_N"] < 0.0:
        raise AssertionError("sphere height force should be negative: attraction toward the plane")

    return {
        "kind": "capacitance_gradient_force_validation",
        "validation_class": True,
        "force_learning": "fixed-voltage electrostatic force is +0.5 V^2 dC/dx along the chosen coordinate",
        "parallel_plate": {
            "eps_r": EPS_R,
            "area_m2": AREA_M2,
            "gap_m": GAP_M,
            "voltage_V": VOLTAGE_V,
            "closed_form": plate,
            "dCdgap_F_per_m": dcdgap,
            "gap_coordinate_force": gap_force,
            "closing_coordinate_force": closing_force,
            "checks": plate_checks,
        },
        "sphere_ground_plane": {
            "radius_m": SPHERE_RADIUS_M,
            "height_m": SPHERE_HEIGHT_M,
            "voltage_V": VOLTAGE_V,
            "step_m": step,
            "capacitance_F": c0,
            "dCdh_F_per_m": dcdh,
            "height_coordinate_force": sphere_force,
            "attractive_force_magnitude_N": -sphere_force["fixed_voltage_force_N"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = add_result_metadata(build_summary(), __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    plate = summary["parallel_plate"]
    sphere = summary["sphere_ground_plane"]
    print("[Capacitance-gradient electrostatic force]")
    print(f"  plate_gap_force_N: {plate['gap_coordinate_force']['fixed_voltage_force_N']:.12g}")
    print(f"  plate_closing_force_N: {plate['closing_coordinate_force']['fixed_voltage_force_N']:.12g}")
    print(f"  plate_pressure_force_N: {plate['closed_form']['force']:.12g}")
    print(f"  sphere_height_force_N: {sphere['height_coordinate_force']['fixed_voltage_force_N']:.12g}")
    print(f"  sphere_attraction_magnitude_N: {sphere['attractive_force_magnitude_N']:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
