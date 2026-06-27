"""Validation-class sampled air-gap Maxwell shear torque.

Machine solvers commonly export air-gap samples of radial and tangential flux
density.  This example validates the post-processing identity

    tau(theta) = Br(theta) Bt(theta) / mu0
    T = r^2 L integral tau(theta) dtheta

against two analytic gates: a uniform field and a sinusoidal harmonic pair.

Run:

    python validation_test/electric_machine/validation_sampled_air_gap_shear_torque.py
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

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    MU0,
    air_gap_shear_torque,
    air_gap_shear_torque_from_angle_samples,
)


OUT_JSON = HERE / "validation_sampled_air_gap_shear_torque_summary.json"

RADIUS_M = 0.04
AXIAL_LENGTH_M = 0.12
SAMPLES = 720


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-9) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def _compact(summary: dict, sample_segments: int = 6) -> dict:
    compact = {key: value for key, value in summary.items() if key != "rows"}
    compact["sample_rows"] = summary["rows"][:sample_segments]
    return compact


def build_summary() -> dict:
    angles = [2.0 * math.pi * index / SAMPLES for index in range(SAMPLES)]

    uniform_br = 0.8
    uniform_bt = 0.1
    uniform = air_gap_shear_torque_from_angle_samples(
        angles,
        [uniform_br] * SAMPLES,
        [uniform_bt] * SAMPLES,
        RADIUS_M,
        axial_length_m=AXIAL_LENGTH_M,
    )
    uniform_expected = air_gap_shear_torque(
        uniform_br,
        uniform_bt,
        RADIUS_M,
        axial_length_m=AXIAL_LENGTH_M,
    )

    harmonic = 3
    phase_rad = math.pi / 4.0
    br0 = 0.9
    bt0 = 0.22
    br = [br0 * math.cos(harmonic * angle) for angle in angles]
    bt = [bt0 * math.cos(harmonic * angle + phase_rad) for angle in angles]
    sinusoidal = air_gap_shear_torque_from_angle_samples(
        angles,
        br,
        bt,
        RADIUS_M,
        axial_length_m=AXIAL_LENGTH_M,
    )
    sinusoidal_average_shear = 0.5 * br0 * bt0 * math.cos(phase_rad) / MU0
    sinusoidal_expected = (
        sinusoidal_average_shear
        * RADIUS_M
        * RADIUS_M
        * AXIAL_LENGTH_M
        * 2.0
        * math.pi
    )

    _assert_close(uniform["torque_Nm"], uniform_expected)
    _assert_close(sinusoidal["average_shear_stress_Pa"], sinusoidal_average_shear, rtol=2.0e-12)
    _assert_close(sinusoidal["torque_Nm"], sinusoidal_expected, rtol=2.0e-12)

    return {
        "kind": "sampled_air_gap_shear_torque_validation",
        "validation_class": True,
        "force_learning": "machine torque from sampled air-gap Br/Bt: T = r^2 L integral Br*Bt/mu0 dtheta",
        "geometry": {
            "radius_m": RADIUS_M,
            "axial_length_m": AXIAL_LENGTH_M,
            "samples": SAMPLES,
        },
        "uniform_gate": {
            "B_radial_T": uniform_br,
            "B_tangential_T": uniform_bt,
            "expected_torque_Nm": uniform_expected,
            "computed": _compact(uniform),
            "abs_error_Nm": abs(uniform["torque_Nm"] - uniform_expected),
        },
        "sinusoidal_gate": {
            "B_radial_amplitude_T": br0,
            "B_tangential_amplitude_T": bt0,
            "harmonic_order": harmonic,
            "phase_rad": phase_rad,
            "expected_average_shear_Pa": sinusoidal_average_shear,
            "expected_torque_Nm": sinusoidal_expected,
            "computed": _compact(sinusoidal),
            "average_shear_abs_error_Pa": abs(
                sinusoidal["average_shear_stress_Pa"] - sinusoidal_average_shear
            ),
            "torque_abs_error_Nm": abs(sinusoidal["torque_Nm"] - sinusoidal_expected),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[Sampled air-gap shear torque]")
    print(f"  uniform_expected_torque_Nm: {summary['uniform_gate']['expected_torque_Nm']:.12g}")
    print(f"  uniform_computed_torque_Nm: {summary['uniform_gate']['computed']['torque_Nm']:.12g}")
    print(f"  sinusoidal_expected_torque_Nm: {summary['sinusoidal_gate']['expected_torque_Nm']:.12g}")
    print(f"  sinusoidal_computed_torque_Nm: {summary['sinusoidal_gate']['computed']['torque_Nm']:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
