"""Validation-class parallel-plate electrostatic force identity.

This example connects three equivalent public textbook views of the same
normal electrostatic force:

    C = eps A / d
    W = 0.5 C V^2
    p = 0.5 eps (V / d)^2
    F = p A = W / d

Run:

    python validation_test/electrostatics/validation_parallel_plate_electrostatic_force.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.electrostatics import (  # noqa: E402
    EPS0,
    parallel_plate_capacitor_energy_force,
)
from radia_mcp.radia_ngsolve.force import electrostatic_traction_summary  # noqa: E402
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_parallel_plate_electrostatic_force_summary.json"

EPS_R = 2.5
AREA_M2 = 2.0e-4
GAP_M = 0.5e-3
VOLTAGE_V = 120.0


def _assert_close(value: float, expected: float, rel: float = 1.0e-12) -> float:
    err = abs(value - expected)
    if err > rel * max(abs(value), abs(expected), 1.0):
        raise AssertionError(f"{value!r} != {expected!r}")
    return err


def build_summary() -> dict:
    out = parallel_plate_capacitor_energy_force(EPS_R, AREA_M2, GAP_M, VOLTAGE_V)
    field = VOLTAGE_V / GAP_M
    eps = EPS0 * EPS_R
    pressure = 0.5 * eps * field * field
    force = pressure * AREA_M2
    traction = electrostatic_traction_summary(
        (0.0, 0.0, field),
        (0.0, 0.0, 1.0),
        area_m2=AREA_M2,
        eps=eps,
    )
    half_gap = parallel_plate_capacitor_energy_force(EPS_R, AREA_M2, 0.5 * GAP_M, VOLTAGE_V)

    checks = {
        "capacitance_error_F": _assert_close(out["C"], eps * AREA_M2 / GAP_M),
        "energy_error_J": _assert_close(out["energy"], 0.5 * out["C"] * VOLTAGE_V * VOLTAGE_V),
        "pressure_error_Pa": _assert_close(out["pressure_Pa"], pressure),
        "force_error_N": _assert_close(out["force"], force),
        "force_energy_gap_error_N": _assert_close(out["force"], out["energy"] / GAP_M),
        "traction_pressure_error_Pa": _assert_close(
            traction["normal_traction_Pa"],
            out["pressure_Pa"],
        ),
        "traction_force_error_N": _assert_close(traction["force_N"][2], out["force"]),
        "half_gap_force_ratio_error": _assert_close(half_gap["force"] / out["force"], 4.0),
    }

    return {
        "kind": "parallel_plate_electrostatic_force",
        "validation_class": True,
        "inputs": {
            "eps_r": EPS_R,
            "area_m2": AREA_M2,
            "gap_m": GAP_M,
            "voltage_V": VOLTAGE_V,
        },
        "parallel_plate": out,
        "maxwell_traction": traction,
        "half_gap_force_N": half_gap["force"],
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = add_result_metadata(build_summary(), __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    plate = summary["parallel_plate"]
    print("[Parallel-plate electrostatic force]")
    print(f"  capacitance: {plate['C']:.12g} F")
    print(f"  field:       {plate['electric_field_V_per_m']:.12g} V/m")
    print(f"  pressure:    {plate['pressure_Pa']:.12g} Pa")
    print(f"  force:       {plate['force']:.12g} N")
    print("[checks]")
    for key, value in summary["checks"].items():
        print(f"  {key}: {value:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
