"""Validation-class radiation force from normal-incidence scattering powers.

For a one-sided wave normally incident on a scatterer with reflectance R and
transmittance T,

    F = (1 + R - T) P_inc / c

Equivalently, with absorptance A = 1 - R - T, the same force is
``(A + 2R) P_inc/c``.  This is the RF port/S-parameter form of the radiation
pressure examples.

Run:

    python validation_test/rf_waveguide/validation_scattering_radiation_force.py
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

from result_metadata import add_result_metadata  # noqa: E402

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    radiation_force_from_power,
    radiation_scattering_force_summary,
)


OUT_JSON = HERE / "validation_scattering_radiation_force_summary.json"
POWER_W = 5.0

CASES = (
    {"label": "perfect_absorber", "R": 0.0, "T": 0.0},
    {"label": "perfect_reflector", "R": 1.0, "T": 0.0},
    {"label": "transparent_through_line", "R": 0.0, "T": 1.0},
    {"label": "lossless_partial_reflector", "R": 0.25, "T": 0.75},
    {"label": "lossy_partial_reflector", "R": 0.5, "T": 0.25},
)


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-18) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    rows = []
    for case in CASES:
        summary = radiation_scattering_force_summary(
            POWER_W,
            reflectance=case["R"],
            transmittance=case["T"],
        )
        expected_factor = 1.0 + case["R"] - case["T"]
        expected_force = expected_factor * POWER_W / C0
        equivalent = radiation_force_from_power(
            POWER_W,
            absorptance=summary["absorptance"],
            reflectance=summary["reflectance"],
        )
        _assert_close(summary["momentum_transfer_factor"], expected_factor)
        _assert_close(summary["force_N"], expected_force)
        _assert_close(summary["force_N"], equivalent)
        rows.append({
            "label": case["label"],
            "summary": summary,
            "expected_factor": expected_factor,
            "expected_force_N": expected_force,
            "force_abs_error_N": abs(summary["force_N"] - expected_force),
        })

    by_label = {row["label"]: row for row in rows}
    return {
        "kind": "scattering_radiation_force_validation",
        "validation_class": True,
        "force_learning": "normal-incidence scattering force: F=(1+R-T)P/c=(A+2R)P/c",
        "incident_power_W": POWER_W,
        "rows": rows,
        "checks": {
            "absorber_force_N": by_label["perfect_absorber"]["summary"]["force_N"],
            "reflector_force_N": by_label["perfect_reflector"]["summary"]["force_N"],
            "transparent_force_N": by_label["transparent_through_line"]["summary"]["force_N"],
            "lossy_partial_factor": by_label["lossy_partial_reflector"]["summary"]["momentum_transfer_factor"],
            "max_force_abs_error_N": max(row["force_abs_error_N"] for row in rows),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary = add_result_metadata(summary, __file__)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[Scattering radiation force]")
    for row in summary["rows"]:
        item = row["summary"]
        print(
            f"  {row['label']}: R={item['reflectance']:.3g}, "
            f"T={item['transmittance']:.3g}, "
            f"A={item['absorptance']:.3g}, "
            f"factor={item['momentum_transfer_factor']:.6g}, "
            f"F={item['force_N']:.6g} N"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
