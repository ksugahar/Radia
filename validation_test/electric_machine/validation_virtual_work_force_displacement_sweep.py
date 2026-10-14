"""Validation-class virtual-work force from displacement energy samples.

Fixed-current magnetic coenergy and fixed-flux stored energy have opposite
force signs:

    fixed current:  F =  dW_co/dx
    fixed flux:     F = -dW/dx

This example keeps the energy model linear so every finite-difference stencil
is exact.  It is a small, readable gate for force sweeps built from repeated
FEM solves at displaced positions.

Run:

    python validation_test/electric_machine/validation_virtual_work_force_displacement_sweep.py
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

from radia_mcp.radia_ngsolve.force import virtual_work_force_summary  # noqa: E402


OUT_JSON = HERE / "validation_virtual_work_force_displacement_sweep_summary.json"

POSITIONS_M = [-2.0e-3, -1.0e-3, 0.0, 1.0e-3, 2.0e-3]
EXPECTED_FORCE_N = 12.5
ENERGY_OFFSET_J = 0.25


def _assert_close(actual: float, expected: float, *, atol: float = 1.0e-12) -> None:
    if abs(actual - expected) > atol:
        raise AssertionError(f"{actual!r} != {expected!r}")


def _strip_rows(summary: dict) -> dict:
    return {key: value for key, value in summary.items() if key != "rows"}


def build_summary() -> dict:
    coenergy = [
        ENERGY_OFFSET_J + EXPECTED_FORCE_N * position
        for position in POSITIONS_M
    ]
    stored_energy = [
        ENERGY_OFFSET_J - EXPECTED_FORCE_N * position
        for position in POSITIONS_M
    ]

    coenergy_summary = virtual_work_force_summary(
        POSITIONS_M,
        coenergy,
        energy_kind="coenergy",
    )
    stored_energy_summary = virtual_work_force_summary(
        POSITIONS_M,
        stored_energy,
        energy_kind="stored_energy",
    )

    for row in coenergy_summary["rows"]:
        _assert_close(row["force_N"], EXPECTED_FORCE_N)
    for row in stored_energy_summary["rows"]:
        _assert_close(row["denergy_dx_N"], -EXPECTED_FORCE_N)
        _assert_close(row["force_N"], EXPECTED_FORCE_N)
    _assert_close(coenergy_summary["force_mean_N"], EXPECTED_FORCE_N)
    _assert_close(stored_energy_summary["force_mean_N"], EXPECTED_FORCE_N)

    return {
        "kind": "virtual_work_force_displacement_sweep",
        "validation_class": True,
        "force_learning": "fixed-current coenergy uses +dWco/dx; stored field energy uses -dW/dx",
        "energy_model": {
            "positions_m": POSITIONS_M,
            "offset_J": ENERGY_OFFSET_J,
            "expected_force_N": EXPECTED_FORCE_N,
            "coenergy_formula": "Wco(x) = offset + F*x",
            "stored_energy_formula": "W(x) = offset - F*x",
        },
        "checks": {
            "coenergy_force_mean_N": coenergy_summary["force_mean_N"],
            "stored_energy_force_mean_N": stored_energy_summary["force_mean_N"],
            "coenergy_force_peak_abs_N": coenergy_summary["force_peak_abs_N"],
            "stored_energy_force_peak_abs_N": stored_energy_summary["force_peak_abs_N"],
        },
        "coenergy_summary": _strip_rows(coenergy_summary),
        "stored_energy_summary": _strip_rows(stored_energy_summary),
        "sample_rows": {
            "coenergy": coenergy_summary["rows"],
            "stored_energy": stored_energy_summary["rows"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[Virtual-work force displacement sweep]")
    print(f"  expected_force_N: {EXPECTED_FORCE_N:.12g}")
    print(f"  coenergy_force_mean_N: {checks['coenergy_force_mean_N']:.12g}")
    print(f"  stored_energy_force_mean_N: {checks['stored_energy_force_mean_N']:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
