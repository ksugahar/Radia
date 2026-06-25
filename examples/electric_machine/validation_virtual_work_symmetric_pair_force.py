"""Validation-class virtual-work force from a symmetric displacement pair.

Many force sweeps use matched displaced geometries at ``x0 - h`` and
``x0 + h``.  This example checks the compact central-pair post-processing
against a five-sample displacement sweep for both fixed-current coenergy and
fixed-flux stored energy sign conventions.

Run:

    python examples/electric_machine/validation_virtual_work_symmetric_pair_force.py
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

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    virtual_work_force_summary,
    virtual_work_symmetric_pair_force_summary,
)


OUT_JSON = HERE / "validation_virtual_work_symmetric_pair_force_summary.json"

DISPLACEMENT_M = 1.0e-3
POSITIONS_M = [-2.0e-3, -1.0e-3, 0.0, 1.0e-3, 2.0e-3]
EXPECTED_FORCE_N = 12.5
STIFFNESS_N_PER_M = 1000.0
ENERGY_OFFSET_J = 0.25


def _coenergy(x):
    return ENERGY_OFFSET_J + EXPECTED_FORCE_N * x + 0.5 * STIFFNESS_N_PER_M * x * x


def _stored_energy(x):
    return ENERGY_OFFSET_J - EXPECTED_FORCE_N * x + 0.5 * STIFFNESS_N_PER_M * x * x


def _strip_rows(summary):
    return {key: value for key, value in summary.items() if key != "rows"}


def _assert_close(actual, expected, atol=1.0e-12):
    if abs(actual - expected) > atol:
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    h = DISPLACEMENT_M
    co_pair = virtual_work_symmetric_pair_force_summary(
        h,
        _coenergy(-h),
        _coenergy(h),
        energy_kind="coenergy",
        energy_center_J=_coenergy(0.0),
    )
    stored_pair = virtual_work_symmetric_pair_force_summary(
        h,
        _stored_energy(-h),
        _stored_energy(h),
        energy_kind="stored_energy",
        energy_center_J=_stored_energy(0.0),
    )
    co_sweep = virtual_work_force_summary(
        POSITIONS_M,
        [_coenergy(x) for x in POSITIONS_M],
        energy_kind="coenergy",
    )
    stored_sweep = virtual_work_force_summary(
        POSITIONS_M,
        [_stored_energy(x) for x in POSITIONS_M],
        energy_kind="stored_energy",
    )

    expected_even = 0.5 * STIFFNESS_N_PER_M * h * h
    _assert_close(co_pair["force_N"], EXPECTED_FORCE_N)
    _assert_close(stored_pair["force_N"], EXPECTED_FORCE_N)
    _assert_close(co_pair["even_energy_residual_J"], expected_even)
    _assert_close(stored_pair["even_energy_residual_J"], expected_even)
    _assert_close(co_sweep["rows"][2]["force_N"], EXPECTED_FORCE_N)
    _assert_close(stored_sweep["rows"][2]["force_N"], EXPECTED_FORCE_N)

    return {
        "kind": "virtual_work_symmetric_pair_force",
        "validation_class": True,
        "force_learning": "matched +/- displacement energies give the center force by a two-point virtual-work central difference",
        "energy_model": {
            "displacement_m": DISPLACEMENT_M,
            "positions_m": POSITIONS_M,
            "expected_force_N": EXPECTED_FORCE_N,
            "stiffness_N_per_m": STIFFNESS_N_PER_M,
            "offset_J": ENERGY_OFFSET_J,
            "coenergy_formula": "Wco(x) = offset + F*x + 0.5*k*x^2",
            "stored_energy_formula": "W(x) = offset - F*x + 0.5*k*x^2",
        },
        "checks": {
            "coenergy_pair_force_N": co_pair["force_N"],
            "stored_energy_pair_force_N": stored_pair["force_N"],
            "coenergy_sweep_center_force_N": co_sweep["rows"][2]["force_N"],
            "stored_energy_sweep_center_force_N": stored_sweep["rows"][2]["force_N"],
            "expected_even_energy_residual_J": expected_even,
            "coenergy_even_energy_residual_J": co_pair["even_energy_residual_J"],
            "stored_energy_even_energy_residual_J": stored_pair["even_energy_residual_J"],
        },
        "coenergy_pair": co_pair,
        "stored_energy_pair": stored_pair,
        "coenergy_sweep_summary": _strip_rows(co_sweep),
        "stored_energy_sweep_summary": _strip_rows(stored_sweep),
        "sweep_center_rows": {
            "coenergy": co_sweep["rows"][2],
            "stored_energy": stored_sweep["rows"][2],
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
    print("[Virtual-work symmetric pair force]")
    print(f"  expected_force_N: {EXPECTED_FORCE_N:.12g}")
    print(f"  coenergy_pair_force_N: {checks['coenergy_pair_force_N']:.12g}")
    print(f"  stored_energy_pair_force_N: {checks['stored_energy_pair_force_N']:.12g}")
    print(f"  even_energy_residual_J: {checks['coenergy_even_energy_residual_J']:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
