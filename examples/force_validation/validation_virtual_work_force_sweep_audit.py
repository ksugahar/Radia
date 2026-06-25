"""Validation-class virtual-work force sweep audit.

Run:

    python examples/force_validation/validation_virtual_work_force_sweep_audit.py

This example samples the separation-dependent coenergy per unit length of two
long line currents,

    W'(d) = -mu0 I1 I2 log(d/d_ref)/(2*pi),

and audits the finite-difference force table against the analytic radial force

    F(d) = dW'/dd = -mu0 I1 I2/(2*pi*d).

The lesson is the sweep form: solver outputs usually arrive as a table of
energy samples, so the validation object should preserve stencils, curvature,
and reference-force errors rather than only one pass/fail number.
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
    virtual_work_force_sweep_audit_summary,
)


OUT_JSON = HERE / "validation_virtual_work_force_sweep_audit_summary.json"


def _line_current_coenergy_and_force(
    separations_m: list[float],
    current1_A: float,
    current2_A: float,
    reference_separation_m: float,
) -> tuple[list[float], list[float]]:
    coefficient = MU0 * current1_A * current2_A / (2.0 * math.pi)
    coenergy = [
        -coefficient * math.log(distance / reference_separation_m)
        for distance in separations_m
    ]
    radial_force = [-coefficient / distance for distance in separations_m]
    return coenergy, radial_force


def build_summary() -> dict[str, object]:
    current1_A = 3.0
    current2_A = 2.0
    separations_m = [0.020 + 0.001 * index for index in range(13)]
    reference_separation_m = separations_m[len(separations_m) // 2]
    coenergy, reference_force = _line_current_coenergy_and_force(
        separations_m,
        current1_A,
        current2_A,
        reference_separation_m,
    )
    sweep = virtual_work_force_sweep_audit_summary(
        separations_m,
        coenergy,
        energy_kind="coenergy",
        reference_force_N=reference_force,
        force_abs_tolerance_N=1.0e-15,
        force_rel_tolerance=1.5e-3,
        comparison_stencils=("central",),
    )

    checks = {
        "n_samples": sweep["n_samples"],
        "reference_checked_count": sweep["reference_checked_count"],
        "status": sweep["status"],
        "max_reference_force_rel_error": sweep["max_reference_force_rel_error"],
        "force_min_N_per_m": sweep["force_min_N"],
        "force_max_N_per_m": sweep["force_max_N"],
        "max_abs_force_gradient_N_per_m2": sweep["max_abs_force_gradient_N_per_m"],
    }

    assert checks["n_samples"] == len(separations_m)
    assert checks["reference_checked_count"] == len(separations_m) - 2
    assert checks["status"] == "ok"
    assert checks["max_reference_force_rel_error"] < 1.5e-3
    assert checks["force_min_N_per_m"] < checks["force_max_N_per_m"] < 0.0
    assert checks["max_abs_force_gradient_N_per_m2"] > 0.0

    return {
        "kind": "virtual_work_force_sweep_audit_validation",
        "validation_class": True,
        "force_learning": (
            "energy/coenergy sweeps should preserve finite-difference stencils, "
            "force-gradient estimates, and reference-force errors"
        ),
        "current1_A": current1_A,
        "current2_A": current2_A,
        "reference_separation_m": reference_separation_m,
        "checks": checks,
        "separations_m": separations_m,
        "coenergy_J_per_m": coenergy,
        "analytic_radial_force_N_per_m": reference_force,
        "sweep": sweep,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[virtual work force sweep audit]")
    print(
        f"  samples={checks['n_samples']} "
        f"checked={checks['reference_checked_count']} status={checks['status']}"
    )
    print(
        f"  force range=[{checks['force_min_N_per_m']:.12g}, "
        f"{checks['force_max_N_per_m']:.12g}] N/m"
    )
    print(
        "  max_reference_force_rel_error="
        f"{checks['max_reference_force_rel_error']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
