"""Validation-class coenergy/torque angle-table consistency check.

Run:

    python examples/force_validation/validation_coenergy_torque_table_consistency.py

This example treats a torque-angle table and a coenergy-angle table as two
views of the same fixed-current virtual-work calculation.  The coenergy includes
a nonperiodic work term, ``T_mean * theta``, so a nonzero mean torque is allowed:

    W'(theta) = T_mean theta + (T_ripple / n) (1 - cos(n theta))
    T(theta) = dW'/dtheta = T_mean + T_ripple sin(n theta)

The validation summary differentiates the coenergy table, compares the selected
central-difference rows against the torque table, and records the integrated
work consistency.
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

from radia_mcp.radia_ngsolve.force import coenergy_torque_table_consistency_summary  # noqa: E402


OUT_JSON = HERE / "validation_coenergy_torque_table_consistency_summary.json"
SAMPLES = 181


def _tables() -> tuple[list[float], list[float], list[float]]:
    mean_torque_Nm = 10.0
    ripple_torque_Nm = 0.5
    ripple_order = 3
    period = 2.0 * math.pi
    angles = [period * index / (SAMPLES - 1) for index in range(SAMPLES)]
    coenergy = [
        mean_torque_Nm * theta
        + (ripple_torque_Nm / ripple_order) * (1.0 - math.cos(ripple_order * theta))
        for theta in angles
    ]
    torque = [
        mean_torque_Nm + ripple_torque_Nm * math.sin(ripple_order * theta)
        for theta in angles
    ]
    return angles, coenergy, torque


def build_summary() -> dict[str, object]:
    angles, coenergy, torque = _tables()
    consistency = coenergy_torque_table_consistency_summary(
        angles,
        coenergy,
        torque,
        periodic=False,
        torque_abs_tolerance_Nm=2.0e-3,
        torque_rel_tolerance=3.0e-4,
        comparison_stencils=("central",),
    )

    checks = {
        "samples": consistency["n_samples"],
        "reference_checked_count": consistency["reference_checked_count"],
        "status": consistency["status"],
        "max_torque_abs_error_Nm": consistency["max_torque_abs_error_Nm"],
        "max_torque_rel_error": consistency["max_torque_rel_error"],
        "coenergy_delta_J": consistency["coenergy_delta_J"],
        "reference_torque_trapezoid_work_J": consistency["reference_torque_trapezoid_work_J"],
        "reference_work_minus_coenergy_delta_J": consistency["reference_work_minus_coenergy_delta_J"],
    }

    assert checks["samples"] == SAMPLES
    assert checks["reference_checked_count"] == SAMPLES - 2
    assert checks["status"] == "ok"
    assert checks["max_torque_abs_error_Nm"] < 2.0e-3
    assert abs(checks["reference_work_minus_coenergy_delta_J"]) < 1.0e-12

    return {
        "kind": "coenergy_torque_table_consistency_validation",
        "validation_class": True,
        "learning_theme": (
            "torque-angle tables and coenergy-angle tables should agree under "
            "the fixed-current virtual-work derivative"
        ),
        "checks": checks,
        "angles_rad": angles,
        "coenergy_J": coenergy,
        "torque_Nm": torque,
        "consistency": consistency,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[coenergy torque table consistency]")
    print(
        f"  samples={checks['samples']} "
        f"checked={checks['reference_checked_count']} status={checks['status']}"
    )
    print(f"  max_torque_abs_error_Nm={checks['max_torque_abs_error_Nm']:.3e}")
    print(f"  max_torque_rel_error={checks['max_torque_rel_error']:.3e}")
    print(
        "  work_error_J="
        f"{checks['reference_work_minus_coenergy_delta_J']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
