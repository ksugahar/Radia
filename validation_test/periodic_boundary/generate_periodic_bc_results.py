"""Generate periodic and anti-periodic manufactured-solution evidence."""

from __future__ import annotations

import json
import math
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import test_periodic_bc as periodic

OUTPUT = Path(__file__).with_name("periodic_bc_results.json")


def main() -> int:
    started = time.perf_counter()
    with ng.TaskManager():
        periodic_error_n8 = float(periodic._periodic_error(8))
        periodic_error_n16 = float(periodic._periodic_error(16))
        periodic_error_n32 = float(periodic._periodic_error(32))
        mesh, field, exact = periodic._antiperiodic_solution()
        antiperiodic_error = float(periodic._l2(field, exact, mesh))
        sign_flip_errors = [
            float(abs(field(mesh(1.0, y)) + field(mesh(0.0, y))))
            for y in (0.3, 0.5, 0.7)
        ]

    h_rate = math.log(periodic_error_n8 / periodic_error_n16) / math.log(2.0)
    observed = {
        "periodic_l2_error_n8": periodic_error_n8,
        "periodic_l2_error_n16": periodic_error_n16,
        "periodic_l2_error_n32": periodic_error_n32,
        "periodic_order3_h_rate": h_rate,
        "antiperiodic_l2_error_n32": antiperiodic_error,
        "antiperiodic_sign_flip_errors": sign_flip_errors,
        "antiperiodic_maximum_sign_flip_error": max(sign_flip_errors),
    }
    thresholds = {
        "maximum_periodic_l2_error_n32": 1e-4,
        "minimum_periodic_order3_h_rate": 3.0,
        "maximum_antiperiodic_l2_error_n32": 1e-4,
        "maximum_antiperiodic_sign_flip_error": 1e-3,
    }
    passed = (
        periodic_error_n32 < thresholds["maximum_periodic_l2_error_n32"]
        and h_rate > thresholds["minimum_periodic_order3_h_rate"]
        and antiperiodic_error < thresholds["maximum_antiperiodic_l2_error_n32"]
        and max(sign_flip_errors)
        < thresholds["maximum_antiperiodic_sign_flip_error"]
    )
    result = {
        "schema": "radia.validation.periodic-boundary.v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "provenance": {
            "command": (
                "python validation_test/periodic_boundary/"
                "generate_periodic_bc_results.py"
            ),
            "python": platform.python_version(),
            "ngsolve": getattr(ng, "__version__", "unknown"),
        },
        "configuration": {
            "element_family": "triangle",
            "polynomial_order": 3,
            "periodic_mesh_subdivisions": [8, 16, 32],
            "antiperiodic_mesh_subdivisions": 32,
            "anti_periodic_phase": -1,
        },
        "runtime_seconds": time.perf_counter() - started,
        "thresholds": thresholds,
        "observed": observed,
        "passed": passed,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
