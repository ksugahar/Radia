"""Generate high-order 3D HEX finite-element convergence evidence."""

from __future__ import annotations

import json
import math
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import numpy as np
import test_hex_highorder_fem as hex_fem

OUTPUT = Path(__file__).with_name("hex_highorder_fem_results.json")


def main() -> int:
    started = time.perf_counter()
    with ng.TaskManager():
        mesh = hex_fem._hex_mesh(8)
        constant_errors = {
            str(order): float(hex_fem._solve_l2(mesh, order))
            for order in (1, 2, 3, 4)
        }
        variable_errors = {
            str(order): float(hex_fem._solve_l2_varcoeff(mesh, order))
            for order in (1, 2, 3, 4)
        }
        order2_error_n4 = float(hex_fem._solve_l2(hex_fem._hex_mesh(4), 2))
        order2_error_n8 = constant_errors["2"]

    order2_h_rate = math.log(order2_error_n4 / order2_error_n8) / math.log(2.0)
    observed = {
        "constant_p_refinement_monotone": all(
            constant_errors[str(order + 1)] < constant_errors[str(order)]
            for order in (1, 2, 3)
        ),
        "constant_order4_l2_error": constant_errors["4"],
        "constant_order2_h_rate": order2_h_rate,
        "variable_p_refinement_monotone": all(
            variable_errors[str(order + 1)] < variable_errors[str(order)]
            for order in (1, 2, 3)
        ),
        "variable_order4_l2_error": variable_errors["4"],
    }
    thresholds = {
        "constant_order4_l2_error_maximum": 1e-6,
        "constant_order2_h_rate_minimum": 2.5,
        "constant_order2_h_rate_maximum": 3.6,
        "variable_order4_l2_error_maximum": 1e-6,
    }
    passed = (
        observed["constant_p_refinement_monotone"]
        and observed["variable_p_refinement_monotone"]
        and observed["constant_order4_l2_error"]
        < thresholds["constant_order4_l2_error_maximum"]
        and thresholds["constant_order2_h_rate_minimum"]
        < observed["constant_order2_h_rate"]
        < thresholds["constant_order2_h_rate_maximum"]
        and observed["variable_order4_l2_error"]
        < thresholds["variable_order4_l2_error_maximum"]
    )
    result = {
        "schema": "radia.validation.hex-highorder-fem.v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "provenance": {
            "command": "python validation_test/cubit/generate_hex_highorder_fem_results.py",
            "python": platform.python_version(),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "numpy": np.__version__,
        },
        "configuration": {
            "element_family": "HEX",
            "mesh_subdivisions": 8,
            "polynomial_orders": [1, 2, 3, 4],
            "constant_coefficient_order2_coarse_subdivisions": 4,
        },
        "runtime_seconds": time.perf_counter() - started,
        "thresholds": thresholds,
        "observed": observed,
        "cases": {
            "constant_coefficient_l2_errors": constant_errors,
            "variable_coefficient_l2_errors": variable_errors,
            "constant_coefficient_order2_l2_error_n4": order2_error_n4,
            "constant_coefficient_order2_l2_error_n8": order2_error_n8,
        },
        "passed": passed,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
