"""Generate nonlinear magnetostatic Newton validation evidence."""

from __future__ import annotations

import json
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import test_nonlinear_magnetostatic_newton as nonlinear

OUTPUT = Path(__file__).with_name("nonlinear_newton_results.json")


def main() -> int:
    started = time.perf_counter()
    with ng.TaskManager():
        default_error = float(nonlinear._newton_error())
        loose_error = float(nonlinear._newton_error(1e-6))
        tight_error = float(nonlinear._newton_error(1e-11))
        saturation_spread = float(nonlinear._saturation_relative_spread())
        residual_history, recovered_root_error = nonlinear._newton_convergence()

    residual_history = [float(value) for value in residual_history]
    drop_ratios = [
        residual_history[index] / residual_history[index + 1]
        for index in range(len(residual_history) - 1)
        if residual_history[index + 1] > 0
    ]
    observed = {
        "default_l2_error": default_error,
        "loose_tolerance_l2_error": loose_error,
        "tight_tolerance_l2_error": tight_error,
        "saturation_relative_spread": saturation_spread,
        "residual_history": residual_history,
        "residual_drop_ratios": drop_ratios,
        "recovered_root_l2_error": float(recovered_root_error),
    }
    thresholds = {
        "maximum_l2_error": 1e-5,
        "minimum_saturation_relative_spread": 0.1,
        "maximum_iterations": 6,
        "maximum_final_residual": 1e-10,
        "minimum_quadratic_drop_ratio": 1e3,
    }
    passed = (
        max(default_error, loose_error, tight_error, recovered_root_error)
        < thresholds["maximum_l2_error"]
        and tight_error <= loose_error * (1 + 1e-6)
        and saturation_spread > thresholds["minimum_saturation_relative_spread"]
        and len(residual_history) <= thresholds["maximum_iterations"]
        and residual_history[-1] < thresholds["maximum_final_residual"]
        and max(drop_ratios) > thresholds["minimum_quadratic_drop_ratio"]
        and drop_ratios[-1] > drop_ratios[0]
    )
    result = {
        "schema": "radia.validation.nonlinear-magnetostatic-newton.v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "provenance": {
            "command": (
                "python validation_test/magnetostatics/"
                "generate_nonlinear_newton_results.py"
            ),
            "python": platform.python_version(),
            "ngsolve": getattr(ng, "__version__", "unknown"),
        },
        "configuration": {
            "mesh_subdivisions": 16,
            "polynomial_order": 3,
            "saturation_alpha": nonlinear.ALPHA,
            "loose_newton_tolerance": 1e-6,
            "tight_newton_tolerance": 1e-11,
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
