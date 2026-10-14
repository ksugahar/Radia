"""Generate representative eddy-current closed-form validation evidence."""

from __future__ import annotations

import cmath
import json
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import numpy as np
import scipy
import test_eddy_ac_resistance as ac_resistance
import test_eddy_cylinder_axial_field as cylinder
import test_eddy_skin_effect as skin
import test_lamination_eddy_loss as lamination
import test_tomega_eddy as tomega

OUTPUT = Path(__file__).with_name("diffusion_closed_form_results.json")


def _relative_error(actual: complex, expected: complex) -> float:
    return float(abs(actual - expected) / abs(expected))


def _round_wire_cases() -> list[dict[str, float]]:
    cases = []
    for a_over_delta in (1.0, 3.0, 5.0):
        fem = float(ac_resistance._fem_rac_over_rdc(a_over_delta))
        exact = float(ac_resistance._bessel_rac_over_rdc(a_over_delta))
        cases.append(
            {
                "a_over_delta": a_over_delta,
                "fem_rac_over_rdc": fem,
                "exact_rac_over_rdc": exact,
                "relative_error": _relative_error(fem, exact),
            }
        )
    return cases


def _cylinder_cases() -> list[dict[str, object]]:
    cases = []
    radius = 0.01
    for a_over_delta in (0.5, 2.0, 5.0):
        mesh, field, k = cylinder._fem1d_profile(radius, a_over_delta)
        errors = []
        for fraction in (0.0, 0.25, 0.5, 0.75):
            fem = complex(field(mesh(fraction * radius)))
            exact = cylinder._bessel_Hz(fraction * radius, radius, k)
            errors.append(_relative_error(fem, exact))
        cases.append(
            {
                "a_over_delta": a_over_delta,
                "sample_radius_fractions": [0.0, 0.25, 0.5, 0.75],
                "maximum_relative_error": max(errors),
            }
        )
    return cases


def _skin_profile_case() -> dict[str, float]:
    mesh, field = skin._solve_skin()
    propagation = (1.0 + 1.0j) / skin.DELTA
    errors = []
    for x_over_delta in np.arange(0.2, 2.51, 0.1):
        x = float(x_over_delta) * skin.DELTA
        exact = skin.H0 * cmath.exp(-propagation * x)
        errors.append(abs(complex(field(mesh(x))) - exact))
    surface_impedance = (
        -(1.0 / skin.SIGMA)
        * skin._grad0(field, mesh)
        / complex(field(mesh(0.0)))
    )
    exact_impedance = (1.0 + 1.0j) / (skin.SIGMA * skin.DELTA)
    return {
        "maximum_absolute_profile_error": float(max(errors)),
        "surface_impedance_relative_error": _relative_error(
            surface_impedance, exact_impedance
        ),
    }


def _lamination_cases() -> list[dict[str, float]]:
    cases = []
    for thickness_m, frequency_hz in ((0.5e-3, 200.0), (0.35e-3, 400.0)):
        fem = float(lamination._fem_loss(thickness_m, frequency_hz))
        exact = float(lamination._exact_loss(thickness_m, frequency_hz))
        cases.append(
            {
                "thickness_m": thickness_m,
                "frequency_hz": frequency_hz,
                "fem_loss_w_per_m3": fem,
                "exact_loss_w_per_m3": exact,
                "relative_error": _relative_error(fem, exact),
            }
        )
    return cases


def _tomega_cases() -> list[dict[str, object]]:
    cases = []
    radius = 0.01
    for a_over_delta in (0.5, 2.0, 5.0):
        mesh, reaction, k = tomega._solve_T_1d(radius, a_over_delta)
        errors = []
        for fraction in (0.0, 0.25, 0.5, 0.75):
            field = complex(reaction(mesh(fraction * radius))) + tomega.H0
            exact = tomega._bessel_Hz(fraction * radius, radius, k)
            errors.append(_relative_error(field, exact))
        cases.append(
            {
                "a_over_delta": a_over_delta,
                "sample_radius_fractions": [0.0, 0.25, 0.5, 0.75],
                "maximum_relative_error": max(errors),
            }
        )
    return cases


def main() -> int:
    started = time.perf_counter()
    with ng.TaskManager():
        round_wire = _round_wire_cases()
        cylinder_cases = _cylinder_cases()
        skin_profile = _skin_profile_case()
        lamination_cases = _lamination_cases()
        tomega_cases = _tomega_cases()

    thresholds = {
        "round_wire_maximum_relative_error": 1e-6,
        "cylinder_maximum_relative_error": 1e-5,
        "skin_maximum_absolute_profile_error": 5e-4,
        "skin_surface_impedance_relative_error": 1e-3,
        "lamination_maximum_relative_error": 1e-3,
        "tomega_maximum_relative_error": 1e-5,
    }
    observed = {
        "round_wire_maximum_relative_error": max(
            case["relative_error"] for case in round_wire
        ),
        "cylinder_maximum_relative_error": max(
            case["maximum_relative_error"] for case in cylinder_cases
        ),
        "skin_maximum_absolute_profile_error": skin_profile[
            "maximum_absolute_profile_error"
        ],
        "skin_surface_impedance_relative_error": skin_profile[
            "surface_impedance_relative_error"
        ],
        "lamination_maximum_relative_error": max(
            case["relative_error"] for case in lamination_cases
        ),
        "tomega_maximum_relative_error": max(
            case["maximum_relative_error"] for case in tomega_cases
        ),
    }
    passed = all(observed[name] < limit for name, limit in thresholds.items())
    result = {
        "schema": "radia.validation.eddy-diffusion-closed-form.v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "provenance": {
            "command": (
                "python validation_test/eddy_current_analytical_validation/"
                "generate_diffusion_closed_form_results.py"
            ),
            "python": platform.python_version(),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "runtime_seconds": time.perf_counter() - started,
        "thresholds": thresholds,
        "observed": observed,
        "cases": {
            "round_wire_ac_resistance": round_wire,
            "cylinder_axial_field": cylinder_cases,
            "skin_profile": skin_profile,
            "lamination_loss": lamination_cases,
            "tomega_reconstruction": tomega_cases,
        },
        "passed": passed,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
