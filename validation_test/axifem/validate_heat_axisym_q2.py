"""Validate axisymmetric thermal Q2 against a Bessel-series cylinder oracle.

This is a validation lane, not a per-commit regression test.  It records the
Q1/Q2 temperature and near-axis-gradient comparison as a durable JSON result.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import tempfile

import numpy as np
from scipy.special import j0, jn_zeros


REPO = Path(__file__).resolve().parents[2]
PANELS = REPO / "src" / "radia" / "panels"
FIXTURES = REPO / "validation_test" / "panels" / "fixtures"
sys.path.insert(0, str(PANELS))
sys.path.insert(0, str(FIXTURES))

from calc_heat_axisym import solve_heat_axisym  # noqa: E402
from generate_heat_cylinder_axisym import main as generate_fixture  # noqa: E402


RADIUS_M = 0.025
HEIGHT_M = 0.025
HEAT_FLUX_W_PER_M2 = 1.0e5
CONDUCTIVITY_W_PER_MK = 46.6
DENSITY_KG_PER_M3 = 7800.0
HEAT_CAPACITY_J_PER_KGK = 467.0
INITIAL_TEMPERATURE_C = 25.0
END_TIME_S = 10.0
TIME_STEP_S = 0.1


def _analytical_temperature(radius_m: float, terms: int = 100) -> float:
    """Constant-flux cylinder solution using roots of J1(lambda)=0."""
    roots = jn_zeros(1, terms)
    diffusivity = CONDUCTIVITY_W_PER_MK / (
        DENSITY_KG_PER_M3 * HEAT_CAPACITY_J_PER_KGK
    )
    particular = HEAT_FLUX_W_PER_M2 / (
        2.0 * CONDUCTIVITY_W_PER_MK * RADIUS_M
    ) * (radius_m**2 - 0.5 * RADIUS_M**2)
    coefficients = -2.0 * HEAT_FLUX_W_PER_M2 * RADIUS_M / (
        CONDUCTIVITY_W_PER_MK * roots**2 * j0(roots)
    )
    transient = np.sum(
        coefficients
        * j0(roots * radius_m / RADIUS_M)
        * np.exp(-diffusivity * roots**2 * END_TIME_S / RADIUS_M**2)
    )
    mean_rise = 2.0 * HEAT_FLUX_W_PER_M2 * END_TIME_S / (
        DENSITY_KG_PER_M3 * HEAT_CAPACITY_J_PER_KGK * RADIUS_M
    )
    return float(INITIAL_TEMPERATURE_C + mean_rise + particular + transient)


def _run_order(order: int, output_dir: Path, reference: dict[str, float]):
    from ngsolve import GridFunction, H1, Mesh, grad

    fixture = FIXTURES / "heat_workpiece_cylinder_R25_H25_axisym.vol"
    result = solve_heat_axisym(
        str(fixture),
        material="steel",
        h_conv=0.0,
        t_ext=INITIAL_TEMPERATURE_C,
        t_initial=INITIAL_TEMPERATURE_C,
        surface_label="outer",
        q_uniform=HEAT_FLUX_W_PER_M2,
        dt=TIME_STEP_S,
        t_end=END_TIME_S,
        linear_solver="sparsecholesky",
        fes_order=order,
        msh_output=str(output_dir / f"heat_axisym_q{order}.msh"),
    )
    mesh = Mesh(result["heat_vol_file"])
    temperature = GridFunction(H1(mesh, order=order))
    temperature.Load(result["T_sol_file"])
    axis_temperature = float(temperature(mesh(0.0, 0.0)))
    surface_temperature = float(temperature(mesh(RADIUS_M, 0.0)))
    axis_gradient = float(grad(temperature)(mesh(0.0, 0.0))[0])
    return {
        "fes_order": order,
        "ndof": int(result["ndof"]),
        "axis_temperature_C": axis_temperature,
        "surface_temperature_C": surface_temperature,
        "axis_temperature_abs_error_C": abs(
            axis_temperature - reference["axis_temperature_C"]
        ),
        "surface_temperature_abs_error_C": abs(
            surface_temperature - reference["surface_temperature_C"]
        ),
        "radial_gradient_at_axis_C_per_m": axis_gradient,
        "reported_min_C": float(result["T_min_C"]),
        "reported_max_C": float(result["T_max_C"]),
        "extrema_sampling": result["T_extrema"],
        "solver_time_s": float(result["t_total_s"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "results_heat_axisym_q2_near_axis_20260903.json"
        ),
    )
    args = parser.parse_args()

    generate_fixture()
    reference = {
        "axis_temperature_C": _analytical_temperature(0.0),
        "surface_temperature_C": _analytical_temperature(RADIUS_M),
        "series_terms": 100,
    }
    with tempfile.TemporaryDirectory(prefix="radia-axisym-heat-") as directory:
        output_dir = Path(directory)
        q1 = _run_order(1, output_dir, reference)
        q2 = _run_order(2, output_dir, reference)

    slope_reduction = abs(q1["radial_gradient_at_axis_C_per_m"]) / max(
        abs(q2["radial_gradient_at_axis_C_per_m"]), 1.0e-15
    )
    checks = {
        "q2_axis_abs_error_below_0_1_C": (
            q2["axis_temperature_abs_error_C"] < 0.1
        ),
        "q2_surface_abs_error_below_0_1_C": (
            q2["surface_temperature_abs_error_C"] < 0.1
        ),
        "q2_axis_gradient_below_1_C_per_m": (
            abs(q2["radial_gradient_at_axis_C_per_m"]) < 1.0
        ),
        "q2_axis_slope_reduction_above_100": slope_reduction > 100.0,
    }
    payload = {
        "schema": "radia.validation.axisym-heat-q2.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "pass" if all(checks.values()) else "fail",
        "runtime": {
            "python": platform.python_version(),
            "ngsolve": importlib.metadata.version("ngsolve"),
            "platform": platform.platform(),
        },
        "problem": {
            "radius_m": RADIUS_M,
            "height_m": HEIGHT_M,
            "heat_flux_W_per_m2": HEAT_FLUX_W_PER_M2,
            "density_kg_per_m3": DENSITY_KG_PER_M3,
            "heat_capacity_J_per_kgK": HEAT_CAPACITY_J_PER_KGK,
            "thermal_conductivity_W_per_mK": CONDUCTIVITY_W_PER_MK,
            "initial_temperature_C": INITIAL_TEMPERATURE_C,
            "end_time_s": END_TIME_S,
            "backward_euler_dt_s": TIME_STEP_S,
            "radial_cells": 9,
            "axial_cells": 9,
        },
        "oracle": {
            "kind": "Bessel series for a constant-flux solid cylinder",
            **reference,
        },
        "results": [q1, q2],
        "q1_to_q2_axis_slope_reduction": slope_reduction,
        "checks": checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, sort_keys=True))
    return 0 if payload["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
