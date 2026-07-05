"""Validation-class axisymmetric-reference -> 3-D force-vector gate.

This script demonstrates the route we want for heavier FEM validation:

1. get a full-revolution axisymmetric force reference,
2. compute or load a 3-D force vector,
3. check axial agreement and transverse cancellation with a reusable gate.

The executable public example uses two coaxial current loops.  The
axisymmetric reference is the exact elliptic-integral mutual-inductance
derivative, while the 3-D side is a direct Biot-Savart line-segment quadrature
around the two loops.  A future NGSolve/Cubit/.vol run can replace only the
3-D force-vector producer and keep the same artifact contract.

Run:

    python validation_test/force_validation/validation_axisymmetric_to_3d_force_gate.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import scipy.special as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.axisymmetric_3d_validation import (  # noqa: E402
    axisymmetric_to_3d_force_gate,
    axisymmetric_to_3d_validation_plan,
)
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_axisymmetric_to_3d_force_gate_summary.json"
MU0 = 4.0e-7 * math.pi
LOOP_RADIUS_M = 1.0
LOOP_SEPARATION_M = 1.5
CURRENT1_A = 1.0
CURRENT2_A = 1.0
N_SEGMENTS = 128


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a: tuple[float, float, float], scale: float) -> tuple[float, float, float]:
    return (scale * a[0], scale * a[1], scale * a[2])


def _mutual_inductance_coaxial_loops(radius_m: float, separation_m: float) -> float:
    m = 4.0 * radius_m * radius_m / (4.0 * radius_m * radius_m + separation_m * separation_m)
    k = math.sqrt(m)
    return MU0 * radius_m * (
        (2.0 / k - k) * sp.ellipk(m) - (2.0 / k) * sp.ellipe(m)
    )


def _axisymmetric_force_reference_N() -> float:
    dz = 1.0e-6 * LOOP_SEPARATION_M
    dmdz = (
        _mutual_inductance_coaxial_loops(LOOP_RADIUS_M, LOOP_SEPARATION_M + dz)
        - _mutual_inductance_coaxial_loops(LOOP_RADIUS_M, LOOP_SEPARATION_M - dz)
    ) / (2.0 * dz)
    return CURRENT1_A * CURRENT2_A * dmdz


def _loop_segments(radius_m: float, z_m: float, n_segments: int) -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    dtheta = 2.0 * math.pi / n_segments
    segments = []
    for index in range(n_segments):
        theta = (index + 0.5) * dtheta
        center = (radius_m * math.cos(theta), radius_m * math.sin(theta), z_m)
        dl = (-radius_m * math.sin(theta) * dtheta, radius_m * math.cos(theta) * dtheta, 0.0)
        segments.append((center, dl))
    return segments


def _biot_savart_loop_force_3d_N(n_segments: int) -> tuple[float, float, float]:
    source = _loop_segments(LOOP_RADIUS_M, 0.0, n_segments)
    target = _loop_segments(LOOP_RADIUS_M, LOOP_SEPARATION_M, n_segments)
    coeff = MU0 * CURRENT1_A / (4.0 * math.pi)
    force = (0.0, 0.0, 0.0)
    for target_center, target_dl in target:
        field = (0.0, 0.0, 0.0)
        for source_center, source_dl in source:
            vector = (
                target_center[0] - source_center[0],
                target_center[1] - source_center[1],
                target_center[2] - source_center[2],
            )
            distance = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
            field = _add(field, _scale(_cross(source_dl, vector), coeff / distance**3))
        force = _add(force, _scale(_cross(target_dl, field), CURRENT2_A))
    return force


def build_summary(n_segments: int = N_SEGMENTS) -> dict[str, object]:
    started = time.perf_counter()
    axisymmetric_reference = _axisymmetric_force_reference_N()
    reference_done = time.perf_counter()
    force_3d = _biot_savart_loop_force_3d_N(n_segments)
    quadrature_done = time.perf_counter()
    full_gate = axisymmetric_to_3d_force_gate(
        axisymmetric_reference,
        force_3d,
        case_id="coaxial_loop_full_revolution",
        axial_axis="z",
        result_basis="full_revolution",
        axial_rtol=1.0e-8,
        transverse_rtol=1.0e-10,
        metadata={"producer": "3d_biot_savart_loop_quadrature", "n_segments": n_segments},
    )
    quarter_sector_force = [0.25 * value for value in force_3d]
    sector_gate = axisymmetric_to_3d_force_gate(
        axisymmetric_reference,
        quarter_sector_force,
        case_id="coaxial_loop_quarter_sector_axial_scaling",
        axial_axis="z",
        result_basis="symmetry_sector",
        sector_angle_deg=90.0,
        axial_rtol=1.0e-8,
        metadata={"producer": "scaled_from_full_3d_quadrature", "sector_angle_deg": 90.0},
    )
    plan = axisymmetric_to_3d_validation_plan("coaxial_loop_force")
    gate_done = time.perf_counter()
    checks = {
        "full_revolution_gate_ok": full_gate["status"] == "ok",
        "sector_axial_scaling_gate_ok": sector_gate["status"] == "ok",
        "transverse_components_checked_for_full_revolution": (
            full_gate["checks"]["transverse_components_cancel"] is True
        ),
        "sector_transverse_is_not_overclaimed": (
            sector_gate["checks"]["transverse_components_cancel"] == "not_checked_for_symmetry_sector"
        ),
        "passed": full_gate["status"] == "ok" and sector_gate["status"] == "ok",
    }
    assert checks["passed"]
    return {
        "kind": "axisymmetric_to_3d_force_gate_validation",
        "validation_class": True,
        "force_learning": (
            "a full 2*pi*r axisymmetric axial force can validate a 3D full-revolution "
            "force vector by axial agreement and transverse cancellation"
        ),
        "case_parameters": {
            "loop_radius_m": LOOP_RADIUS_M,
            "loop_separation_m": LOOP_SEPARATION_M,
            "current1_A": CURRENT1_A,
            "current2_A": CURRENT2_A,
            "n_segments": n_segments,
        },
        "checks": checks,
        "run_duration_s": round(gate_done - started, 6),
        "timing_breakdown_s": {
            "axisymmetric_reference": round(reference_done - started, 6),
            "three_d_biot_savart_quadrature": round(quadrature_done - reference_done, 6),
            "gate_and_plan_build": round(gate_done - quadrature_done, 6),
        },
        "plan": plan,
        "gates": [full_gate, sector_gate],
    }


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    parser.add_argument("--segments", type=int, default=N_SEGMENTS)
    args = parser.parse_args()
    if args.segments < 16:
        raise ValueError("--segments must be >= 16")

    summary = add_result_metadata(_json_clean(build_summary(args.segments)), __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    full_gate = summary["gates"][0]
    print("[axisymmetric -> 3D force gate]")
    print(f"  segments={args.segments}")
    print(f"  axisymmetric Fz={full_gate['axisymmetric_reference']['axial_force_N']:.12e} N")
    print(f"  3D F={full_gate['three_d_result']['full_revolution_force_vector_N']} N")
    print(f"  axial_rel_error={full_gate['errors']['axial_rel_error']:.3e}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
