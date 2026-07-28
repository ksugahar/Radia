"""Executable Kelvin-transform validation for static open boundaries."""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Mapping, Sequence
from typing import Any

from .fem_bem_coupling import (
    kelvin_dtn_eigenvalue,
    kelvin_twosphere_shell_dipole,
)


SCHEMA = "radia.kelvin-open-boundary-validation.v1"


def _sha256(value: object) -> str:
    raw = json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def _finite_positive(value: object, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return parsed


def run_kelvin_open_boundary_validation(request: Mapping[str, Any]) -> dict[str, Any]:
    """Run 2-D mode and genuine 3-D two-sphere Kelvin validations."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    if str(request.get("method", "kelvin_transform")).lower() != "kelvin_transform":
        raise ValueError("magnetostatic retirement evidence requires kelvin_transform")
    if bool(request.get("pml", False)):
        raise ValueError("PML is not accepted for the static open-boundary policy")
    if str(request.get("physics_regime", "magnetostatic")).lower() not in {
        "magnetostatic",
        "static_poisson",
    }:
        raise ValueError("Kelvin validation is static-only and cannot validate wave radiation")
    if str(request.get("wave_boundary_inference", "forbidden")).lower() != "forbidden":
        raise ValueError("wave-boundary inference from static Kelvin evidence is forbidden")

    degrees_raw = request.get("degrees_2d", [1, 2, 3])
    if not isinstance(degrees_raw, Sequence) or isinstance(degrees_raw, (str, bytes)):
        raise ValueError("degrees_2d must be a list")
    degrees = [int(value) for value in degrees_raw]
    if not degrees or len(degrees) != len(set(degrees)) or any(value < 1 or value > 8 for value in degrees):
        raise ValueError("degrees_2d must contain unique integers in [1, 8]")
    order_2d = int(request.get("order_2d", max(degrees)))
    if order_2d < max(degrees):
        raise ValueError("order_2d must be at least the maximum requested degree")
    maxh_2d = _finite_positive(request.get("maxh_2d", 0.35), "maxh_2d")
    max_mode_error = _finite_positive(
        request.get("max_mode_relative_error", 0.005), "max_mode_relative_error"
    )
    max_3d_error = _finite_positive(
        request.get("max_three_dimensional_relative_error", 0.025),
        "max_three_dimensional_relative_error",
    )
    max_residual = _finite_positive(
        request.get("max_free_residual_inf", 1.0e-10), "max_free_residual_inf"
    )

    normalized_request = {
        "method": "kelvin_transform",
        "physics_regime": str(request.get("physics_regime", "magnetostatic")).lower(),
        "wave_boundary_inference": "forbidden",
        "degrees_2d": degrees,
        "order_2d": order_2d,
        "maxh_2d": maxh_2d,
        "max_mode_relative_error": max_mode_error,
        "three_dimensional": {
            "inner_radius": _finite_positive(request.get("inner_radius", 0.5), "inner_radius"),
            "outer_radius": _finite_positive(request.get("outer_radius", 1.0), "outer_radius"),
            "offset": _finite_positive(request.get("offset", 3.0), "offset"),
            "maxh": _finite_positive(request.get("maxh_3d", 0.18), "maxh_3d"),
            "order": int(request.get("order_3d", 2)),
            "curve_order": int(request.get("curve_order_3d", 3)),
        },
        "max_three_dimensional_relative_error": max_3d_error,
        "max_free_residual_inf": max_residual,
    }
    three = normalized_request["three_dimensional"]
    if three["inner_radius"] >= three["outer_radius"]:
        raise ValueError("inner_radius must be smaller than outer_radius")
    if three["offset"] <= 2.0 * three["outer_radius"]:
        raise ValueError("offset must keep the physical and Kelvin spheres disjoint")
    if three["order"] < 1 or three["curve_order"] < 1:
        raise ValueError("3-D field and curve orders must be positive")

    started = time.perf_counter()
    mode_rows = [
        kelvin_dtn_eigenvalue(
            R=three["outer_radius"],
            degree=degree,
            maxh=maxh_2d,
            order=order_2d,
            dim=2,
        )
        for degree in degrees
    ]
    modes_done = time.perf_counter()
    three_dimensional = kelvin_twosphere_shell_dipole(
        R_inner=three["inner_radius"],
        R_outer=three["outer_radius"],
        offset=three["offset"],
        maxh=three["maxh"],
        order=three["order"],
        curve_order=three["curve_order"],
    )
    solved = time.perf_counter()
    finite_modes = all(
        math.isfinite(float(row["rel_err"]))
        and math.isfinite(float(row["lam"]))
        and math.isfinite(float(row["lam_exact"]))
        for row in mode_rows
    )
    checks = {
        "static_kelvin_policy_without_pml": True,
        "two_dimensional_modes_are_finite": finite_modes,
        "two_dimensional_modes_match_exact_dtn": finite_modes
        and max(float(row["rel_err"]) for row in mode_rows) <= max_mode_error,
        "three_dimensional_two_sphere_periodic_solve": three_dimensional["ndof"] > 0
        and three_dimensional["mesh"]["volume_elements"] > 0,
        "three_dimensional_dipole_matches_exact_exterior": math.isfinite(
            float(three_dimensional["rel_err"])
        )
        and float(three_dimensional["rel_err"]) <= max_3d_error,
        "three_dimensional_free_residual_is_bounded": float(
            three_dimensional["free_residual_inf"]
        )
        <= max_residual,
        "mesh_and_operator_are_content_addressed": all(
            len(str(three_dimensional[name])) == 64
            for name in ("mesh_sha256", "operator_sha256")
        ),
        "periodic_map_direction_is_explicit": three_dimensional["mesh"]["periodic_map"]
        == "kelvin_ext_to_kelvin_int_inverse_translation",
    }
    result_core = {
        "schema": SCHEMA,
        "request_sha256": _sha256(normalized_request),
        "request": normalized_request,
        "mode_rows_2d": mode_rows,
        "three_dimensional": three_dimensional,
        "checks": checks,
        "validated_capabilities": ["kelvin_open_boundary", "open_boundary"],
        "wave_boundary_policy": "not_applicable_do_not_infer_from_static_kelvin",
    }
    result_sha = _sha256(result_core)
    expected_sha = str(request.get("expected_result_sha256", "")).lower()
    identity_matches = not expected_sha or expected_sha == result_sha
    checks["expected_result_identity_matches"] = identity_matches
    passed = all(checks.values())
    return {
        **result_core,
        "result_sha256": result_sha,
        "status": "verified" if passed else "needs_attention",
        "pass": passed,
        "issues": [name for name, ok in checks.items() if not ok],
        "solver_launched": True,
        "owned_worker_required": True,
        "timing_seconds": {
            "two_dimensional_modes": modes_done - started,
            "three_dimensional_solve": solved - modes_done,
            "validation_and_hash": time.perf_counter() - solved,
            "total": time.perf_counter() - started,
        },
    }
