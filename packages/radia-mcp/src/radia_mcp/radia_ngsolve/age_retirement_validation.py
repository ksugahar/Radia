"""Self-contained no-remesh AGE retirement validation."""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from netgen.geom2d import SplineGeometry

from .age_periodic_motion import solve_age_periodic_motion


SCHEMA = "radia.age-retirement-validation.v1"


def _sha(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def _validation_vol() -> str:
    geometry = SplineGeometry()
    geometry.AddCircle((0, 0), 0.2, leftdomain=0, rightdomain=1, bc="rotor_inner")
    geometry.AddCircle((0, 0), 0.8, leftdomain=1, rightdomain=0, bc="rotor_ring")
    geometry.AddCircle((0, 0), 1.0, leftdomain=0, rightdomain=2, bc="stator_ring")
    geometry.AddCircle((0, 0), 1.4, leftdomain=2, rightdomain=0, bc="outer")
    geometry.SetMaterial(1, "rotor")
    geometry.SetMaterial(2, "stator")
    mesh = geometry.GenerateMesh(maxh=0.18)
    Path(r"C:\temp").mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="radia-age-", dir=r"C:\temp") as folder:
        path = Path(folder) / "generated_age_validation.vol"
        mesh.Save(str(path))
        return path.read_text(encoding="utf-8")


def run_age_retirement_validation(request: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    if str(request.get("method", "ngsolve_age_phase_only")) != "ngsolve_age_phase_only":
        raise ValueError("method must be ngsolve_age_phase_only")
    samples = int(request.get("angle_samples", 8))
    if samples < 4 or samples > 64:
        raise ValueError("angle_samples must be in [4, 64]")
    sector = {
        "slots": 12,
        "poles": 4,
        "sector_count": 4,
        "sector_angle_deg": 90.0,
        "boundary": str(request.get("boundary", "anti-periodic")),
        "boundary_phase": float(request.get("boundary_phase", -1.0)),
    }
    period = 2.0 * math.pi
    angles = [period * index / samples for index in range(samples)]
    solve_request = {
        "vol_text": _validation_vol(),
        "source_name": "generated_age_validation.vol",
        "airgap": {
            "inner_radius_m": 0.8,
            "outer_radius_m": 1.0,
            "rotor_ring": "rotor_ring",
            "stator_ring": "stator_ring",
            "rotor_inner": "rotor_inner",
            "outer": "outer",
            "rotor_material": "rotor",
            "stator_material": "stator",
            "harmonics": [1],
        },
        "materials": {
            "rotor": {"relative_permeability": 1.0, "conductivity_s_per_m": 0.0},
            "stator": {"relative_permeability": 1.0, "conductivity_s_per_m": 0.0},
        },
        "periodic_sector": sector,
        "excitation": {"1": {"rotor_amplitude": 1.0, "stator_amplitude": 0.8}},
        "rotor_angles_rad": angles,
        "axial_length_m": 0.05,
        "frequency_hz": 0.0,
        "element_order": 2,
    }
    solved = solve_age_periodic_motion(solve_request)
    summary = solved["torque_summary"]
    checks = {
        "generated_vol2d_solved": solved["status"] == "solved" and solved["mesh_contract"]["dimension"] == 2,
        "periodic_sector_is_anti_periodic": solved["periodicity_contract"]["boundary"] == "anti-periodic" and solved["periodicity_contract"]["boundary_phase"] == -1.0,
        "mesh_operator_and_factorization_reused": solved["mesh_reused_all_angles"] is True and solved["operator_reused_all_angles"] is True and solved["factorization_reused_all_angles"] is True,
        "rotation_is_harmonic_phase_only": solved["rotation_method"] == "rotor_harmonic_phase_only_no_remesh",
        "torque_closes_over_period": summary["closure_relative_error"] <= 1.0e-8,
        "torque_phase_sign_reverses": summary["phase_sign_reversal_observed"] is True,
        "execution_is_content_addressed": all(
            len(str(solved[key])) == 64
            for key in (
                "operator_sha256",
                "age_factorization_sha256",
                "angle_grid_sha256",
                "excitation_sha256",
                "torque_output_sha256",
            )
        ),
    }
    core = {
        "schema": SCHEMA,
        "request": {"method": "ngsolve_age_phase_only", "angle_samples": samples, **sector},
        "mesh_contract": solved["mesh_contract"],
        "periodicity_contract": solved["periodicity_contract"],
        "operator_sha256": solved["operator_sha256"],
        "age_factorization_sha256": solved["age_factorization_sha256"],
        "angle_grid_sha256": solved["angle_grid_sha256"],
        "excitation_sha256": solved["excitation_sha256"],
        "torque_output_sha256": solved["torque_output_sha256"],
        "torque_summary": summary,
        "checks": checks,
        "validated_capabilities": ["moving_band", "periodic_boundary"],
        "solver_launched": True,
        "source_solver_launched": False,
        "timing_seconds": solved["timing_breakdown_s"],
    }
    result_sha = _sha(core)
    expected = str(request.get("expected_result_sha256", "")).lower()
    checks["expected_result_identity_matches"] = not expected or expected == result_sha
    passed = all(checks.values())
    return {
        **core,
        "result_sha256": result_sha,
        "status": "verified" if passed else "needs_attention",
        "pass": passed,
        "issues": [name for name, ok in checks.items() if not ok],
    }
