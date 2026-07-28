"""Executable planar and axisymmetric harmonic magnetic validation."""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Mapping
from typing import Any

import ngsolve as ng
import numpy as np
from netgen.geom2d import unit_square


SCHEMA = "radia.harmonic-magnetic-validation.v1"


def _sha256(value: object) -> str:
    raw = json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def _csr_sha256(matrix: Any) -> str:
    digest = hashlib.sha256()
    for array in matrix.CSR():
        contiguous = np.ascontiguousarray(np.asarray(array))
        digest.update(str(contiguous.dtype).encode("ascii"))
        digest.update(np.asarray(contiguous.shape, dtype=np.int64).tobytes())
        digest.update(contiguous.tobytes())
    return digest.hexdigest()


def _positive(value: object, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return parsed


def _solve_case(
    formulation: str,
    *,
    frequency_hz: float,
    maxh: float,
    order: int,
    reluctivity: float,
    conductivity: float,
) -> dict[str, Any]:
    mesh = ng.Mesh(unit_square.GenerateMesh(maxh=maxh))
    fes = ng.H1(
        mesh,
        order=order,
        complex=True,
        dirichlet="left|right|top|bottom",
    )
    u, v = fes.TnT()
    x, y = ng.x, ng.y
    omega = 2.0 * math.pi * frequency_hz
    amplitude = 1.0 + 0.5j
    exact = amplitude * ng.sin(math.pi * x) * ng.sin(math.pi * y)
    exact_dx = amplitude * math.pi * ng.cos(math.pi * x) * ng.sin(math.pi * y)
    form = ng.BilinearForm(fes)
    if formulation == "planar_az":
        source = (
            2.0 * math.pi**2 * reluctivity * exact
            + 1j * omega * conductivity * exact
        )
        form += reluctivity * ng.grad(u) * ng.grad(v) * ng.dx
        form += 1j * omega * conductivity * u * v * ng.dx
        loss_weight = 1.0
        physical_domain = "unit_square_planar"
    elif formulation == "axisymmetric_psi_r_aphi":
        radius = 1.0 + x
        source = reluctivity * (
            2.0 * math.pi**2 * exact / radius + exact_dx / radius**2
        ) + 1j * omega * conductivity * exact / radius
        form += reluctivity / radius * ng.grad(u) * ng.grad(v) * ng.dx
        form += 1j * omega * conductivity / radius * u * v * ng.dx
        loss_weight = 1.0 / radius
        physical_domain = "annulus_cross_section_r_1_to_2_z_0_to_1"
    else:
        raise ValueError("formulation must be planar_az or axisymmetric_psi_r_aphi")
    rhs = ng.LinearForm(fes)
    rhs += source * v * ng.dx
    form.Assemble()
    rhs.Assemble()
    solution = ng.GridFunction(fes)
    solution.vec.data = form.mat.Inverse(
        fes.FreeDofs(), inverse="umfpack"
    ) * rhs.vec

    error_sq = float(ng.Integrate(ng.Norm(solution - exact) ** 2, mesh))
    reference_sq = float(ng.Integrate(ng.Norm(exact) ** 2, mesh))
    residual = rhs.vec.CreateVector()
    residual.data = rhs.vec - form.mat * solution.vec
    free = fes.FreeDofs()
    residual_inf = max(
        (abs(complex(residual[index])) for index in range(fes.ndof) if free[index]),
        default=0.0,
    )
    rhs_inf = max(
        (abs(complex(rhs.vec[index])) for index in range(fes.ndof) if free[index]),
        default=0.0,
    )
    loss = 0.5 * omega**2 * conductivity * float(
        ng.Integrate(ng.Norm(solution) ** 2 * loss_weight, mesh)
    )
    center = complex(solution(mesh(0.5, 0.5)))
    mesh_identity = {
        "vertices": mesh.nv,
        "volume_elements": mesh.ne,
        "maxh": maxh,
        "order": order,
        "complex_space": True,
        "physical_domain": physical_domain,
    }
    core = {
        "formulation": formulation,
        "frequency_hz": frequency_hz,
        "omega_rad_s": omega,
        "ndof": fes.ndof,
        "relative_l2_error": math.sqrt(error_sq / reference_sq),
        "relative_free_residual_inf": residual_inf / max(rhs_inf, 1.0),
        "ohmic_loss": loss,
        "center_solution": [center.real, center.imag],
        "mesh": mesh_identity,
        "mesh_sha256": _sha256(mesh_identity),
        "operator_sha256": _csr_sha256(form.mat),
    }
    return {**core, "case_sha256": _sha256(core)}


def run_harmonic_magnetic_validation(request: Mapping[str, Any]) -> dict[str, Any]:
    """Run dimension-distinct complex magnetic-diffusion manufactured solves."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    if str(request.get("method", "complex_magnetic_diffusion")).lower() != (
        "complex_magnetic_diffusion"
    ):
        raise ValueError("method must be complex_magnetic_diffusion")
    planar_frequency = _positive(
        request.get("planar_frequency_hz", 1.0), "planar_frequency_hz"
    )
    axisymmetric_frequency = _positive(
        request.get("axisymmetric_frequency_hz", 50000.0),
        "axisymmetric_frequency_hz",
    )
    maxh = _positive(request.get("maxh", 0.12), "maxh")
    order = int(request.get("order", 2))
    if order not in {1, 2, 3}:
        raise ValueError("order must be 1, 2, or 3")
    max_error = _positive(request.get("max_relative_l2_error", 0.002), "max_relative_l2_error")
    max_residual = _positive(
        request.get("max_relative_free_residual_inf", 1.0e-10),
        "max_relative_free_residual_inf",
    )
    normalized = {
        "method": "complex_magnetic_diffusion",
        "planar_frequency_hz": planar_frequency,
        "axisymmetric_frequency_hz": axisymmetric_frequency,
        "maxh": maxh,
        "order": order,
        "reluctivity": _positive(request.get("reluctivity", 2.3), "reluctivity"),
        "conductivity": _positive(request.get("conductivity", 4.0), "conductivity"),
        "max_relative_l2_error": max_error,
        "max_relative_free_residual_inf": max_residual,
    }
    started = time.perf_counter()
    common = {
        "maxh": maxh,
        "order": order,
        "reluctivity": normalized["reluctivity"],
        "conductivity": normalized["conductivity"],
    }
    planar = _solve_case(
        "planar_az", frequency_hz=planar_frequency, **common
    )
    planar_done = time.perf_counter()
    axisymmetric = _solve_case(
        "axisymmetric_psi_r_aphi",
        frequency_hz=axisymmetric_frequency,
        **common,
    )
    solved = time.perf_counter()
    rows = [planar, axisymmetric]
    checks = {
        "both_dimension_specific_complex_solves_executed": all(row["ndof"] > 0 for row in rows),
        "manufactured_solutions_match": all(row["relative_l2_error"] <= max_error for row in rows),
        "relative_free_residuals_are_bounded": all(
            row["relative_free_residual_inf"] <= max_residual for row in rows
        ),
        "ohmic_losses_are_positive_and_finite": all(
            math.isfinite(row["ohmic_loss"]) and row["ohmic_loss"] > 0.0 for row in rows
        ),
        "complex_phase_is_preserved": all(
            abs(row["center_solution"][1]) > 0.1 for row in rows
        ),
        "mesh_operator_and_case_are_content_addressed": all(
            all(len(row[key]) == 64 for key in ("mesh_sha256", "operator_sha256", "case_sha256"))
            for row in rows
        ),
    }
    core = {
        "schema": SCHEMA,
        "request": normalized,
        "request_sha256": _sha256(normalized),
        "cases": rows,
        "checks": checks,
        "validated_capabilities": ["axisymmetric_ac_magnetic", "planar_ac_magnetic"],
    }
    result_sha = _sha256(core)
    expected_sha = str(request.get("expected_result_sha256", "")).lower()
    checks["expected_result_identity_matches"] = not expected_sha or expected_sha == result_sha
    passed = all(checks.values())
    return {
        **core,
        "result_sha256": result_sha,
        "status": "verified" if passed else "needs_attention",
        "pass": passed,
        "issues": [name for name, ok in checks.items() if not ok],
        "solver_launched": True,
        "owned_worker_required": True,
        "timing_seconds": {
            "planar_solve": planar_done - started,
            "axisymmetric_solve": solved - planar_done,
            "validation_and_hash": time.perf_counter() - solved,
            "total": time.perf_counter() - started,
        },
    }
