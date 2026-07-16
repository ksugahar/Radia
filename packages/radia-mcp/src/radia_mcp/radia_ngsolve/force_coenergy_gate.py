"""Virtual-work consistency gate for displacement-force sweeps."""
from __future__ import annotations

import math


def force_coenergy_displacement_gate(
    positions_m,
    coenergy_j,
    forces_along_displacement_n,
    *,
    energy_kind: str = "constant_current_coenergy",
    max_central_relative_error: float = 0.02,
    min_sample_count: int = 5,
    artifact_identity: dict | None = None,
):
    """Compare direct force with the central derivative of magnetic coenergy.

    The caller must project the direct force onto the increasing displacement
    coordinate before calling this gate.  Endpoints are reported using one-sided
    differences but are not part of the acceptance metric.
    """
    x = [float(value) for value in positions_m]
    w = [float(value) for value in coenergy_j]
    force = [float(value) for value in forces_along_displacement_n]
    if not (len(x) == len(w) == len(force)):
        raise ValueError("positions, coenergy, and force must have the same length")
    if min_sample_count < 5:
        raise ValueError("min_sample_count must be >= 5")
    if max_central_relative_error < 0.0:
        raise ValueError("max_central_relative_error must be >= 0")

    identity_present = isinstance(artifact_identity, dict)
    force_snapshot_ok = True
    mesh_family_ok = True
    if artifact_identity is not None and not identity_present:
        force_snapshot_ok = False
        mesh_family_ok = False
    elif identity_present:
        direct = artifact_identity.get("direct_force_snapshot")
        derivative = artifact_identity.get("coenergy_derivative_snapshot")
        if not isinstance(direct, dict) or not isinstance(derivative, dict):
            force_snapshot_ok = False
        else:
            direct_step = str(direct.get("load_step_id", ""))
            derivative_step = str(derivative.get("load_step_id", ""))
            try:
                direct_time = float(direct["time_s"])
                derivative_time = float(derivative["time_s"])
            except (KeyError, TypeError, ValueError):
                direct_time = math.nan
                derivative_time = math.nan
            force_snapshot_ok = (
                bool(direct_step)
                and direct_step == derivative_step
                and math.isfinite(direct_time)
                and math.isfinite(derivative_time)
                and direct_time == derivative_time
            )
        generations = artifact_identity.get("coenergy_mesh_family_generations")
        mesh_family_ok = (
            isinstance(generations, list)
            and len(generations) == len(x)
            and all(isinstance(value, str) and bool(value) for value in generations)
            and len(set(generations)) == 1
        )

    finite = all(math.isfinite(value) for value in x + w + force)
    increasing = finite and all(right > left for left, right in zip(x, x[1:]))
    rows = []
    central_errors = []
    if finite and increasing and len(x) >= 2:
        for index in range(len(x)):
            if index == 0:
                derivative = (w[1] - w[0]) / (x[1] - x[0])
                stencil = "forward"
            elif index == len(x) - 1:
                derivative = (w[-1] - w[-2]) / (x[-1] - x[-2])
                stencil = "backward"
            else:
                derivative = (w[index + 1] - w[index - 1]) / (x[index + 1] - x[index - 1])
                stencil = "central"
            scale = max(abs(force[index]), abs(derivative), 1.0e-30)
            relative_error = abs(derivative - force[index]) / scale
            if stencil == "central":
                central_errors.append(relative_error)
            rows.append(
                {
                    "index": index,
                    "position_m": x[index],
                    "direct_force_N": force[index],
                    "coenergy_derivative_force_N": derivative,
                    "stencil": stencil,
                    "relative_error": relative_error,
                }
            )

    max_error = max(central_errors) if central_errors else math.inf
    checks = {
        "sample_count_sufficient": len(x) >= min_sample_count,
        "all_finite": finite,
        "positions_strictly_increase": increasing,
        "constant_current_coenergy_recorded": energy_kind == "constant_current_coenergy",
        "coenergy_nontrivial": finite and bool(w) and max(w) > min(w),
        "central_rows_available": len(central_errors) >= 3,
        "central_virtual_work_matches_direct_force": max_error <= max_central_relative_error,
        "force_and_coenergy_share_load_step_snapshot": force_snapshot_ok,
        "coenergy_stencil_uses_one_mesh_family_generation": mesh_family_ok,
    }
    return {
        "policy": "force_coenergy_displacement_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(x),
        "central_sample_count": len(central_errors),
        "max_central_relative_error": max_error,
        "mean_central_relative_error": (
            sum(central_errors) / len(central_errors) if central_errors else None
        ),
        "endpoint_errors_are_diagnostic_only": True,
        "checks": checks,
        "warnings": [] if identity_present else ["artifact_identity_not_recorded"],
        "rows": rows,
        "lesson": (
            "At fixed current, direct force projected onto the displacement axis "
            "must match dW'/dx. Use central differences for the acceptance gate; "
            "one-sided endpoint errors are diagnostics only."
        ),
    }
