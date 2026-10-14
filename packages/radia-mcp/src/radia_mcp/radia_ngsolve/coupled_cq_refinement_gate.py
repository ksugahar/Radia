"""Solver-neutral refinement gate for coupled FEM/BEM convolution quadrature."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from .cq_urn import bdf_delta


def _numpy():
    import numpy as np

    return np


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (positive and parsed <= 0.0):
        requirement = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {requirement}")
    return parsed


def _array(value: object, name: str, count: int):
    np = _numpy()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = np.asarray(value, dtype=float).ravel()
    if result.size != count or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain {count} finite values")
    return result


def coupled_cq_refinement_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Validate CQ symbols, contour balance, refinement trend, and replay."""
    np = _numpy()
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    contract = summary.get("model_contract")
    if not isinstance(contract, Mapping):
        raise ValueError("model_contract must be an object")
    runs = summary.get("runs")
    if (
        not isinstance(runs, Sequence)
        or isinstance(runs, (str, bytes))
        or len(runs) != 5
        or any(not isinstance(run, Mapping) for run in runs)
    ):
        raise ValueError("runs must contain exactly five objects")
    comparison = summary.get("comparison")
    if not isinstance(comparison, Mapping):
        raise ValueError("comparison must be an object")

    expected_roles = {
        "bdf1_coarse",
        "bdf2_coarse",
        "bdf2_medium",
        "bdf2_fine",
        "bdf2_medium_replay",
    }
    roles = set()
    symbol_errors = []
    contour_errors = []
    min_real_nodes = []
    residuals = []
    imaginary_ratios = []
    mesh_rows = []
    for index, run in enumerate(runs):
        prefix = f"runs[{index}]"
        role = str(run.get("name") or "").strip()
        method = str(run.get("method") or "").strip().lower()
        if "bdf1" in method:
            bdf_method = "bdf1"
        elif "bdf2" in method:
            bdf_method = "bdf2"
        else:
            raise ValueError(f"{prefix}.method must identify BDF1 or BDF2")
        count = int(_finite(run.get("num_time"), f"{prefix}.num_time", positive=True))
        step = _finite(run.get("time_step_s"), f"{prefix}.time_step_s", positive=True)
        radius = _finite(run.get("cq_radius"), f"{prefix}.cq_radius", positive=True)
        if radius >= 1.0:
            raise ValueError(f"{prefix}.cq_radius must be below one")
        real = _array(run.get("cq_laplace_parameter_real"), f"{prefix}.cq_laplace_parameter_real", count)
        imag = _array(run.get("cq_laplace_parameter_imag"), f"{prefix}.cq_laplace_parameter_imag", count)
        actual = real + 1j * imag
        zeta = radius * np.exp(-2j * np.pi * np.arange(count) / count)
        expected = bdf_delta(zeta, bdf_method) / step
        symbol_error = float(
            np.max(np.abs(actual - expected)) / max(float(np.max(np.abs(expected))), 1.0e-300)
        )
        symbol_errors.append(symbol_error)
        contour_errors.append(abs(radius**count - 1.0e-4))
        min_real_nodes.append(float(np.min(real)))
        residuals.append(_finite(run.get("max_relative_residual"), f"{prefix}.max_relative_residual"))
        interior_scale = _finite(
            run.get("max_abs_interior_pressure"),
            f"{prefix}.max_abs_interior_pressure",
            positive=True,
        )
        exterior_scale = _finite(
            run.get("max_abs_exterior_pressure"),
            f"{prefix}.max_abs_exterior_pressure",
            positive=True,
        )
        imaginary_ratios.extend([
            _finite(
                run.get("max_imag_interior_before_real"),
                f"{prefix}.max_imag_interior_before_real",
            )
            / max(interior_scale, 1.0e-300),
            _finite(
                run.get("max_imag_exterior_before_real"),
                f"{prefix}.max_imag_exterior_before_real",
            )
            / max(exterior_scale, 1.0e-300),
        ])
        mesh_rows.append((
            int(run.get("num_volume_nodes", 0)),
            int(run.get("num_boundary_nodes", 0)),
            int(run.get("num_interior_only_nodes", 0)),
        ))
        roles.add(role)

    bdf1_error = _finite(
        comparison.get("bdf1_coarse_to_bdf2_fine_relative_error"),
        "comparison.bdf1_coarse_to_bdf2_fine_relative_error",
        positive=True,
    )
    bdf2_coarse_error = _finite(
        comparison.get("bdf2_coarse_to_fine_relative_error"),
        "comparison.bdf2_coarse_to_fine_relative_error",
        positive=True,
    )
    bdf2_medium_error = _finite(
        comparison.get("bdf2_medium_to_fine_relative_error"),
        "comparison.bdf2_medium_to_fine_relative_error",
        positive=True,
    )
    refinement_ratio = _finite(
        comparison.get("bdf2_refinement_error_ratio"),
        "comparison.bdf2_refinement_error_ratio",
        positive=True,
    )
    replay_error = _finite(
        comparison.get("bdf2_medium_replay_relative_error"),
        "comparison.bdf2_medium_replay_relative_error",
    )
    initial_alias = _finite(
        comparison.get("initial_exterior_relative_amplitude"),
        "comparison.initial_exterior_relative_amplitude",
    )
    alias_target = _finite(
        comparison.get("contour_alias_target"),
        "comparison.contour_alias_target",
        positive=True,
    )
    trusted_fraction = _finite(
        comparison.get("trusted_window_fraction"),
        "comparison.trusted_window_fraction",
        positive=True,
    )

    checks = {
        "p1_johnson_nedelec_cq_contract": contract
        == {
            "volume_element": "P1 tetrahedron",
            "boundary_element": "P1 triangle",
            "coupling": "JohnsonNedelecCalderon",
            "time_method_family": "Lubich CQ",
            "absorbing_boundary": "none",
        },
        "five_required_refinement_roles": roles == expected_roles,
        "same_nontrivial_first_order_mesh": len(set(mesh_rows)) == 1
        and mesh_rows[0][0] > mesh_rows[0][1] > 0
        and mesh_rows[0][2] > 0,
        "bdf_symbols_match_saved_laplace_nodes": max(symbol_errors) <= 1.0e-14,
        "cq_laplace_nodes_are_in_open_right_half_plane": min(min_real_nodes) > 0.0,
        "contour_alias_target_is_balanced": alias_target == 1.0e-4
        and max(contour_errors) <= 1.0e-12 * alias_target,
        "coupled_solve_residuals_are_small": max(residuals) <= 1.0e-10,
        "inverse_cq_reconstruction_is_real": max(imaginary_ratios) <= 1.0e-8,
        "trusted_window_is_declared": trusted_fraction == 0.75
        and _finite(comparison.get("trusted_window_end_s"), "comparison.trusted_window_end_s", positive=True)
        > 0.0,
        "bdf2_refinement_reduces_trusted_window_error": bdf2_medium_error
        < 0.9 * bdf2_coarse_error
        and refinement_ratio < 0.9,
        "bdf2_coarse_improves_on_bdf1_control": bdf2_coarse_error < bdf1_error,
        "medium_replay_is_deterministic": replay_error <= 1.0e-12,
        "initial_alias_is_bounded_by_contour_target": initial_alias <= 10.0 * alias_target,
    }
    return {
        "policy": "coupled_cq_refinement_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "maximum_bdf_symbol_relative_error": max(symbol_errors),
            "minimum_real_laplace_parameter": min(min_real_nodes),
            "maximum_coupled_residual": max(residuals),
            "maximum_imaginary_reconstruction_ratio": max(imaginary_ratios),
            "bdf1_coarse_error": bdf1_error,
            "bdf2_coarse_error": bdf2_coarse_error,
            "bdf2_medium_error": bdf2_medium_error,
            "bdf2_refinement_error_ratio": refinement_ratio,
            "bdf2_medium_replay_error": replay_error,
            "initial_alias_relative_amplitude": initial_alias,
        },
        "lesson": (
            "CQ validation must bind the BDF symbol and contour radius before "
            "interpreting a time trace. Compare refinement on a declared trusted "
            "window; the rho^-n inverse transform can amplify round-off near the "
            "end of the DFT window even when every Laplace solve has a tiny residual."
        ),
    }
