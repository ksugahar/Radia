"""Production validation for mapped/non-affine HEX BDM2 HDiv-VIM.

This is intentionally a validation_test workload.  Its two dense Gram builds,
linear/nonlinear solves, and full/reduced IMA comparisons run on mdx or hibino,
not in the lightweight CI test suite.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from ngsolve.meshes import MakeStructured3DMesh

import radia
from radia.vim import FieldFromSolution, MagnetizationSource, Solve
from radia.vim import _vim as vim_core

HALF_WIDTH = 0.02
FIELD_POINTS = np.asarray(
    ((0.030, 0.005, 0.004), (0.040, -0.004, 0.003), (0.0, 0.035, 0.002)),
    dtype=float,
)


def _physical_map(x: float, y: float, z: float) -> tuple[float, float, float]:
    scale2 = HALF_WIDTH**2
    return (
        x * (1.0 + 0.20 * y * z / scale2),
        y + 0.08 * x * x * z / scale2,
        z + 0.06 * x * x * y / scale2,
    )


def _mesh(*, half: bool) -> ng.Mesh:
    x_min = 0.0 if half else -HALF_WIDTH
    nx = 1 if half else 2
    return MakeStructured3DMesh(
        hexes=True,
        nx=nx,
        ny=2,
        nz=2,
        mapping=lambda x, y, z: _physical_map(
            x_min + (HALF_WIDTH - x_min) * x,
            -HALF_WIDTH + 2.0 * HALF_WIDTH * y,
            -HALF_WIDTH + 2.0 * HALF_WIDTH * z,
        ),
    )


def _relative(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1.0e-300))


def _field_relative(a: np.ndarray, b: np.ndarray) -> float:
    denominator = np.maximum(np.linalg.norm(b, axis=1), 1.0e-300)
    return float(np.max(np.linalg.norm(a - b, axis=1) / denominator))


def _average_gridfunction(gf: ng.GridFunction, mesh: ng.Mesh) -> list[float]:
    volume = float(ng.Integrate(1.0, mesh))
    return [float(ng.Integrate(gf[c], mesh) / volume) for c in range(3)]


def _operator_sweep(mesh: ng.Mesh, rules: list[tuple[int, int]]) -> dict:
    fes = ng.HDiv(mesh, order=2)
    basis = vim_core._charge_basis_hex(fes, cob_quad=3)
    mass = sp.csr_matrix(basis["M_mass"]).toarray()
    excitation = ng.GridFunction(fes)
    excitation.Set(ng.CF((0.0, 0.0, 1000.0)))
    rhs = mass @ excitation.vec.FV().NumPy()
    operators: dict[tuple[int, int], np.ndarray] = {}
    rows = []
    for outer, inner in rules:
        started = time.perf_counter()
        _, gram, _, _ = vim_core._build_charge_gram_hex(
            fes,
            glout_n=outer,
            glin_n=inner,
            eps=1.0e-14,
            leafsize=4096,
        )
        _, gram, _ = vim_core._configure_cpp_operator(
            basis["B"], gram, basis["M_mass"], basis["M_mass_ngsolve"]
        )
        operator = np.asarray(gram.demag_matrix().ToDense(), dtype=float)
        operator = 0.5 * (operator + operator.T)
        operators[(outer, inner)] = operator
        eigenvalues = la.eigh(operator, mass, eigvals_only=True)
        coefficients = la.solve(mass / 99.0 + operator, rhs, assume_a="sym")
        excitation.vec.FV().NumPy()[:] = coefficients
        rows.append(
            {
                "outer_order": outer,
                "inner_order": inner,
                "build_wall_s": time.perf_counter() - started,
                "minimum_generalized_eigenvalue": float(eigenvalues.min()),
                "maximum_generalized_eigenvalue": float(eigenvalues.max()),
                "eigenvalues_outside_physical_interval": int(
                    np.count_nonzero(
                        (eigenvalues < -1.0e-10) | (eigenvalues > 1.0 + 5.0e-4)
                    )
                ),
                "material_average_magnetization": _average_gridfunction(
                    excitation, mesh
                ),
                "coefficients": coefficients,
                "stats": dict(gram.stats()),
            }
        )

    candidate = rows[0]
    reference = rows[-1]
    candidate_key = (candidate["outer_order"], candidate["inner_order"])
    reference_key = (reference["outer_order"], reference["inner_order"])
    difference = operators[candidate_key] - operators[reference_key]
    spectral_difference = la.eigh(difference, mass, eigvals_only=True)
    dc = candidate.pop("coefficients")
    dr = reference.pop("coefficients")
    for row in rows[1:-1]:
        row.pop("coefficients")
    coefficient_difference = dc - dr
    comparison = {
        "candidate": list(candidate_key),
        "reference": list(reference_key),
        "operator_frobenius_relative": _relative(
            operators[candidate_key], operators[reference_key]
        ),
        "operator_mass_spectral_difference": float(np.max(np.abs(spectral_difference))),
        "material_solution_euclidean_relative": _relative(dc, dr),
        "material_solution_mass_relative": float(
            np.sqrt(coefficient_difference @ mass @ coefficient_difference)
            / max(np.sqrt(dr @ mass @ dr), 1.0e-300)
        ),
        "material_average_relative": _relative(
            np.asarray(candidate["material_average_magnetization"]),
            np.asarray(reference["material_average_magnetization"]),
        ),
    }
    return {"hdiv_dofs": int(fes.ndof), "rules": rows, "comparison": comparison}


def _solve_summary(result: dict) -> dict:
    return {
        "converged": bool(result["last_solve_converged"]),
        "iterations": int(result["iters"]),
        "final_relative_residual": float(result["last_solve_final_relative_residual"]),
        "linear_solver": str(result["linear_solver"]),
        "M_avg": np.asarray(result["M_avg"], dtype=float).tolist(),
        "timing": dict(result.get("timing", {})),
    }


def _solve_and_ima() -> dict:
    full_mesh = _mesh(half=False)
    half_mesh = _mesh(half=True)
    solve_options = {
        "mu_r": 100.0,
        "H_ext": ng.CF((0.0, 0.0, 1000.0)),
        "order": 2,
        "gram_eps": 1.0e-14,
        "leaf": 4096,
        "tol": 1.0e-9,
        "maxit": 2000,
    }
    full = Solve(full_mesh, **solve_options)
    reduced = Solve(half_mesh, image="+x", **solve_options)
    full_field = np.asarray(FieldFromSolution(full, FIELD_POINTS, algorithm="direct"))
    reduced_field = np.asarray(
        FieldFromSolution(reduced, FIELD_POINTS, algorithm="direct")
    )

    mu0 = 4.0e-7 * np.pi
    linear_bh = np.asarray(
        ((0.0, 0.0), (1.0e3, mu0 * 100.0e3), (1.0e5, mu0 * 100.0e5)),
        dtype=float,
    )
    nonlinear = Solve(
        _mesh(half=False),
        bh_table=linear_bh,
        H_ext=ng.CF((0.0, 0.0, 1000.0)),
        order=2,
        gram_eps=1.0e-14,
        leaf=4096,
        tol=1.0e-9,
        nl_maxit=30,
    )
    return {
        "linear_full": _solve_summary(full),
        "linear_reduced": _solve_summary(reduced),
        "nonlinear_full": _solve_summary(nonlinear),
        "ima_solve_comparison": {
            "magnetization_relative": _relative(
                np.asarray(full["M_avg"], dtype=float),
                np.asarray(reduced["M_avg"], dtype=float),
            ),
            "field_relative": _field_relative(full_field, reduced_field),
            "field_maximum_absolute": float(np.max(np.abs(full_field - reduced_field))),
            "field_scale": float(np.max(np.abs(full_field))),
            "field_full": full_field.tolist(),
            "field_reduced": reduced_field.tolist(),
        },
    }


def _prescribed_source_roundoff() -> dict:
    magnetization = (0.0, 2.0e4, 1.0e5)
    full = MagnetizationSource(_mesh(half=False), magnetization, order=2)
    reduced = MagnetizationSource(_mesh(half=True), magnetization, order=2, image="+x")
    full_field = np.asarray(full.Field(FIELD_POINTS, algorithm="direct"))
    reduced_field = np.asarray(reduced.Field(FIELD_POINTS, algorithm="direct"))
    scale = float(np.max(np.abs(full_field)))
    absolute = float(np.max(np.abs(full_field - reduced_field)))
    return {
        "magnetization": list(magnetization),
        "full_stats": full.stats,
        "reduced_stats": reduced.stats,
        "field_full": full_field.tolist(),
        "field_reduced": reduced_field.tolist(),
        "field_maximum_absolute": absolute,
        "field_scale": scale,
        "field_relative": _field_relative(full_field, reduced_field),
        "field_error_in_machine_eps": float(
            absolute / max(scale * np.finfo(float).eps, 1.0e-300)
        ),
    }


def run(rules: list[tuple[int, int]]) -> dict:
    started = time.perf_counter()
    operator_mesh = _mesh(half=False)
    with ng.TaskManager():
        operator = _operator_sweep(operator_mesh, rules)
        solves = _solve_and_ima()
        prescribed = _prescribed_source_roundoff()
    comparison = operator["comparison"]
    checks = {
        "all_quadrature_rules_have_physical_spectrum": all(
            row["eigenvalues_outside_physical_interval"] == 0
            for row in operator["rules"]
        ),
        "default_quadrature_has_sub_per_mille_material_response_error": (
            comparison["material_solution_mass_relative"] <= 1.0e-3
        ),
        "linear_full_converged": solves["linear_full"]["converged"],
        "linear_reduced_converged": solves["linear_reduced"]["converged"],
        "nonlinear_energy_newton_converged": (
            solves["nonlinear_full"]["converged"]
            and solves["nonlinear_full"]["linear_solver"] == "energy-newton-cpp"
        ),
        "ima_independent_solves_match": (
            solves["ima_solve_comparison"]["magnetization_relative"] <= 1.0e-10
            and solves["ima_solve_comparison"]["field_relative"] <= 1.0e-10
        ),
        "ima_prescribed_source_field_matches_within_ten_eps": (
            prescribed["field_error_in_machine_eps"] <= 10.0
        ),
    }
    return {
        "schema": "radia.validation.mapped-hex-bdm2-production.v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "machine": platform.node(),
        "versions": {
            "radia": getattr(radia, "__version__", "unknown"),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "python": platform.python_version(),
        },
        "geometry": {
            "topology": "HEX",
            "hdiv_family": "BDM2",
            "geometry_order": int(operator_mesh.GetCurveOrder()),
            "full_cells": 8,
            "reduced_cells": 4,
            "map": "reflection-symmetric non-affine trilinear HEX map",
        },
        "operator_quadrature": operator,
        "solves": solves,
        "prescribed_source_roundoff": prescribed,
        "checks": checks,
        "pass": all(checks.values()),
        "timing_s": {"total": time.perf_counter() - started},
    }


def _parse_rule(value: str) -> tuple[int, int]:
    try:
        outer, inner = value.split(":", maxsplit=1)
        return int(outer), int(inner)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("quadrature rule must be OUTER:INNER") from exc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rule",
        type=_parse_rule,
        action="append",
        default=None,
        help="outer:inner tensor/Duffy order; repeat for a convergence pair",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("mapped_hex_bdm2_production_summary.json"),
    )
    args = parser.parse_args()
    rules = args.rule or [(9, 12), (10, 16)]
    if len(rules) < 2:
        parser.error("at least two --rule values are required")
    summary = run(rules)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "pass": summary["pass"],
                "checks": summary["checks"],
                "timing_s": summary["timing_s"],
            },
            indent=2,
        ),
        flush=True,
    )
    if not summary["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
