"""Diagnose the native HDiv material-mass operator on an ESRF #3 response mesh.

This explicit investigation tool checks the native identity
``A_W x = W x + N x`` before the nonlinear response path is accepted.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import ngsolve as ng
import numpy as np

from radia.vim import DemagOperator


def _relative_norm(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(float(np.linalg.norm(expected)), 1.0))


def _mass_matrix(fes, coefficient):
    trial, test = fes.TnT()
    form = ng.BilinearForm(fes)
    form += coefficient * trial * test * ng.dx
    form.Assemble()
    return form.mat


def _apply(matrix, values: np.ndarray) -> np.ndarray:
    source = matrix.CreateRowVector()
    source.FV().NumPy()[:] = values
    result = matrix.CreateColVector()
    matrix.Mult(source, result)
    return np.ascontiguousarray(result.FV().NumPy(), dtype=np.float64)


def _solve_recovery(gram, *, inv_chi: float, rhs: np.ndarray,
                    expected: np.ndarray, tol: float = 1.0e-9,
                    maxit: int = 1000) -> dict[str, object]:
    try:
        solved = gram.solve_configured_linear_material_mass_riesz(
            float(inv_chi), np.ascontiguousarray(rhs, dtype=np.float64),
            tol=float(tol), maxit=int(maxit), symmetric=True,
        )
    except RuntimeError as error:
        return {"ok": False, "error": str(error)}
    recovered = np.asarray(solved["m"], dtype=float)
    timings = dict(solved.get("timings", {}))
    converged = bool(timings.get("last_solve_converged", 0.0) >= 0.5)
    return {
        "ok": True,
        "iterations": int(solved["iters"]),
        "converged_before_limit": converged and int(solved["iters"]) <= int(maxit),
        "maximum_iterations": int(maxit),
        "relative_residual_tolerance": float(tol),
        "recovery_relative_error": _relative_norm(recovered, expected),
        "timings": timings,
    }


def run(mesh_path: Path, *, order: int, gram_eps: float, seed: int,
        threads: int | None = None, exact_dense_memory_mb: int | None = None) -> dict[str, object]:
    if threads is not None:
        ng.SetNumThreads(int(threads))
    mesh = ng.Mesh(str(mesh_path))
    materials = tuple(mesh.GetMaterials())
    if materials != ("iron",):
        raise ValueError(f"expected iron-only HDiv response mesh, got {materials!r}")
    fes = ng.HDiv(mesh, order=int(order))
    alpha = 1.0 / 999.0
    rng = np.random.default_rng(seed)
    x = np.ascontiguousarray(rng.standard_normal(fes.ndof), dtype=np.float64)
    y = np.ascontiguousarray(rng.standard_normal(fes.ndof), dtype=np.float64)

    with ng.TaskManager():
        operator = DemagOperator(fes, eps=float(gram_eps), curve_order=2, leafsize=16, eta=2.0)
        geometry = _mass_matrix(fes, ng.CoefficientFunction(1.0))
        weighted = _mass_matrix(fes, ng.CoefficientFunction(alpha))
        geometry_x, geometry_y = _apply(geometry, x), _apply(geometry, y)
        weighted_x, weighted_y = _apply(weighted, x), _apply(weighted, y)
        native_geometry_x = np.asarray(operator._G.apply_configured_geometry_mass(x))
        native_geometry_y = np.asarray(operator._G.apply_configured_geometry_mass(y))
        charges_x = np.ascontiguousarray(operator._B @ x, dtype=np.float64)
        charges_y = np.ascontiguousarray(operator._B @ y, dtype=np.float64)
        gram_x = np.asarray(operator._G.matvec_sym(charges_x))
        gram_y = np.asarray(operator._G.matvec_sym(charges_y))
        demag_x = np.asarray(operator._G.apply_configured_demag(x, True))
        demag_y = np.asarray(operator._G.apply_configured_demag(y, True))
        operator._G.configure_mass_matrix_ngsolve(weighted)
        operator_x = np.asarray(operator._G.apply_configured_linear_material_operator(1.0, x))
        operator_y = np.asarray(operator._G.apply_configured_linear_material_operator(1.0, y))
        native_weighted_x = operator_x - demag_x
        native_weighted_y = operator_y - demag_y
        # This is the exact first PCG direction for x0=0.  If p^T A p is not
        # positive while M^{-1} itself recovers correctly, the defect is in
        # the compressed Gram's PSD contract rather than in material assembly
        # or PARDISO's mass factor.
        riesz_direction = np.asarray(operator._G.apply_configured_mass_riesz(operator_x))
        riesz_applied = np.asarray(
            operator._G.apply_configured_linear_material_operator(1.0, riesz_direction)
        )
        riesz_charges = np.ascontiguousarray(operator._B @ riesz_direction, dtype=np.float64)
        hmatrix_riesz_gram_quadratic = float(
            riesz_charges @ np.asarray(operator._G.matvec_sym(riesz_charges))
        )
        raw_riesz_gram_quadratic = float(
            operator._G.raw_symmetric_quadratic_form(riesz_charges)
        )
        hmatrix_material_recovery = _solve_recovery(
            operator._G, inv_chi=1.0, rhs=operator_x, expected=x)
        operator._G.restore_geometry_mass_matrix()
        hmatrix_geometry_operator_x = np.asarray(
            operator._G.apply_configured_linear_material_operator(alpha, x)
        )
        operator._G.configure_mass_matrix_ngsolve(weighted)
        exact_dense = None
        exact_dense_material_recovery = None
        if exact_dense_memory_mb is not None:
            exact_dense = dict(operator._G.build_exact_dense_normalized_gram(
                int(exact_dense_memory_mb)))
            dense_riesz_gram_quadratic = float(
                riesz_charges @ np.asarray(operator._G.matvec_sym(riesz_charges)))
            dense_operator_x = np.asarray(
                operator._G.apply_configured_linear_material_operator(1.0, x)
            )
            exact_dense_material_recovery = _solve_recovery(
                operator._G, inv_chi=1.0, rhs=dense_operator_x, expected=x)
        else:
            dense_riesz_gram_quadratic = None
        original_charge = operator._B.tocsr()
        empty_indptr = np.zeros(original_charge.shape[0] + 1, dtype=np.int32)
        operator._G.configure_charge_map(
            empty_indptr, np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.float64), int(fes.ndof),
        )
        weighted_preconditioner_recovery = _solve_recovery(
            operator._G, inv_chi=1.0, rhs=weighted_x, expected=x)
        operator._G.restore_geometry_mass_matrix()
        geometry_preconditioner_recovery = _solve_recovery(
            operator._G, inv_chi=alpha, rhs=alpha * geometry_x, expected=x)
        operator._G.configure_charge_map(
            np.ascontiguousarray(original_charge.indptr, dtype=np.int32),
            np.ascontiguousarray(original_charge.indices, dtype=np.int32),
            np.ascontiguousarray(original_charge.data, dtype=np.float64),
            int(fes.ndof),
        )
        geometry_operator_x = np.asarray(
            operator._G.apply_configured_linear_material_operator(alpha, x)
        )
        geometry_recovery = _solve_recovery(
            operator._G, inv_chi=alpha, rhs=geometry_operator_x, expected=x)
        hmatrix_stats = dict(operator._G.stats())

    report = {
        "schema": "radia.validation.hdiv-material-mass-diagnostic.v1",
        "mesh": str(mesh_path),
        "materials": list(materials),
        "elements": int(mesh.ne),
        "vertices": int(mesh.nv),
        "ndof": int(fes.ndof),
        "order": int(order),
        "threads": None if threads is None else int(threads),
        "uniform_alpha": alpha,
        "geometry_native_relative_error": _relative_norm(native_geometry_x, geometry_x),
        "weighted_ngsolve_scaling_relative_error": _relative_norm(weighted_x, alpha * geometry_x),
        "weighted_native_relative_error": _relative_norm(native_weighted_x, weighted_x),
        "weighted_native_scaling_relative_error": _relative_norm(native_weighted_x, alpha * native_geometry_x),
        "gram_symmetry_relative_error": abs(float(charges_x @ gram_y - charges_y @ gram_x)) / max(abs(float(charges_x @ gram_y)), abs(float(charges_y @ gram_x)), 1.0),
        "demag_native_sparse_composition_relative_error": _relative_norm(demag_x, operator._B.T @ gram_x),
        "operator_composition_relative_error": _relative_norm(operator_x, weighted_x + demag_x),
        "uniform_operator_path_relative_error": _relative_norm(
            hmatrix_geometry_operator_x, operator_x),
        "operator_symmetry_relative_error": abs(float(x @ operator_y - y @ operator_x)) / max(abs(float(x @ operator_y)), abs(float(y @ operator_x)), 1.0),
        "operator_quadratic_form": float(x @ operator_x),
        "first_mass_riesz_direction_quadratic_form": float(riesz_direction @ riesz_applied),
        "first_mass_riesz_direction_norm": float(np.linalg.norm(riesz_direction)),
        "first_mass_riesz_hmatrix_gram_quadratic_form": hmatrix_riesz_gram_quadratic,
        "first_mass_riesz_raw_gram_quadratic_form": raw_riesz_gram_quadratic,
        "first_mass_riesz_exact_dense_gram_quadratic_form": dense_riesz_gram_quadratic,
        "exact_dense_gram": exact_dense,
        "hmatrix_material_recovery": hmatrix_material_recovery,
        "exact_dense_material_recovery": exact_dense_material_recovery,
        "weighted_quadratic_form": float(x @ native_weighted_x),
        "demag_quadratic_form": float(x @ demag_x),
        "second_vector_weighted_native_relative_error": _relative_norm(native_weighted_y, weighted_y),
        "second_vector_geometry_native_relative_error": _relative_norm(native_geometry_y, geometry_y),
        "geometry_mass_riesz_recovery": geometry_recovery,
        "weighted_preconditioner_only_recovery": weighted_preconditioner_recovery,
        "geometry_preconditioner_only_recovery": geometry_preconditioner_recovery,
        "hmatrix": hmatrix_stats,
    }
    tolerance = 5.0e-12
    failed = {name: value for name, value in report.items() if name.endswith("relative_error") and float(value) > tolerance}
    if report["operator_quadratic_form"] <= 0.0:
        failed["operator_quadratic_form"] = report["operator_quadratic_form"]
    if report["weighted_quadratic_form"] <= 0.0:
        failed["weighted_quadratic_form"] = report["weighted_quadratic_form"]
    for name in (
        "geometry_mass_riesz_recovery",
        "weighted_preconditioner_only_recovery",
        "geometry_preconditioner_only_recovery",
    ):
        recovery = report[name]
        if not recovery["ok"]:
            failed[name] = recovery["error"]
        elif not recovery["converged_before_limit"]:
            failed[name] = "did not converge before iteration limit"
        elif recovery["recovery_relative_error"] > 3.0e-7:
            failed[name] = recovery["recovery_relative_error"]
    if exact_dense_memory_mb is None:
        recovery = report["hmatrix_material_recovery"]
        if not recovery["ok"]:
            failed["hmatrix_material_recovery"] = recovery["error"]
        elif recovery["recovery_relative_error"] > 3.0e-7:
            failed["hmatrix_material_recovery"] = recovery["recovery_relative_error"]
        if report["first_mass_riesz_direction_quadratic_form"] <= 0.0:
            failed["first_mass_riesz_direction_quadratic_form"] = (
                report["first_mass_riesz_direction_quadratic_form"])
    else:
        recovery = report["exact_dense_material_recovery"]
        if not recovery or not recovery["ok"]:
            failed["exact_dense_material_recovery"] = (
                None if recovery is None else recovery["error"])
        elif not recovery["converged_before_limit"]:
            failed["exact_dense_material_recovery"] = "did not converge before iteration limit"
        elif recovery["recovery_relative_error"] > 3.0e-7:
            failed["exact_dense_material_recovery"] = recovery["recovery_relative_error"]
        dense_q = report["first_mass_riesz_exact_dense_gram_quadratic_form"]
        raw_q = report["first_mass_riesz_raw_gram_quadratic_form"]
        if dense_q is None or abs(dense_q - raw_q) > 1.0e-8 * max(abs(raw_q), 1.0):
            failed["exact_dense_raw_gram_agreement"] = (dense_q, raw_q)
    if failed:
        raise RuntimeError(f"native material-mass contract failed: {failed}; report={report}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--gram-eps", type=float, default=1.0e-7)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--exact-dense-memory-mb", type=int)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    encoded = json.dumps(
        run(options.mesh, order=options.order, gram_eps=options.gram_eps,
            seed=options.seed, threads=options.threads,
            exact_dense_memory_mb=options.exact_dense_memory_mb),
        indent=2,
    ) + "\n"
    if options.output is None:
        print(encoded, end="")
    else:
        options.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()
