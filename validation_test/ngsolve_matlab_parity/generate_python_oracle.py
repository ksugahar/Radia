"""Generate public-NGSolve Python results for the 100-case MATLAB oracle."""

from __future__ import annotations

import importlib.metadata
import json
import platform
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from netgen.geom2d import unit_square
from netgen.occ import unit_cube
from ngsolve import (
    BilinearForm,
    CF,
    H1,
    HCurl,
    HDiv,
    InnerProduct,
    Mesh,
    TaskManager,
    curl,
    div,
    dx,
    grad,
)
from scipy.io import savemat
from scipy.sparse import coo_matrix

from case_catalog import DOF_LIMIT, MESH_SPECS, build_case_catalog


SCHEMA = "radia.ngsolve-python-matlab-100-case-oracle.v1"


def _generate_mesh(spec: dict[str, Any], output_directory: Path) -> dict[str, Any]:
    generator = unit_square if spec["generator"] == "unit_square" else unit_cube
    path = output_directory / f"{spec['id']}.vol"
    netgen_mesh = generator.GenerateMesh(maxh=float(spec["maxh"]))
    netgen_mesh.Save(str(path))
    mesh = Mesh(str(path))
    return {
        **spec,
        "path": str(path.resolve()),
        "vertices": int(mesh.nv),
        "elements": int(mesh.ne),
    }


def _space(mesh: Mesh, spec: dict[str, Any]):
    options: dict[str, Any] = {"order": int(spec["order"])}
    if spec["dirichlet"]:
        options["dirichlet"] = spec["dirichlet"]
    if spec["space"] == "h1":
        return H1(mesh, **options)
    if spec["space"] == "hcurl":
        return HCurl(mesh, nograds=True, **options)
    if spec["space"] == "hdiv":
        return HDiv(mesh, **options)
    raise ValueError(f"unsupported space {spec['space']!r}")


def _assemble(mesh: Mesh, spec: dict[str, Any]):
    space = _space(mesh, spec)
    trial, test = space.TnT()
    weight = CF(float(spec["weight"]))
    form = BilinearForm(space, symmetric=True)
    if spec["space"] == "h1" and spec["form"] == "mass":
        form += weight * trial * test * dx
    elif spec["space"] == "h1" and spec["form"] == "stiffness":
        form += weight * InnerProduct(grad(trial), grad(test)) * dx
    elif spec["space"] == "hcurl" and spec["form"] == "mass":
        form += weight * InnerProduct(trial, test) * dx
    elif spec["space"] == "hcurl" and spec["form"] == "curlcurl":
        form += weight * InnerProduct(curl(trial), curl(test)) * dx
    elif spec["space"] == "hdiv" and spec["form"] == "mass":
        form += weight * InnerProduct(trial, test) * dx
    elif spec["space"] == "hdiv" and spec["form"] == "divdiv":
        form += weight * div(trial) * div(test) * dx
    else:
        raise ValueError(f"invalid operator specification: {spec}")
    with TaskManager():
        form.Assemble()
    return space, form


def _sparse_matrix(matrix) -> coo_matrix:
    rows, cols, values = matrix.COO()
    result = coo_matrix((values, (rows, cols)), shape=matrix.shape)
    result.eliminate_zeros()
    return result


def _values(size: int, case_number: int, phase: float) -> np.ndarray:
    indices = np.arange(1, size + 1, dtype=float)
    return (
        np.sin(indices * (0.071 + 0.0007 * case_number) + phase)
        + 0.25 * np.cos(indices * (0.037 + 0.0003 * case_number) - phase)
    )


def _vector(matrix, values: np.ndarray):
    result = matrix.CreateColVector()
    result.FV().NumPy()[:] = values
    return result


def _case_oracle(mesh: Mesh, spec: dict[str, Any], case_number: int):
    started = time.perf_counter()
    space, form = _assemble(mesh, spec)
    matrix = form.mat
    ndof = int(space.ndof)
    if ndof >= DOF_LIMIT:
        raise RuntimeError(
            f"{spec['case_id']} has {ndof} DoFs, exceeding the {DOF_LIMIT} limit"
        )
    free_dofs = space.FreeDofs()
    free_mask = np.fromiter(
        (bool(free_dofs[index]) for index in range(ndof)),
        dtype=bool,
        count=ndof,
    )
    input_values = _values(ndof, case_number, 0.2)
    input_vector = _vector(matrix, input_values)
    matvec = matrix.CreateColVector()
    with TaskManager():
        matvec.data = matrix * input_vector

    result: dict[str, Any] = {
        "matrix": _sparse_matrix(matrix),
        "input": input_values,
        "matvec": matvec.FV().NumPy().copy(),
        "free_dofs": free_mask,
        "ndof": float(ndof),
        "nnz": float(_sparse_matrix(matrix).nnz),
    }
    residual_norm = 0.0
    if spec["solve"]:
        rhs_values = _values(ndof, case_number, 0.7)
        rhs = _vector(matrix, rhs_values)
        solution = matrix.CreateColVector()
        with TaskManager():
            solution.data = matrix.Inverse(free_dofs) * rhs
        residual = matrix.CreateColVector()
        with TaskManager():
            residual.data = matrix * solution - rhs
        residual_values = residual.FV().NumPy()
        residual_norm = float(np.linalg.norm(residual_values[free_mask]))
        result["rhs"] = rhs_values
        result["solution"] = solution.FV().NumPy().copy()
    result["residual_norm"] = residual_norm
    result["duration_s"] = time.perf_counter() - started
    return result


def generate(output_directory: Path) -> tuple[Path, Path]:
    started = time.perf_counter()
    output_directory.mkdir(parents=True, exist_ok=True)
    mesh_rows = [_generate_mesh(spec, output_directory) for spec in MESH_SPECS]
    mesh_by_id = {row["id"]: row for row in mesh_rows}
    meshes = {key: Mesh(row["path"]) for key, row in mesh_by_id.items()}
    cases = build_case_catalog()
    mat_values: dict[str, Any] = {"fixture_schema": SCHEMA}
    actual_cases: list[dict[str, Any]] = []

    for case_number, case in enumerate(cases, start=1):
        oracle = _case_oracle(meshes[case["mesh_id"]], case, case_number)
        key = case["oracle_key"]
        for name in ("matrix", "input", "matvec", "free_dofs", "ndof", "nnz"):
            mat_values[f"{key}_{name}"] = oracle[name]
        if case["solve"]:
            mat_values[f"{key}_rhs"] = oracle["rhs"]
            mat_values[f"{key}_solution"] = oracle["solution"]
        actual_cases.append(
            {
                **case,
                "mesh_path": mesh_by_id[case["mesh_id"]]["path"],
                "python_ndof": int(oracle["ndof"]),
                "python_nnz": int(oracle["nnz"]),
                "python_residual_norm": oracle["residual_norm"],
                "python_duration_s": oracle["duration_s"],
            }
        )

    manifest = {
        "schema": SCHEMA,
        "generated_utc": datetime.now(UTC).isoformat(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": importlib.metadata.version("scipy"),
        "ngsolve_version": importlib.metadata.version("ngsolve"),
        "netgen_mesher_version": importlib.metadata.version("netgen-mesher"),
        "dof_limit": DOF_LIMIT,
        "case_count": len(actual_cases),
        "solve_case_count": sum(bool(row["solve"]) for row in actual_cases),
        "dimension_counts": {
            "two_d": sum(row["dimension"] == 2 for row in actual_cases),
            "three_d": sum(row["dimension"] == 3 for row in actual_cases),
        },
        "maximum_dofs": max(row["python_ndof"] for row in actual_cases),
        "meshes": mesh_rows,
        "cases": actual_cases,
        "python_total_duration_s": time.perf_counter() - started,
    }
    if manifest["case_count"] != 100 or manifest["dimension_counts"]["three_d"] == 0:
        raise RuntimeError("the oracle catalog must contain exactly 100 cases including 3D")

    mat_path = output_directory / "python_oracle.mat"
    manifest_path = output_directory / "python_oracle_manifest.json"
    savemat(mat_path, mat_values, do_compression=False, oned_as="column")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest_path, mat_path


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: generate_python_oracle.py OUTPUT_DIRECTORY")
    manifest, oracle = generate(Path(sys.argv[1]))
    print(json.dumps({"manifest": str(manifest), "oracle": str(oracle)}))


if __name__ == "__main__":
    main()

