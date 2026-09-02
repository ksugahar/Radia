"""Generate the 500-case breadth and 20-case scale Python oracles."""

from __future__ import annotations

import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import ngsolve as ng
import numpy as np
from scipy.io import savemat
from scipy.sparse import coo_matrix

from extended_catalog import (
    DOF_LIMIT,
    SCALE_DOF_FLOOR,
    build_breadth_catalog,
    build_scale_catalog,
    SCALE_MESH_IDS,
    SMALL_MESH_IDS,
)


SCHEMA = "radia.ngsolve-python-matlab-extended-oracle.v1"


def _space(mesh: ng.Mesh, spec: dict[str, Any]):
    options: dict[str, Any] = {
        "order": int(spec["order"]),
        "complex": bool(spec["weight_imag"]),
    }
    if spec["dirichlet"]:
        options["dirichlet"] = spec["dirichlet"]
    if spec["space"] == "h1":
        return ng.H1(mesh, **options)
    if spec["space"] == "hcurl":
        return ng.HCurl(mesh, nograds=True, **options)
    if spec["space"] == "hdiv":
        return ng.HDiv(mesh, **options)
    raise ValueError(f"unsupported space {spec['space']!r}")


def _coefficient(spec: dict[str, Any]):
    weight = complex(spec["weight_real"], spec["weight_imag"])
    if weight.imag == 0.0:
        weight = weight.real
    if spec["coefficient_kind"] == "constant":
        return ng.CF(weight)
    if spec["coefficient_kind"] == "affine_x":
        return weight + 0.2 * ng.x
    raise ValueError(f"unsupported coefficient kind {spec['coefficient_kind']!r}")


def _assemble(mesh: ng.Mesh, spec: dict[str, Any]):
    space = _space(mesh, spec)
    trial, test = space.TnT()
    coefficient = _coefficient(spec)
    form = ng.BilinearForm(space, symmetric=True)
    if spec["boundary"]:
        form += coefficient * trial * test * ng.ds
    elif spec["space"] == "h1" and spec["form"] == "mass":
        form += coefficient * trial * test * ng.dx
    elif spec["space"] == "h1" and spec["form"] == "stiffness":
        form += coefficient * ng.InnerProduct(ng.grad(trial), ng.grad(test)) * ng.dx
    elif spec["space"] == "hcurl" and spec["form"] == "mass":
        form += coefficient * ng.InnerProduct(trial, test) * ng.dx
    elif spec["space"] == "hcurl" and spec["form"] == "curlcurl":
        form += coefficient * ng.InnerProduct(ng.curl(trial), ng.curl(test)) * ng.dx
    elif spec["space"] == "hdiv" and spec["form"] == "mass":
        form += coefficient * ng.InnerProduct(trial, test) * ng.dx
    elif spec["space"] == "hdiv" and spec["form"] == "divdiv":
        form += coefficient * ng.div(trial) * ng.div(test) * ng.dx
    else:
        raise ValueError(f"invalid operator specification: {spec}")
    with ng.TaskManager():
        form.Assemble()
    return space, form


def _sparse_matrix(matrix) -> coo_matrix:
    rows, cols, values = matrix.COO()
    result = coo_matrix((values, (rows, cols)), shape=matrix.shape)
    result.eliminate_zeros()
    return result


def _values(size: int, case_number: int, phase: float, is_complex=False):
    indices = np.arange(1, size + 1, dtype=float)
    real = (
        np.sin(indices * (0.047 + 0.00011 * case_number) + phase)
        + 0.2 * np.cos(indices * (0.031 + 0.00007 * case_number) - phase)
    )
    if not is_complex:
        return real
    imaginary = 0.3 * np.cos(indices * 0.019 + phase)
    return real + 1j * imaginary


def _vector(matrix, values):
    result = matrix.CreateColVector()
    result.FV().NumPy()[:] = values
    return result


def _free_mask(space, ndof: int) -> np.ndarray:
    free = space.FreeDofs()
    return np.fromiter((bool(free[index]) for index in range(ndof)),
                       dtype=bool, count=ndof)


def _breadth_oracle(mesh: ng.Mesh, spec: dict[str, Any], number: int):
    started = time.perf_counter()
    space, form = _assemble(mesh, spec)
    matrix = form.mat
    sparse = _sparse_matrix(matrix)
    ndof = int(space.ndof)
    if ndof >= DOF_LIMIT:
        raise RuntimeError(f"{spec['case_id']} exceeds the DoF limit")
    free = _free_mask(space, ndof)
    is_complex = bool(spec["weight_imag"])
    input_values = _values(ndof, number, 0.2, is_complex)
    input_vector = _vector(matrix, input_values)
    output = matrix.CreateColVector()
    with ng.TaskManager():
        output.data = matrix * input_vector
    result: dict[str, Any] = {
        "matrix": sparse,
        "input": input_values,
        "matvec": output.FV().NumPy().copy(),
        "free_dofs": free,
        "ndof": float(ndof),
        "nnz": float(sparse.nnz),
    }
    solve_performed = bool(spec["solve"] and np.any(free))
    if solve_performed:
        rhs_values = _values(ndof, number, 0.7)
        rhs = _vector(matrix, rhs_values)
        solution = matrix.CreateColVector()
        with ng.TaskManager():
            solution.data = matrix.Inverse(space.FreeDofs()) * rhs
        result["rhs"] = rhs_values
        result["solution"] = solution.FV().NumPy().copy()
    result["solve_performed"] = solve_performed
    result["duration_s"] = time.perf_counter() - started
    return result


def _scale_oracle(mesh: ng.Mesh, spec: dict[str, Any], number: int):
    started = time.perf_counter()
    space, form = _assemble(mesh, spec)
    matrix = form.mat
    ndof = int(space.ndof)
    if not SCALE_DOF_FLOOR <= ndof < DOF_LIMIT:
        raise RuntimeError(
            f"{spec['case_id']} has {ndof} DoFs; expected "
            f"[{SCALE_DOF_FLOOR}, {DOF_LIMIT})"
        )
    input_values = _values(ndof, number, 0.35)
    input_vector = _vector(matrix, input_values)
    output = matrix.CreateColVector()
    with ng.TaskManager():
        output.data = matrix * input_vector
    output_values = output.FV().NumPy().copy()
    return {
        "input": input_values,
        "matvec": output_values,
        "ndof": float(ndof),
        "input_norm": float(np.linalg.norm(input_values)),
        "matvec_norm": float(np.linalg.norm(output_values)),
        "energy": float(np.dot(input_values, output_values)),
        "duration_s": time.perf_counter() - started,
    }


def _generate_mesh_files(output_directory: Path, mesh_ids) -> list[dict[str, Any]]:
    writer = Path(__file__).with_name("write_extended_mesh.py")
    rows: list[dict[str, Any]] = []
    for mesh_id in mesh_ids:
        diagnostics = ""
        completed_code = -1
        mesh = None
        path = None
        for attempt in range(1, 4):
            candidate = output_directory / f"{mesh_id}_attempt{attempt}.vol"
            if candidate.is_file() and candidate.stat().st_size > 0:
                try:
                    existing = ng.Mesh(str(candidate))
                except Exception:
                    pass
                else:
                    if existing.ne > 0 and existing.nv > 0:
                        mesh = existing
                        path = candidate
                        completed_code = -3
                        break
            try:
                completed = subprocess.run(
                    [sys.executable, str(writer), mesh_id, str(candidate)],
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=65,
                )
                completed_code = completed.returncode
                diagnostics = (
                    f"exit={completed.returncode}; stdout={completed.stdout!r}; "
                    f"stderr={completed.stderr!r}"
                )
            except subprocess.TimeoutExpired as exc:
                completed_code = -2
                diagnostics = f"writer timeout: {exc}"
            if candidate.is_file() and candidate.stat().st_size > 0:
                try:
                    mesh = ng.Mesh(str(candidate))
                except Exception as exc:
                    diagnostics += f"; reload={exc!r}"
                else:
                    if mesh.ne > 0 and mesh.nv > 0:
                        path = candidate
                        break
            diagnostics += f"; attempt={attempt}"
        if mesh is None or path is None:
            raise RuntimeError(f"could not generate {mesh_id}: {diagnostics}")
        if mesh.ne <= 0 or mesh.nv <= 0:
            raise RuntimeError(f"generated mesh {mesh_id} is empty")
        rows.append({
            "id": mesh_id,
            "path": str(path.resolve()),
            "dimension": int(mesh.dim),
            "vertices": int(mesh.nv),
            "elements": int(mesh.ne),
            "writer_exit_code": int(completed_code),
        })
    return rows


def generate(output_directory: Path) -> tuple[Path, Path]:
    started = time.perf_counter()
    output_directory.mkdir(parents=True, exist_ok=True)
    small_rows = _generate_mesh_files(output_directory, SMALL_MESH_IDS)
    scale_rows = _generate_mesh_files(output_directory, SCALE_MESH_IDS)
    small_meshes = {row["id"]: ng.Mesh(row["path"]) for row in small_rows}
    scale_meshes = {row["id"]: ng.Mesh(row["path"]) for row in scale_rows}
    breadth_cases = build_breadth_catalog(small_rows)
    scale_cases = build_scale_catalog(scale_rows)
    mat_values: dict[str, Any] = {"fixture_schema": SCHEMA}

    realized_breadth: list[dict[str, Any]] = []
    for number, case in enumerate(breadth_cases, start=1):
        oracle = _breadth_oracle(small_meshes[case["mesh_id"]], case, number)
        key = case["oracle_key"]
        for name in ("matrix", "input", "matvec", "free_dofs", "ndof", "nnz"):
            mat_values[f"{key}_{name}"] = oracle[name]
        if oracle["solve_performed"]:
            mat_values[f"{key}_rhs"] = oracle["rhs"]
            mat_values[f"{key}_solution"] = oracle["solution"]
        realized_breadth.append({
            **case,
            "solve": oracle["solve_performed"],
            "python_ndof": int(oracle["ndof"]),
            "python_nnz": int(oracle["nnz"]),
            "python_duration_s": oracle["duration_s"],
        })

    realized_scale: list[dict[str, Any]] = []
    for number, case in enumerate(scale_cases, start=1):
        oracle = _scale_oracle(scale_meshes[case["mesh_id"]], case, number)
        key = case["oracle_key"]
        for name in ("input", "matvec", "ndof", "input_norm",
                     "matvec_norm", "energy"):
            mat_values[f"{key}_{name}"] = oracle[name]
        realized_scale.append({
            **case,
            "python_ndof": int(oracle["ndof"]),
            "python_duration_s": oracle["duration_s"],
        })

    manifest = {
        "schema": SCHEMA,
        "generated_utc": datetime.now(UTC).isoformat(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": importlib.metadata.version("scipy"),
        "ngsolve_version": importlib.metadata.version("ngsolve"),
        "netgen_mesher_version": importlib.metadata.version("netgen-mesher"),
        "dof_limit": DOF_LIMIT,
        "scale_dof_floor": SCALE_DOF_FLOOR,
        "small_meshes": small_rows,
        "scale_meshes": scale_rows,
        "breadth_case_count": len(realized_breadth),
        "scale_case_count": len(realized_scale),
        "breadth_cases": realized_breadth,
        "scale_cases": realized_scale,
        "maximum_breadth_dofs": max(row["python_ndof"] for row in realized_breadth),
        "minimum_scale_dofs": min(row["python_ndof"] for row in realized_scale),
        "maximum_scale_dofs": max(row["python_ndof"] for row in realized_scale),
        "python_total_duration_s": time.perf_counter() - started,
    }
    if len(realized_breadth) != 500 or len(realized_scale) != 20:
        raise RuntimeError("extended oracle catalog size changed")

    mat_path = output_directory / "python_extended_oracle.mat"
    manifest_path = output_directory / "python_extended_manifest.json"
    savemat(mat_path, mat_values, do_compression=False, oned_as="column")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest_path, mat_path


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: generate_extended_oracle.py OUTPUT_DIRECTORY")
    manifest, oracle = generate(Path(sys.argv[1]))
    print(json.dumps({"manifest": str(manifest), "oracle": str(oracle)}))


if __name__ == "__main__":
    main()
