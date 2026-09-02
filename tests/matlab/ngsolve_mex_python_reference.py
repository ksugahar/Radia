"""Generate the Python/pybind NGSolve oracle for the MATLAB MEX gateway.

The fixture deliberately uses public NGSolve Python objects, not Radia's MEX
implementation.  MATLAB consumes the same Netgen ``.vol`` mesh and compares
native MEX results for finite-element spaces, coefficient evaluation,
GridFunction values, assembled sparse matrices, and a native linear solve.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from ngsolve import (
    BilinearForm,
    CF,
    CoefficientFunction,
    GridFunction,
    H1,
    HCurl,
    HDiv,
    InnerProduct,
    LinearForm,
    Mesh,
    curl,
    div,
    dx,
    grad,
    x,
    y,
    z,
)
from scipy.io import savemat
from scipy.sparse import coo_matrix


def _sparse_matrix(matrix) -> coo_matrix:
    """Return an NGSolve matrix as an explicitly shaped SciPy COO matrix.

    SciPy is used only to serialize the Python oracle into a MATLAB-readable
    MAT fixture.  Production MATLAB calculation stays in native NGSolve/MEX
    handles or MATLAB's own sparse representation.
    """

    rows, cols, values = matrix.COO()
    sparse = coo_matrix((values, (rows, cols)), shape=matrix.shape)
    sparse.eliminate_zeros()
    return sparse


def _space(mesh: Mesh, name: str, order: int):
    if name == "h1":
        return H1(mesh, order=order)
    if name == "hcurl":
        return HCurl(mesh, order=order, nograds=True)
    if name == "hdiv":
        return HDiv(mesh, order=order)
    raise ValueError(f"unsupported NGSolve space {name!r}")


def _assemble(mesh: Mesh, name: str, order: int, form: str):
    space = _space(mesh, name, order)
    trial, test = space.TnT()
    bilinear = BilinearForm(space, symmetric=True)
    if form == "mass":
        bilinear += InnerProduct(trial, test) * dx if name != "h1" else trial * test * dx
    elif form == "stiffness" and name == "h1":
        bilinear += InnerProduct(grad(trial), grad(test)) * dx
    elif form == "curlcurl" and name == "hcurl":
        bilinear += InnerProduct(curl(trial), curl(test)) * dx
    elif form == "divdiv" and name == "hdiv":
        bilinear += div(trial) * div(test) * dx
    else:
        raise ValueError(f"invalid form {form!r} for {name!r}")
    bilinear.Assemble()
    return space, bilinear


def _evaluate(coefficient: CoefficientFunction, mesh: Mesh, points: np.ndarray) -> np.ndarray:
    values = [coefficient(mesh(float(px), float(py), float(pz))) for px, py, pz in points]
    return np.asarray(values, dtype=float)


def build_reference(mesh_path: Path) -> dict[str, object]:
    mesh = Mesh(str(mesh_path))
    points = np.array(
        [[0.10, 0.10, 0.10], [0.20, 0.10, 0.10], [0.10, 0.20, 0.10]],
        dtype=float,
    )
    cases = (
        ("h1_mass", "h1", 2, "mass"),
        ("h1_stiffness", "h1", 2, "stiffness"),
        ("hcurl_mass", "hcurl", 2, "mass"),
        ("hcurl_curlcurl", "hcurl", 1, "curlcurl"),
        ("hdiv_mass", "hdiv", 2, "mass"),
        ("hdiv_divdiv", "hdiv", 1, "divdiv"),
    )
    result: dict[str, object] = {
        "fixture_schema": "radia.ngsolve-matlab-mex-parity.v1",
        "points": points,
        "mesh_dimension": mesh.dim,
        "mesh_vertices": mesh.nv,
        "mesh_elements": mesh.ne,
        "coordinates": _evaluate(CF((x, y, z)), mesh, points),
        "constant_vector": _evaluate(CF((1.25, -0.5, 0.75)), mesh, points),
    }

    assembled: dict[str, tuple[object, object]] = {}
    for prefix, space_name, order, form in cases:
        space, bilinear = _assemble(mesh, space_name, order, form)
        sparse = _sparse_matrix(bilinear.mat)
        result[f"{prefix}_matrix"] = sparse
        result[f"{prefix}_ndof"] = float(space.ndof)
        result[f"{prefix}_nnz"] = float(sparse.nnz)
        assembled[prefix] = (space, bilinear)

    h1_space, h1_mass = assembled["h1_mass"]
    grid = GridFunction(h1_space)
    grid_values = np.arange(1, h1_space.ndof + 1, dtype=float) / h1_space.ndof
    grid.vec.FV().NumPy()[:] = grid_values
    result["h1_grid_values"] = grid_values
    result["h1_grid_evaluation"] = _evaluate(grid, mesh, points)

    source = LinearForm(h1_space)
    source += 2.5 * h1_space.TestFunction() * dx
    source.Assemble()
    result["h1_constant_rhs"] = source.vec.FV().NumPy().copy()

    rhs_values = np.arange(1, h1_space.ndof + 1, dtype=float)
    rhs = h1_mass.mat.CreateColVector()
    rhs.FV().NumPy()[:] = rhs_values
    solution = h1_mass.mat.CreateColVector()
    solution.data = h1_mass.mat.Inverse(h1_space.FreeDofs()) * rhs
    result["h1_solver_rhs"] = rhs_values
    result["h1_solver_solution"] = solution.FV().NumPy().copy()
    residual = h1_mass.mat.CreateColVector()
    residual.data = h1_mass.mat * solution - rhs
    result["h1_solver_residual"] = residual.FV().NumPy().copy()
    return result


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: ngsolve_mex_python_reference.py INPUT.vol OUTPUT.mat"
        )
    savemat(sys.argv[2], build_reference(Path(sys.argv[1])), do_compression=False, oned_as="column")


if __name__ == "__main__":
    main()
