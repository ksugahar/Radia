"""Fast public-pipeline regression for flat HEX/WEDGE RT2 HDiv-VIM."""

import numpy as np
import pytest
import scipy.linalg as sla
import scipy.sparse as sp

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import (  # noqa: E402
    ChargeGram, FieldFromSolution, MagnetizationSource, Solve,
)


def _cube(kind):
    mapping = lambda x, y, z: (  # noqa: E731
        0.02 * (x - 0.5), 0.02 * (y - 0.5), 0.02 * (z - 0.5))
    kwargs = {"hexes": True} if kind == "hex" else {"prism": True}
    try:
        mesh = MakeStructured3DMesh(
            nx=1, ny=1, nz=1, mapping=mapping, **kwargs)
    except TypeError:
        pytest.skip("this NGSolve build cannot generate a structured prism mesh")
    expected_vertices = 8 if kind == "hex" else 6
    if {len(el.vertices) for el in mesh.Elements(ng.VOL)} != {expected_vertices}:
        pytest.skip(f"structured {kind} generator returned a different topology")
    return mesh


def _dense_demag(B, gram):
    B = sp.csr_matrix(B)
    operator = np.empty((B.shape[1], B.shape[1]))
    for column in range(B.shape[1]):
        basis = np.zeros(B.shape[1])
        basis[column] = 1.0
        charge = np.ascontiguousarray(B @ basis, dtype=np.float64)
        operator[:, column] = B.T @ np.asarray(gram.matvec(charge))
    return 0.5 * (operator + operator.T)


@pytest.mark.parametrize(
    ("kind", "expected_fes_dofs", "expected_charge_dofs"),
    [("hex", 108, 81), ("wedge", 141, 96)],
)
def test_rt2_charge_solve_and_field_pipeline(kind, expected_fes_dofs,
                                             expected_charge_dofs):
    mesh = _cube(kind)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        B, gram, mass = ChargeGram(fes, eps=1e-14, leafsize=256)
        result = Solve(
            mesh, mu_r=100.0, H_ext=ng.CF((0, 0, 1000.0)),
            order=2, tol=1e-9)

    assert fes.ndof == expected_fes_dofs
    assert B.shape == (expected_charge_dofs, expected_fes_dofs)
    assert mass.shape == (expected_fes_dofs, expected_fes_dofs)
    assert gram.ndof() == expected_charge_dofs
    assert gram.hex_state_check()["ctor"] == gram.hex_state_check()["now"]
    assert result["order"] == 2
    assert abs(result["demag"] - 1.0 / 3.0) < 5e-3
    assert result["iters"] < 50

    field = np.asarray(FieldFromSolution(
        result, np.asarray([[0.04, 0.0, 0.0], [0.0, 0.0, 0.04]]),
        algorithm="direct"))
    assert field.shape == (2, 3)
    assert np.isfinite(field).all()
    assert np.max(np.abs(field[:, 2])) > 1.0


@pytest.mark.parametrize("kind", ["hex", "wedge"])
def test_rt2_prescribed_magnetization_uses_the_native_field_path(kind):
    mesh = _cube(kind)
    with ng.TaskManager():
        source = MagnetizationSource(mesh, (0.0, 0.0, 1.0e5), order=2)

    assert source.stats["order"] == 2
    assert source.stats["hmatrix_built"] is False
    assert source.stats["projection_relative_residual"] < 1e-12
    field = np.asarray(source.Field(
        [[0.04, 0.0, 0.0], [0.0, 0.0, 0.04]], algorithm="direct"))
    assert field.shape == (2, 3)
    assert np.isfinite(field).all()
    assert np.max(np.abs(field[:, 2])) > 1.0


@pytest.mark.parametrize("kind", ["hex", "wedge"])
def test_rt2_demag_spectrum_is_physical(kind):
    mesh = _cube(kind)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        B, gram, mass = ChargeGram(fes, eps=1e-14, leafsize=256)
        operator = _dense_demag(B, gram)
    eigenvalues = sla.eigh(
        operator, sp.csr_matrix(mass).toarray(), eigvals_only=True)
    assert eigenvalues.min() > -1e-10
    assert eigenvalues.max() <= 1.0 + 2e-5
