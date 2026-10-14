"""Production validation for the public tetrahedral BDM1/BDM2 HDiv-VIM path.

The test intentionally uses only ``radia.vim`` APIs.  Flat/Curve(2) TET BDM2
``Solve``, ``ChargeGram``, and ``DemagOperator`` are public production routes.
HEX/WEDGE BDM2 is locked by its topology-specific tests.  The TET BDM2 route includes IMA,
Curve(2), and persistent C++ field evaluation and is locked by the contract tests.
"""

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

from radia.vim import DemagOperator, Solve  # noqa: E402


def _cube(h=0.9):
    return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=h))


def test_public_rt2_uniform_demag_is_order_invariant():
    """Uniform cube magnetization gives the analytic one-third factor at BDM1 and BDM2."""
    mesh = _cube()
    with ng.TaskManager():
        values = {
            p: DemagOperator(ng.HDiv(mesh, order=p), eps=1e-7).DemagFactor(ng.CF((0, 0, 1)))
            for p in (1, 2)
        }
    assert all(0.31 < value < 0.345 for value in values.values()), values
    assert abs(values[2] - 1.0 / 3.0) < abs(values[1] - 1.0 / 3.0), values
    assert abs(values[2] - values[1]) < 2e-3, values


def test_public_rt2_nonuniform_operator_improves_with_order():
    """The nonuniform M=(0,0,z) energy moves toward the high-order limit at BDM2."""
    mesh = _cube()
    with ng.TaskManager():
        values = {
            p: DemagOperator(ng.HDiv(mesh, order=p), eps=1e-7).DemagFactor(ng.CF((0, 0, ng.z)))
            for p in (1, 2)
        }
    assert abs(values[2] - 0.410) < abs(values[1] - 0.410), values
    assert abs(values[2] - values[1]) < 3e-3, values


def test_public_rt2_linear_and_nonlinear_material_solve():
    """BDM2 uses the all-C++ mass-Riesz CG and energy-Newton material paths."""
    mesh = _cube(h=1.2)
    applied = ng.CF((0, 0, 1000.0))
    mu0 = 4e-7 * np.pi
    linear_bh = [[0.0, 0.0], [1e3, mu0 * 100.0 * 1e3], [1e5, mu0 * 100.0 * 1e5]]
    with ng.TaskManager():
        linear = Solve(mesh, mu_r=100.0, H_ext=applied, order=2, gram_eps=1e-7)
        nonlinear = Solve(
            mesh,
            bh_table=linear_bh,
            H_ext=applied,
            order=2,
            gram_eps=1e-7,
            nl_maxit=30,
        )
    assert linear["order"] == nonlinear["order"] == 2
    assert linear["linear_solver"] == "mass-riesz-cg"
    assert nonlinear["nonlinear"] is True
    assert nonlinear["iters"] < 30
    assert abs(linear["demag"] - 1.0 / 3.0) < 1e-3
    rel = np.linalg.norm(np.asarray(nonlinear["M_avg"]) - np.asarray(linear["M_avg"])) / np.linalg.norm(
        linear["M_avg"]
    )
    assert rel < 1e-5, rel
