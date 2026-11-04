"""Device AMG, AMS and BDDC solves against sparse Cholesky on identical inputs.

Guards: H1 AMG with CG, real and complex curl-curl with AMS, the complex SIBC form with a boundary
term, and BDDC at p=2,3 for H1 and HCurl each reach the sparse Cholesky solution of the same form,
with the Dirichlet lift kept and only permitted payload families crossing. A wrong Galerkin
product, gradient, interpolation, harmonic extension or Schur complement breaks convergence or the
agreement; a wrong true residual raises.
Assumes: the CPU direct solve is the reference on identical NGSolve inputs; this is parity of one
discretization, not accuracy. SOLUTION_REL and ITERATION_LIMIT are PROPOSED detectors, not adopted
criteria. Optional GPU lane: skips visibly unless RADIA_FEM_NATIVE_TEST_DIR is set.
"""
import os

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

SOLUTION_REL = 1e-6
ITERATION_LIMIT = 120
TOLERANCE = 1e-10


@pytest.fixture(scope="module")
def route():
    if not os.environ.get("RADIA_FEM_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_FEM_NATIVE_TEST_DIR")
    import radia.fem_gpu as module
    return module


def _cube(n, hexes=False):
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=hexes, nx=n, ny=n, nz=n,
                                mapping=lambda x, y, z: (x + 0.1 * y, y + 0.1 * z, z))


def _rhs(fes, kind):
    v = fes.TestFunction()
    form = ng.LinearForm(fes)
    if kind == "H1":
        form += (1 + ng.x * ng.y) * v * ng.dx
    else:
        form += ng.InnerProduct(ng.CF((1 + ng.y, ng.z, ng.x * ng.y)), v) * ng.dx
    form.Assemble()
    return np.array(form.vec.FV().NumPy(), copy=True)


def _check(route, solve, fes, coefficients, rhs, *, lift=None, bonus=(0, 0, 0), cycles=None):
    expected = route.cpu_solve(fes, coefficients, rhs, lift=lift, preconditioner="direct", bonus=bonus)
    actual = solve(fes, coefficients, rhs, Backend="gpu", lift=lift, bonus=bonus, tolerance=TOLERANCE, cycles=cycles)
    assert actual["backend"] == "gpu"
    difference = float(np.linalg.norm(actual["solution"] - expected["solution"]) / np.linalg.norm(expected["solution"]))
    assert difference <= SOLUTION_REL, f"solution differs by {difference}, limit {SOLUTION_REL}; iterations {actual['iterations']}"
    assert 0 < actual["iterations"] <= ITERATION_LIMIT, f"iterations {actual['iterations']}, limit {ITERATION_LIMIT}"
    assert actual["true_residual"] <= 1e-6, f"true residual {actual['true_residual']}"
    names = {row["payload"] for row in actual["ledger"]}
    assert not any(name.startswith(("matrix", "operator")) for name in names), names
    assert not any(row["direction"] == "D2H" and row["payload"].startswith(("mesh.", "table.", "coef.", "rhs.", "lift."))
                   for row in actual["ledger"]), actual["ledger"]
    return actual


def test_h1_amg_with_lift(route):
    mesh = _cube(10)
    fes = ng.H1(mesh, order=1, dirichlet="left|right")
    coefficients = route.Coefficients(stiff=2.0, mass=0.5, stiff_aux=2.0, mass_aux=0.5)
    lift = ng.GridFunction(fes)
    lift.Set(1 + ng.x, definedon=mesh.Boundaries("left|right"))
    result = _check(route, route.solve_h1, fes, coefficients, _rhs(fes, "H1"), lift=np.array(lift.vec.FV().NumPy(), copy=True))
    assert len(result["hierarchy"][0]["rows"]) >= 2, result["hierarchy"]


def test_h1_amg_neumann_block_pins_the_kernel(route):
    # A coarse level of a pure-Neumann block is singular; the dense inverse must stay finite.
    mesh = _cube(3)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    coefficients = route.Coefficients(stiff=1.0, stiff_aux=1.0)
    _check(route, route.solve_h1, fes, coefficients, _rhs(fes, "H1"))


def test_real_curl_curl_ams(route):
    mesh = _cube(8)
    fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet=".*")
    coefficients = route.Coefficients(stiff=1.0, mass=10.0, stiff_aux=1.0, mass_aux=10.0)
    result = _check(route, route.solve_hcurl, fes, coefficients, _rhs(fes, "HCurl"))
    assert len(result["hierarchy"]) == 4, result["hierarchy"]


def test_complex_curl_curl_with_sibc_term(route):
    mesh = _cube(8)
    fes = ng.HCurl(mesh, order=1, nograds=True, complex=True, dirichlet="left|right")
    robin = 40.0
    coefficients = route.Coefficients(stiff=1.0, mass=1e-3, stiff_aux=1.0, mass_aux=1e-3 + 20.0, trace_im={"back": robin},
                                      trace_aux={"back": robin})
    _check(route, route.solve_hcurl, fes, coefficients, _rhs(fes, "HCurl") + 0j, bonus=(4, 0, 0))


@pytest.mark.parametrize("order", [2, 3])
def test_h1_bddc(route, order):
    mesh = _cube(4)
    fes = ng.H1(mesh, order=order, dirichlet="left", wb_withedges=False)
    coefficients = route.Coefficients(stiff=3.0, stiff_aux=3.0)
    _check(route, route.solve_h1, fes, coefficients, _rhs(fes, "H1"))


@pytest.mark.parametrize("order", [2, 3])
def test_complex_hcurl_bddc(route, order):
    mesh = _cube(3)
    fes = ng.HCurl(mesh, order=order, complex=True, dirichlet=".*")
    coefficients = route.Coefficients(stiff=1.0, mass_im=30.0, stiff_aux=1.0, mass_aux=30.0)
    _check(route, route.solve_hcurl, fes, coefficients, _rhs(fes, "HCurl") + 0j, bonus=(2, 2, 2), cycles=3)


def test_real_right_hand_side_is_complex_checked(route):
    mesh = _cube(2)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    coefficients = route.Coefficients(stiff=1.0, stiff_aux=1.0)
    system = route.DeviceSystem(fes, coefficients)
    system.setup("amg")
    with pytest.raises(ValueError, match="complex but the system is real"):
        system.solve(np.ones(fes.ndof) * 1j)
    with pytest.raises(ValueError, match="must be a finite vector"):
        system.solve(np.ones(fes.ndof + 1))
    with pytest.raises(RuntimeError, match="right-hand side is zero"):
        system.solve(np.zeros(fes.ndof))


def test_true_residual_limit_raises(route, monkeypatch):
    mesh = _cube(4)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    coefficients = route.Coefficients(stiff=1.0, stiff_aux=1.0)
    monkeypatch.setattr(route, "RELATIVE_LIMIT", 1e-30)
    with pytest.raises(RuntimeError, match="true relative residual"):
        route.solve_h1(fes, coefficients, _rhs(fes, "H1"), Backend="gpu", tolerance=1e-6)


def test_bddc_rejects_lowest_order(route):
    mesh = _cube(2)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    system = route.DeviceSystem(fes, route.Coefficients(stiff=1.0, stiff_aux=1.0))
    with pytest.raises(ValueError, match="BDDC applies to p>1"):
        system.setup("bddc")
