"""The direct HCurl backend is native SparseCholesky, with true-residual gates."""
import pytest
ng = pytest.importorskip("ngsolve")
from radia import vector_potential_solver as vps


def _small_curl_curl_system(mass=1.0):
    """curl-curl plus ``mass`` times the identity on the unit cube, order 1.

    With ``mass=1`` the system is well conditioned and two exact
    factorisations agree to round-off; the gauged production systems use a
    1e-6 mass and are compared on their residuals, not their solutions.
    """
    from netgen.occ import unit_cube
    from ngsolve import (HCurl, BilinearForm, LinearForm, CoefficientFunction,
                         curl, InnerProduct, dx)
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.25))
    fes = HCurl(mesh, order=1)
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += InnerProduct(curl(u), curl(v)) * dx + mass * InnerProduct(u, v) * dx
    f = LinearForm(fes)
    f += InnerProduct(CoefficientFunction((0.0, 0.0, 1.0)), v) * dx
    with ng.TaskManager():
        a.Assemble()
        f.Assemble()
    return fes, a, f


@pytest.mark.parametrize("mass", [1.0, 1e-6])
def test_sparsecholesky_solves_well_conditioned_and_gauged_systems(mass):
    fes, a, f = _small_curl_curl_system(mass)
    solution, residual = f.vec.CreateVector(), f.vec.CreateVector()
    with ng.TaskManager():
        solution.data = a.mat.Inverse(fes.FreeDofs(), inverse=vps.direct_inverse_type()) * f.vec
        residual.data = f.vec - a.mat * solution
    assert residual.Norm() / f.vec.Norm() <= vps.LINEAR_RELATIVE_RESIDUAL_LIMIT


def test_direct_selector_has_no_pardiso_fallback(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "ngsolve.solvers.mkl_pardiso", None)
    assert vps.direct_inverse_type() == "sparsecholesky"
    assert vps.direct_inverse_type("sparsecholesky") == "sparsecholesky"
    with pytest.raises(ValueError):
        vps.direct_inverse_type("pardiso")
