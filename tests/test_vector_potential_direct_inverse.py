"""radia's direct HCurl factorisation: METIS ordering through NGSolve's own hooks.

The pip NGSolve PARDISO wrapper forces the minimum-degree ordering and ignores
the ``params`` it is given; radia registers its own SPD subclass that switches
to METIS after construction.  These tests pin that the registered type is what
``solve_linear`` / the Picard loop use, that the ordering really is METIS on
the constructed object, and that the solution is the shipped wrapper's.
"""
import pytest

ng = pytest.importorskip("ngsolve")

from radia import vector_potential_solver as vps  # noqa: E402


def _wrapper_has_pardiso():
    import ngsolve.solvers.mkl_pardiso as mp
    return mp._pardiso is not None


pytestmark = pytest.mark.skipif(
    not _wrapper_has_pardiso(),
    reason="ngsolve.solvers.mkl_pardiso found no MKL PARDISO in this environment")


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


def test_registered_type_is_spd_pardiso_with_metis_ordering():
    name = vps.direct_inverse_type()
    assert name == vps.DIRECT_INVERSE_TYPE
    assert vps.direct_inverse_type() == name          # idempotent
    fes, a, f = _small_curl_curl_system()
    with ng.TaskManager():
        inverse = a.mat.Inverse(fes.FreeDofs(), inverse=name)
    assert inverse.is_spd
    assert inverse._params[1] == 2, "ordering must be METIS nested dissection"
    with ng.TaskManager():
        shipped = a.mat.Inverse(fes.FreeDofs(), inverse="pardisospd")
    assert shipped._params[1] == 0, "the shipped wrapper still forces minimum degree"


def test_metis_ordering_returns_the_shipped_solution():
    fes, a, f = _small_curl_curl_system()
    x_metis, x_shipped, residual, difference = (f.vec.CreateVector() for _ in range(4))
    with ng.TaskManager():
        x_metis.data = a.mat.Inverse(fes.FreeDofs(), inverse=vps.direct_inverse_type()) * f.vec
        x_shipped.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec
        residual.data = f.vec - a.mat * x_metis
    difference.data = x_metis - x_shipped
    assert difference.Norm() <= 1.0e-10 * x_shipped.Norm()
    assert residual.Norm() <= 1.0e-10 * f.vec.Norm()


def test_gauged_system_residual_matches_the_shipped_wrapper():
    """The 1e-6-gauged production-like system is ill conditioned: two exact
    factorisations differ at 1e-7 relative and both leave a residual near
    5e-8 (measured for mass 1, 1e-3, 1e-6: residuals 5e-14, 5e-11, 5e-8 for
    every variant).  The ordering must not make that worse, and the residual
    must clear radia's own contract limit."""
    fes, a, f = _small_curl_curl_system(mass=1.0e-6)
    residuals = {}
    for name in (vps.direct_inverse_type(), "pardiso"):
        x, residual = f.vec.CreateVector(), f.vec.CreateVector()
        with ng.TaskManager():
            x.data = a.mat.Inverse(fes.FreeDofs(), inverse=name) * f.vec
            residual.data = f.vec - a.mat * x
        residuals[name] = residual.Norm() / f.vec.Norm()
    metis, shipped = residuals[vps.direct_inverse_type()], residuals["pardiso"]
    assert metis <= 3.0 * shipped
    assert metis <= vps.LINEAR_RELATIVE_RESIDUAL_LIMIT


def test_direct_paths_ask_for_the_registered_type():
    import inspect
    source = inspect.getsource(vps.VectorPotentialSolver.solve_linear)
    assert "inverse=direct_inverse_type()" in source
    assert "a.mat.Inverse(fes.FreeDofs())" not in source
    picard = inspect.getsource(vps.VectorPotentialSolver.solve_nonlinear)
    assert "inverse=direct_inverse_type()" in picard
    assert "inverse='pardisospd'" not in picard
