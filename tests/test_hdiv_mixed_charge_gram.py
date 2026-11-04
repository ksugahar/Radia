"""Typed charge Gram of connected tet, pyramid and hex meshes against analytic demag identities.

Guards: ``_build_charge_gram_mixed`` on a unit cube split into six pyramids, on one hex, one
pyramid and ten tets, and on a 2x2x2 grid with real far pairs, BDM1 and BDM2: a uniformly
magnetized cube has the demag factor 1/3 per axis by symmetry; the factors of a body with
warped pyramids sum to 1; eig(M^-1 B^T G B) stays in [0, 1]; G is symmetric; the CPU solver
reaches the dense equation of the same operators. Loud refusals of order 3, curved
geometry and images. Dropping the typed near rule, the rational pyramid charge or the face
orientation moves the factors far beyond DEMAG_REL. The energies of the rational charge on one
pyramid match a reduction that shares no rule with the Gram. The mass matches a long rule on a mesh
mapped onto a ball, whose cells the default rule of u*v*dx underintegrates.
Assumes: NGSolve owns the HDiv basis, orientation and Piola maps; the values guard code,
not the accuracy of a physical model, and DEMAG_REL, RATIONAL_REL and SOLVE_REL are PROPOSED detectors.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia.vim import _solve as solve
from radia.vim import _vim as V
from tests.test_hdiv_connected_meshes import connected_mesh
from validation_test.hmatrix_cuda.hdiv_mixed_meshes import build_mesh, mixed_grid, oriented, pyramid_cube, to_ball, transition_cube
from validation_test.hmatrix_cuda.rational_pyramid_reference import energies

DEMAG_REL = 1e-5
RATIONAL_REL = 5e-4
EIG_ROUNDOFF = 1e-12
SOLVE_REL = 1e-7
MASS_REL = 1e-5
_CUBES = dict(pyramids=pyramid_cube, transition=transition_cube, grid=lambda: mixed_grid((2, 2, 2)))


def _assemble(mesh, order):
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        B, G, M, _ = V._build_charge_gram_mixed(fes, build_hmatrix=False)
    n = B.shape[0]
    gram = np.array([[G.entry(i, j) for j in range(n)] for i in range(n)])
    return fes, B.toarray(), gram, M.toarray()


def _uniform(fes, direction, M):
    with ng.TaskManager():
        source = ng.LinearForm(fes)
        source += ng.InnerProduct(ng.CF(direction), fes.TestFunction())*ng.dx(bonus_intorder=8)
        source.Assemble()
    return np.linalg.solve(M, source.vec.FV().NumPy().copy())


def _factors(fes, B, G, M):
    factors = []
    for direction in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        u = _uniform(fes, direction, M)
        q = B @ u
        factors.append(float(q @ G @ q/(u @ M @ u)))
    return np.array(factors)


@pytest.mark.parametrize("name,order", [("pyramids", 1), ("pyramids", 2), ("transition", 1), ("transition", 2), ("grid", 1)])
def test_mixed_charge_gram_reproduces_the_cube(request, name, order):
    fes, B, G, M = _assemble(_CUBES[name](), order)
    factors = _factors(fes, B, G, M)
    eigenvalues = np.linalg.eigvals(np.linalg.solve(M, B.T @ G @ B)).real
    asymmetry = float(np.abs(G - G.T).max()/np.abs(G).max())
    request.node.user_properties += [("dofs", fes.ndof), ("factors", factors.tolist())]
    assert np.abs(factors - 1/3).max() <= DEMAG_REL, f"{name} order={order} factors={factors} expected 1/3"
    assert eigenvalues.min() >= -EIG_ROUNDOFF and eigenvalues.max() <= 1 + EIG_ROUNDOFF, (
        f"{name} order={order} eig(M^-1 N) in [{eigenvalues.min()}, {eigenvalues.max()}], expected [0, 1]")
    assert asymmetry == 0, f"{name} order={order} relative asymmetry {asymmetry}"


@pytest.mark.parametrize("order", [1, 2])
def test_mixed_charge_gram_trace_on_warped_pyramids(order):
    fes, B, G, M = _assemble(pyramid_cube(.3), order)
    factors = _factors(fes, B, G, M)
    assert abs(factors.sum() - 1) <= DEMAG_REL, f"order={order} factors={factors} sum={factors.sum()} expected 1"


def test_rational_charge_energies_match_their_reduction():
    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1)]
    mesh = build_mesh(points, [oriented(points, [0, 1, 2, 3, 4])])
    assert [tuple(mesh[v].point) for v in mesh[ng.ElementId(ng.VOL, 0)].vertices] == points
    _, _, G, _ = _assemble(mesh, 2)
    reference = energies()
    actual = dict(constant=4*np.pi*G[0, 0], rational=4*np.pi*G[4, 4], mixed=4*np.pi*G[0, 4])
    for name, value in actual.items():
        assert abs(value - reference[name]) <= RATIONAL_REL*reference[name], f"{name} energy {value}, expected {reference[name]}"


@pytest.mark.parametrize("order", [1, 2])
def test_mixed_mass_rule_is_long_enough_on_warped_cells(order):
    mesh = mixed_grid((2, 2, 2), lengths=(2., 2., 2.), origin=(-1., -1., -1.), mapping=to_ball(1.))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        M = V._build_charge_gram_mixed(fes, build_hmatrix=False)[2].toarray()
        u, v = fes.TnT()
        form = ng.BilinearForm(fes)
        form += u*v*ng.dx(bonus_intorder=16)
        form.Assemble()
    rows, cols, values = form.mat.COO()
    expected = np.zeros_like(M)
    expected[np.asarray(rows), np.asarray(cols)] = np.asarray(values)
    relative = float(np.linalg.norm(M - expected)/np.linalg.norm(expected))
    assert relative <= MASS_REL, f"order={order} relative={relative:.3e} limit={MASS_REL:.0e}"


@pytest.mark.parametrize("order", [1, 2])
def test_mixed_cpu_solver_reaches_the_dense_equation(order):
    mesh = connected_mesh("tet-pyramid-hex")
    fes, B, G, M = _assemble(mesh, order)
    applied = (0., 0., 1e4)
    mu_r = 11.
    with ng.TaskManager():
        result = solve.hdiv_demag_solve(mesh, mu_r=mu_r, H_ext=ng.CF(applied), order=order, tol=1e-11)
        source = ng.LinearForm(fes)
        source += ng.InnerProduct(ng.CF(applied), fes.TestFunction())*ng.dx(bonus_intorder=8)
        source.Assemble()
    x = np.linalg.solve(B.T @ G @ B + M/(mu_r - 1), source.vec.FV().NumPy().copy())
    expected = _mean(mesh, fes, x)
    relative = np.linalg.norm(result["M_avg"] - expected)/np.linalg.norm(expected)
    assert result["ndof"] == fes.ndof and relative <= SOLVE_REL, (
        f"order={order} M_avg={result['M_avg']} expected={expected} relative={relative:.3e}")


def _mean(mesh, fes, x):
    field = ng.GridFunction(fes, novisual=True)
    field.vec.FV().NumPy()[:] = x
    volume = ng.Integrate(ng.CF(1), mesh, order=4)
    with ng.TaskManager():
        return np.array(ng.Integrate(field, mesh, order=8))/volume


def test_mixed_charge_gram_refuses_unsupported_requests():
    mesh = connected_mesh("tet-pyramid-hex")
    with ng.TaskManager():
        with pytest.raises(NotImplementedError, match="order 1 or 2"):
            V._charge_basis_mixed(ng.HDiv(mesh, order=3))
        fes = ng.HDiv(mesh, order=1)
        with pytest.raises(ValueError, match="image_masks"):
            V.build_charge_gram(fes, image_masks=[1], image_signs=[1.])
        with pytest.raises(NotImplementedError, match="images? symmetry"):
            solve.hdiv_demag_solve(mesh, mu_r=11., H_ext=ng.CF((0., 0., 1.)), order=1, image="x")
        curved = connected_mesh("tet-pyramid-hex")
        curved.Curve(2)
        with pytest.raises(NotImplementedError, match="straight undeformed"):
            V._charge_basis_mixed(ng.HDiv(curved, order=1))
