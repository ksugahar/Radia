"""Typed charge Gram of connected tet, pyramid and hex meshes against analytic demag identities.

Guards: ``_build_charge_gram_mixed`` on a unit cube split into six pyramids, on one hex, one
pyramid and ten tets, and on a 2x2x2 grid with real far pairs, BDM1 and BDM2: a uniformly
magnetized cube has the demag factor 1/3 per axis by symmetry; the factors of a body with
warped pyramids sum to 1; eig(M^-1 B^T G B) stays in [0, 1]; G is symmetric; the CPU solver
reaches the dense equation of the same operators. Loud refusals of order 3, curved
geometry and images. Dropping the typed near rule, the rational pyramid charge or the face
orientation moves the factors far beyond DEMAG_REL. The energies of the rational charge on one
pyramid match a reduction that shares no rule with the Gram. The mass matches a long rule on a mesh
mapped onto a ball, whose cells the default rule of u*v*dx underintegrates. Both routes refuse a cell
whose map changes sign, the bound behind that refusal holds against the element map, and every mesh
family of the record stays one-signed. On warped cells the tangent closed form with its bounded remainder stays
near the radial, fan and level rules at count 0, which integrate the same entries with no closed form of the source;
dropping the remainder or the rational Taylor term moves the operator beyond WARPED_INNER_REL.
Assumes: NGSolve owns the HDiv basis, orientation and Piola maps; the values guard code, not the accuracy of a
physical model, and DEMAG_REL, RATIONAL_REL, SOLVE_REL and WARPED_INNER_REL are PROPOSED detectors.
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
WARPED_INNER_REL = 1e-2
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


@pytest.mark.parametrize("name,order", [("warped", 1), ("warped", 2), ("ball", 1)])
def test_tangent_rule_stays_near_the_rules_without_a_closed_form(monkeypatch, name, order):
    mesh = pyramid_cube(.2) if name == "warped" else mixed_grid((2, 2, 2), lengths=(2., 2., 2.), origin=(-1., -1., -1.), mapping=to_ball(1.))
    rules = list(V._mixed_gram_rules(order))
    assert rules[6] > 0, f"production tangent count {rules[6]}, expected a tangent rule"
    _, B, G, M = _assemble(mesh, order)
    rules[1], rules[6] = (6 if order == 1 else 14), 0
    monkeypatch.setattr(V, "_mixed_gram_rules", lambda p: tuple(rules))
    _, B0, G0, _ = _assemble(mesh, order)
    assert np.array_equal(B, B0)
    difference = B.T @ (G - G0) @ B
    distance = float(np.abs(np.linalg.eigvals(np.linalg.solve(M, .5*(difference + difference.T)))).max())
    assert distance <= WARPED_INNER_REL, f"{name} order={order} operator distance {distance:.3e}, limit {WARPED_INNER_REL:.0e}"


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


_TWISTED_PYRAMID = [(0, 0, 0), (1, 0, 0), (1, 1, 2.2), (0, 1, 0), (.5, .5, .5)]
_TWISTED_HEX = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (.2, .2, .2), (0, 1, 1)]


def _one_cell(points):
    return build_mesh(points, [list(range(len(points)))])


@pytest.mark.parametrize("points", [_TWISTED_PYRAMID, _TWISTED_HEX], ids=["pyramid", "hex"])
@pytest.mark.parametrize("route", ["cpu", "device"])
def test_mixed_route_refuses_a_cell_whose_map_changes_sign(points, route):
    with ng.TaskManager():
        fes = ng.HDiv(_one_cell(points), order=1)
        with pytest.raises(ValueError, match="changes sign"):
            if route == "cpu":
                V._build_charge_gram_mixed(fes, build_hmatrix=False)
            else:
                V._prepare_gpu_mixed_charge_inputs(fes, backend="gpu")


@pytest.mark.parametrize("name", ["warped", "twisted", "hex"])
def test_mixed_jacobian_extremes_bound_the_element_map(name):
    points = dict(warped=[(0, 0, 0), (1, 0, 0), (1.2, 1, .3), (0, 1, 0), (.5, .5, .8)], twisted=_TWISTED_PYRAMID, hex=_TWISTED_HEX)[name]
    mesh = _one_cell(points)
    element = next(iter(mesh.Elements(ng.VOL)))
    jacobian = ng.specialcf.JacobianMatrix(3)
    values = np.asarray(jacobian(mesh.GetTrafo(element)(ng.IntegrationRule(element.type, 14)))).reshape(-1, 3, 3)
    sampled = np.linalg.det(values)
    low, high = V._mixed_jacobian_extremes(np.array([v.point for v in mesh.vertices], dtype=float), np.array([[v.nr for v in element.vertices]]))
    # The element's own orientation sets one overall sign; the bounds must then hold with the same sign.
    sign = np.sign(sampled[np.argmax(np.abs(sampled))])*np.sign(low[0] if abs(low[0]) > abs(high[0]) else high[0])
    bound_low, bound_high = sorted((sign*low[0], sign*high[0]))
    assert bound_low <= sampled.min() + 1e-12 and sampled.max() <= bound_high + 1e-12, (
        f"{name}: sampled [{sampled.min()}, {sampled.max()}] outside bound [{bound_low}, {bound_high}]")


@pytest.mark.parametrize("name,mesh", [("cube", pyramid_cube), ("warped", lambda: pyramid_cube(.3)), ("transition", transition_cube),
                                       ("grid", lambda: mixed_grid((3, 3, 3))), ("connected", lambda: connected_mesh("tet-pyramid-hex")),
                                       ("ball", lambda: mixed_grid((6, 6, 6), lengths=(2., 2., 2.), origin=(-1., -1., -1.), mapping=to_ball(1.)))])
def test_mixed_meshes_of_the_record_are_one_signed(name, mesh):
    V._require_one_signed_mixed_jacobians(mesh())


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
