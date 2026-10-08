"""Single-solve linear total-A contracts on an independently prescribed current."""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia.p1_linear import MU0, solve_p1_linear
from radia.p1_newton import solve_p1_newton

POINTS = np.array([[0.0, 0.0, 0.3], [0.4, 0.0, 0.0], [0.0, 0.4, 0.0]])


def make_mesh(boundary="outer"):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    iron = Box(Pnt(-.2, -.2, -.2), Pnt(.2, .2, .2)).mat("iron")
    coil = Box(Pnt(.5, -.1, -.1), Pnt(.7, .1, .1)).mat("coil")
    outer = Box(Pnt(-1, -1, -1), Pnt(1, 1, 1)).mat("air")
    outer.faces.name = boundary
    return ng.Mesh(OCCGeometry(Glue([iron, coil, outer-iron-coil])).GenerateMesh(maxh=.4))


@pytest.fixture(scope="module")
def linear_mesh():
    return make_mesh()


def current():
    # curl(0, 0, psi), psi=(x-.5)(.7-x)(.01-y^2); zero normal coil current.
    return 1e7 * ng.CF((-2*ng.y*(ng.x-.5)*(.7-ng.x),
                        -(1.2-2*ng.x)*(.01-ng.y*ng.y), 0))


def solve(mesh, **kwargs):
    options = dict(current_cf=current(), current_materials="coil",
                   mu_r_dict={"iron": 500.}, cg_tolerance=1e-10,
                   observation_points=POINTS)
    options.update(kwargs)
    return solve_p1_linear(mesh, **options)


def test_one_pcg_call_and_independent_matrix_equation(linear_mesh, monkeypatch):
    import radia.sparsesolv_ngsolve as native

    real = native.NativePCG
    calls = []

    class CountedPCG:
        def __init__(self, *args):
            self.pcg = real(*args)

        def Solve(self, *args):
            calls.append(1)
            return self.pcg.Solve(*args)

    monkeypatch.setattr(native, "NativePCG", CountedPCG)
    result = solve(linear_mesh)
    assert len(calls) == result["stats"]["linear_solves"] == 1
    assert result["stats"]["final_residual_relative"] <= 1e-10*(1+1e-9)
    fes = result["fes"]
    u, v = fes.TnT()
    nu = linear_mesh.MaterialCF({"iron": 1/(MU0*500)}, default=1/MU0)
    # Independent assembly uses two distinct material coefficients. A
    # permutation between Netgen element order and native element order
    # would change this residual (iron=500, air/coil=1).
    form = ng.BilinearForm(fes, symmetric=False)
    form += nu*ng.curl(u)*ng.curl(v)*ng.dx
    rhs = ng.LinearForm(fes)
    rhs += ng.InnerProduct(current(), v)*ng.dx("coil")
    with ng.TaskManager():
        form.Assemble()
        rhs.Assemble()
        residual = rhs.vec.CreateVector()
        residual.data = rhs.vec-form.mat*result["A"].vec
    free = np.array(fes.FreeDofs(), dtype=bool)
    relative = np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(rhs.vec.FV().NumPy()[free])
    assert relative <= 1e-9
    phase_sum = sum(result["stats"]["phases_s"].values())
    assert 0 < phase_sum <= result["stats"]["total_s"]


def test_matches_tight_newton_and_independent_gauged_direct(linear_mesh):
    fast = solve(linear_mesh)
    table = [[0., 0.], [10., MU0*500*10], [1e5, MU0*500*1e5]]
    baseline = solve_p1_newton(linear_mesh, table, current_cf=current(), current_materials="coil",
                               cg_tolerance=1e-10, newton_tolerance=1e-10,
                               observation_points=POINTS)
    np.testing.assert_allclose(fast["observation_B_T"], baseline["observation_B_T"],
                               rtol=1e-6, atol=1e-10)
    # An independent NGSolve assembly with a small mass gauge as an SPD control.
    fes = fast["fes"]
    u, v = fes.TnT()
    nu = linear_mesh.MaterialCF({"iron": 1/(MU0*500)}, default=1/MU0)
    matrix = ng.BilinearForm(fes, symmetric=True)
    matrix += nu*ng.curl(u)*ng.curl(v)*ng.dx + 1e-8/MU0*u*v*ng.dx
    rhs = ng.LinearForm(fes)
    rhs += current()*v*ng.dx("coil")
    a = ng.GridFunction(fes)
    with ng.TaskManager():
        matrix.Assemble()
        rhs.Assemble()
        a.vec.data = matrix.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*rhs.vec
    b = ng.curl(a)
    direct = np.asarray([list(b(linear_mesh(*p))) for p in POINTS])
    assert np.linalg.norm(direct-fast["observation_B_T"])/np.linalg.norm(direct) <= 1e-6
    with ng.TaskManager():
        e_fast = ng.Integrate(.5*nu*fast["B_cf"]*fast["B_cf"], linear_mesh)
        e_direct = ng.Integrate(.5*nu*b*b, linear_mesh)
    assert abs(e_fast-e_direct)/e_direct <= 1e-6


def test_zero_source_and_requested_precision(linear_mesh):
    out = solve(linear_mesh, current_cf=ng.CF((0., 0., 0.)), mixed_precision=False)
    assert out["stats"]["cg_iterations"] == 0
    assert out["stats"]["linear_solves"] == 1
    assert out["stats"]["final_residual_relative"] == 0
    assert not out["stats"]["mixed_precision"]
    assert np.all(out["A"].vec.FV().NumPy() == 0)


def test_rejects_incompatible_current(linear_mesh):
    with pytest.raises(ValueError, match="orthogonal"):
        solve(linear_mesh, current_cf=ng.CF((1e6, 0, 0)))


def test_rejects_taskmanager_and_nonconvergence(linear_mesh):
    with ng.TaskManager(), pytest.raises(RuntimeError, match="outside ngsolve.TaskManager"):
        solve(linear_mesh)
    with pytest.raises(RuntimeError, match="residual"):
        solve(linear_mesh, cg_max_iterations=1)


@pytest.mark.parametrize("kwargs", [
    {"mu_r_dict": {"iron": 0}}, {"mu_r_dict": {"iron": float("nan")}},
    {"mu_r_dict": {"missing": 2}}, {"cg_tolerance": 0},
    {"cg_max_iterations": 1.5}, {"cg_max_iterations": True},
    {"current_materials": []}, {"current_materials": "missing"},
    {"current_cf": ng.CF((1+1j, 0, 0))}, {"current_cf": ng.CF(1)},
    {"observation_points": [[0, 0, float("nan")]]}, {"mixed_precision": "auto"},
])
def test_invalid_inputs(linear_mesh, kwargs):
    with pytest.raises(ValueError):
        solve(linear_mesh, **kwargs)


def test_thread_agreement(linear_mesh):
    ng.SetNumThreads(1)
    try:
        serial = solve(linear_mesh)
        ng.SetNumThreads(8)
        parallel = solve(linear_mesh)
        np.testing.assert_allclose(serial["observation_B_T"], parallel["observation_B_T"],
                                   rtol=1e-6, atol=1e-10)
    finally:
        ng.SetNumThreads(1)


def test_curved_mesh_rejected():
    mesh = make_mesh()
    mesh.Curve(2)
    with pytest.raises(ValueError, match="straight"):
        solve(mesh)


def test_nondefault_dirichlet_label():
    mesh = make_mesh(boundary="wall")
    fast = solve(mesh, dirichlet="wall")
    table = [[0., 0.], [10., MU0*500*10], [1e5, MU0*500*1e5]]
    control = solve_p1_newton(mesh, table, current_cf=current(), current_materials="coil",
                              dirichlet="wall", cg_tolerance=1e-10, newton_tolerance=1e-10,
                              observation_points=POINTS)
    np.testing.assert_allclose(fast["observation_B_T"], control["observation_B_T"],
                               rtol=1e-6, atol=1e-10)
