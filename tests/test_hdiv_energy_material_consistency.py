"""Real NGSolve weak-form derivative checks without a Radia native binary."""
import ast
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "energy_material_under_test", ROOT / "src/radia/vim/_nonlinear.py")
material_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(material_module)


def _table():
    h = np.array([0., 10., 100., 1000., 10000., 100000.])
    return np.column_stack((h, 4e-7*np.pi*(h + 1e6*h/(h+300))))


def test_inverse_curve_energy_and_tangent_are_derivatives():
    table = _table()
    fields, energy, maximum = material_module._bh_inverse_funcs(*table.T)
    for m in (maximum*.1, maximum*.5, maximum*1.01):
        step = max(m*1e-6, 1e-4)
        secant, tangent = fields(np.array([m]))
        derivative = (energy(np.array([m+step]))-energy(np.array([m-step])))/(2*step)
        hplus = fields(np.array([m+step]))[0]*(m+step)
        hminus = fields(np.array([m-step]))[0]*(m-step)
        np.testing.assert_allclose(derivative, secant*m, rtol=2e-7)
        np.testing.assert_allclose((hplus-hminus)/(2*step), tangent, rtol=2e-7)
    secant, tangent = fields(np.array([0.]))
    assert secant[0] > 0
    np.testing.assert_allclose(secant, tangent)
    assert energy(np.array([0.]))[0] == 0


@pytest.mark.parametrize("kind", ["tet", "hex", "wedge"])
@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("curved", [False, True])
def test_discrete_material_energy_gradient_and_hessian(kind, order, curved):
    kwargs = {"hexes": kind == "hex", "prism": kind == "wedge"}
    mesh = MakeStructured3DMesh(nx=1, ny=1, nz=1, **kwargs)
    if curved:
        deformation = ng.GridFunction(ng.VectorH1(mesh, order=2))
        deformation.Set(ng.CF((.05*ng.x*ng.z, .04*ng.x*ng.y, .03*ng.y*ng.z)))
        mesh.SetDeformation(deformation)
    fes = ng.HDiv(mesh, order=order)
    field, direction = ng.GridFunction(fes), ng.GridFunction(fes)
    with ng.TaskManager():
        field.Set(ng.CF((2e5+1e5*ng.x, 1e5*ng.y, 1e5*ng.z)))
        direction.Set(ng.CF((1e5*ng.y, 1e5*ng.x, 1e5*(1+ng.z))))
        law = material_module._EnergyMaterialQuadrature(
            fes, {name: _table() for name in mesh.GetMaterials()})
        m = field.vec.FV().NumPy().copy()
        d = direction.vec.FV().NumPy().copy()
        u, v = fes.TnT()

        def gradient(values):
            law.update(values)
            f = ng.LinearForm(fes)
            f += law.field*v*law.measure
            f.Assemble()
            return f.vec.FV().NumPy().copy()

        g = gradient(m)
        hessian = ng.BilinearForm(fes)
        hessian += law.bilinear_integrator(ng.InnerProduct(law.tangent*u, v))
        hessian.Assemble()
        result = direction.vec.CreateVector()
        hessian.mat.Mult(direction.vec, result)
        hd = result.FV().NumPy().copy()
        step = 1e-5
        de = (law.update(m+step*d)-law.update(m-step*d))/(2*step)
        np.testing.assert_allclose(de, g@d, rtol=2e-6)
        dg = (gradient(m+step*d)-gradient(m-step*d))/(2*step)
        np.testing.assert_allclose(dg, hd, rtol=3e-6, atol=1e-7)


def test_distinct_materials_and_zero_field_tangent():
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    first = Box(Pnt(0, 0, 0), Pnt(1, 1, 1)).mat("first")
    second = Box(Pnt(1, 0, 0), Pnt(2, 1, 1)).mat("second")
    mesh = ng.Mesh(OCCGeometry(Glue([first, second])).GenerateMesh(maxh=1))
    fes = ng.HDiv(mesh, order=2)
    table = _table()
    other = table.copy()
    other[:, 0] *= 2
    with ng.TaskManager():
        law = material_module._EnergyMaterialQuadrature(fes, {"first": table, "second": other})
        assert law.update(np.zeros(fes.ndof)) == 0
        assert set(law.regions) == {0, 1}
        for i, values in enumerate((table, other)):
            slope = material_module._bh_inverse_funcs(*values.T)[0](np.array([0.]))[1][0]
            tensor = law.tensor_samples.vec.FV().NumPy().reshape(9, law.npoints).T
            np.testing.assert_allclose(tensor[law.regions == i],
                                       np.tile((slope*np.eye(3)).ravel(), ((law.regions == i).sum(), 1)))
        g = ng.GridFunction(fes)
        g.Set(ng.CF((2e5, 1e5*ng.x, 0.)))
        m = g.vec.FV().NumPy().copy()
        law.update(m)
        load = ng.LinearForm(fes)
        load += law.field*fes.TestFunction()*law.measure
        load.Assemble()
        h = 1e-5
        de = (law.update((1+h)*m)-law.update((1-h)*m))/(2*h)
        np.testing.assert_allclose(de, load.vec.FV().NumPy()@m, rtol=2e-6)


def _line_search(energy):
    source = ROOT / "src/radia/vim/_solve.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    outer = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name == "_solve_nonlinear_energy_cpp")
    inner = next(n for n in outer.body if isinstance(n, ast.FunctionDef)
                 and n.name == "_line_search")
    namespace = {"np": np, "_energy": energy}
    exec(compile(ast.Module(body=[inner], type_ignores=[]), str(source), "exec"), namespace)
    return namespace["_line_search"]


@pytest.mark.parametrize("value", [2., float("nan"), float("inf")])
def test_armijo_exhaustion_never_returns_a_step(value):
    search = _line_search(lambda m, rhs: value)
    step, energy, count = search(np.ones(2), -np.ones(2), None, 1., 1.)
    assert step is None and energy is None and count == 34


@pytest.mark.parametrize("decrement", [0., -1., float("nan")])
def test_armijo_requires_a_finite_descent_direction(decrement):
    search = _line_search(lambda m, rhs: 0.)
    assert search(np.ones(2), -np.ones(2), None, 1., decrement) == (None, None, 0)


def test_armijo_accepts_an_evaluated_step():
    search = _line_search(lambda m, rhs: float(m@m)/2)
    assert search(np.ones(2), -np.ones(2), None, 1., 2.) == (1., 0., 0)


def _outer_solver(mass, *, reject_steps=False):
    """Real production Newton orchestration, with a dense SPD demag test double."""
    source = ROOT / "src/radia/vim/_solve.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    outer = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name == "_solve_nonlinear_energy_cpp")
    outer.body = [n for n in outer.body if not isinstance(n, ast.ImportFrom)]
    captured, calls = [], []

    class Demag:
        def apply_configured_demag(self, v, configured):
            return .005*mass@v

    def solve_native(operator, matrix, n, factor, rhs, tolerance, limit, **kwargs):
        import scipy.sparse as sparse
        rows, cols, values = matrix.COO()
        w = sparse.coo_matrix((values, (rows, cols)), shape=mass.shape).toarray()
        calls.append(np.linalg.norm(rhs))
        answer = np.linalg.solve(w+.005*mass, rhs)
        if reject_steps:
            answer *= -1
        return {"m": answer, "iters": 1}

    namespace = dict(np=np, ng=ng, _EnergyMaterialQuadrature=material_module._EnergyMaterialQuadrature,
                     _geometry_mass_apply=lambda op, h: mass@h,
                     _h_solve_mass_riesz=solve_native, _h_solve_auto_prec=solve_native,
                     _capture_cpp_solve_timings=lambda result: None,
                     _capture_nonlinear_solve_stats=lambda stats: captured.append(dict(stats)),
                     _f64=lambda v: np.asarray(v, dtype=float))
    exec(compile(ast.Module(body=[outer], type_ignores=[]), str(source), "exec"), namespace)
    return namespace[outer.name], Demag(), captured, calls


@pytest.mark.parametrize("mode", ["linear", "m0", "continuation", "zero", "reject"])
def test_production_newton_requires_a_residual_root(mode):
    import scipy.sparse as sparse
    mesh = MakeStructured3DMesh(nx=1, ny=1, nz=1, hexes=True)
    fes = ng.HDiv(mesh, order=2)
    with ng.TaskManager():
        u, v = fes.TnT()
        a = ng.BilinearForm(fes)
        a += ng.InnerProduct(u, v)*ng.dx
        a.Assemble()
        rows, cols, values = a.mat.COO()
        mass = sparse.coo_matrix((values, (rows, cols)), shape=(fes.ndof, fes.ndof)).toarray()
        solve, demag, captured, calls = _outer_solver(mass, reject_steps=mode == "reject")
        source = ng.GridFunction(fes)
        source.Set(ng.CF((0., 0., 0.)) if mode == "zero"
                   else ng.CF((1000.+100.*ng.x, 0., 100.*ng.z)))
        h = source.vec.FV().NumPy().copy()
        table = _table()
        if mode == "linear":
            table = np.array([[0., 0.], [1e5, 4e-7*np.pi*101.*1e5]])
        options = dict(reuse_tangent_steps=2, inner_preconditioner="mass-riesz")
        if mode in ("m0", "reject"):
            options["m0"] = np.zeros(fes.ndof)
        if mode == "continuation":
            options["continuation_steps"] = 3
        if mode == "reject":
            with pytest.raises(RuntimeError, match="Armijo line search failed") as exc:
                solve(mesh, fes, table, demag, fes.ndof, h, 1e-10, 500, 80, 1e-8, **options)
            assert captured[-1]["nonlinear_converged_final_stage"] is False
            assert captured[-1]["nonlinear_line_search_exhausted"] is True
            assert captured[-1]["nonlinear_convergence_mode"] == "line-search-exhausted"
            assert exc.value.nonlinear_stats == captured[-1]
            assert captured[-1]["nonlinear_rejected_newton_step_norm"] > 0
            assert captured[-1]["nonlinear_rejected_relative_newton_step"] is None
            return
        m, _ = solve(mesh, fes, table, demag, fes.ndof, h, 1e-10, 500, 80, 1e-8, **options)
        assert captured[-1]["nonlinear_final_relative_residual"] <= 1e-8
        assert captured[-1]["nonlinear_convergence_mode"] == "tolerance"
        assert captured[-1]["nonlinear_line_search_exhausted"] is False
        law = material_module._EnergyMaterialQuadrature(fes, table)
        law.update(m)
        load = ng.LinearForm(fes)
        load += law.field*v*law.measure
        load.Assemble()
        residual = load.vec.FV().NumPy()+.005*mass@m-mass@h
        assert np.linalg.norm(residual) <= 1e-8*max(np.linalg.norm(mass@h), 1.)
        if mode == "zero":
            assert calls == []
            np.testing.assert_array_equal(m, 0.)
        if mode == "linear":
            np.testing.assert_allclose(m, h/(.01+.005), rtol=1e-7, atol=1e-7)
