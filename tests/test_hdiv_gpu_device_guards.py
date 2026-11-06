"""Fail-loud behaviour and known roots of the device material and solver.

Guards: on device-built tables, ``_ResidentBHMaterial`` matches an analytic cubic
law evaluated by NGSolve on the same rule (energy, gradient, tangent action, below
and above the cap, and finite differences); refuses NaN and overflowing states and
a tangent of a stale linearization; reuses its tables without another upload. The
shared Newton solve recovers a root manufactured from the analytic law and the
device operator, and the iteration limits of the linear and Newton solves raise and
leave the operator usable. The operator action of BDM1 and BDM2 matches the CPU
equation and repeats bitwise, with the same device ownership, over 100 applications.
Assumes: the CPU exact dense Gram, the charge map and NGSolve are the references
on identical inputs. Limits are those of the retired resident checks
(1e-11/1e-12 material, 1e-6 differences, 1e-10 root, 1e-11 residual). Optional GPU
lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane
loads the binding through the provenance-checked loader.
"""

import os
from pathlib import Path

import numpy as np
import pytest

CASES = [("tet", 1), ("tet", 2), ("hex", 1)]
INV_CHI = 2.
# H(M) = M + |M|^2 M up to the cap 1.5, then a straight barrier of slope 7.
LAW = dict(knots=np.array([0., 1., 2.]), coefficients=np.array([1., 0., 1., 0., 1., 3., 4., 2.]),
           integral_base=np.array([0., .75, 6.]), cap=1.5, barrier_slope=7., last_field=4.875,
           initial_reluctivity=1.)
TABLE = np.column_stack([[0., 100., 1000., 10000., 100000.], 4e-7*np.pi*np.array([0., 200., 1500., 12000., 110000.])])


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


def _case(tmp_path, topology, order, cuda):
    import ngsolve as ng
    from netgen.csg import unit_cube
    from radia import hdiv_gpu
    from radia.vim._vim import build_charge_gram
    if topology == "tet":
        path = tmp_path/f"cube-{order}.vol"
        unit_cube.GenerateMesh(maxh=.5).Save(str(path))
    else:
        path = Path(__file__).parent/"fixtures/hdiv/hex_pair.vol"
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        fes = ng.HDiv(mesh, order=order)
        setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=None, compression="dense")
        B, gram, _ = build_charge_gram(fes, curve_order=0, curve_gauss=8, gram_backend="exact-dense",
                                       exact_dense_memory_mb=256)
        demag = np.column_stack([B.T @ np.asarray(gram.matvec_sym(B[:, k].toarray().ravel())) for k in range(fes.ndof)])
    return fes, setup, demag


def _analytic(fes, state, direction):
    """Energy, gradient and tangent action of LAW, by NGSolve on the material rule."""
    import ngsolve as ng
    from radia.vim._nonlinear import _EnergyMaterialQuadrature
    measure = _EnergyMaterialQuadrature(fes, TABLE).measure
    field = ng.GridFunction(fes)
    field.vec.FV().NumPy()[:] = state
    amplitude = ng.sqrt(ng.InnerProduct(field, field))
    extra = ng.IfPos(amplitude - 1.5, amplitude - 1.5, 0.)
    secant = ng.IfPos(amplitude - 1.5, (4.875 + 7*extra)/amplitude, 1 + amplitude*amplitude)
    weight = ng.IfPos(amplitude - 1.5, (7 - secant)/(amplitude*amplitude), 2.)
    energy = ng.IfPos(amplitude - 1.5, 2.390625 + 4.875*extra + 3.5*extra*extra,
                      .5*amplitude*amplitude + .25*amplitude**4)
    u, v = fes.TnT()
    load = ng.LinearForm(fes)
    load += secant*field*v*measure
    load.Assemble()
    tangent = ng.BilinearForm(fes)
    tangent += (secant*ng.InnerProduct(u, v) + weight*ng.InnerProduct(field, u)*ng.InnerProduct(field, v))*measure
    tangent.Assemble()
    action = tangent.mat.CreateRowVector()
    source = action.CreateVector()
    source.FV().NumPy()[:] = direction
    action.data = tangent.mat*source
    return float(ng.Integrate(energy*measure, fes.mesh)), load.vec.FV().NumPy().copy(), action.FV().NumPy().copy()


@pytest.mark.parametrize("topology,order", CASES)
@pytest.mark.parametrize("scale", [.05, 2.], ids=["below-cap", "above-cap"])
def test_device_material_matches_analytic_law(tmp_path, request, cuda, topology, order, scale):
    import ngsolve as ng
    fes, setup, _ = _case(tmp_path, topology, order, cuda)
    nd = setup["dofs"]
    rng = np.random.default_rng(order)
    state, direction = scale*rng.standard_normal(nd), rng.standard_normal(nd)
    material = cuda._ResidentBHMaterial.from_flat_tables(**setup["evaluation_tables"], dofs=nd, law=LAW, backend="gpu")
    try:
        device = material.evaluate(state, direction)
        with ng.TaskManager():
            energy, gradient, action = _analytic(fes, state, direction)
        errors = dict(energy=abs(device["energy"] - energy)/abs(energy),
                      gradient=np.linalg.norm(device["gradient"] - gradient)/np.linalg.norm(gradient),
                      tangent=np.linalg.norm(device["tangent_action"] - action)/np.linalg.norm(action))
        step = 1e-6*max(1., scale)
        plus, minus = material.evaluate(state + step*direction, direction), material.evaluate(state - step*direction, direction)
        difference = (plus["gradient"] - minus["gradient"])/(2*step)
        # The barrier kink is not differentiable, so differences run below the cap only.
        below = scale < 1
        errors["difference"] = float(np.linalg.norm(difference - device["tangent_action"])/np.linalg.norm(action)) if below else 0.
        request.node.user_properties += [("dofs", nd)] + [(k, float(v)) for k, v in errors.items()]
        assert errors["energy"] <= 1e-12 and errors["gradient"] <= 1e-11 and errors["tangent"] <= 1e-11, errors
        assert errors["difference"] <= 1e-6, errors
    finally:
        material.close()


def test_device_material_refuses_invalid_state_and_stale_tangent(tmp_path, cuda):
    _, setup, _ = _case(tmp_path, "tet", 1, cuda)
    nd = setup["dofs"]
    rng = np.random.default_rng(27)
    state, direction = .05*rng.standard_normal(nd), rng.standard_normal(nd)
    material = cuda._ResidentBHMaterial.from_flat_tables(**setup["evaluation_tables"], dofs=nd, law=LAW, backend="gpu")
    try:
        first = material.evaluate(state, direction)
        uploaded = material.stats()["setup_upload_bytes"]
        for _ in range(5):
            again = material.evaluate(state, direction)
            assert np.array_equal(again["gradient"], first["gradient"])
            assert np.array_equal(again["tangent_action"], first["tangent_action"])
        assert material.stats()["setup_upload_bytes"] == uploaded and material.stats()["evaluations"] == 6
        for invalid in (np.full(nd, np.nan), np.full(nd, 1e308)):
            with pytest.raises((ValueError, RuntimeError)):
                material.evaluate(invalid, direction)
        assert np.array_equal(material.evaluate(state, direction)["gradient"], first["gradient"])
        assert np.array_equal(material.evaluate(np.zeros(nd), direction)["gradient"], np.zeros(nd))
        epoch = material.linearize(state)["linearization_epoch"]
        np.testing.assert_allclose(material.apply_tangent(direction, epoch), first["tangent_action"], rtol=1e-11, atol=1e-13)
        material.evaluate(state + .01*direction, np.empty(0))
        with pytest.raises((ValueError, RuntimeError), match="epoch"):
            material.apply_tangent(direction, epoch)
        fresh = material.linearize(state)["linearization_epoch"]
        assert fresh != epoch
        np.testing.assert_allclose(material.apply_tangent(direction, fresh), first["tangent_action"], rtol=1e-11, atol=1e-13)
    finally:
        material.close()
    assert material.stats()["closed"] and material.stats()["device_payload_bound_bytes"] == 0
    with pytest.raises((ValueError, RuntimeError)):
        material.evaluate(state, direction)


@pytest.mark.parametrize("topology,order", CASES)
def test_device_newton_recovers_a_manufactured_root(tmp_path, request, cuda, topology, order):
    import ngsolve as ng
    fes, setup, _ = _case(tmp_path, topology, order, cuda)
    nd = setup["dofs"]
    root = .05*np.random.default_rng(order).standard_normal(nd)
    with ng.TaskManager():
        _, gradient, _ = _analytic(fes, root, np.zeros(nd))
    material = cuda._ResidentBHMaterial.from_flat_tables(**setup["evaluation_tables"], dofs=nd, law=LAW, backend="gpu")
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    try:
        rhs = operator.apply(root) + gradient
        solve = lambda iterations: operator.solve_nonlinear_material(
            material, rhs, np.zeros(nd), 1e-11, 5e-14, iterations, 10000, "reference-diagonal", 32, .01)
        result = solve(30)
        with ng.TaskManager():
            _, found, _ = _analytic(fes, result["x"], np.zeros(nd))
        error = np.linalg.norm(result["x"] - root)/np.linalg.norm(root)
        residual = np.linalg.norm(operator.apply(result["x"]) + found - rhs)/np.linalg.norm(rhs)
        request.node.user_properties += [("dofs", nd), ("root_error", float(error)), ("residual", float(residual))]
        assert error <= 1e-10 and residual <= 1e-11, f"root={error:.3e} residual={residual:.3e}"
        with pytest.raises(RuntimeError, match="Newton"):
            solve(1)
        np.testing.assert_allclose(solve(30)["x"], result["x"], rtol=1e-9, atol=1e-12)
    finally:
        operator.close()
        material.close()


@pytest.mark.parametrize("order", [1, 2])
def test_device_linear_solve_limits_and_operator_ownership(tmp_path, request, cuda, order):
    _, setup, demag = _case(tmp_path, "tet", order, cuda)
    nd = setup["dofs"]
    rng = np.random.default_rng(3)
    rhs = setup["rhs_operator"] @ np.array([0., 0., 1e5])
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    try:
        with pytest.raises(RuntimeError, match=r"iterations=1.*requested=.*actual="):
            operator.solve_linear(mass, rhs, INV_CHI, 1e-12, 1, "mass-diagonal")
        with pytest.raises(ValueError, match="preconditioner"):
            operator.solve_linear(mass, rhs, INV_CHI, 1e-12, 100, "automatic")
        with pytest.raises(ValueError, match="dimension"):
            operator.solve_linear(mass, rhs[:-1], INV_CHI, 1e-12, 100, "mass-diagonal")
        assert operator.solve_linear(mass, rhs, INV_CHI, 1e-11, 4000, "mass-diagonal")["relative_residual"] <= 1e-11
        owned = operator.stats()["owned_device_bytes"]
        x = rng.standard_normal(nd)
        expected = demag @ x
        first = operator.apply(x)
        relative = np.linalg.norm(first - expected)/np.linalg.norm(expected)
        request.node.user_properties.append(("relative", float(relative)))
        assert relative <= 1e-8, f"order={order} relative={relative:.3e} limit=1e-8"
        for _ in range(100):
            assert np.array_equal(operator.apply(x), first)
            assert operator.stats()["owned_device_bytes"] == owned
        with pytest.raises(ValueError, match="dimension"):
            operator.apply(np.ones(nd + 1))
    finally:
        operator.close()
        mass.close()
    assert operator.stats()["closed"] and operator.stats()["owned_device_bytes"] == 0
    with pytest.raises(RuntimeError, match="closed"):
        operator.apply(np.ones(nd))
