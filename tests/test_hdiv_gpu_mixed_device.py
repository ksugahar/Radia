"""Connected tet, pyramid and hex device route against the CPU typed charge Gram.

Guards: on a unit cube of six pyramids, the same with warped bases, a hex, a pyramid and ten
tets, and a 2x2x2 grid with real far pairs, BDM1 and BDM2, the device entries, the energy
operator B^T G B (device B from the per-class moment tables), the reference mass, the linear
and nonlinear solves and the three element-type evaluation groups agree with the CPU
``_build_charge_gram_mixed`` route, the staged material quadrature and the CPU solver on one
mesh. Dropping a host type, the rational pyramid charge, the near rule, an evaluation group or
the pyramid Jacobian moves these far beyond the detectors.
Assumes: the CPU typed route is the reference on identical NGSolve inputs; this is parity of
one discretization, not independent accuracy (``test_hdiv_mixed_charge_gram`` checks the
analytic cube). The REL values are PROPOSED detectors, not adopted criteria; MASS_REL leaves
room for the pyramid map of unpatched Netgen 6.2.2608, which NGSolve scales by 1 - 1e-10 and the
device lattice does not. Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is
set; a configured lane loads the binding through the provenance-checked loader and fails if
it cannot.
"""

import os

import numpy as np
import pytest
from scipy.sparse import csc_matrix, diags

ENTRY_REL = 1e-10
OPERATOR_REL = 1e-9
MASS_REL = 1e-9
TANGENT_REL = 1e-12
LINEAR_REL = 1e-7
MEAN_REL = 1e-7
INV_CHI = 0.1
PAIR_SAMPLE = 20000
ACA = dict(eps=1e-10, eta=2.0, leaf_size=8, max_rank=200, rank_cap="dense")
# The kernel expansion is a detector setting here, not an adopted profile: dropping the charge scale
# or crossing the Chebyshev axes moves the solution by order one, far beyond EXPANSION_REL (PROPOSED).
EXPANSION = dict(ACA, eps=1e-5, eta=1.5, leaf_size=64, expansion_order=6)
EXPANSION_REL = 5e-3
STOP_ACA = dict(ACA, eps=1e-5, leaf_size=32)
STOP_REL = 2e-5
BH_TABLE = np.array([[0., 0.], [50., .3], [100., .6], [200., 1.], [500., 1.3], [1e3, 1.45], [5e3, 1.65],
                     [2e4, 1.85], [1e5, 2.1], [1e6, 3.24]])
_GEOMETRY = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents",
             "tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule", "far_factor", "family",
             "hex_rules", "hex_bands", "host_types")
_NAMES = ("pyramids", "warped", "transition", "grid")
_CASES = [(name, order) for name in _NAMES for order in (1, 2) if (name, order) != ("grid", 2)]


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.fixture(scope="module")
def meshes(tmp_path_factory):
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import mixed_grid, pyramid_cube, to_ball, transition_cube
    root = tmp_path_factory.mktemp("mixed")
    target = {}
    for name, mesh in (("pyramids", pyramid_cube()), ("warped", pyramid_cube(.3)), ("transition", transition_cube()),
                       ("grid", mixed_grid((2, 2, 2))),
                       ("ball", mixed_grid((2, 2, 2), lengths=(2., 2., 2.), origin=(-1., -1., -1.), mapping=to_ball(1.)))):
        target[name] = root / (name + ".vol")
        mesh.ngmesh.Save(str(target[name]))
    return target


def _case(path, order, cuda, compression="dense"):
    import ngsolve as ng
    from radia import hdiv_gpu
    from radia.vim._vim import _build_charge_gram_mixed
    aca = None if compression == "dense" else ACA
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=aca, compression=compression)
        B, gram, _, _ = _build_charge_gram_mixed(ng.HDiv(ng.Mesh(str(path)), order=order), build_hmatrix=False)
    assert setup["family"] == "mixed" and setup["dofs"] == B.shape[1]
    assert len(setup["charge_gram"]["charge_host"]) == B.shape[0]
    return setup, B, gram


def _dense(gram, n):
    return np.array([[gram.entry(i, j) for j in range(n)] for i in range(n)])


def _material_mass(path, order, cuda):
    staging = cuda._MaterialStagingContext(str(path), order)
    try:
        payload = staging.quadrature()
    finally:
        staging.close()
    E = csc_matrix((payload["evaluation_values"], payload["evaluation_row_indices"], payload["evaluation_column_offsets"]),
                   shape=(payload["evaluation_rows"], payload["dofs"])).tocsr()
    return E.T @ diags(np.tile(np.asarray(payload["weights"]), 3)) @ E


@pytest.mark.parametrize("name,order", _CASES)
def test_mixed_device_entries_match_cpu_charge_gram(meshes, request, cuda, name, order):
    setup, _, gram = _case(meshes[name], order, cuda)
    device = {key: setup["charge_gram"][key] for key in _GEOMETRY}
    n = len(device["charge_host"])
    rng = np.random.default_rng(20261104)
    pairs = np.unique(np.r_[np.c_[np.arange(n), np.arange(n)], rng.integers(0, n, size=(PAIR_SAMPLE, 2))], axis=0)
    rows, cols = (np.ascontiguousarray(pairs[:, k], dtype=np.uint64) for k in (0, 1))
    actual = np.concatenate([cuda._diagnostic_flat_charge_entries(**device, rows=rows[s:s + 65536],
                                                                  cols=cols[s:s + 65536], backend="gpu")
                             for s in range(0, len(rows), 65536)])
    expected = np.array([gram.entry(int(i), int(j)) for i, j in pairs])
    diagonal = np.array([gram.entry(i, i) for i in range(n)])
    relative = np.abs(actual - expected)/np.sqrt(np.abs(diagonal[pairs[:, 0]]*diagonal[pairs[:, 1]]))
    worst = int(np.argmax(relative))
    request.node.user_properties += [("charges", n), ("pairs", len(pairs)), ("max_relative", float(relative.max()))]
    assert relative.max() <= ENTRY_REL, (
        f"{name} order={order} pair={pairs[worst].tolist()} actual={actual[worst]!r} "
        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={ENTRY_REL:.0e}")


@pytest.mark.parametrize("name,order", _CASES)
def test_mixed_device_operator_and_mass_match_cpu(meshes, request, cuda, name, order):
    from radia import hdiv_gpu
    setup, B, gram = _case(meshes[name], order, cuda)
    nd = setup["dofs"]
    M = _material_mass(meshes[name], order, cuda)
    G = _dense(gram, B.shape[0])
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        assert hdiv_gpu._checked_device_assembly(operator.assembly(), setup["charge_gram"])["route"] == "device-mixed"
        assert hdiv_gpu._checked_material_assembly(mass.assembly(), setup["evaluation_tables"],
                                                   mass=True)["route"] == "device-mixed-mass"
        rng = np.random.default_rng(order)
        worst = dict(operator=0., mass=0.)
        for _ in range(4):
            x = rng.standard_normal(nd)
            expected = B.T @ (G @ (B @ x))
            actual = operator.apply(np.ascontiguousarray(x))
            worst["operator"] = max(worst["operator"], np.linalg.norm(actual - expected)/np.linalg.norm(expected))
            y = mass.reference_mass_apply(np.ascontiguousarray(x))["y"]
            worst["mass"] = max(worst["mass"], np.linalg.norm(y - M @ x)/np.linalg.norm(M @ x))
        request.node.user_properties += [("dofs", nd)] + [("max_relative_" + k, float(v)) for k, v in worst.items()]
        assert worst["operator"] <= OPERATOR_REL and worst["mass"] <= MASS_REL, (
            f"{name} order={order} {worst} limits operator={OPERATOR_REL:.0e} mass={MASS_REL:.0e}")
    finally:
        mass.close()
        operator.close()


@pytest.mark.parametrize("name,order", [("ball", 1), ("ball", 2), ("transition", 2), ("grid", 1)])
def test_mixed_frozen_tangent_matches_the_matrix_free_action(meshes, cuda, name, order):
    from radia import hdiv_gpu
    from radia.vim._nonlinear import _bh_inverse_payload
    with_tables = hdiv_gpu._prepare_device_components(meshes[name], order, cuda, aca=ACA, compression="aca")
    law = _bh_inverse_payload(BH_TABLE[:, 0], BH_TABLE[:, 1])
    native_law = dict(law)
    for key in ("knots", "coefficients", "integral_base"):
        native_law[key] = np.ascontiguousarray(native_law[key], dtype=np.float64).ravel()
    nd = with_tables["dofs"]
    material = cuda._ResidentBHMaterial.from_flat_tables(**with_tables["evaluation_tables"], dofs=nd, law=native_law, backend="gpu")
    try:
        rng = np.random.default_rng(order)
        state = 1e3*rng.standard_normal(nd)
        direction = rng.standard_normal(nd)
        expected = material.evaluate(np.ascontiguousarray(state), np.ascontiguousarray(direction))["tangent_action"]
        epoch = material.linearize(np.ascontiguousarray(state))["linearization_epoch"]
        actual = material.apply_tangent(np.ascontiguousarray(direction), epoch)
        relative = float(np.linalg.norm(actual - expected)/np.linalg.norm(expected))
        assert relative <= TANGENT_REL, f"{name} order={order} relative={relative:.3e} limit={TANGENT_REL:.0e}"
    finally:
        material.close()


@pytest.mark.parametrize("name,order", [("ball", 1), ("ball", 2), ("transition", 2)])
def test_mixed_cpu_mass_takes_the_device_rule(meshes, cuda, name, order):
    import ngsolve as ng
    from radia.vim._vim import _build_charge_gram_mixed
    with ng.TaskManager():
        _, _, cpu, _ = _build_charge_gram_mixed(ng.HDiv(ng.Mesh(str(meshes[name])), order=order), build_hmatrix=False)
    device = _material_mass(meshes[name], order, cuda)
    relative = float(np.linalg.norm((cpu - device).toarray())/np.linalg.norm(device.toarray()))
    assert relative <= MASS_REL, f"{name} order={order} relative={relative:.3e} limit={MASS_REL:.0e}"


@pytest.mark.parametrize("name,order,compression,preconditioner", [
    ("pyramids", 2, "aca", "mass-diagonal"), ("warped", 1, "dense", "mass-diagonal"), ("transition", 1, "aca", "mass-cholesky"),
    ("transition", 2, "dense", "mass-cholesky"), ("grid", 1, "aca", "mass-diagonal")])
def test_mixed_device_linear_solve_matches_cpu_equation(meshes, request, cuda, name, order, compression, preconditioner):
    setup, B, gram = _case(meshes[name], order, cuda, compression)
    nd = setup["dofs"]
    M = _material_mass(meshes[name], order, cuda).toarray()
    rhs = setup["rhs_operator"] @ np.array([0., 0., 1e5])
    expected = np.linalg.solve(B.T @ (_dense(gram, B.shape[0]) @ B.toarray()) + INV_CHI*M, rhs)
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        result = operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 4000, preconditioner)
        relative = np.linalg.norm(result["x"] - expected)/np.linalg.norm(expected)
        request.node.user_properties += [("dofs", nd), ("relative", float(relative)), ("iterations", result["iterations"])]
        assert result["relative_residual"] <= 1e-11 and operator.assembly()["compression"] == compression
        assert relative <= LINEAR_REL, (f"{name} order={order} {compression} {preconditioner} relative={relative:.3e} "
                                        f"limit={LINEAR_REL:.0e} iterations={result['iterations']}")
    finally:
        mass.close()
        operator.close()


@pytest.mark.parametrize("name,order", [("transition", 1), ("transition", 2), ("warped", 1)])
def test_mixed_device_nonlinear_solve_matches_the_cpu_solver(meshes, request, cuda, name, order):
    import ngsolve as ng
    from radia.vim import _solve
    from radia.vim._nonlinear import _bh_inverse_payload
    setup, _, _ = _case(meshes[name], order, cuda, "aca")
    nd = setup["dofs"]
    applied = np.array([0., 0., 2e4])
    law = _bh_inverse_payload(BH_TABLE[:, 0], BH_TABLE[:, 1])
    law = {key: np.ascontiguousarray(value, dtype=np.float64).ravel() if key in ("knots", "coefficients", "integral_base") else value
           for key, value in law.items()}
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    material = cuda._ResidentBHMaterial.from_flat_tables(**setup["evaluation_tables"], dofs=nd, law=law, backend="gpu")
    try:
        result = operator.solve_nonlinear_material(material, np.ascontiguousarray(setup["rhs_operator"] @ applied), np.zeros(nd),
                                                   1e-10, 5e-13, 60, 10000, "reference-diagonal", 32, .01)
    finally:
        material.close()
        operator.close()
    device = setup["mean_operator"] @ result["x"]
    with ng.TaskManager():
        cpu = _solve.hdiv_demag_solve(ng.Mesh(str(meshes[name])), bh_table=BH_TABLE.tolist(), H_ext=ng.CF(tuple(applied)),
                                      order=order, tol=1e-11, nl_tol=1e-9)
    relative = np.linalg.norm(device - cpu["M_avg"])/np.linalg.norm(cpu["M_avg"])
    request.node.user_properties += [("dofs", nd), ("relative", float(relative))]
    assert relative <= MEAN_REL, f"{name} order={order} device={device} cpu={cpu['M_avg']} relative={relative:.3e}"


def _solved(path, order, cuda, aca, compression):
    from radia import hdiv_gpu
    setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=aca, compression=compression)
    nd = setup["dofs"]
    rhs = np.ascontiguousarray(setup["rhs_operator"] @ np.array([0., 0., 1e5]))
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        assembly = hdiv_gpu._checked_device_assembly(operator.assembly(), setup["charge_gram"])
        return operator.solve_linear(mass, rhs, INV_CHI, 1e-11, 6000, "mass-diagonal"), assembly
    finally:
        mass.close()
        operator.close()


@pytest.mark.parametrize("order", (1, 2))
def test_mixed_kernel_expansion_stays_near_the_exact_device_gram(tmp_path, request, cuda, order):
    """Guards: far blocks built from the kernel at Chebyshev nodes solve to nearly the dense device solution.

    Assumes: the dense device route is the reference on identical inputs; this is not an accuracy claim.
    """
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import mixed_grid
    path = tmp_path / "grid.vol"
    mixed_grid((4, 4, 4)).ngmesh.Save(str(path))
    exact, _ = _solved(path, order, cuda, None, "dense")
    result, assembly = _solved(path, order, cuda, EXPANSION, "hca")
    counters = assembly["hmatrix"]
    assert counters["expansion_order"] == 6 and counters["lowrank_leaves"] > 0, counters
    assert counters["moment_coefficients"] > 0 and counters["kernel_rank_sum"] > 0, counters
    relative = np.linalg.norm(result["x"] - exact["x"])/np.linalg.norm(exact["x"])
    request.node.user_properties += [("relative", float(relative)), ("rank_sum", counters["rank_sum"])]
    assert relative <= EXPANSION_REL, f"order={order} relative={relative:.3e} limit={EXPANSION_REL:.0e}"


def test_aca_stop_verifies_rows_before_a_leaf_stops(tmp_path, request, cuda):
    """Guards: a small last term does not stop a leaf until probe rows agree with its factors.

    Assumes: the dense device route is the reference on identical inputs. Without the probes the distance is six times
    the one with them on this mesh, because a few leaves stop early by hundreds of eps. The limit is a proposed value.
    """
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import mixed_grid
    path = tmp_path / "grid.vol"
    mixed_grid((6, 6, 6)).ngmesh.Save(str(path))
    exact, _ = _solved(path, 1, cuda, None, "dense")
    result, assembly = _solved(path, 1, cuda, STOP_ACA, "aca")
    counters = assembly["hmatrix"]
    assert counters["lowrank_leaves"] > 100, counters
    relative = np.linalg.norm(result["x"] - exact["x"])/np.linalg.norm(exact["x"])
    request.node.user_properties += [("relative", float(relative)), ("rank_sum", counters["rank_sum"])]
    assert relative <= STOP_REL, f"relative={relative:.3e} limit={STOP_REL:.0e} lowrank={counters['lowrank_leaves']}"


def test_tangent_rule_is_refused_outside_mixed_meshes(tmp_path, cuda):
    """Guards: a pure hex mesh keeps its own near rules and never takes the tangent rule of the mixed family."""
    from radia import hdiv_gpu
    from ngsolve.meshes import MakeStructured3DMesh
    path = tmp_path / "hex.vol"
    MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1).ngmesh.Save(str(path))
    setup = hdiv_gpu._prepare_device_components(path, 1, cuda, aca=None, compression="dense")
    assert setup["family"] == "hex"
    gram = dict(setup["charge_gram"], hex_rules=np.vstack([setup["charge_gram"]["hex_rules"], [10., .5, 1.]]))
    with pytest.raises(ValueError, match="belongs to mixed meshes"):
        cuda._ResidentHDiv.from_flat_charges(**gram, dofs=setup["dofs"], backend="gpu")


def test_kernel_expansion_is_refused_outside_mixed_meshes(tmp_path, cuda):
    """Guards: the route is not silently replaced by another compression on a family it does not cover."""
    from radia import hdiv_gpu
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import tet_cube
    path = tmp_path / "tets.vol"
    tet_cube().ngmesh.Save(str(path))
    setup = hdiv_gpu._prepare_device_components(path, 1, cuda, aca=EXPANSION, compression="hca")
    assert setup["family"] != "mixed"
    with pytest.raises(ValueError, match="belongs to mixed meshes"):
        cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=setup["dofs"], backend="gpu")


def test_mixed_device_refuses_unsupported_requests(meshes, cuda):
    from radia import hdiv_gpu
    with pytest.raises(NotImplementedError, match="order 1 or 2 on straight geometry"):
        hdiv_gpu._prepare_device_components(meshes["transition"], 3, cuda, aca=ACA, compression="aca")
    with pytest.raises(NotImplementedError, match="order 1 or 2 on straight geometry"):
        hdiv_gpu._prepare_device_components(meshes["transition"], 1, cuda, aca=ACA, compression="aca", curve_order=2)
