"""Curved (P2) TET device route against the CPU curved equation.

Guards: on an NGSolve ``Curve(2)`` ball, BDM1 and BDM2, the device-built
curved charge Gram entries, the energy operator B^T G B (device B from the
NGSolve moment tables and the per-host curved Gram solve), the reference
mass and the linear material solve agree with the CPU curved route:
``_curved_tet_charge_inputs`` (``_charge_basis_curved`` B and the native
curved ChargeGram) and the NGSolve mass on the material rule. Dropping the
curved Jacobian (affine map on P2 nodes), the |det J| weight of the curved
Gram mass, or the curved touching/Duffy branch moves these far beyond the
detectors. A P2 ``.vol`` stored on disk is solved as stored, and declaring it
linear is refused.
Assumes: the CPU curved route is the reference on identical NGSolve inputs;
this is parity of one discretization, not independent accuracy. The REL
values are PROPOSED detectors, not adopted criteria. Optional GPU lane: skips
visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane loads the
binding through the provenance-checked loader and fails if it cannot.
"""

import os

import numpy as np
import pytest

ENTRY_REL = 1e-10
OPERATOR_REL = 1e-9
MASS_REL = 1e-12
LINEAR_REL = 1e-8
INV_CHI = 0.1
PAIR_SAMPLE = 20000
_GEOMETRY = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents",
             "tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule", "far_factor", "family",
             "cell_numbers", "face_numbers", "curve_gl", "curve_gw")


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.fixture(scope="module")
def ball(tmp_path_factory):
    from netgen.occ import OCCGeometry, Pnt, Sphere
    import ngsolve as ng
    path = tmp_path_factory.mktemp("curved") / "ball.vol"
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.04))
    with ng.TaskManager():
        mesh.Curve(2)  # curved where the CAD lives; the device prepare uses the stored order
    mesh.ngmesh.Save(str(path))
    return path


@pytest.fixture(scope="module")
def saved_ball(ball, tmp_path_factory):
    import ngsolve as ng
    from radia.mesh_curve import load_vol, vol_geometry_kind
    path = tmp_path_factory.mktemp("curved-saved") / "ball-p2.vol"
    lines = ball.read_text().splitlines()
    path.write_text("\n".join(lines[:lines.index("endmesh") + 1]) + "\n")  # CAD-free stored P2, as a Cubit export
    assert vol_geometry_kind(path) is None and load_vol(path).GetCurveOrder() == 2
    return path


def _case(path, order, cuda, compression="dense"):
    import ngsolve as ng
    from radia import hdiv_gpu
    from radia.mesh_curve import load_vol
    from radia.vim._vim import _curved_tet_charge_inputs
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=None, compression=compression,
                                                    curve_order=2)
        mesh = load_vol(path)
        assert mesh.GetCurveOrder() == 2
        fes = ng.HDiv(mesh, order=order)
        cbk, kw = _curved_tet_charge_inputs(fes, max(2, 3*order), eps=1e-10, leafsize=64, eta=2., far_quad=3,
                                            ho_far_factor=2., curve_gauss=8, image_masks=[], image_signs=[],
                                            materialize_mass=False)
        # The device mass is E^T W E on the material rule IntegrationRule(TET, 2*max(6, 2p+2)).
        u, v = fes.TnT()
        rule = {ng.TET: ng.IntegrationRule(ng.TET, 2*max(6, 2*order + 2))}
        mass = ng.BilinearForm(fes)
        mass += u*v*ng.dx(intrules=rule)
        mass.Assemble()
    rows, cols, values = mass.mat.COO()
    from scipy.sparse import coo_matrix
    M = coo_matrix((values, (rows, cols)), shape=(fes.ndof, fes.ndof)).tocsr()
    assert setup["dofs"] == fes.ndof and cbk["B"].shape == (len(kw["charge_host"]), fes.ndof)
    return setup, cbk["B"], kw, M


def _cpu_gram(kw):
    from radia.vim._vim import _rp
    gram = _rp._ChargeGramHMatrix(**kw, build=False)
    gram.build_exact_dense_normalized_gram(maximum_memory_mb=256)
    return gram


@pytest.mark.parametrize("order", [1, 2])
def test_curved_device_entries_match_cpu_charge_gram(ball, request, cuda, order):
    setup, _, kw, _ = _case(ball, order, cuda)
    gram = _cpu_gram(kw)
    device = {name: setup["charge_gram"][name] for name in _GEOMETRY}
    n = len(device["charge_host"])
    numbers = [device["cell_numbers"].reshape(-1, 4), device["face_numbers"].reshape(-1, 3)]
    corners = [set(numbers[k][h]) for k, h in zip(device["charge_kind"], device["charge_host"])]
    rng = np.random.default_rng(20261103)
    pairs = np.unique(np.r_[np.c_[np.arange(n), np.arange(n)], rng.integers(0, n, size=(PAIR_SAMPLE, 2))], axis=0)
    touching = np.array([bool(corners[i] & corners[j]) for i, j in pairs])
    # Self, touching and separated pairs must all be present, or a branch is unobserved.
    assert touching.any() and (~touching).any() and (pairs[:, 0] == pairs[:, 1]).any()
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
        f"order={order} pair={pairs[worst].tolist()} touching={touching[worst]} actual={actual[worst]!r} "
        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={ENTRY_REL:.0e}")


@pytest.mark.parametrize("source", ["ball", "saved_ball"])
@pytest.mark.parametrize("order", [1, 2])
def test_curved_device_operator_and_mass_match_cpu(source, request, cuda, order):
    setup, B, kw, M = _case(request.getfixturevalue(source), order, cuda)
    gram = _cpu_gram(kw)
    nd = setup["dofs"]
    assert setup["charge_gram"]["family"] == "curved-tet" and setup["evaluation_tables"]["tables"][0]["vertices"].shape[1] == 30
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        from radia import hdiv_gpu
        assert hdiv_gpu._checked_device_assembly(operator.assembly(), setup["charge_gram"])["route"] == "device-curved-tet"
        assert hdiv_gpu._checked_material_assembly(mass.assembly(), setup["evaluation_tables"],
                                                   mass=True)["route"] == "device-curved-tet-mass"
        rng = np.random.default_rng(order)
        worst = dict(operator=0., mass=0.)
        for _ in range(4):
            x = rng.standard_normal(nd)
            expected = B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B @ x)))
            actual = operator.apply(np.ascontiguousarray(x))
            worst["operator"] = max(worst["operator"], np.linalg.norm(actual - expected)/np.linalg.norm(expected))
            y = mass.reference_mass_apply(np.ascontiguousarray(x))["y"]
            worst["mass"] = max(worst["mass"], np.linalg.norm(y - M @ x)/np.linalg.norm(M @ x))
        request.node.user_properties += [("dofs", nd)] + [("max_relative_" + k, float(v)) for k, v in worst.items()]
        assert worst["operator"] <= OPERATOR_REL and worst["mass"] <= MASS_REL, (
            f"order={order} {worst} limits operator={OPERATOR_REL:.0e} mass={MASS_REL:.0e}")
    finally:
        mass.close()
        operator.close()


def test_saved_curved_mesh_declared_linear_is_refused(saved_ball, cuda):
    import ngsolve as ng
    from radia import hdiv_gpu
    with ng.TaskManager(), pytest.raises(ValueError, match="affine or P2-curved"):
        hdiv_gpu._prepare_device_components(saved_ball, 1, cuda, aca=None, compression="dense", curve_order=1)


@pytest.mark.parametrize("order,compression", [(1, "aca"), (2, "dense"), (2, "aca")])
def test_curved_device_linear_solve_matches_cpu_equation(ball, request, cuda, order, compression):
    setup, B, kw, M = _case(ball, order, cuda, compression)
    gram = _cpu_gram(kw)
    nd = setup["dofs"]
    columns = [B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B[:, k].toarray().ravel()))) for k in range(nd)]
    rhs = setup["rhs_operator"] @ np.array([0., 0., 1e5])
    expected = np.linalg.solve(np.column_stack(columns) + INV_CHI*M.toarray(), rhs)
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        result = operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 4000, "mass-diagonal")
        relative = np.linalg.norm(result["x"] - expected)/np.linalg.norm(expected)
        request.node.user_properties += [("dofs", nd), ("relative", float(relative)), ("iterations", result["iterations"])]
        assert result["relative_residual"] <= 1e-11 and operator.assembly()["compression"] == compression
        assert relative <= LINEAR_REL, (f"order={order} {compression} relative={relative:.3e} limit={LINEAR_REL:.0e} "
                                        f"iterations={result['iterations']}")
    finally:
        mass.close()
        operator.close()


def test_native_staging_refuses_a_stored_order_below_the_request(ball, saved_ball, cuda, tmp_path):
    """The native staging reopens the .vol itself: it uses the stored order and never curves.

    A linear .vol requested at order 2 is refused before any space or table is built, with or
    without CAD in the file; the stored P2 ball, with and without its CAD, is accepted.
    """
    from netgen.occ import OCCGeometry, Pnt, Sphere
    linear = tmp_path / "ball-linear.vol"
    OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.04).Save(str(linear))
    lines = linear.read_text().splitlines()
    linear_nocad = tmp_path / "ball-linear-nocad.vol"
    linear_nocad.write_text("\n".join(lines[:lines.index("endmesh") + 1]) + "\n")
    for path in (linear, linear_nocad):
        with pytest.raises(ValueError, match="stored at geometry order 2; it stores order 1"):
            cuda._MaterialStagingContext(str(path), 1, 2)
    for path in (ball, saved_ball):
        staging = cuda._MaterialStagingContext(str(path), 1, 2)
        assert staging.flat_tables()["dofs"] > 0
