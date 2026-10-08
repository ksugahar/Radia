"""Small independent basis-change and source-work checks on a closed surface."""
import numpy as np
import pytest
import json
from pathlib import Path
import ngsolve as ng
import netgen.meshing as nm

from radia.bem_loop_work import _loop_work_row, _p0_single_layer, MU_0
from radia.bem_sibc_solver import ScalarBIESIBCSolver


def test_genus1_study_exercises_production_power_gate():
    """The historical low-level study bypassed the CLI power rejection.

    Its stored self-authored ring powers would be rejected by the production
    gate; the corrected three-level powers pass the identical tolerance.
    """
    from radia.workpiece_surface import _check_sibc_reaction_power
    result = Path(__file__).resolve().parents[1] / "validation_test" / "induction_heating" / "results" / "scalar_loop_work_ring.json"
    data = json.loads(result.read_text())
    omega = 2*np.pi*data["model"]["frequency_hz"]
    for level in data["values"]:
        for case in level["cases"]:
            old, new = case["old"], case["new"]
            with pytest.raises(RuntimeError, match="power balance failed"):
                _check_sibc_reaction_power(old["power"], -.5*omega*old["reaction"][1])
            assert _check_sibc_reaction_power(new["power"], -.5*omega*new["reaction"][1]) < .1


@pytest.fixture(scope="module", params=["dense", "fmm"])
def surface(request):
    points = np.array([[0., 0., 0.], [1., 0., 0.],
                       [0., 1., 0.], [0., 0., 1.]])
    triangles = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    mesh = nm.Mesh(dim=3)
    face = mesh.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    vertices = [mesh.Add(nm.MeshPoint(nm.Pnt(*p))) for p in points]
    for t in triangles:
        mesh.Add(nm.Element2D(face, [vertices[i] for i in t]))
    mesh = ng.Mesh(mesh)
    p = points[triangles]
    cross = np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])
    areas = np.linalg.norm(cross, axis=1)/2
    normals = cross/(2*areas[:, None])
    gradients = np.array([np.cross(normals, p[:, (k+2) % 3]-p[:, (k+1) % 3])
                          /(2*areas[:, None]) for k in range(3)]).transpose(1, 0, 2)
    with ng.TaskManager():
        bem = ScalarBIESIBCSolver(mesh, order=1, use_intree_bem=True,
                                 assemble_dense=True, intree_geom_order=1,
                                 intree_singular_n_q=10,
                                 intree_regular_quad_degree=11,
                                 loop_work_backend=request.param)
    return bem, points, triangles, areas, normals, gradients


def incident(points):
    return .5*MU_0*np.c_[-points[:, 1], points[:, 0], np.zeros(len(points))]


@pytest.mark.parametrize("zero_mean", [False, True])
def test_complete_work_changes_covariantly_with_generic_lift(surface, zero_mean):
    bem, p, t, areas, normals, gradients = surface
    g = np.array([.3, -.2, .4, .7])
    q = np.array([.1, -.2, .3, -.1])
    w = np.array([.2, .1, -.3, .4])
    mass_row = bem.M.sum(axis=1)
    if zero_mean:
        w -= mass_row @ w/mass_row.sum()
    z = .03+.02j
    phi_inc = -p[:, 2].astype(complex)
    B = .5*bem.M-bem.DL
    shifted_q = q-np.linalg.solve(bem.SL, B @ w)
    with ng.TaskManager():
        old = _loop_work_row(bem, p, t, t, np.arange(4), gradients, normals,
                             areas, g, z, incident, q, phi_inc, quadrature_bonus=12)
        new = _loop_work_row(bem, p, t, t, np.arange(4), gradients, normals,
                             areas, g+w, z, incident, shifted_q, phi_inc,
                             quadrature_bonus=12)
        S = _p0_single_layer(bem, areas, quadrature_bonus=12)
    # Independently build ordinary exterior energy, not the mixed BIE matrix.
    C = np.zeros((4, 4, 3))
    for ti, triangle in enumerate(t):
        C[ti, triangle] = np.cross(normals[ti], -gradients[ti])
    with ng.TaskManager():
        Q = (sum(C[:, :, c].T @ S @ C[:, :, c] for c in range(3))
             + B.T @ np.linalg.solve(bem.SL, B))
    np.testing.assert_allclose(new[0]-old[0], z*bem.K @ w, atol=2e-14, rtol=1e-12)
    np.testing.assert_allclose(new[2]-old[2], MU_0*(w @ Q), atol=2e-15, rtol=1e-7)
    np.testing.assert_allclose(new[1]-old[1], 2*(old[0] @ w)+z*(w @ bem.K @ w),
                               atol=2e-14, rtol=1e-12)
    np.testing.assert_allclose(new[3]-old[3], 2*(old[2] @ w)+MU_0*(w @ Q @ w),
                               atol=2e-15, rtol=1e-7)
    # Closed-surface Stokes identity gives the incident work of the added lift.
    incident_w = MU_0*np.sum(areas*w[t].mean(axis=1)*normals[:, 2])
    np.testing.assert_allclose(new[4]-old[4],
                               incident_w-(new[2]-old[2]) @ phi_inc,
                               atol=2e-15, rtol=1e-7)
    for result in (old, new):
        assert result[5]["magnetic_work_symmetry"] < 1e-6
        assert result[5]["gradient_work_residual"] < 1e-12
        assert result[5]["conormal_flux_relative_mismatch"] < 1e-12
    # Every ordinary and gauge row transforms, including a nonzero mean lift.
    omega = 300.
    ordinary = B+z/(1j*omega*MU_0)*(bem.SL @ bem.M_inv @ bem.K)
    alpha_old = bem.SL @ (bem.M_inv @ old[0]/(1j*omega*MU_0)-q)
    alpha_new = bem.SL @ (bem.M_inv @ new[0]/(1j*omega*MU_0)-shifted_q)
    np.testing.assert_allclose(alpha_new-alpha_old, ordinary @ w, atol=1e-13)
    transform = np.eye(6, dtype=complex)
    transform[:4, 4] = w
    original = np.zeros((6, 6), dtype=complex)
    original[:4, :4] = ordinary
    original[:4, 4] = alpha_old
    original[:4, 5] = original[5, :4] = mass_row
    original[4, :4] = old[0]+1j*omega*old[2]
    original[4, 4] = old[1]+1j*omega*old[3]
    changed = original @ transform
    np.testing.assert_allclose(changed[:4, 4], alpha_new, atol=1e-13)
    assert changed[5, 4] == pytest.approx(mass_row @ w, abs=1e-14)
    np.testing.assert_array_equal(changed[:, 5], original[:, 5])


def test_stale_geometry_cache_is_rejected(surface):
    bem, _, _, areas, _, _ = surface
    previous = getattr(bem, "_loop_work_geometry_identity", None)
    try:
        bem._loop_work_geometry_identity = "stale"
        with pytest.raises(ValueError, match="geometry changed"):
            _p0_single_layer(bem, areas)
    finally:
        if previous is None:
            del bem._loop_work_geometry_identity
        else:
            bem._loop_work_geometry_identity = previous


def test_negative_exterior_work_is_rejected(surface, monkeypatch):
    from radia import bem_loop_work as module
    bem, p, t, areas, normals, gradients = surface
    with ng.TaskManager():
        positive = _p0_single_layer(bem, areas, quadrature_bonus=12)
    monkeypatch.setattr(module, "_p0_single_layer", lambda *a, **k: -positive)
    with pytest.raises(RuntimeError, match="magnetic work must be positive"):
        _loop_work_row(bem, p, t, t, np.arange(4), gradients, normals, areas,
                       np.array([.3, -.2, .4, .7]), .03+.02j, incident,
                       np.zeros(4), -p[:, 2])


@pytest.mark.parametrize("curved", [True, False])
def test_unsupported_geometry_rejected_before_assembly(curved):
    from types import SimpleNamespace
    from radia.bem_loop_extension import solve_loop_extended
    mesh = SimpleNamespace(GetCurveOrder=lambda: 2 if curved else 1,
                           deformation=None if curved else object())
    with pytest.raises(ValueError, match="undeformed flat"):
        solve_loop_extended(SimpleNamespace(mesh=mesh), np.zeros(4), .03+.02j,
                             300., incident)


def test_native_loop_operator_preserves_permutation_and_complex_products(surface):
    from ngsolve.bem import LaplaceSL
    from radia.bem_loop_work import _P0SingleLayer
    bem, _, _, _, _, _ = surface
    with ng.TaskManager():
        space=ng.SurfaceL2(bem.mesh,order=0,dual_mapping=False)
        u,v=space.TnT();d=ng.ds(bonus_intorder=12)
        native=LaplaceSL(u*d,use_fmm=False)*v*d
        rows,cols,values=native.mat.COO()
        from scipy.sparse import coo_matrix
        dense=coo_matrix((values,(rows,cols)),shape=(space.ndof,space.ndof)).toarray()
        permutation=np.arange(space.ndof)[::-1]
        expected=dense[np.ix_(permutation,permutation)]
        operator=_P0SingleLayer(native,space,permutation,dense.nbytes)
        vectors=np.c_[np.arange(space.ndof)+.2,np.arange(space.ndof)**2-.7].astype(complex)
        vectors+=1j*vectors[::-1]
        np.testing.assert_allclose(operator@vectors,expected@vectors,rtol=1e-13,atol=1e-15)
        np.testing.assert_allclose(operator.T@vectors,expected.T@vectors,rtol=1e-13,atol=1e-15)
        np.testing.assert_allclose(vectors.T@operator,vectors.T@expected,rtol=1e-13,atol=1e-15)
        assert operator.dense_array_bytes==0


def test_loop_conormal_diagnostic_rejects_nonmanifold_edge_counts(surface):
    bem,p,t,areas,normals,gradients=surface
    duplicate=t.copy();duplicate[1]=duplicate[0]
    with ng.TaskManager(),pytest.raises(ValueError,match='closed manifold'):
        _loop_work_row(bem,p,duplicate,t,np.arange(4),gradients,normals,areas,
            np.array([.3,-.2,.4,.7]),.03+.02j,incident,np.zeros(4),-p[:,2])


def test_loop_fmm_failure_is_explicit_without_fallback(surface, monkeypatch):
    from ngsolve import bem as native_bem
    bem, _, _, areas, _, _ = surface
    monkeypatch.setattr(bem, "loop_work_backend", "fmm")
    calls = []

    def unavailable(*args, **kwargs):
        calls.append(kwargs)
        raise AttributeError("FMM unavailable in this build")

    monkeypatch.setattr(native_bem, "LaplaceSL", unavailable)
    with ng.TaskManager(), pytest.raises(RuntimeError, match="--wp-loop-work-backend dense"):
        _p0_single_layer(bem, areas, quadrature_bonus=13)
    assert len(calls) == 1 and calls[0]["use_fmm"] is True


def test_loop_route_has_declared_threshold_and_rejects_invalid_option():
    from types import SimpleNamespace
    from radia.bem_loop_work import _loop_work_route
    bem = SimpleNamespace(loop_work_backend="auto")
    assert _loop_work_route(bem, 511) == "dense"
    assert _loop_work_route(bem, 512) == "fmm"
    bem.loop_work_backend = "invalid"
    with pytest.raises(ValueError, match="auto, dense, or fmm"):
        _loop_work_route(bem, 512)
