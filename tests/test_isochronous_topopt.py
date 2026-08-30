"""Adjoint-gradient and reciprocity locks for radia.isochronous_topopt.

Stage-1 gate of docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md: the
per-element density adjoint gradient must match central finite differences to
the 1e-6 class, and the dipole-array reciprocity load must reproduce the
independent C++ analytic charge evaluator.  Bands are set from measured values
with 2-3 decades of margin (promoted research run 2026-07-28:
directional 8.1e-10, per-element
worst 3.9e-7 at the FD noise floor, reciprocity 1.1e-10 at bonus_intorder=10).
"""
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")

from netgen.occ import OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import HDiv, InnerProduct, Mesh, SetNumThreads, TaskManager  # noqa: E402

from radia.ffag_topopt import FFAGCyclicDensityMap  # noqa: E402
from radia.isochronous_topopt import (  # noqa: E402
    CHI_MIN, MU0, DensityAdjointVIM, HelmholtzFilter, HeavisideProjection,
    _accept_deep_restoration, _solve_minimax_lp_update,
    demag_field_from_solution, demag_field_evaluator,
    density_gradient_from_s_gradient,
    density_discreteness, density_to_s, field_functional_load,
    gradient_pair_points, iron_only_verification_ready,
    iron_only_mesh, optimize_density, orbit_arc_points, uniform_field_load,
    verify_design_iron_only, track_sector_orbit, sector_linear_optics,
    combined_function_exit_metrics,
    combined_function_exit_metrics_from_field_response,
    combined_function_linear_optics,
    combined_function_transfer_map,
    combined_function_transfer_map_from_field_response,
    design_combined_function_transfer_map_target,
    design_achromatic_gradient_profile,
    straightened_bend_validation,
    isochronous_increment_targets, isochronous_total_field_bands,
    isochronous_profile_metrics, restore_projected_volume,
    transfer_map_reachability,
)


def test_combined_function_sector_map_and_analytic_gradient_chain():
    rho=2.3;theta=0.41
    optics=combined_function_linear_optics(
        [1.0/rho],[0.0],[rho*theta],gradient_jacobian=[[1.0]])
    np.testing.assert_allclose(
        optics.dispersion[-1],
        [rho*(1.0-np.cos(theta)),np.sin(theta)],rtol=2e-14,atol=2e-14)
    expected=np.array([[np.cos(theta),rho*np.sin(theta)],
                       [-np.sin(theta)/rho,np.cos(theta)]])
    np.testing.assert_allclose(optics.radial_matrix,expected,rtol=2e-14,atol=2e-14)
    assert optics.radial_stable and not optics.vertical_stable

    step=2e-6
    plus=combined_function_linear_optics([1/rho],[step],[rho*theta])
    minus=combined_function_linear_optics([1/rho],[-step],[rho*theta])
    fd=(plus.dispersion[-1]-minus.dispersion[-1])/(2*step)
    np.testing.assert_allclose(optics.endpoint_jacobian[:,0],fd,
                               rtol=2e-9,atol=2e-11)


def test_combined_function_exit_metrics_match_sector_and_score_drift():
    length = 8.0
    bend_angle = np.pi / 6.0
    curvature = bend_angle / length
    radius = 1.0 / curvature
    score_drift = 0.3

    metrics = combined_function_exit_metrics(
        [curvature], [0.0], [length],
        reference_curvature=curvature,
        downstream_drift=score_drift,
    )

    assert metrics.x0_m == pytest.approx(0.0, abs=1.0e-14)
    assert metrics.psi0_rad == pytest.approx(0.0, abs=1.0e-14)
    assert metrics.eta_m == pytest.approx(
        radius * (1.0 - np.cos(bend_angle))
        + score_drift * np.sin(bend_angle),
        abs=2.0e-14,
    )
    assert metrics.eta_prime_rad == pytest.approx(
        np.sin(bend_angle), abs=2.0e-14)


def test_positive_curvature_error_bends_central_orbit_opposite_dispersion():
    metrics = combined_function_exit_metrics(
        [0.08], [0.0], [0.2], reference_curvature=0.05)
    assert metrics.psi0_rad < 0.0
    assert metrics.x0_m < 0.0
    assert metrics.eta_prime_rad > 0.0
    assert metrics.eta_m > 0.0


def test_combined_function_four_response_jacobian_is_analytic():
    curvature = np.array([0.05, 0.07, 0.06])
    gradient = np.array([-0.01, 0.02, -0.015])
    lengths = np.array([0.7, 1.1, 0.9])
    reference = np.full(3, 0.06)
    curvature_jacobian = np.array([
        [0.02, -0.01],
        [-0.03, 0.015],
        [0.01, 0.025],
    ])
    gradient_jacobian = np.array([
        [0.08, -0.04],
        [-0.02, 0.06],
        [0.03, 0.01],
    ])
    metrics = combined_function_exit_metrics(
        curvature, gradient, lengths,
        reference_curvature=reference,
        downstream_drift=0.3,
        curvature_jacobian=curvature_jacobian,
        gradient_jacobian=gradient_jacobian,
    )

    # Finite differences are used only here as a regression check of the
    # production Frechet derivative, never by the optimizer.
    step = 1.0e-6
    finite_difference = np.empty_like(metrics.response_jacobian)
    for parameter in range(curvature_jacobian.shape[1]):
        plus = combined_function_exit_metrics(
            curvature + step * curvature_jacobian[:, parameter],
            gradient + step * gradient_jacobian[:, parameter],
            lengths,
            reference_curvature=reference,
            downstream_drift=0.3,
        )
        minus = combined_function_exit_metrics(
            curvature - step * curvature_jacobian[:, parameter],
            gradient - step * gradient_jacobian[:, parameter],
            lengths,
            reference_curvature=reference,
            downstream_drift=0.3,
        )
        finite_difference[:, parameter] = (
            plus.response - minus.response) / (2.0 * step)

    np.testing.assert_allclose(
        metrics.response_jacobian,
        finite_difference,
        rtol=2.0e-8,
        atol=2.0e-10,
    )


def test_hdiv_field_response_fuses_into_four_optics_rows_with_explicit_signs():
    lengths = np.array([0.4, 0.6])
    rigidity = 2.0
    bz = np.array([-0.10, -0.14])
    gradient = np.array([3.0, -2.0])
    field_jacobian = np.array([
        [0.2, -0.1],
        [0.1, 0.3],
        [1.0, -2.0],
        [-0.5, 0.4],
    ])
    fused = combined_function_exit_metrics_from_field_response(
        np.r_[bz, gradient], lengths, rigidity,
        reference_curvature=0.06,
        downstream_drift=0.2,
        field_response_jacobian=field_jacobian,
        curvature_sign=-1.0,
        gradient_sign=-1.0,
    )
    direct = combined_function_exit_metrics(
        -bz / rigidity, -gradient / rigidity, lengths,
        reference_curvature=0.06,
        downstream_drift=0.2,
        curvature_jacobian=-field_jacobian[:2] / rigidity,
        gradient_jacobian=-field_jacobian[2:] / rigidity,
    )
    np.testing.assert_allclose(fused.response, direct.response)
    np.testing.assert_allclose(
        fused.response_jacobian, direct.response_jacobian)


def test_combined_function_transfer_map_matches_uniform_sector():
    radius = 2.3
    angle = 0.41
    length = radius * angle
    result = combined_function_transfer_map(
        [1.0 / radius], [0.0], [length])
    matrix = result.matrix
    expected_radial = np.array([
        [np.cos(angle), radius * np.sin(angle)],
        [-np.sin(angle) / radius, np.cos(angle)],
    ])
    np.testing.assert_allclose(matrix[:2, :2], expected_radial,
                               rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(
        matrix[:2, 5],
        [radius * (1.0 - np.cos(angle)), np.sin(angle)],
        rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(
        matrix[2:4, 2:4], [[1.0, length], [0.0, 1.0]],
        rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(
        matrix[4, [0, 1, 5]],
        [np.sin(angle),
         radius * (1.0 - np.cos(angle)),
         radius * (angle - np.sin(angle))],
        rtol=2e-14, atol=2e-14)
    assert matrix.flags.c_contiguous
    assert result.response.shape == (13,)
    assert result.response_jacobian.shape == (13, 0)


def test_transfer_map_field_chain_uses_forward_ad_jacobian():
    lengths = np.array([0.4, 0.6, 0.5])
    rigidity = 2.0
    field = np.array([-0.10, -0.14, -0.11, 3.0, -2.0, 0.5])
    field_jacobian = np.array([
        [0.2, -0.1], [0.1, 0.3], [-0.05, 0.12],
        [1.0, -2.0], [-0.5, 0.4], [0.3, 0.8],
    ])
    length_jacobian = np.array([
        [0.03, -0.02], [-0.01, 0.04], [0.02, 0.01],
    ])
    differentiated = combined_function_transfer_map_from_field_response(
        field, lengths, rigidity,
        field_response_jacobian=field_jacobian,
        segment_length_jacobian=length_jacobian,
        curvature_sign=-1.0, gradient_sign=-1.0)
    assert differentiated.derivative_backend==(
        "forward-mode-expm-frechet-ad")

    # Finite differences are a regression oracle only; production topology
    # sensitivities are forward-mode AD tangents. The matrix exponential is
    # differentiated by its exact Frechet derivative primitive.
    step = 1.0e-6
    finite_difference = np.empty_like(differentiated.response_jacobian)
    for parameter in range(field_jacobian.shape[1]):
        plus = combined_function_transfer_map_from_field_response(
            field + step * field_jacobian[:, parameter],
            lengths + step * length_jacobian[:, parameter], rigidity,
            curvature_sign=-1.0, gradient_sign=-1.0)
        minus = combined_function_transfer_map_from_field_response(
            field - step * field_jacobian[:, parameter],
            lengths - step * length_jacobian[:, parameter], rigidity,
            curvature_sign=-1.0, gradient_sign=-1.0)
        finite_difference[:, parameter] = (
            plus.response - minus.response) / (2.0 * step)
    np.testing.assert_allclose(
        differentiated.response_jacobian, finite_difference,
        rtol=3.0e-8, atol=3.0e-10)


def test_transfer_map_reachability_rejects_uncontrolled_matrix_entry():
    failed = transfer_map_reachability(
        [0.0, 0.0], [[1.0, 0.0], [0.0, 0.0]],
        [1.0, 1.0], [1.0, 1.0], acceptance_ratio=0.5)
    assert failed.numerical_rank == 1
    np.testing.assert_allclose(failed.predicted_response, [1.0, 0.0])
    np.testing.assert_allclose(failed.residual, [0.0, -1.0])
    assert not failed.reachable

    passed = transfer_map_reachability(
        [0.0, 0.0], [[1.0], [0.0]], [1.0, 0.0], [1.0, 1.0])
    assert passed.reachable


def test_transfer_map_reachability_aca_qr_tsvd_matches_dense_oracle():
    row = np.array([1.0, -0.4, 0.2, 0.7, -0.3])
    left = np.array([1.0, 0.5, -0.25, 0.75, 0.1, -0.6])
    second = np.array([-0.2, 0.8, 0.3, -0.5, 0.4])
    matrix = np.outer(left, row) + np.outer(
        np.array([0.3, -0.1, 0.9, 0.2, -0.7, 0.4]), second)
    target_step = np.array([0.2, -0.15, 0.1, 0.05, -0.08])
    target = matrix @ target_step
    dense = transfer_map_reachability(
        np.zeros(6), matrix, target, np.ones(6),
        relative_tolerance=1e-11, method="dense")
    compressed = transfer_map_reachability(
        np.zeros(6), matrix, target, np.ones(6),
        relative_tolerance=1e-9, method="aca_qr_tsvd",
        aca_tolerance=1e-13)
    assert dense.numerical_rank == compressed.numerical_rank == 2
    assert compressed.factorization_method == "aca_qr_tsvd"
    assert compressed.aca_rank == 2
    assert compressed.reachable
    np.testing.assert_allclose(
        compressed.predicted_response, dense.predicted_response,
        rtol=2e-11, atol=2e-12)
    np.testing.assert_allclose(
        compressed.parameter_step, dense.parameter_step,
        rtol=2e-11, atol=2e-12)


def test_transfer_map_reachability_completes_aca_pivot_breakdown(monkeypatch):
    import radia.stream_function as stream_function

    matrix = np.array([
        [4.0, 0.0, 1.0], [0.0, 2.0, -0.5],
        [1.0, -1.0, 0.25], [0.5, 0.75, -0.0625],
    ])
    exact_u, exact_s, exact_vh = np.linalg.svd(matrix, full_matrices=False)

    def rank_one_seed(*_args, **_kwargs):
        return stream_function.StreamTSVD(
            U=exact_u[:, :1], S=exact_s[:1], V=exact_vh[:1].T,
            k_aca=1, method="aca_qr_tsvd")

    monkeypatch.setattr(stream_function, "aca_tsvd", rank_one_seed)
    target = matrix@np.array([0.2, -0.1, 0.05])
    completed = transfer_map_reachability(
        np.zeros(4), matrix, target, np.ones(4),
        relative_tolerance=1.0e-11, method="aca_qr_tsvd",
        aca_tolerance=1.0e-13)
    dense = transfer_map_reachability(
        np.zeros(4), matrix, target, np.ones(4),
        relative_tolerance=1.0e-11, method="dense")

    assert completed.aca_seed_rank == 1
    assert completed.aca_residual_completion_rank == 2
    assert completed.aca_rank == completed.numerical_rank == 3
    assert completed.factorization_relative_error < 1.0e-14
    np.testing.assert_allclose(
        completed.predicted_response, dense.predicted_response,
        rtol=2.0e-13, atol=2.0e-14)


def test_ideal_optics_calculates_realisable_achromatic_transfer_target():
    radial = np.array([[0.0, -3.0], [1.0 / 3.0, 0.0]])
    vertical = np.array([[0.0, 4.0], [-0.25, 0.0]])
    design = design_combined_function_transfer_map_target(
        length=8.0, bend_angle=np.pi / 6.0,
        radial_matrix=radial, vertical_matrix=vertical,
        n_segments=16, normalized_gradient_limit=12.0)
    np.testing.assert_allclose(
        design.transfer_map.matrix[:2, :2], radial,
        rtol=0.0, atol=2e-10)
    np.testing.assert_allclose(
        design.transfer_map.matrix[2:4, 2:4], vertical,
        rtol=0.0, atol=2e-10)
    np.testing.assert_allclose(
        design.transfer_map.matrix[:2, 5], 0.0,
        rtol=0.0, atol=2e-11)
    np.testing.assert_allclose(
        design.transfer_map.matrix[4, :2], 0.0,
        rtol=0.0, atol=2e-11)
    assert design.matrix[4, 5] > 0.0
    assert design.maximum_scaled_residual < 1e-9
    assert design.transfer_map.optics.radial_stable
    assert design.transfer_map.optics.vertical_stable


def test_analytic_achromatic_gradient_design_has_zero_endpoint_eta():
    design=design_achromatic_gradient_profile(
        length=4.0,bend_angle=np.pi/6,n_segments=4,
        normalized_gradient_limit=5.0)
    optics=design.optics
    assert optics.dispersion[0,0] == 0.0
    assert abs(optics.dispersion[-1,0]) < 1e-10
    assert optics.radial_stable and optics.vertical_stable
    assert np.all(np.abs(design.normalized_gradient)<=5.0+1e-12)
    np.testing.assert_allclose(
        design.curvature@design.segment_lengths,np.pi/6,rtol=0,atol=1e-14)


def test_analytic_achromatic_gradient_design_supports_reachable_sign_branch():
    design=design_achromatic_gradient_profile(
        length=4.0,bend_angle=np.pi/6,n_segments=4,
        normalized_gradient_limit=5.0,
        initial_normalized_gradient=[-2.3,2.85,.5,1.9],
        gradient_sign_pattern=[-1,1,1,1])
    assert abs(design.optics.dispersion[-1,0])<1e-10
    assert design.optics.radial_stable and design.optics.vertical_stable
    assert design.normalized_gradient[0]<=0.0
    assert np.all(design.normalized_gradient[1:]>=0.0)


def test_iron_only_mesh_extracts_hex_and_names_shape_boundaries():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    keep=np.array([True,False])
    extracted=iron_only_mesh(
        mesh,keep,boundary_classifier=lambda center,normal:
        "pole" if normal[0]>0.9 else "fixed")
    assert extracted.ne==1
    assert set(extracted.GetBoundaries())=={"fixed","pole"}
    assert len(tuple(extracted.Elements(ng.VOL))[0].vertices)==8


def test_iron_only_mesh_can_handoff_active_hexes_as_conforming_tets():
    from ngsolve.meshes import MakeStructured3DMesh
    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    extracted=iron_only_mesh(
        mesh,np.array([True,False]),tetrahedralize_hex=True,
        boundary_classifier=lambda center,normal:
        "pole" if normal[0]>0.9 else "fixed")
    elements=tuple(extracted.Elements(ng.VOL))
    assert len(elements)==6
    assert {len(element.vertices) for element in elements}=={4}
    assert set(extracted.GetBoundaries())=={"fixed","pole"}
    np.testing.assert_allclose(
        ng.Integrate(1.0,extracted),0.5,rtol=0.0,atol=2e-14)


def test_heaviside_projection_chain_and_grayness_gate():
    projection = HeavisideProjection(beta=4.0, eta=0.45)
    rho = np.array([0.05, 0.3, 0.45, 0.7, 0.95])
    gradient = np.array([1.0, -2.0, 0.5, 0.3, -0.7])
    step = 1e-7
    direction = np.array([0.2, -0.1, 0.3, -0.2, 0.1])
    analytic = projection.chain(rho, gradient)@direction
    fd = gradient@(projection.apply(rho+step*direction)
                   - projection.apply(rho-step*direction))/(2*step)
    np.testing.assert_allclose(analytic, fd, rtol=2e-8, atol=2e-10)
    projected = projection.apply(rho)
    assert np.all(np.diff(projected)>0) and np.all((projected>=0)&(projected<=1))
    ready, metrics = iron_only_verification_ready([0., 0.02, .98, 1.])
    assert ready and metrics["intermediate_fraction"] == 0.0
    ready, metrics = iron_only_verification_ready(np.full(20, .5))
    assert not ready and metrics == density_discreteness(np.full(20, .5))


@pytest.fixture(scope="module")
def problem():
    SetNumThreads(4)
    with TaskManager():
        mesh = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.45))
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        rng = np.random.default_rng(1234)
        s0 = np.exp(rng.uniform(np.log(1e-2), 0.0, prob.n_el))
        theta = np.linspace(0.0, 2.0 * pi, 6, endpoint=False)
        orbit = np.stack([1.6 * np.cos(theta), 1.6 * np.sin(theta),
                          0.3 * np.ones_like(theta)], axis=1)
        w_orbit = np.array([1.0, -0.7, 1.3, 0.5, -1.1, 0.9])
        pts, wts = gradient_pair_points(orbit, w_orbit, delta=0.1, axis=0)
        f_state = uniform_field_load(fes, (0.0, 0.0, 1.0))
        f_adj = field_functional_load(fes, pts, wts, axis=2, scale=MU0,
                                      bonus_intorder=10)
        base = prob.objective_and_gradient(s0, f_state, f_adj)
    return SimpleNamespace(mesh=mesh, fes=fes, prob=prob, s0=s0,
                           dipole_points=pts, dipole_weights=wts,
                           f_state=f_state, f_adj=f_adj, base=base, rng=rng)


@pytest.fixture(scope="module")
def broken_problem(problem):
    with TaskManager():
        fes = ng.HDiv(problem.mesh, order=1, discontinuous=True)
        prob = DensityAdjointVIM(
            fes, eps=1e-7, internal_interfaces=True)
    return SimpleNamespace(mesh=problem.mesh, fes=fes, prob=prob)


def _objective_at(problem, s):
    with TaskManager():
        gf, _ = problem.prob.solve(s, problem.f_state)
    return float(InnerProduct(problem.f_adj.vec, gf.vec))


def test_adjoint_gradient_directional_fd(problem):
    """All-element check: grad . ds == central FD along a random direction.

    h = 1e-4: the FD signal |2 h grad.ds| must sit well above the CG solution
    noise (tol * J ~ 1e-19); at h = 1e-5 this direction's derivative is
    cancellation-small and the comparison hits the solver noise floor.
    """
    d = np.random.default_rng(7).standard_normal(problem.prob.n_el)
    ds = problem.s0 * d              # multiplicative direction keeps s > 0
    predicted = float(problem.base.gradient @ ds)
    h = 1e-4
    fd = (_objective_at(problem, problem.s0 + h * ds)
          - _objective_at(problem, problem.s0 - h * ds)) / (2.0 * h)
    rel = abs(fd - predicted) / abs(predicted)
    assert rel < 1e-6, (fd, predicted, rel)


def test_adjoint_gradient_per_element_fd(problem):
    """Per-element FD on the strongest-gradient element plus random picks."""
    g = problem.base.gradient
    picks = {int(np.abs(g).argmax())}
    rng = np.random.default_rng(21)
    while len(picks) < 4:
        picks.add(int(rng.integers(0, problem.prob.n_el)))
    for e in sorted(picks):
        h = 1e-4 * problem.s0[e]
        sp = problem.s0.copy(); sp[e] += h
        sm = problem.s0.copy(); sm[e] -= h
        fd = (_objective_at(problem, sp) - _objective_at(problem, sm)) / (2.0 * h)
        rel = abs(fd - g[e]) / max(abs(g[e]), abs(fd))
        assert rel < 1e-5, (e, g[e], fd, rel)


def test_reciprocity_matches_cpp_field_evaluator(problem):
    """f_adj^T m == MU0 * sum_i w_i H_d[m]_z(x_i) via the independent C++ path."""
    with TaskManager():
        H_d = demag_field_from_solution(problem.prob.demag, problem.base.gfM,
                                        problem.dipole_points)
    independent = float(MU0 * np.dot(problem.dipole_weights, H_d[:, 2]))
    rel = abs(problem.base.objective - independent) / abs(independent)
    assert rel < 1e-7, (problem.base.objective, independent, rel)


def test_reusable_demag_field_evaluator_matches_batch_path(problem):
    points = problem.dipole_points[:3]
    reusable = demag_field_evaluator(problem.prob.demag, problem.base.gfM)
    batch = demag_field_from_solution(problem.prob.demag, problem.base.gfM, points)
    np.testing.assert_allclose(reusable(points), batch, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(reusable(points[0]), batch[0], rtol=0.0, atol=0.0)


def test_native_flat_tet_field_functional_rows_match_exact_batch(broken_problem):
    """Sparse native observation rows reuse the exact analytic TET/TRI kernel."""
    rng=np.random.default_rng(20260730)
    coefficients=rng.normal(size=broken_problem.fes.ndof)
    gf=ng.GridFunction(broken_problem.fes)
    gf.vec.FV().NumPy()[:]=coefficients
    points=np.array([[1.45,0.13,0.17],[-1.37,0.22,-0.19],
                     [0.18,1.51,0.23]],dtype=float)
    weights=rng.normal(size=(3,len(points),3))
    with TaskManager():
        rows=np.asarray(
            broken_problem.prob.demag._G.configured_field_functional_rows(
                points,weights),dtype=float)
        field=demag_field_from_solution(
            broken_problem.prob.demag,gf,points)
    expected=np.einsum("rpc,pc->r",weights,field)
    assert rows.shape==(3,broken_problem.fes.ndof)
    assert rows.flags.c_contiguous
    np.testing.assert_allclose(rows@coefficients,expected,
                               rtol=2e-13,atol=2e-13)


def test_native_affine_hex_bdm1_field_functional_rows_match_exact_batch():
    """Affine HEX BDM1 rows and direct fields share the analytic TET/TRI source."""
    from ngsolve.meshes import MakeStructured3DMesh

    mesh=MakeStructured3DMesh(
        hexes=True,nx=2,ny=1,nz=1,
        mapping=lambda x,y,z:(1.3*x-.2*y,.7*y+.1*z,.45*z))
    fes=ng.HDiv(mesh,order=1,discontinuous=True)
    points=np.array([
        [.41,.27,.451],
        [1.05,.18,.49],
        [-.08,.33,.22],
        [.62,-.07,.12],
    ],dtype=float)
    rng=np.random.default_rng(20260803)
    weights=rng.normal(size=(3,len(points),3))
    coefficients=rng.normal(size=fes.ndof)
    with TaskManager():
        problem=DensityAdjointVIM(
            fes,eps=1e-7,internal_interfaces=True)
        rows=np.asarray(
            problem.demag._G.configured_field_functional_rows(points,weights))
        evaluator=problem.demag._G.create_field_evaluator(coefficients)
        field=np.asarray(evaluator.field(points,"direct"))/(4.0*np.pi)
    expected=np.einsum("rpc,pc->r",weights,field)
    np.testing.assert_allclose(rows@coefficients,expected,
                               rtol=3e-12,atol=3e-12)
    assert dict(evaluator.stats())["source_representation"]=="analytic-tet"

    # Independent NGSolve/Piola reciprocity is reliable once targets are well
    # separated from the body.  It checks the affine trilinear-to-physical
    # cubic conversion, not merely agreement between two users of that source.
    far_points=np.array([[.2,.1,2.8],[1.1,.6,3.2],[-.7,.4,2.5]])
    far_weights=rng.normal(size=(len(far_points),3))
    with TaskManager():
        far_row=np.asarray(problem.demag._G.configured_field_functional_rows(
            far_points,far_weights[None,:,:]))[0]
        reciprocal=np.zeros(fes.ndof)
        for axis in range(3):
            reciprocal+=field_functional_load(
                fes,far_points,far_weights[:,axis],axis=axis,scale=1.0,
                bonus_intorder=14).vec.FV().NumPy()
    np.testing.assert_allclose(far_row,reciprocal,rtol=3e-9,atol=3e-11)


def test_native_affine_hex_bdm2_field_functional_rows_match_exact_batch():
    """Affine HEX BDM2 rows retain every Q2 volume and facet charge mode."""
    from ngsolve.meshes import MakeStructured3DMesh

    mesh=MakeStructured3DMesh(
        hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(1.1*x-.15*y,.8*y+.08*z,.55*z))
    fes=ng.HDiv(mesh,order=2,discontinuous=True)
    points=np.array([
        [.31,.17,.73],
        [1.04,.42,.61],
        [-.22,.29,.36],
    ],dtype=float)
    rng=np.random.default_rng(20260830)
    weights=rng.normal(size=(3,len(points),3))
    coefficients=rng.normal(size=fes.ndof)
    with TaskManager():
        problem=DensityAdjointVIM(
            fes,eps=1e-12,internal_interfaces=True)
        rows=np.asarray(
            problem.demag._G.configured_field_functional_rows(points,weights))
        evaluator=problem.demag._G.create_field_evaluator(coefficients)
        field=np.asarray(evaluator.field(points,"direct"))/(4.0*np.pi)
    expected=np.einsum("rpc,pc->r",weights,field)
    assert rows.shape==(3,fes.ndof)
    assert rows.flags.c_contiguous
    np.testing.assert_allclose(rows@coefficients,expected,
                               rtol=2e-11,atol=2e-11)
    assert dict(evaluator.stats())["source_representation"]=="analytic-tet"

    far_points=np.array([[.2,.1,2.8],[1.1,.6,3.2],[-.7,.4,2.5]])
    far_weights=rng.normal(size=(len(far_points),3))
    with TaskManager():
        far_row=np.asarray(problem.demag._G.configured_field_functional_rows(
            far_points,far_weights[None,:,:]))[0]
        reciprocal=np.zeros(fes.ndof)
        for axis in range(3):
            reciprocal+=field_functional_load(
                fes,far_points,far_weights[:,axis],axis=axis,scale=1.0,
                bonus_intorder=20).vec.FV().NumPy()
    np.testing.assert_allclose(far_row,reciprocal,rtol=3e-8,atol=3e-10)


def test_native_flat_tet_field_rows_are_finite_on_coplanar_panel_extension():
    """A point in a face plane but outside the face has a finite field.

    A topology removal can place a new pole end face at the same longitudinal
    station as an orbit sample.  The quadratic surface-charge kernel must use
    the analytic coplanar exterior limit instead of forming ``0/0`` from its
    normal inverse-cube moment.
    """
    from ngsolve.meshes import MakeStructured3DMesh

    parent=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    mesh=iron_only_mesh(
        parent,np.array([True,False]),tetrahedralize_hex=True)
    fes=ng.HDiv(mesh,order=1)
    point=np.array([[0.25,-0.25,0.0]],dtype=float)
    weights=np.array([[[0.0,0.0,1.0]],[[0.3,-0.2,0.4]]],dtype=float)
    epsilon=1.0e-7
    with TaskManager():
        prob=DensityAdjointVIM(fes,eps=1e-7)
        rows=np.asarray(
            prob.demag._G.configured_field_functional_rows(point,weights))
        plus=np.asarray(prob.demag._G.configured_field_functional_rows(
            point+np.array([[0.0,0.0,epsilon]]),weights))
        minus=np.asarray(prob.demag._G.configured_field_functional_rows(
            point-np.array([[0.0,0.0,epsilon]]),weights))
    assert rows.flags.c_contiguous and np.all(np.isfinite(rows))
    np.testing.assert_allclose(
        rows,0.5*(plus+minus),rtol=2e-8,atol=2e-10)


def test_warm_start_reproduces_and_saves_iterations(problem):
    """Warm-started re-solve at the same design matches and converges fast."""
    with TaskManager():
        again = problem.prob.objective_and_gradient(
            problem.s0, problem.f_state, problem.f_adj, warm=problem.base)
    rel = abs(again.objective - problem.base.objective) / abs(problem.base.objective)
    assert rel < 1e-9
    assert again.state_iterations <= max(3, problem.base.state_iterations // 5)


def test_native_jacobi_solver_satisfies_true_residual(problem):
    """Native convergence is checked against b-Ax, not its CG recurrence."""
    with TaskManager():
        gf, _ = problem.prob.solve(
            problem.s0, problem.f_state, tol=1e-10,
            solver="native-jacobi")
        _, operator, _ = problem.prob._system(problem.s0)
        residual = problem.f_state.vec.CreateVector()
        residual.data = problem.f_state.vec - operator * gf.vec
        relative = ng.Norm(residual) / ng.Norm(problem.f_state.vec)
    assert relative < 1.1e-10


def test_native_batched_multi_rhs_is_row_major_and_true_residual(problem):
    """Batched H-leaf traversal changes execution, not the solution."""
    with TaskManager():
        mass, operator, _ = problem.prob._system(problem.s0)
        gram = problem.prob.demag._G
        gram.configure_mass_matrix_ngsolve(mass.mat)
        rhs = np.ascontiguousarray(np.stack([
            problem.f_state.vec.FV().NumPy(),
            problem.f_adj.vec.FV().NumPy(),
        ]), dtype=float)
        result = gram.solve_configured_linear_material_auto_prec_many(
            1.0, rhs, 1e-10, 20000)
        solutions = np.asarray(result["m"], dtype=float)
        assert solutions.shape == rhs.shape and solutions.flags.c_contiguous
        assert len(result["iters"]) == 2
        assert int(result["coarse_dim"]) == 0
        assert int(result["recycle_dim"]) == 0
        for row, load in zip(solutions, [problem.f_state, problem.f_adj]):
            gf = ng.GridFunction(problem.fes)
            gf.vec.FV().NumPy()[:] = row
            residual = load.vec.CreateVector()
            residual.data = load.vec - operator * gf.vec
            relative = ng.Norm(residual) / ng.Norm(load.vec)
            assert relative < 1.1e-10
        reference, _ = problem.prob.solve(
            problem.s0, problem.f_state, tol=1e-10,
            solver="native-jacobi")
        reference_array = np.asarray(reference.vec.FV().NumPy(), dtype=float)
        relative_solution = np.linalg.norm(solutions[0] - reference_array) / np.linalg.norm(reference_array)
        assert relative_solution < 1e-9

        # The measured CPU default stays cluster-free, but the explicit
        # cluster-tree contraction/deflation surface remains numerically
        # locked for accelerator or future GPU backends.
        clustered = gram.solve_configured_linear_material_auto_prec_many(
            1.0, rhs, 1e-10, 20000, 8, 2, 1)
        assert int(clustered["coarse_dim"]) == 2
        assert int(clustered["recycle_dim"]) == 1
        for row, load in zip(np.asarray(clustered["m"]),
                             [problem.f_state, problem.f_adj]):
            gf = ng.GridFunction(problem.fes)
            gf.vec.FV().NumPy()[:] = row
            residual = load.vec.CreateVector()
            residual.data = load.vec - operator * gf.vec
            assert ng.Norm(residual) / ng.Norm(load.vec) < 1.1e-10

        # The multi-RHS mass-Riesz path retains independent CG recurrences
        # while sharing the expensive row-major operator/preconditioner
        # traversals and preserving the same true-residual contract.
        mass_batch = gram.solve_configured_linear_material_auto_prec_many(
            1.0, rhs, 1e-10, 20000, mass_riesz=True)
        assert np.asarray(mass_batch["m"]).flags.c_contiguous
        for row, load in zip(np.asarray(mass_batch["m"]),
                             [problem.f_state, problem.f_adj]):
            gf = ng.GridFunction(problem.fes)
            gf.vec.FV().NumPy()[:] = row
            residual = load.vec.CreateVector()
            residual.data = load.vec - operator * gf.vec
            assert ng.Norm(residual) / ng.Norm(load.vec) < 1.1e-10


def test_dipole_point_inside_mesh_raises(problem):
    with pytest.raises(ValueError, match="INSIDE"):
        field_functional_load(problem.fes, [[0.0, 0.0, 0.0]], [1.0])


def test_wrong_design_vector_raises(problem):
    with pytest.raises(ValueError, match="entries"):
        problem.prob.solve(np.ones(3), problem.f_state)
    with pytest.raises(ValueError, match="positive"):
        problem.prob.solve(np.zeros(problem.prob.n_el), problem.f_state)


def test_broken_hdiv_interface_charge_preserves_uniform_field(
        problem, broken_problem):
    """All internal jump rows cancel for a globally uniform magnetization."""
    with TaskManager():
        gf = ng.GridFunction(broken_problem.fes)
        gf.Set(ng.CoefficientFunction((0.0, 0.0, 1.0)))
        charges = np.asarray(
            broken_problem.prob.demag._B @ gf.vec.FV().NumPy(), dtype=float)
        face_charges = charges[problem.mesh.ne:].reshape(
            problem.mesh.nfacet, 3)
        owners = {facet.nr: 0 for facet in problem.mesh.facets}
        for el in problem.mesh.Elements(ng.VOL):
            for facet in el.facets:
                owners[facet.nr] += 1
        internal = np.array(
            [owners[facet.nr] == 2 for facet in problem.mesh.facets])
        assert np.max(np.abs(face_charges[internal])) < 1e-10
        reference = problem.prob.demag.DemagFactor(
            ng.CoefficientFunction((0.0, 0.0, 1.0)))
        measured = broken_problem.prob.demag.DemagFactor(
            ng.CoefficientFunction((0.0, 0.0, 1.0)))
    # The interface-charge MECHANISM is pinned by the per-face 1e-10
    # assert above; this end-to-end line compares two DIFFERENT solver
    # pipelines (conforming vs broken spaces), which agree to solver
    # tolerance, not machine epsilon.  The sigma-normalized charge-Gram
    # storage (133537e09) legitimately reshuffled the roundoff path and
    # moved the gap to 1.48e-9 (deterministic, identical digits on the
    # CI runner and LAB), so the old 1e-9 was calibration, not physics:
    # 5e-9 absolute (~1.5e-8 relative) still fails loudly for any real
    # interface-charge defect, which shows up orders of magnitude higher.
    assert abs(measured - reference) < 5e-9


def test_internal_interface_charge_rejects_conforming_hdiv(problem):
    with TaskManager(), pytest.raises(ValueError, match="discontinuous=True"):
        DensityAdjointVIM(
            problem.fes, eps=1e-7, internal_interfaces=True)


def test_broken_rt0_is_reduced_complete_topology_space(problem):
    """RT0 keeps one cell-divergence and one normal-jump mode per facet."""
    with TaskManager():
        fes = ng.HDiv(problem.mesh, order=0, discontinuous=True)
        reduced = DensityAdjointVIM(
            fes, eps=1e-7, internal_interfaces=True)
        assert fes.ndof == 4 * problem.mesh.ne
        assert reduced.demag._B.shape == (
            problem.mesh.ne + problem.mesh.nfacet, fes.ndof)
        factor = reduced.demag.DemagFactor(
            ng.CoefficientFunction((0.0, 0.0, 1.0)))
    assert 0.32 < factor < 0.35


def test_density_mapping_floor_and_chain_rule():
    chi_iron = 100.0
    s = density_to_s(np.array([0.0, 1.0, 0.5]), chi_iron)
    assert s[0] == pytest.approx(1.0 / CHI_MIN)
    assert s[1] == pytest.approx(1.0 / chi_iron)
    rho = np.array([0.1, 0.4, 0.9])
    grad_s = np.array([2.0, -1.0, 0.5])
    h = 1e-7
    fd = (grad_s * (density_to_s(rho + h, chi_iron)
                    - density_to_s(rho - h, chi_iron)) / (2.0 * h))
    chain = density_gradient_from_s_gradient(rho, grad_s, chi_iron)
    assert np.allclose(chain, fd, rtol=1e-5)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        density_to_s(np.array([1.5]), chi_iron)


def test_gradient_pair_points_contract():
    pts, wts = gradient_pair_points([[1.0, 2.0, 3.0]], [4.0], delta=0.2, axis=1)
    assert pts.shape == (2, 3) and wts.shape == (2,)
    np.testing.assert_allclose(pts[0], [1.0, 2.1, 3.0])
    np.testing.assert_allclose(pts[1], [1.0, 1.9, 3.0])
    np.testing.assert_allclose(wts, [20.0, -20.0])


def test_element_volumes_sum_to_mesh_volume(problem):
    with TaskManager():
        direct = float(ng.Integrate(ng.CoefficientFunction(1.0), problem.mesh))
    total = float(problem.prob.element_volumes.sum())
    assert abs(total - direct) / direct < 1e-12
    deficit = (4.0 * pi / 3.0 - total) / (4.0 * pi / 3.0)
    assert 0.0 < deficit < 0.12    # measured 7.0% straight-facet deficit at maxh=0.45


# ---------------------------------------------------------------- Stage 2
def test_helmholtz_filter_transpose_and_invariants(problem):
    with TaskManager():
        filt = HelmholtzFilter(problem.mesh, radius=0.15)
        n_el = problem.prob.n_el
        rng = np.random.default_rng(3)
        rho = rng.uniform(0.2, 0.8, n_el)
        g = rng.standard_normal(n_el)
        # constants are exact fixed points (pure-Neumann Helmholtz)
        np.testing.assert_allclose(filt.apply(np.ones(n_el)), 1.0, atol=1e-9)
        # chain == transpose of apply: FD of g . apply(rho) per component
        chain = filt.chain(g)
        h = 1e-6
        for e in [0, n_el // 2, n_el - 1]:
            rp = rho.copy(); rp[e] += h
            rm = rho.copy(); rm[e] -= h
            fd = (g @ filt.apply(rp) - g @ filt.apply(rm)) / (2.0 * h)
            assert abs(fd - chain[e]) <= 1e-8 * max(1.0, abs(chain[e]))
        # smoothing: a spike loses amplitude; the P1 Helmholtz realization
        # undershoots slightly (measured -1.0e-3 here) -- the design loop
        # clips the filtered density with the piecewise-exact chain rule.
        spike = np.zeros(n_el); spike[n_el // 3] = 1.0
        smoothed = filt.apply(spike)
        assert smoothed.max() < 0.9 and smoothed.min() > -5e-3


def test_native_linearize_matches_ngsolve_reference(problem):
    with TaskManager():
        cpts, _ = orbit_arc_points(1.4, -0.3, 6)
        f_con = field_functional_load(problem.fes, cpts,
                                      np.full(len(cpts), 1.0 / len(cpts)),
                                      axis=2, scale=MU0, bonus_intorder=10)
        s = problem.s0
        lin = problem.prob.linearize(s, problem.f_state,
                                     [problem.f_adj, f_con])
        single = problem.prob.objective_and_gradient(
            s, problem.f_state, f_con, solver="ngsolve-cg")
    assert lin.values[0] == pytest.approx(problem.base.objective, rel=1e-10)
    np.testing.assert_allclose(lin.jacobians[0], problem.base.gradient,
                               rtol=1e-8, atol=1e-18)
    assert lin.values[1] == pytest.approx(single.objective, rel=1e-10)
    np.testing.assert_allclose(lin.jacobians[1], single.gradient,
                               rtol=1e-8, atol=1e-18)


def test_optimize_density_constrained_monotone(problem):
    with TaskManager():
        filt = HelmholtzFilter(problem.mesh, radius=0.12)
        cpts, _ = orbit_arc_points(1.4, -0.3, 8)
        f_con = field_functional_load(problem.fes, cpts,
                                      np.full(len(cpts), 1.0 / len(cpts)),
                                      axis=2, scale=MU0, bonus_intorder=10)
        rho0 = np.full(problem.prob.n_el, 0.5)
        lin0 = problem.prob.linearize(
            density_to_s(filt.apply(rho0), 100.0), problem.f_state,
            [problem.f_adj, f_con])
        target = float(lin0.values[1])
        result = optimize_density(
            problem.prob, problem.f_state, problem.f_adj, [f_con], [target],
            chi_iron=100.0, volume_fraction=0.5, density_filter=filt,
            move_limit=0.1, max_iterations=6)
    hist = result.history
    assert len(hist) >= 3
    objectives = [h["objective"] for h in hist]
    assert objectives[-1] > float(lin0.values[0])          # real ascent
    assert all(b >= a * (1.0 - 1e-6) for a, b in zip(objectives,
                                                     objectives[1:]))
    volume_max = 0.5 * problem.prob.element_volumes.sum()
    for h in hist:
        assert max(np.array(h["violation"]) / np.array(h["band"])) <= 1.25 + 1e-9
        assert h["volume"] <= volume_max * (1.0 + 1e-9)
    assert result.density.min() >= 0.0 and result.density.max() <= 1.0
    assert result.solves >= len(hist) + 1


def test_optimize_density_unconstrained_runs(problem):
    with TaskManager():
        result = optimize_density(
            problem.prob, problem.f_state, problem.f_adj,
            chi_iron=100.0, volume_fraction=0.5, max_iterations=3)
    objectives = [h["objective"] for h in result.history]
    assert len(objectives) == 3
    assert all(b >= a * (1.0 - 1e-6) for a, b in zip(objectives,
                                                     objectives[1:]))


def test_optimize_density_validation_errors(problem):
    with pytest.raises(ValueError, match="targets"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         [problem.f_adj], [], chi_iron=100.0,
                         volume_fraction=0.5)
    with pytest.raises(ValueError, match="volume_fraction"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         chi_iron=100.0, volume_fraction=0.0)
    with pytest.raises(ValueError, match="move_min"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         chi_iron=100.0, volume_fraction=0.5,
                         move_limit=-0.1)
    with pytest.raises(ValueError, match="volume budget"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         chi_iron=100.0, volume_fraction=0.5,
                         initial_density=np.ones(problem.prob.n_el))


def test_optimize_density_deep_restoration_prioritizes_feasibility():
    """An infeasible functional start must improve even when J must decrease."""
    chi_iron = 100.0

    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0 / np.asarray(s, dtype=float)
            rho = (chi - CHI_MIN) / (chi_iron - CHI_MIN)
            # Objective and constraint both increase with rho[0].  Reaching
            # the target therefore requires a deliberate objective decrease.
            grad_s = np.zeros(2)
            grad_s[0] = -(chi[0] ** 2) / (chi_iron - CHI_MIN)
            return SimpleNamespace(
                values=np.array([rho[0], rho[0]]),
                jacobians=np.stack([grad_s, grad_s]),
                gfM=None, gfLambdas=(None, None), state_iterations=0,
                adjoint_iterations=(0, 0))

    result = optimize_density(
        LinearProblem(), object(), object(), [object()], [0.2],
        chi_iron=chi_iron, volume_fraction=1.0,
        initial_density=np.array([0.9, 0.1]), band_floor=0.01,
        move_limit=0.2, max_iterations=2)
    assert len(result.history) == 2
    assert all(entry["band_mode"] == "deep" for entry in result.history)
    assert result.history[1]["objective"] < result.history[0]["objective"]
    assert result.history[1]["violation"][0] < result.history[0]["violation"][0]


def test_deep_restoration_reduces_worst_band_normalized_violation():
    """Restoration must not trade a worse max row for a smaller L1 total."""
    chi_iron = 100.0

    class UnevenLinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0 / np.asarray(s, dtype=float)
            rho = (chi - CHI_MIN) / (chi_iron - CHI_MIN)
            grad0 = np.array([
                -(chi[0] ** 2) / (chi_iron - CHI_MIN), 0.0])
            grad1 = np.array([
                0.0, -(chi[1] ** 2) / (chi_iron - CHI_MIN)])
            return SimpleNamespace(
                values=np.array([rho[1], rho[0], rho[1]]),
                jacobians=np.stack([grad1, grad0, grad1]),
                gfM=None, gfLambdas=(None, None, None),
                state_iterations=0, adjoint_iterations=(0, 0, 0))

    result = optimize_density(
        UnevenLinearProblem(), object(), object(), [object(), object()],
        [0.2, 0.2], chi_iron=chi_iron, volume_fraction=1.0,
        initial_density=np.array([0.9, 0.5]), band_floor=[0.01, 0.01],
        move_limit=0.2, max_iterations=2)
    merits = [entry["max_violation_over_band"] for entry in result.history]
    assert len(merits) == 2
    assert merits[1] < merits[0] < 70.0


def test_deep_restoration_allows_active_worst_row_to_change():
    band = np.ones(2)
    # The first row rises by 40 %, but remains below the improved maximum.
    # This is legitimate Chebyshev progress and must not be blocked by an
    # individual-row trust guard.
    assert _accept_deep_restoration(
        np.array([5.0, -8.0]), np.array([7.0, -7.5]), band, True)
    # A smaller L1 total cannot compensate for a worse maximum.
    assert not _accept_deep_restoration(
        np.array([5.0, -8.0]), np.array([1.0, -8.1]), band, True)
    assert not _accept_deep_restoration(
        np.array([5.0, -8.0]), np.array([4.0, -7.0]), band, False)


def test_deep_restoration_minimax_lp_uses_common_epigraph_cap():
    rho = np.array([0.5, 0.5, 0.37])
    gradients = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0]])
    violation = np.array([5.0, -8.0])
    update = _solve_minimax_lp_update(
        rho, gradients, violation, np.ones(2), np.ones(3), 3.0,
        move_limit=0.2)
    linearized = violation + gradients @ update.delta
    np.testing.assert_allclose(np.max(np.abs(linearized)), 6.0, atol=1e-10)
    np.testing.assert_allclose(update.predicted_objective, 6.0, atol=2e-8)
    assert update.predicted_objective < np.max(np.abs(violation))
    # The proximal second LP leaves the response-null density untouched.
    np.testing.assert_allclose(update.density[2], rho[2], atol=1e-12)


def test_optimize_density_applies_projection_and_reports_discreteness():
    chi_iron = 100.0

    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0/np.asarray(s, dtype=float)
            material_rho = (chi-CHI_MIN)/(chi_iron-CHI_MIN)
            grad_s = np.array([-(chi[0]**2)/(chi_iron-CHI_MIN), 0.0])
            return SimpleNamespace(values=np.array([material_rho[0]]),
                jacobians=grad_s[None,:],gfM=None,gfLambdas=(None,),
                state_iterations=0,adjoint_iterations=(0,))

    checkpoints=[]
    result = optimize_density(LinearProblem(),object(),object(),
        chi_iron=chi_iron,volume_fraction=.5,
        initial_density=np.array([.4,.6]),density_projection=HeavisideProjection(2.),
        move_limit=.1,max_iterations=1,
        checkpoint_callback=lambda entry,rho:checkpoints.append((entry,rho)))
    assert len(result.history)==1 and result.density[0]>.4
    assert "intermediate_fraction" in result.history[0]
    assert len(checkpoints)==1
    np.testing.assert_allclose(checkpoints[0][1],result.density)


def test_optimize_density_uses_periodic_reduced_variables_and_full_checkpoints():
    chi_iron = 100.0
    density_map = FFAGCyclicDensityMap(
        np.array([0, 0]), ((0, 1),), 1,
        ("periodic_min", "periodic_max"))

    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0/np.asarray(s, dtype=float)
            material_rho = (chi-CHI_MIN)/(chi_iron-CHI_MIN)
            weights = np.array([1.0, 2.0])
            gradient_s = -weights*chi**2/(chi_iron-CHI_MIN)
            return SimpleNamespace(
                values=np.array([weights @ material_rho]),
                jacobians=gradient_s[None, :], gfM=None,
                gfLambdas=(None,), state_iterations=0,
                adjoint_iterations=(0,))

    checkpoints = []
    result = optimize_density(
        LinearProblem(), object(), object(), chi_iron=chi_iron,
        volume_fraction=0.5, design_map=density_map,
        initial_density=np.array([0.4, 0.4]), move_limit=0.1,
        max_iterations=1,
        checkpoint_callback=lambda entry, rho: checkpoints.append(rho))

    assert result.density.shape == (2,)
    np.testing.assert_allclose(result.density[0], result.density[1])
    np.testing.assert_allclose(checkpoints[0], result.density)
    assert result.history[0]["design_volume"] == pytest.approx(
        np.sum(result.density))
    with pytest.raises(ValueError, match="not equal"):
        optimize_density(
            LinearProblem(), object(), object(), chi_iron=chi_iron,
            volume_fraction=0.5, design_map=density_map,
            initial_density=np.array([0.4, 0.5]), max_iterations=1)


def test_optimize_density_rejects_projected_volume_over_budget():
    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            raise AssertionError("an infeasible projected design must fail before solve")

    with pytest.raises(ValueError, match="projected iron volume budget"):
        optimize_density(
            LinearProblem(), object(), object(), chi_iron=100.0,
            volume_fraction=0.5, initial_density=np.full(2, 0.5),
            density_projection=HeavisideProjection(beta=4.0, eta=0.3),
            max_iterations=1)


def test_restore_projected_volume_repairs_beta_continuation_start():
    rho = np.array([0.2, 0.45, 0.7, 0.9])
    volumes = np.array([1.0, 2.0, 1.5, 0.5])
    projection = HeavisideProjection(beta=8.0, eta=0.4)
    before = float(volumes @ projection.apply(rho))
    limit = 0.5 * float(volumes.sum())
    assert before > limit

    feasible, info = restore_projected_volume(
        rho, volumes, 0.5, density_projection=projection)
    after = float(volumes @ projection.apply(feasible))
    assert info["changed"]
    assert info["shift"] > 0.0
    assert after <= limit + 1e-12
    assert abs(after-limit) <= 1e-10
    np.testing.assert_allclose(
        feasible, np.clip(rho-info["shift"], 0.0, 1.0), atol=0.0)
    # The scalar correction preserves the topology ranking.
    assert np.all(np.diff(feasible) >= 0.0)


def test_restore_projected_volume_leaves_feasible_design_unchanged():
    rho = np.array([0.1, 0.2, 0.3])
    feasible, info = restore_projected_volume(
        rho, np.ones(3), 0.5,
        density_projection=HeavisideProjection(beta=4.0))
    np.testing.assert_array_equal(feasible, rho)
    assert not info["changed"]
    assert info["iterations"] == 0


def test_optimize_density_retries_nonconverged_trial_with_smaller_move():
    chi_iron = 100.0

    class RetryProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0 / np.asarray(s, dtype=float)
            rho = (chi - CHI_MIN) / (chi_iron - CHI_MIN)
            if rho[0] > 0.575:
                raise RuntimeError("synthetic CG did not converge")
            grad_s = np.array([
                -(chi[0]**2) / (chi_iron - CHI_MIN), 0.0])
            return SimpleNamespace(
                values=np.array([rho[0]]), jacobians=grad_s[None, :],
                gfM=None, gfLambdas=(None,), state_iterations=0,
                adjoint_iterations=(0,))

    result = optimize_density(
        RetryProblem(), object(), object(), chi_iron=chi_iron,
        volume_fraction=0.5, initial_density=np.array([0.5, 0.5]),
        move_limit=0.1, move_min=0.01, max_iterations=1)
    np.testing.assert_allclose(result.density, [0.55, 0.45], atol=1e-12)
    assert result.solves == 3
    assert result.history[0]["trials"] == 2


# ---------------------------------------------------------------- Stage 3
def test_iron_only_extraction_and_ersatz_band(problem):
    """Hemisphere extraction: exact volume, sane demag physics (surface
    orientation lock), and the embedded-vs-exact-void ersatz sign.

    Bands from measured values (maxh=0.45 ball, 86 kept tets):
    <Mz> extracted 2.284 (a flipped surface orientation runs away to 18+),
    embedded 0/1 <Mz> 2.089, ersatz -8.6 % (charge clamp under-magnetizes
    the embedded iron; -8.3/-8.8 % on the maxh .35/.22 research meshes).
    """
    with TaskManager():
        vols = problem.prob.element_volumes
        cz = np.asarray(ng.Integrate(ng.z, problem.mesh, element_wise=True),
                        float) / vols
        keep = cz > 0.0
        hemi = iron_only_mesh(problem.mesh, keep)
        V_kept = float(vols[keep].sum())
        V_new = float(ng.Integrate(ng.CoefficientFunction(1.0), hemi))
        assert abs(V_new - V_kept) / V_kept < 1e-12
        assert hemi.ne == int(keep.sum())
        prob_h = DensityAdjointVIM(ng.HDiv(hemi, order=1), eps=1e-7)
        gf, _ = prob_h.solve(np.full(prob_h.n_el, 1.0 / 100.0),
                             uniform_field_load(prob_h.fes, (0, 0, 1.0)))
        Mz_extracted = float(ng.Integrate(gf[2], hemi)) / V_new
        s_bin = np.where(keep, 1.0 / 100.0, 1.0 / CHI_MIN)
        gf_b, _ = problem.prob.solve(s_bin, problem.f_state)
        vals = np.asarray(ng.Integrate(gf_b[2], problem.mesh,
                                       element_wise=True), float)
        Mz_embedded = float(vals[keep].sum() / V_kept)
    assert 2.2 < Mz_extracted < 2.4, Mz_extracted
    ersatz = (Mz_embedded - Mz_extracted) / Mz_extracted
    assert -0.15 < ersatz < -0.03, ersatz


def test_verify_design_iron_only_protocol(problem):
    """Matched-0/1 verification bands from the promoted protocol function.

    Measured on this mesh: band +4.4 % for an external-field functional
    (values -1.77e-7 embedded vs -1.85e-7 iron-only)."""
    with TaskManager():
        vols = problem.prob.element_volumes
        cz = np.asarray(ng.Integrate(ng.z, problem.mesh, element_wise=True),
                        float) / vols
        density = (cz > 0.0).astype(float)
        cpts, _ = orbit_arc_points(1.5, 0.4, 6)

        def state_builder(fes):
            return uniform_field_load(fes, (0.0, 0.0, 1.0))

        def functional_builder(fes):
            return field_functional_load(fes, cpts, np.full(6, 1.0 / 6.0),
                                         axis=2, scale=MU0, bonus_intorder=10)

        ver = verify_design_iron_only(
            problem.prob, density, state_builder, [functional_builder],
            chi_iron=100.0, gram_kwargs=dict(eps=1e-7))
    assert ver.keep.sum() == int((cz > 0.0).sum())
    assert ver.values_embedded.shape == ver.values_iron_only.shape == (1,)
    assert ver.values_embedded[0] < 0.0 and ver.values_iron_only[0] < 0.0
    assert 0.0 < ver.bands[0] < 0.12, ver.bands
    assert ver.embedded_iterations > 0 and ver.iron_only_iterations > 0
    assert ver.embedded_solution.space is problem.fes
    assert ver.iron_solution.space is ver.iron_problem.fes


def test_broken_hdiv_binary_topology_matches_exact_void(broken_problem):
    """Facet jumps remove the conforming-density interface model error."""
    with TaskManager():
        vols = broken_problem.prob.element_volumes
        cz = np.asarray(
            ng.Integrate(ng.z, broken_problem.mesh, element_wise=True),
            dtype=float) / vols
        density = (cz > 0.0).astype(float)
        cpts, _ = orbit_arc_points(1.5, 0.4, 6)

        def state_builder(fes):
            return uniform_field_load(fes, (0.0, 0.0, 1.0))

        def functional_builder(fes):
            return field_functional_load(
                fes, cpts, np.full(6, 1.0 / 6.0), axis=2,
                scale=MU0, bonus_intorder=10)

        ver = verify_design_iron_only(
            broken_problem.prob, density, state_builder,
            [functional_builder], chi_iron=100.0,
            gram_kwargs=dict(eps=1e-7), linear_solver="native")
    assert ver.iron_problem.demag.internal_interfaces
    assert abs(ver.bands[0]) < 2e-5, ver.bands


def test_iron_only_mesh_guards(problem):
    n_el = problem.prob.n_el
    with pytest.raises(ValueError, match="proper non-empty subset"):
        iron_only_mesh(problem.mesh, np.ones(n_el, bool))
    with pytest.raises(ValueError, match="mask has"):
        iron_only_mesh(problem.mesh, np.ones(3, bool))
    with pytest.raises(ValueError, match="shape"):
        verify_design_iron_only(problem.prob, np.ones(3),
                                lambda fes: None, [], chi_iron=100.0)


def test_density_penalty_mapping():
    chi_iron = 1000.0
    rho = np.array([0.0, 0.5, 1.0])
    s1 = density_to_s(rho, chi_iron, penalty=1.0)
    s3 = density_to_s(rho, chi_iron, penalty=3.0)
    # endpoints penalty-invariant; interior penalized (chi smaller -> s larger)
    assert s3[0] == s1[0] and s3[2] == s1[2]
    assert s3[1] > s1[1]
    # chain rule vs FD at penalty=3
    rho = np.array([0.2, 0.5, 0.9])
    grad_s = np.array([1.0, -2.0, 0.5])
    h = 1e-7
    fd = grad_s * (density_to_s(rho + h, chi_iron, penalty=3.0)
                   - density_to_s(rho - h, chi_iron, penalty=3.0)) / (2 * h)
    chain = density_gradient_from_s_gradient(rho, grad_s, chi_iron,
                                             penalty=3.0)
    np.testing.assert_allclose(chain, fd, rtol=1e-5)
    with pytest.raises(ValueError, match="penalty"):
        density_to_s(rho, chi_iron, penalty=0.5)


def test_orbit_arc_and_direction_pairs():
    pts, radial = orbit_arc_points(2.0, 0.5, 4)
    assert pts.shape == (4, 3) and radial.shape == (4, 3)
    np.testing.assert_allclose(np.linalg.norm(radial, axis=1), 1.0)
    np.testing.assert_allclose(pts[0], [2.0, 0.0, 0.5], atol=1e-14)
    # full circle omits the duplicate endpoint
    assert not np.allclose(pts[-1], pts[0])
    ppts, wts = gradient_pair_points(pts, np.ones(4), delta=0.2,
                                     direction=radial)
    assert ppts.shape == (8, 3)
    np.testing.assert_allclose(ppts[0], [2.1, 0.0, 0.5], atol=1e-14)
    np.testing.assert_allclose(ppts[4], [1.9, 0.0, 0.5], atol=1e-14)
    np.testing.assert_allclose(wts[:4], 5.0)
    np.testing.assert_allclose(wts[4:], -5.0)
    with pytest.raises(ValueError, match="zero direction"):
        gradient_pair_points(pts, np.ones(4), delta=0.2,
                             direction=np.zeros((4, 3)))


def test_sector_orbit_and_linear_optics_uniform_field():
    radius, field_strength = 1.2, 2.0
    span = (0.0, pi/3.0)

    def field(_point):
        return np.array([0.0, 0.0, field_strength])

    inverse_rigidity = -1.0/(radius*field_strength)
    orbit = track_sector_orbit(
        field, radius, span, inverse_rigidity, n_steps=360)
    assert orbit.exit_radius == pytest.approx(radius, rel=2e-11)
    assert orbit.exit_angle == pytest.approx(span[1], abs=2e-11)
    optics = sector_linear_optics(
        field, radius, span, inverse_rigidity, n_steps=360)
    assert optics.radial_determinant == pytest.approx(1.0, rel=2e-6)
    assert optics.vertical_determinant == pytest.approx(1.0, rel=2e-6)


def test_straightened_bend_validation_uniform_combined_function():
    length=4.0;rigidity=2.0;angle=pi/6.0
    stations=np.linspace(0.0,length,81)
    bz=np.full(stations.shape,rigidity*angle/length)
    result=straightened_bend_validation(
        stations,bz,np.zeros_like(stations),rigidity)
    assert result.bend_angle == pytest.approx(angle,rel=2e-15)
    np.testing.assert_allclose(result.position[-1,0],0.0,atol=2e-14)
    expected_chord=2.0*(length/angle)*np.sin(0.5*angle)
    assert result.position[-1,1] == pytest.approx(expected_chord,rel=2e-14)
    assert result.optics.dispersion[0,0] == 0.0
    assert result.optics.dispersion[0,1] == 0.0
    assert result.optics.dispersion[-1,0] == pytest.approx(
        (length/angle)*(1.0-np.cos(angle)),rel=2e-14)


def test_isochronous_profile_metrics_exact_gamma_law():
    radii = np.linspace(0.1, 0.4, 5)
    gamma = 1.0/np.sqrt(1.0-(radii/0.8)**2)
    exact = isochronous_profile_metrics(radii, 1.7*gamma, gamma)
    assert exact["max_abs_field_error"] < 3e-16
    assert exact["max_abs_period_error"] < 3e-16


def test_isochronous_increment_targets_include_fixed_coil_field():
    gamma = np.array([1.0, 1.1, 1.25])
    targets = isochronous_increment_targets(
        gamma, reference_increment=0.2, external_bz=1.0,
        reference_index=1)
    total = 1.0 + targets
    np.testing.assert_allclose(total/total[1], gamma/gamma[1], rtol=1e-15)
    assert not np.allclose(targets, 0.2*gamma/gamma[1])
    np.testing.assert_allclose(
        isochronous_total_field_bands(targets, 1.0, 5e-3),
        5e-3*np.abs(total))
