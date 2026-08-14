import numpy as np

from radia.accelerator_magnet_topopt import (
    CoilBuilderHDivSource,
    CoilHDivTotalField,
    MultiMomentumTransferMatrixObjective,
    PlanarDesignOrbit,
    PlanarTransferMatrixObjective,
    build_planar_orbit_field_response_matrix,
    optimize_hdiv_mmm_magnet_from_transfer_matrix,
    planar_orbit_field_observations,
    run_transfer_matrix_material_inverse_pipeline,
    solve_transfer_matrix_field_correction,
)
from radia.isochronous_topopt import (
    combined_function_transfer_map_from_field_response,
)


def _one_segment_arc(*, radius=10.0, angle=0.1, rigidity=1.5):
    positions=np.array([
        [0.0,0.0,2.0],
        [radius*np.sin(angle),radius*(1.0-np.cos(angle)),2.0],
    ])
    tangents=np.array([
        [1.0,0.0,0.0],
        [np.cos(angle),np.sin(angle),0.0],
    ])
    return PlanarDesignOrbit(
        positions,tangents,magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0,0.0,1.0]))


def test_planar_transfer_matrix_objective_uses_orbit_curvature_and_full_map():
    orbit=_one_segment_arc()
    bend_field=float(orbit.magnetic_rigidity*orbit.signed_curvature[0])
    raw=np.array([bend_field,0.17])
    transfer=combined_function_transfer_map_from_field_response(
        raw,orbit.segment_lengths,orbit.magnetic_rigidity,
        response_entries=tuple(
            (row,column) for row in range(6) for column in range(6)))
    objective=PlanarTransferMatrixObjective(
        orbit,transfer.matrix,transfer_matrix_band=1e-6,
        bend_field_band=1e-7)

    np.testing.assert_allclose(objective.required_bend_field,[bend_field])
    np.testing.assert_allclose(
        objective.transform(raw),objective.response_target,atol=1e-14)
    jacobian=objective.transform_jacobian(raw)
    assert jacobian.shape==(37,2)
    np.testing.assert_array_equal(jacobian[0],[1.0,0.0])


def test_transfer_matrix_inverse_produces_field_target_before_material_step():
    orbit=_one_segment_arc(radius=9.0,angle=0.08,rigidity=1.8)
    bend=float(orbit.magnetic_rigidity*orbit.signed_curvature[0])
    current=np.array([bend,0.15])
    wanted=np.array([bend,0.17])
    target_map=combined_function_transfer_map_from_field_response(
        wanted,orbit.segment_lengths,orbit.magnetic_rigidity,
        response_entries=tuple(
            (row,column) for row in range(6) for column in range(6)))
    objective=PlanarTransferMatrixObjective(
        orbit,target_map.matrix,transfer_matrix_band=2e-5,
        bend_field_band=2e-6)

    correction=solve_transfer_matrix_field_correction(
        objective,current,relative_tolerance=1e-10,
        line_search_steps=8)

    assert correction.numerical_rank==2
    assert correction.step_scale>0.0
    assert correction.nonlinear_max_band_ratio<(
        correction.current_max_band_ratio)
    assert abs(correction.target_field_response[1]-wanted[1])<(
        abs(current[1]-wanted[1]))
    np.testing.assert_allclose(
        correction.target_field_response,
        correction.current_field_response+correction.field_correction)
    assert np.all(correction.field_response_band>0.0)
    assert correction.status==(
        "dense optics TSVD-Chebyshev solve to orbit-field target")
    assert correction.derivative_backend=="forward-mode-expm-frechet-ad"
    assert correction.field_to_design_jacobian.shape==(37,2)


def test_field_to_map_to_aca_qr_tsvd_material_pipeline_is_auditable():
    orbit=_one_segment_arc(radius=9.0,angle=0.08,rigidity=1.8)
    bend=float(orbit.magnetic_rigidity*orbit.signed_curvature[0])
    current=np.array([bend,0.15])
    wanted=np.array([bend,0.17])
    target_map=combined_function_transfer_map_from_field_response(
        wanted,orbit.segment_lengths,orbit.magnetic_rigidity,
        response_entries=tuple(
            (row,column) for row in range(6) for column in range(6)))
    objective=PlanarTransferMatrixObjective(
        orbit,target_map.matrix,transfer_matrix_band=2e-5,
        bend_field_band=2e-6)

    # Candidate 10 supplies the desired gradient change. Candidate 11 is a
    # small opposite bend correction. The production factorization sees both
    # columns through ACA -> thin QR -> TSVD and the binary LP chooses the
    # useful whole element.
    candidate_delta=np.column_stack((
        np.array([0.0,0.02]),
        np.array([-1.0e-6,0.0]),
    ))
    result=run_transfer_matrix_material_inverse_pipeline(
        objective,current,candidate_elements=np.array([10,11]),
        candidate_field_response_delta=candidate_delta,
        candidate_volumes=np.ones(2),volume_budget=2.0,
        candidate_material_active=np.zeros(2,dtype=bool),
        maximum_changed_elements=2,
        field_inverse_relative_tolerance=1e-10,
        material_relative_tolerance=1e-10)

    assert result.stage_order==(
        "magnetic-field-distribution",
        "forward-ad-transfer-matrix",
        "target-transfer-matrix-difference",
        "tsvd-minimax-field-correction",
        "aca-thin-qr-tsvd-material-inverse")
    np.testing.assert_allclose(result.field_distribution,current)
    np.testing.assert_allclose(
        result.transfer_matrix_difference,
        target_map.matrix[None]-result.realized_transfer_matrices)
    np.testing.assert_allclose(
        result.normalized_transfer_matrix_difference,
        result.transfer_matrix_difference/2e-5)
    automatic=result.automatic_differentiation
    assert automatic.backend=="forward-mode-expm-frechet-ad"
    assert automatic.full_jacobian.shape==(37,2)
    np.testing.assert_allclose(
        automatic.directional_jacobian,automatic.full_jacobian)
    np.testing.assert_allclose(
        result.field_correction.target_field_response,
        current+result.field_correction.field_correction)
    selection=result.material_selection
    assert selection.aca_rank>=1
    assert selection.numerical_rank>=1
    np.testing.assert_array_equal(selection.selected_elements,[10])
    assert selection.predicted_max_band_ratio<1.0
    assert result.status==(
        "all-candidate band-normalized TSVD plus binary LP")


def test_transfer_matrix_inverse_uses_minimax_direction_in_tsvd_space():
    class CompetingResponseObjective:
        raw_field_response_size=1
        response_target=np.array([2.0,1.8])
        response_band=np.ones(2)

        @staticmethod
        def transform(field):
            value=float(np.asarray(field).reshape(-1)[0])
            return np.array([value,-10.0*value])

        @staticmethod
        def transform_jacobian(field):
            return np.array([[1.0],[-10.0]])

    correction=solve_transfer_matrix_field_correction(
        CompetingResponseObjective(),np.array([0.0]),
        relative_tolerance=1e-12,line_search_steps=8)

    # Ordinary weighted least squares points in the negative direction and
    # increases the active infinity-norm error.  The TSVD-retained Chebyshev
    # solve chooses the positive minimax direction instead.
    assert correction.field_correction[0]>0.0
    assert correction.step_scale==1.0
    assert correction.nonlinear_max_band_ratio<(
        correction.current_max_band_ratio)


def test_planar_orbit_observations_measure_binormal_field_and_normal_gradient():
    orbit=_one_segment_arc()
    points,weights=planar_orbit_field_observations(
        orbit,gradient_offset=0.02)
    tangent=np.sum(orbit.tangents,axis=0)
    tangent/=np.linalg.norm(tangent)
    normal=np.cross(orbit.bend_axis,tangent)
    field=(2.0+3.0*(points@normal))[:,None]*orbit.bend_axis[None,:]
    response=np.einsum("rpc,pc->r",weights,field)
    assert points.shape==(3,3)
    assert weights.shape==(2,3,3)
    np.testing.assert_allclose(response[1],3.0,rtol=0.0,atol=2e-14)


def test_coilbuilder_source_owns_hdiv_rhs_incident_rows_and_total_field():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.biot_savart import h_segments_batch
    from radia.coil_builder import CoilBuilder
    from radia.vim._vim import build_charge_gram

    coil = (CoilBuilder(current=1200.0)
            .set_start([0.2, 0.0, 2.0])
            .set_cross_section(0.01, 0.01)
            .add_arc(radius=0.2, arc_angle=360.0))
    source = CoilBuilderHDivSource.from_coilbuilders(coil, n_arc=24)
    points = np.array([[0.1, 0.1, 0.5], [0.2, 0.2, 0.6]])
    segments, current = coil.to_wire_segments(n_arc=24)
    np.testing.assert_allclose(
        source.h_field(points),
        h_segments_batch(np.asarray(segments), points, current=current),
        rtol=2e-14, atol=1e-14)
    assert source.segment_count == 24

    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    fes = ng.HDiv(mesh, order=1, discontinuous=True)
    with ng.TaskManager():
        rhs = source.assemble_hdiv_rhs(fes, bonus_intorder=3)
        _, gram, _ = build_charge_gram(
            fes, eps=1e-10, leafsize=256, eta=2.0,
            internal_interfaces=True)
    assert rhs.shape == (fes.ndof,)
    assert np.linalg.norm(rhs) > 0.0

    orbit = _one_segment_arc(radius=8.0, angle=0.08, rigidity=2.0)
    objective = MultiMomentumTransferMatrixObjective(
        (orbit,), np.asarray([np.eye(6)]), 1.0, 1.0)
    incident = source.incident_orbit_field_response(
        objective, gradient_offset=0.02)
    observation_points, weights = planar_orbit_field_observations(
        orbit, gradient_offset=0.02)
    np.testing.assert_allclose(
        incident,
        np.einsum("rpc,pc->r", weights, source.b_field(observation_points)),
        rtol=2e-14, atol=1e-15)

    total = CoilHDivTotalField(
        source, gram, np.zeros(fes.ndof), source_scale=1.7,
        hdiv_order=1)
    np.testing.assert_allclose(
        total.b_field(points), 1.7 * source.b_field(points),
        rtol=2e-14, atol=1e-15)


def test_orbit_and_transfer_matrix_create_target_whole_hex_magnet():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.topology_optimization import solve_hdiv_mmm_active_elements
    from radia.vim._vim import build_charge_gram

    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1,discontinuous=True)
    with ng.TaskManager():
        _,gram,mass=build_charge_gram(
            fes,eps=1e-10,leafsize=256,eta=2.0,
            internal_interfaces=True)
    rng=np.random.default_rng(20260811)
    rhs=np.asarray(mass@rng.normal(size=fes.ndof))
    initial=np.array([True,False])
    target_active=np.ones(2,dtype=bool)
    zero_response=np.zeros((1,fes.ndof))
    initial_state=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=zero_response,active_elements=initial,
        solve_tolerance=1e-11)[0]
    target_state=solve_hdiv_mmm_active_elements(
        charge_gram=gram,fes=fes,inv_chi=.2,rhs=rhs,
        response_matrix=zero_response,active_elements=target_active,
        solve_tolerance=1e-11)[0]

    orbit=_one_segment_arc(radius=8.0,angle=0.08,rigidity=2.0)
    with ng.TaskManager():
        native_rows=build_planar_orbit_field_response_matrix(
            gram,orbit,gradient_offset=0.025,field_scale=1.0)
    assert native_rows.shape==(2,fes.ndof)
    assert native_rows.flags.c_contiguous
    bend_field=float(orbit.magnetic_rigidity*orbit.signed_curvature[0])
    target_gradient=0.11
    states=np.vstack((initial_state,target_state))
    field_response=np.vstack((
        np.linalg.lstsq(states,np.array([0.0,bend_field]),rcond=None)[0],
        np.linalg.lstsq(states,np.array([0.0,target_gradient]),rcond=None)[0],
    ))
    target_map=combined_function_transfer_map_from_field_response(
        [bend_field,target_gradient],orbit.segment_lengths,
        orbit.magnetic_rigidity,
        response_entries=tuple(
            (row,column) for row in range(6) for column in range(6)))
    volumes=np.asarray(ng.Integrate(1.0,mesh,element_wise=True))
    objective=PlanarTransferMatrixObjective(
        orbit,target_map.matrix,transfer_matrix_band=2e-7,
        bend_field_band=2e-8)
    field_correction=solve_transfer_matrix_field_correction(
        objective,field_response@initial_state,
        relative_tolerance=1e-10)

    result=optimize_hdiv_mmm_magnet_from_transfer_matrix(
        orbit,target_map.matrix,transfer_matrix_band=2e-7,
        bend_field_band=2e-8,charge_gram=gram,fes=fes,inv_chi=.2,
        rhs=rhs,field_response_matrix=field_response,
        field_correction=field_correction,
        active_elements=initial,element_volumes=volumes,
        volume_max=float(np.sum(volumes))+1e-14,
        fixed_active_elements=initial,maximum_batch_elements=1,
        graph_front_proposal_limit=0,max_iterations=1,
        solve_tolerance=1e-11)

    assert result.converged
    assert result.field_correction is field_correction
    np.testing.assert_array_equal(result.active_elements,target_active)
    np.testing.assert_allclose(
        result.realized_transfer_matrix,target_map.matrix,atol=2e-12)
    assert result.orbit_field_max_band_ratio<1e-3
    assert result.transfer_matrix_max_band_ratio<1e-3
    assert result.target_symplectic_residual<1e-12
    assert result.realized_symplectic_residual<1e-12
    assert result.topology.valid
