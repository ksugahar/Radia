import numpy as np
import pytest

from radia.accelerator_magnet_topopt import (
    MultiMomentumTransferMatrixObjective,
    build_multi_orbit_field_response_matrix,
    optimize_hdiv_mmm_magnet_from_transfer_matrices,
    solve_transfer_matrix_field_correction,
)
from radia.ffag_topopt import (
    FFAGSoftEdgeCellSpec,
    build_ffag_cell_target_family,
    enge_fringe_integrals,
    magnetic_rigidity_from_kinetic_energy,
    recover_periodic_planar_closed_orbit,
    recover_periodic_planar_closed_orbit_native,
)
from radia.isochronous_topopt import (
    combined_function_transfer_map_from_field_response,
)


def _one_segment_arc(*, radius=10.0, angle=0.1, rigidity=1.5):
    from radia.accelerator_magnet_topopt import PlanarDesignOrbit

    positions = np.array([
        [0.0, 0.0, 2.0],
        [radius * np.sin(angle), radius * (1.0 - np.cos(angle)), 2.0],
    ])
    tangents = np.array([
        [1.0, 0.0, 0.0],
        [np.cos(angle), np.sin(angle), 0.0],
    ])
    return PlanarDesignOrbit(
        positions, tangents, magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0, 0.0, 1.0]))


def test_ffag_validation_restart_replays_accepted_binary_history(tmp_path):
    import importlib.util
    import json
    from pathlib import Path

    runner_path=(Path(__file__).parents[1]/"validation_test"/"ffag_topopt"/
                 "validation_ffag_full_field_c_yoke.py")
    spec=importlib.util.spec_from_file_location("ffag_restart_runner",runner_path)
    runner=importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(runner)
    saved={
        "schema":"radia.ffag-fixed-one-pass-c-yoke/v1",
        "mesh":{"elements":4,"dofs":12,"final_active_elements":3},
        "coil":{"optimized_source_scale":2.25},
        "optics":{"energies_mev":[31.0,250.0]},
        "optimization":{
            "final_max_band_ratio":3.5,
            "history":[
                {"added_elements":[1,2],"removed_elements":[]},
                {"added_elements":[3],"removed_elements":[1]},
            ],
        },
    }
    result_path=tmp_path/"restart.json"
    result_path.write_text(json.dumps(saved),encoding="utf-8")

    active,metadata=runner._restart_active_elements(
        result_path,fixed_active=np.array([True,False,False,False]),
        mesh_elements=4,dofs=12,energies=[31.0,250.0])

    np.testing.assert_array_equal(active,[True,False,True,True])
    assert metadata["source_final_active_elements"]==3
    assert metadata["source_final_max_band_ratio"]==3.5
    assert metadata["source_scale"]==2.25
    assert metadata["active_set_contract"]=="accepted-history-replay"
    saved["mesh"]["final_active_elements"]=2
    result_path.write_text(json.dumps(saved),encoding="utf-8")
    with pytest.raises(RuntimeError,match="saved final count"):
        runner._restart_active_elements(
            result_path,fixed_active=np.array([True,False,False,False]),
            mesh_elements=4,dofs=12,energies=[31.0,250.0])


def test_ffag_validation_restart_prefers_explicit_final_active_set(tmp_path):
    import importlib.util
    import json
    from pathlib import Path

    runner_path=(Path(__file__).parents[1]/"validation_test"/"ffag_topopt"/
                 "validation_ffag_full_field_c_yoke.py")
    spec=importlib.util.spec_from_file_location(
        "ffag_explicit_restart_runner",runner_path)
    runner=importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(runner)
    saved={
        "schema":"radia.ffag-fixed-one-pass-c-yoke/v1",
        "mesh":{
            "elements":5,
            "dofs":15,
            "final_active_elements":4,
            "final_active_element_ids":[0,1,3,4],
        },
        "coil":{"optimized_source_scale":2.25},
        "optics":{"energies_mev":[31.0,250.0]},
        "optimization":{
            "final_max_band_ratio":3.5,
            # This history is intentionally relative to an earlier restart and
            # therefore cannot reconstruct the complete incumbent by itself.
            "history":[{"added_elements":[4],"removed_elements":[2]}],
        },
    }
    result_path=tmp_path/"restart-explicit.json"
    result_path.write_text(json.dumps(saved),encoding="utf-8")

    active,metadata=runner._restart_active_elements(
        result_path,fixed_active=np.array([True,True,False,False,False]),
        mesh_elements=5,dofs=15,energies=[31.0,250.0])

    np.testing.assert_array_equal(active,[True,True,False,True,True])
    assert metadata["active_set_contract"]=="explicit-final-elements"
    saved["mesh"]["final_active_element_ids"]=[0,1,1,4]
    result_path.write_text(json.dumps(saved),encoding="utf-8")
    with pytest.raises(RuntimeError,match="explicit active set"):
        runner._restart_active_elements(
            result_path,fixed_active=np.array([True,True,False,False,False]),
            mesh_elements=5,dofs=15,energies=[31.0,250.0])


def test_ffag_validation_restart_chains_delta_histories(tmp_path):
    import importlib.util
    import json
    from pathlib import Path

    runner_path=(Path(__file__).parents[1]/"validation_test"/"ffag_topopt"/
                 "validation_ffag_full_field_c_yoke.py")
    spec=importlib.util.spec_from_file_location(
        "ffag_chained_restart_runner",runner_path)
    runner=importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(runner)
    common={
        "schema":"radia.ffag-fixed-one-pass-c-yoke/v1",
        "coil":{"optimized_source_scale":2.25},
        "optics":{"energies_mev":[31.0,250.0]},
    }
    parent={**common,
        "mesh":{"elements":5,"dofs":15,"final_active_elements":3},
        "restart":None,
        "optimization":{
            "final_max_band_ratio":4.0,
            "history":[{"added_elements":[2],"removed_elements":[]}],
        },
    }
    parent_path=tmp_path/"parent.json"
    parent_path.write_text(json.dumps(parent),encoding="utf-8")
    child={**common,
        "mesh":{"elements":5,"dofs":15,"final_active_elements":3},
        "restart":{"path":str(parent_path)},
        "optimization":{
            "final_max_band_ratio":3.5,
            "history":[{"added_elements":[4],"removed_elements":[2]}],
        },
    }
    child_path=tmp_path/"child.json"
    child_path.write_text(json.dumps(child),encoding="utf-8")

    active,metadata=runner._restart_active_elements(
        child_path,fixed_active=np.array([True,True,False,False,False]),
        mesh_elements=5,dofs=15,energies=[31.0,250.0])

    np.testing.assert_array_equal(active,[True,True,False,False,True])
    assert metadata["active_set_contract"]=="parent-restart-plus-history"


def test_relativistic_proton_rigidity_and_tanh_enge_integrals():
    np.testing.assert_allclose(
        magnetic_rigidity_from_kinetic_energy([31.0, 250.0]),
        [0.8111430402581485, 2.4321284387255044], rtol=2e-14)
    coordinate = np.linspace(-12.0, 12.0, 48001)
    field = 0.5 * (1.0 - np.tanh(coordinate))
    result = enge_fringe_integrals(
        coordinate, field, body_field_t=1.0, full_gap_m=1.0)
    assert abs(result.effective_boundary_m) < 2e-14
    np.testing.assert_allclose(result.i1, np.pi**2 / 24.0, atol=2e-8)
    np.testing.assert_allclose(result.i2, 0.5, atol=2e-10)
    assert abs(result.equal_integral_residual) < 1e-14


def test_bell_abell_soft_edge_family_is_periodic_and_hits_its_exact_map():
    spec = FFAGSoftEdgeCellSpec.bell_abell()
    family = build_ffag_cell_target_family(
        [31.0, 140.0, 250.0], spec=spec, n_segments=96,
        transfer_matrix_band=2e-3, bend_field_band=2e-3)
    raw = np.concatenate([
        reference.field_response for reference in family.references])
    np.testing.assert_allclose(
        [reference.bend_angle_rad for reference in family.references],
        spec.cell_bend_angle_rad, atol=8e-15)
    assert max(reference.periodic_position_residual_m
               for reference in family.references) < 8e-14
    assert max(reference.periodic_tangent_residual
               for reference in family.references) < 8e-14
    np.testing.assert_allclose(
        family.objective.transform(raw),
        family.objective.response_target, atol=3e-12)
    assert family.objective.raw_field_response_size == 6 * 96
    assert np.all(np.diff([
        reference.transverse_offset_m
        for reference in family.references]) > 0.0)


def test_multi_momentum_map_jacobian_matches_regression_difference():
    orbits = (
        _one_segment_arc(radius=8.0, angle=0.07, rigidity=1.1),
        _one_segment_arc(radius=9.0, angle=0.09, rigidity=2.0),
    )
    raw = np.array([
        orbits[0].magnetic_rigidity * orbits[0].signed_curvature[0], 0.17,
        orbits[1].magnetic_rigidity * orbits[1].signed_curvature[0], -0.08,
    ])
    matrices = np.asarray([
        combined_function_transfer_map_from_field_response(
            raw[2*index:2*index+2], orbit.segment_lengths,
            orbit.magnetic_rigidity,
            response_entries=((0, 0), (0, 5), (2, 2))).matrix
        for index, orbit in enumerate(orbits)])
    objective = MultiMomentumTransferMatrixObjective(
        orbits, matrices, 1e-3, 1e-3,
        response_entries=((0, 0), (0, 5), (2, 2)))
    analytic = objective.transform_jacobian(raw)
    regression = np.empty_like(analytic)
    step = 2.0e-7
    for column in range(raw.size):
        delta = np.zeros_like(raw)
        delta[column] = step
        regression[:, column] = (
            objective.transform(raw + delta)
            - objective.transform(raw - delta)) / (2.0 * step)
    np.testing.assert_allclose(analytic, regression, rtol=2e-7, atol=2e-9)


def test_full_field_periodic_orbit_recovers_uniform_field_circle():
    class UniformField:
        def __init__(self, bending_field):
            self.bending_field = float(bending_field)

        def b_field(self, points):
            values = np.asarray(points, dtype=float)
            return np.broadcast_to(
                [0.0, 0.0, self.bending_field], values.shape)

    rigidity = 1.7
    radius = 3.2
    angle = 0.25
    bending_field = rigidity / radius
    result = recover_periodic_planar_closed_orbit(
        UniformField(bending_field), magnetic_rigidity=rigidity,
        cell_angle_rad=angle, initial_radius_m=1.03 * radius,
        initial_incidence_angle_rad=0.02, n_segments=32,
        gradient_offset=0.002)

    np.testing.assert_allclose(result.entrance_radius_m, radius, atol=3e-10)
    np.testing.assert_allclose(result.path_length_m, radius * angle,
                               atol=3e-10)
    np.testing.assert_allclose(
        result.orbit.path_length_stations,
        np.linspace(0.0, radius * angle, 33),
        atol=3e-10,
    )
    np.testing.assert_allclose(
        result.orbit.signed_curvature, 1.0 / radius, atol=2e-11
    )
    assert abs(result.entrance_incidence_angle_rad) < 2e-10
    assert result.periodic_position_residual_m < 3e-10
    assert result.periodic_tangent_residual < 3e-10
    assert result.vertical_position_residual_m == 0.0
    assert result.vertical_tangent_residual == 0.0
    np.testing.assert_allclose(result.field_response[:32], bending_field)
    np.testing.assert_allclose(result.field_response[32:], 0.0, atol=1e-13)
    assert result.transfer.matrix.shape == (6, 6)


def test_native_periodic_orbit_recovery_matches_uniform_field_circle():
    class UniformField:
        def __init__(self, bending_field):
            self.bending_field = float(bending_field)

        def b_field(self, points):
            values = np.asarray(points, dtype=float)
            return np.broadcast_to(
                [0.0, 0.0, self.bending_field], values.shape)

    rigidity = 1.7
    radius = 3.2
    angle = 0.25
    bending_field = rigidity / radius
    result = recover_periodic_planar_closed_orbit_native(
        UniformField(bending_field),
        constant_field_t=(0.0, 0.0, bending_field),
        magnetic_rigidity=rigidity, cell_angle_rad=angle,
        initial_radius_m=1.03 * radius,
        initial_incidence_angle_rad=0.02, n_segments=32,
        gradient_offset=0.002, tracking_step_m=5.0e-4)

    np.testing.assert_allclose(result.entrance_radius_m, radius, atol=3e-9)
    np.testing.assert_allclose(
        result.path_length_m, radius * angle, atol=3e-9)
    np.testing.assert_allclose(
        result.orbit.signed_curvature, 1.0 / radius, atol=2e-12)
    assert abs(result.entrance_incidence_angle_rad) < 2e-9
    assert result.periodic_position_residual_m < 3e-9
    assert result.periodic_tangent_residual < 3e-9
    assert result.vertical_position_residual_m == 0.0
    assert result.vertical_tangent_residual == 0.0
    np.testing.assert_allclose(result.field_response[:32], bending_field)
    np.testing.assert_allclose(result.field_response[32:], 0.0, atol=1e-13)


def test_two_rigidity_bdm1_hex_topology_reaches_both_cell_maps():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.accelerator_magnet_topopt import PlanarDesignOrbit
    from radia.isochronous_topopt import MU0, uniform_field_load
    from radia.topology_optimization import solve_hdiv_mmm_active_elements
    from radia.vim._vim import build_charge_gram

    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=1, nz=1)
    fes = ng.HDiv(mesh, order=1, discontinuous=True)
    with ng.TaskManager():
        _, gram, _ = build_charge_gram(
            fes, eps=1e-10, leafsize=256, eta=2.0,
            internal_interfaces=True)
    with ng.TaskManager():
        source = uniform_field_load(fes, (0.0, 0.0, 1.0e5))
    rhs = np.asarray(source.vec.FV().NumPy(), dtype=float).copy()
    initial = np.array([True, False])
    target_active = np.ones(2, dtype=bool)
    zero_response = np.zeros((1, fes.ndof))
    initial_state = solve_hdiv_mmm_active_elements(
        charge_gram=gram, fes=fes, inv_chi=.2, rhs=rhs,
        response_matrix=zero_response, active_elements=initial,
        solve_tolerance=1e-11)[0]
    target_state = solve_hdiv_mmm_active_elements(
        charge_gram=gram, fes=fes, inv_chi=.2, rhs=rhs,
        response_matrix=zero_response, active_elements=target_active,
        solve_tolerance=1e-11)[0]

    provisional_orbits = (
        _one_segment_arc(radius=7.0, angle=0.06, rigidity=1.3),
        _one_segment_arc(radius=9.0, angle=0.08, rigidity=2.1),
    )
    provisional_maps = np.asarray([np.eye(6), np.eye(6)])
    provisional_objective = MultiMomentumTransferMatrixObjective(
        provisional_orbits, provisional_maps, 1.0, 1.0)
    with ng.TaskManager():
        native_rows = build_multi_orbit_field_response_matrix(
            gram, provisional_objective, gradient_offset=(0.02, 0.025),
            field_scale=MU0)
    incident = np.array([MU0 * 1.0e5, 0.0, MU0 * 1.0e5, 0.0])
    target_raw = native_rows @ target_state + incident
    orbits = tuple(PlanarDesignOrbit(
        orbit.positions, orbit.tangents,
        magnetic_rigidity=float(
            target_raw[2 * index] / orbit.signed_curvature[0]),
        bend_axis=orbit.bend_axis)
        for index, orbit in enumerate(provisional_orbits))
    assert all(orbit.magnetic_rigidity > 0.0 for orbit in orbits)
    maps = np.asarray([
        combined_function_transfer_map_from_field_response(
            target_raw[2*index:2*index+2], orbit.segment_lengths,
            orbit.magnetic_rigidity,
            response_entries=tuple(
                (row, column) for row in range(6) for column in range(6)
            )).matrix
        for index, orbit in enumerate(orbits)])
    assert native_rows.shape == (4, fes.ndof)
    assert native_rows.flags.c_contiguous
    initial_raw = native_rows @ initial_state + incident
    bend_change = np.max(np.abs(
        target_raw[[0, 2]] - initial_raw[[0, 2]]))
    map_change = max(np.max(np.abs(
        maps[index] - combined_function_transfer_map_from_field_response(
            initial_raw[2*index:2*index+2], orbit.segment_lengths,
            orbit.magnetic_rigidity,
            response_entries=tuple(
                (row, column) for row in range(6) for column in range(6)
            )).matrix)) for index, orbit in enumerate(orbits))
    assert bend_change > 0.0 and map_change > 0.0
    bend_band = 0.2 * bend_change
    map_band = 0.2 * map_change
    volumes = np.asarray(ng.Integrate(1.0, mesh, element_wise=True))
    objective = MultiMomentumTransferMatrixObjective(
        orbits, maps, map_band, bend_band)
    field_correction = solve_transfer_matrix_field_correction(
        objective, initial_raw, relative_tolerance=1e-10)

    result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
        orbits, maps, transfer_matrix_band=map_band,
        bend_field_band=bend_band, charge_gram=gram, fes=fes, inv_chi=.2,
        rhs=rhs, field_response_matrix=native_rows,
        incident_field_response=incident,
        field_correction=field_correction,
        active_elements=initial, element_volumes=volumes,
        volume_max=float(np.sum(volumes)) + 1e-14,
        fixed_active_elements=initial, maximum_batch_elements=1,
        graph_front_proposal_limit=0, max_iterations=1,
        solve_tolerance=1e-11)

    assert result.converged
    assert result.field_correction is field_correction
    np.testing.assert_array_equal(result.active_elements, target_active)
    np.testing.assert_allclose(result.realized_transfer_matrices, maps,
                               atol=2e-12)
    assert np.max(result.orbit_field_max_band_ratios) < 1e-7
    assert np.max(result.transfer_matrix_max_band_ratios) < 1e-7
    assert np.max(result.realized_symplectic_residuals) < 1e-12
    assert result.topology.valid


def test_full_field_outer_loop_tracks_only_before_and_after_binary_batch(
        monkeypatch):
    from types import SimpleNamespace

    import radia.ffag_topopt as ffag
    import radia.topology_optimization as topopt
    from radia.accelerator_magnet_topopt import (
        CoilBuilderHDivSource,
        MultiMomentumAcceleratorMagnetTopologyResult,
    )
    from radia.ffag_topopt import FullFieldClosedOrbit

    family=build_ffag_cell_target_family(
        [31.0,40.0],n_segments=16,
        transfer_matrix_band=2e-4,bend_field_band=2e-4,
        response_entries=((0,0),(0,5),(2,2)))
    target_raw=np.concatenate([
        reference.field_response for reference in family.references])
    current_raw=target_raw.copy()
    current_raw[16]+=0.02

    def recovered(raw):
        values=[];offset=0
        for reference,objective in zip(
                family.references,family.objective.objectives):
            count=objective.raw_field_response_size
            field=np.asarray(raw[offset:offset+count],dtype=float)
            transfer=combined_function_transfer_map_from_field_response(
                field,reference.orbit.segment_lengths,
                reference.orbit.magnetic_rigidity,
                response_entries=objective.response_entries)
            values.append(FullFieldClosedOrbit(
                reference.magnetic_rigidity_tm,reference.orbit,
                float(np.sum(reference.orbit.segment_lengths)),
                float(np.linalg.norm(reference.orbit.positions[0,:2])),
                0.0,0.0,0.0,0.0,0.0,1,field,transfer))
            offset+=count
        return tuple(values)

    recovered_calls=[]
    def fake_recover(*args,**kwargs):
        recovered_calls.append(kwargs.get("initial_references"))
        return recovered(current_raw if len(recovered_calls)==1 else target_raw)

    correction_calls=[]
    active_initial=np.array([True,False])
    active_candidate=np.array([True,True])
    dummy_topology=SimpleNamespace(valid=True)
    def fake_optimize(orbits,matrices,**kwargs):
        correction=kwargs.get("field_correction")
        assert correction is not None
        correction_calls.append(correction)
        generation=topopt.HDivMMMGenerationResult(
            active_candidate.copy(),np.ones(2),
            correction.target_field_response.copy(),tuple(),False,1.0,
            correction.target_field_response.copy(),"one mocked field batch")
        objective=MultiMomentumTransferMatrixObjective(
            tuple(orbits),np.asarray(matrices),
            kwargs["transfer_matrix_band"],kwargs["bend_field_band"],
            kwargs["response_entries"])
        split=objective.split_raw_response(correction.target_field_response)
        return MultiMomentumAcceleratorMagnetTopologyResult(
            objective,generation,split,
            np.asarray(matrices),np.ones(len(orbits)),
            np.ones(len(orbits)),np.zeros(len(orbits)),
            np.zeros(len(orbits)),dummy_topology,correction)

    segments=np.array([
        [[0.0,0.0,0.0],[1.0,0.0,0.0]],
        [[1.0,0.0,0.0],[0.0,0.0,0.0]],
    ])
    source=CoilBuilderHDivSource(((segments,1.0),))
    monkeypatch.setattr(
        CoilBuilderHDivSource,"assemble_hdiv_rhs",
        lambda self,fes:np.zeros(2))
    monkeypatch.setattr(
        CoilBuilderHDivSource,"incident_orbit_field_response",
        lambda self,objective,gradient_offset:current_raw.copy())
    monkeypatch.setattr(
        topopt,"solve_hdiv_mmm_active_elements",
        lambda **kwargs:(np.zeros(2),np.zeros(1),0))
    monkeypatch.setattr(
        topopt,"ngsolve_growth_topology",
        lambda mesh,active:dummy_topology)
    monkeypatch.setattr(ffag,"CoilHDivTotalField",lambda *args,**kwargs:object())
    monkeypatch.setattr(ffag,"recover_ffag_closed_orbit_family",fake_recover)
    monkeypatch.setattr(
        ffag,"build_multi_orbit_field_response_matrix",
        lambda gram,objective,gradient_offset:np.zeros((
            objective.raw_field_response_size,2)))
    monkeypatch.setattr(
        ffag,"optimize_hdiv_mmm_magnet_from_transfer_matrices",fake_optimize)
    fes=SimpleNamespace(ndof=2,mesh=object())

    result=ffag.optimize_ffag_hdiv_mmm_from_transfer_matrices(
        family,source=source,charge_gram=object(),fes=fes,inv_chi=0.1,
        active_elements=active_initial,element_volumes=np.ones(2),
        volume_max=2.0,optimize_source_scale=False,orbit_segments=16,
        max_outer_iterations=1,inner_iterations=1,
        outer_initial_material_move_fraction=1.0,
        outer_trust_region_trials=1)

    assert result.converged
    np.testing.assert_array_equal(result.active_elements,active_candidate)
    assert len(correction_calls)==1
    assert correction_calls[0].step_scale>0.0
    # One realized-orbit recovery creates the field target and one verifies
    # the accepted binary batch.  Candidate screening performs no tracking.
    assert len(recovered_calls)==2


def test_fixed_design_orbit_path_never_runs_periodic_orbit_recovery(
        monkeypatch):
    from types import SimpleNamespace

    import radia.ffag_topopt as ffag
    import radia.topology_optimization as topopt
    from radia.accelerator_magnet_topopt import (
        CoilBuilderHDivSource,
        MultiMomentumAcceleratorMagnetTopologyResult,
    )

    family=build_ffag_cell_target_family(
        [31.0,40.0],n_segments=16,
        transfer_matrix_band=2e-4,bend_field_band=2e-4,
        response_entries=((0,0),(0,5),(2,2)))
    target_raw=np.concatenate([
        reference.field_response for reference in family.references])
    current_raw=target_raw.copy()
    current_raw[16]+=0.02
    active_initial=np.array([True,False])
    active_candidate=np.array([True,True])
    dummy_topology=SimpleNamespace(valid=True)
    correction_calls=[]
    active_calls=[]
    material_iteration_calls=[]

    def fake_optimize(orbits,matrices,**kwargs):
        correction=kwargs["field_correction"]
        correction_calls.append(correction)
        active_calls.append(np.asarray(kwargs["active_elements"]).copy())
        material_iteration_calls.append(kwargs["max_iterations"])
        ratio=(1.0e12 if len(correction_calls)==1 else
               (2.0 if len(correction_calls)==2 else 1.0))
        generation=topopt.HDivMMMGenerationResult(
            active_candidate.copy(),np.ones(2),
            correction.target_field_response.copy(),tuple(),False,1.0,
            correction.target_field_response.copy(),"one fixed-orbit batch")
        objective=MultiMomentumTransferMatrixObjective(
            tuple(orbits),np.asarray(matrices),
            kwargs["transfer_matrix_band"],kwargs["bend_field_band"],
            kwargs["response_entries"])
        split=objective.split_raw_response(correction.target_field_response)
        return MultiMomentumAcceleratorMagnetTopologyResult(
            objective,generation,split,np.asarray(matrices),
            np.full(len(orbits),ratio),np.zeros(len(orbits)),
            np.zeros(len(orbits)),np.zeros(len(orbits)),
            dummy_topology,correction)

    segments=np.array([
        [[0.0,0.0,0.0],[1.0,0.0,0.0]],
        [[1.0,0.0,0.0],[0.0,0.0,0.0]],
    ])
    source=CoilBuilderHDivSource(((segments,1.0),))
    monkeypatch.setattr(
        CoilBuilderHDivSource,"assemble_hdiv_rhs",
        lambda self,fes:np.zeros(2))
    monkeypatch.setattr(
        CoilBuilderHDivSource,"incident_orbit_field_response",
        lambda self,objective,gradient_offset:current_raw.copy())
    monkeypatch.setattr(
        topopt,"solve_hdiv_mmm_active_elements",
        lambda **kwargs:(np.zeros(2),np.zeros(1),0))
    monkeypatch.setattr(
        ffag,"build_multi_orbit_field_response_matrix",
        lambda gram,objective,gradient_offset:np.zeros((
            objective.raw_field_response_size,2)))
    monkeypatch.setattr(
        ffag,"optimize_hdiv_mmm_magnet_from_transfer_matrices",fake_optimize)
    monkeypatch.setattr(
        ffag,"recover_ffag_closed_orbit_family",
        lambda *args,**kwargs:(_ for _ in ()).throw(
            AssertionError("fixed one-pass path must not recover a closed orbit")))
    fes=SimpleNamespace(ndof=2,mesh=object())

    result=ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
        family,source=source,charge_gram=object(),fes=fes,inv_chi=0.1,
        active_elements=active_initial,element_volumes=np.ones(2),
        volume_max=2.0,optimize_source_scale=False,
        max_optics_iterations=2,material_iterations_per_optics=1,
        map_trust_region_trials=2,initial_material_move_fraction=1.0,
        maximum_material_move_fraction=1.0)

    np.testing.assert_array_equal(result.active_elements,active_candidate)
    assert result.field_correction is correction_calls[2]
    assert result.source_scale==1.0
    assert len(result.optics_history)==2
    np.testing.assert_array_equal(active_calls[0],active_initial)
    np.testing.assert_array_equal(active_calls[1],active_initial)
    np.testing.assert_array_equal(active_calls[2],active_candidate)
    assert material_iteration_calls==[1,1,1]
    assert [trial.accepted for trial in result.map_trust_history]==[
        False,True,True]
    assert [trial.material_move_fraction
            for trial in result.map_trust_history]==[1.0,0.5,1.0]
    assert result.stop_reason=="fixed one-pass transfer bands reached"

    common=dict(
        source=source,charge_gram=object(),fes=fes,inv_chi=0.1,
        active_elements=active_initial,element_volumes=np.ones(2),
        volume_max=2.0,optimize_source_scale=False)
    with pytest.raises(TypeError,match="material_iterations_per_optics"):
        ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
            family,max_iterations=1,**common)
    with pytest.raises(ValueError,match="iteration counts must be positive"):
        ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
            family,max_optics_iterations=0,**common)
    with pytest.raises(ValueError,match="map_ratio_tolerance"):
        ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
            family,map_ratio_tolerance=np.inf,**common)
    with pytest.raises(ValueError,match="beam width and depth"):
        ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
            family,direct_map_oracle_exact_beam_width=2,**common)
    with pytest.raises(ValueError,match="graph-front proposal limit"):
        ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
            family,direct_map_oracle_graph_front_proposal_limit=-1,**common)

    reject_calls=[]

    def fake_reject(orbits,matrices,**kwargs):
        correction=kwargs["field_correction"]
        incumbent=np.asarray(kwargs["active_elements"],dtype=bool).copy()
        is_baseline=kwargs["max_iterations"]==0
        returned_active=(incumbent if is_baseline else
                         active_candidate.copy())
        ratio=5.0 if is_baseline else 1.0e12
        reject_calls.append((kwargs["max_iterations"],returned_active.copy()))
        generation=topopt.HDivMMMGenerationResult(
            returned_active,np.ones(2),
            correction.target_field_response.copy(),tuple(),False,1.0,
            correction.target_field_response.copy(),
            "baseline" if is_baseline else "rejected proposal")
        objective=MultiMomentumTransferMatrixObjective(
            tuple(orbits),np.asarray(matrices),
            kwargs["transfer_matrix_band"],kwargs["bend_field_band"],
            kwargs["response_entries"])
        split=objective.split_raw_response(correction.target_field_response)
        return MultiMomentumAcceleratorMagnetTopologyResult(
            objective,generation,split,np.asarray(matrices),
            np.full(len(orbits),ratio),np.zeros(len(orbits)),
            np.zeros(len(orbits)),np.zeros(len(orbits)),
            dummy_topology,correction)

    monkeypatch.setattr(
        ffag,"optimize_hdiv_mmm_magnet_from_transfer_matrices",fake_reject)
    rejected=ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
        family,source=source,charge_gram=object(),fes=fes,inv_chi=0.1,
        active_elements=active_initial,element_volumes=np.ones(2),
        volume_max=2.0,optimize_source_scale=False,
        max_optics_iterations=1,material_iterations_per_optics=1,
        map_trust_region_trials=2,initial_material_move_fraction=1.0,
        maximum_material_move_fraction=1.0)

    np.testing.assert_array_equal(rejected.active_elements,active_initial)
    assert rejected.stop_reason=="map-level trust-region proposals rejected"
    assert rejected.optics_history==tuple()
    assert [trial.accepted for trial in rejected.map_trust_history]==[
        False,False]
    assert [trial.material_move_fraction
            for trial in rejected.map_trust_history]==[1.0,0.5]
    assert [call[0] for call in reject_calls]==[1,1,0]

    oracle_calls=[]

    def fake_oracle(orbits,matrices,**kwargs):
        correction=kwargs.get("field_correction")
        incumbent=np.asarray(kwargs["active_elements"],dtype=bool).copy()
        direct=correction is None
        returned_active=active_candidate.copy() if direct else incumbent
        ratio=10.0 if direct else 1.0e12
        oracle_calls.append((
            "direct-map-jacobian" if direct else "field-target",
            kwargs.get("exact_beam_width"),
            kwargs.get("exact_beam_depth"),
            kwargs.get("graph_front_proposal_limit")))
        raw=(current_raw.copy() if correction is None else
             correction.target_field_response.copy())
        generation=topopt.HDivMMMGenerationResult(
            returned_active,np.ones(2),raw.copy(),tuple(),False,1.0,
            raw.copy(),"direct oracle" if direct else "rejected field target")
        objective=MultiMomentumTransferMatrixObjective(
            tuple(orbits),np.asarray(matrices),
            kwargs["transfer_matrix_band"],kwargs["bend_field_band"],
            kwargs["response_entries"])
        split=objective.split_raw_response(raw)
        return MultiMomentumAcceleratorMagnetTopologyResult(
            objective,generation,split,np.asarray(matrices),
            np.full(len(orbits),ratio),np.zeros(len(orbits)),
            np.zeros(len(orbits)),np.zeros(len(orbits)),
            dummy_topology,correction)

    monkeypatch.setattr(
        ffag,"optimize_hdiv_mmm_magnet_from_transfer_matrices",fake_oracle)
    oracle=ffag.optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
        family,source=source,charge_gram=object(),fes=fes,inv_chi=0.1,
        active_elements=active_initial,element_volumes=np.ones(2),
        volume_max=2.0,optimize_source_scale=False,
        max_optics_iterations=1,material_iterations_per_optics=1,
        map_trust_region_trials=1,initial_material_move_fraction=1.0,
        maximum_material_move_fraction=1.0,
        direct_map_oracle_fallback=True,
        direct_map_oracle_exact_beam_width=4,
        direct_map_oracle_exact_beam_depth=2,
        direct_map_oracle_graph_front_proposal_limit=3)

    np.testing.assert_array_equal(oracle.active_elements,active_candidate)
    assert oracle_calls==[("field-target",None,None,None),
                          ("direct-map-jacobian",4,2,3)]
    assert [trial.proposal_model for trial in oracle.map_trust_history]==[
        "field-target","direct-map-jacobian"]
    assert [trial.accepted for trial in oracle.map_trust_history]==[
        False,True]
    assert oracle.topology_result.field_correction is None
