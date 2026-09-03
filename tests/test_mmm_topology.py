from __future__ import annotations

import numpy as np
import pytest

import radia.mmm_topology as mmm
from radia.sheet_metal_optimization import (
    ShapeModelEvaluation,
    TopologyPreservingShapeResult,
    TopologyPreservingShapeState,
)
from radia.topology_optimization import HDivMMMGenerationResult


def _generation(*, converged=False):
    return HDivMMMGenerationResult(
        active_elements=np.asarray([True, False, True]),
        state=np.asarray([1.0, 2.0]),
        response=np.asarray([0.5]),
        history=(),
        converged=converged,
        stop_reason="test-stage-complete",
    )


def _shape_state():
    return TopologyPreservingShapeState(
        mesh=object(),
        model=object(),
        reference_parameters=np.asarray([0.0]),
        parameters=np.asarray([0.0]),
        evaluation=ShapeModelEvaluation(
            objective=2.0, response=np.asarray([0.25])),
    )


def test_default_policy_names_binary_lego_then_gettrafo_without_sculpt():
    policy = mmm.DEFAULT_MMM_TOPOLOGY_POLICY
    assert policy.name == "MMM-topology"
    assert policy.lego_stage == "aca-qr-tsvd-exact-schur-binary-lego"
    assert policy.shape_stage == "ngsolve-gettrafo-full-resolve"
    assert policy.primary_cad_mesher == "coreform-cubit"
    assert policy.fallback_mesher is None
    assert policy.sculpt_in_core_loop is False


def test_two_stage_facade_runs_lego_before_shape_and_checks_final_bands(
        monkeypatch):
    calls = []
    generation = _generation()
    initial = _shape_state()
    final_state = TopologyPreservingShapeState(
        initial.mesh, initial.model, initial.reference_parameters,
        np.asarray([0.1]),
        ShapeModelEvaluation(0.5, np.asarray([0.02])),
    )
    shape = TopologyPreservingShapeResult(final_state, (), True)

    def fake_lego(**options):
        calls.append(("lego", options))
        return generation

    def build_shape(observed):
        calls.append(("build-shape", observed))
        return initial

    def fake_shape(observed, **options):
        calls.append(("gettrafo", observed, options))
        return shape

    monkeypatch.setattr(mmm, "grow_hdiv_mmm_by_superposition", fake_lego)
    monkeypatch.setattr(mmm, "optimize_topology_preserving_shape", fake_shape)
    result = mmm.optimize_mmm_topology(
        lego_options={"marker": 17},
        build_shape_state=build_shape,
        linearize_shape_step="linearize",
        deformation_factory="deform",
        rebuild_shape_model="rebuild",
        evaluate_shape_model="evaluate",
        move_limit=0.05,
        shape_options={"max_iterations": 3},
        final_acceptance=lambda lego, smooth: bool(
            lego is generation
            and np.max(np.abs(smooth.state.evaluation.response)) <= 0.1),
    )

    assert [row[0] for row in calls] == ["lego", "build-shape", "gettrafo"]
    assert calls[0][1] == {"marker": 17}
    assert calls[2][2]["max_iterations"] == 3
    assert calls[2][2]["move_limit"] == pytest.approx(0.05)
    assert result.generation is generation
    assert result.shape is shape
    assert result.final_target_accepted is True
    np.testing.assert_array_equal(result.active_elements, [True, False, True])
    assert result.final_evaluation is final_state.evaluation


def test_two_stage_facade_does_not_claim_acceptance_without_band_checker(
        monkeypatch):
    generation = _generation(converged=True)
    state = _shape_state()
    shape = TopologyPreservingShapeResult(state, (), True)
    monkeypatch.setattr(
        mmm, "grow_hdiv_mmm_by_superposition", lambda **kwargs: generation)
    monkeypatch.setattr(
        mmm, "optimize_topology_preserving_shape",
        lambda initial, **kwargs: shape)

    result = mmm.optimize_mmm_topology(
        lego_options={},
        build_shape_state=lambda observed: state,
        linearize_shape_step=object(),
        deformation_factory=object(),
        rebuild_shape_model=object(),
        evaluate_shape_model=object(),
        move_limit=0.1,
    )
    assert result.final_target_accepted is None


def test_two_stage_facade_rejects_sculpt_core_policy_and_reserved_overrides():
    with pytest.raises(ValueError, match="Sculpt"):
        mmm.optimize_mmm_topology(
            lego_options={},
            build_shape_state=object(),
            linearize_shape_step=object(),
            deformation_factory=object(),
            rebuild_shape_model=object(),
            evaluate_shape_model=object(),
            move_limit=0.1,
            policy=mmm.MMMTopologyPolicy(sculpt_in_core_loop=True),
        )

    with pytest.raises(TypeError, match="cannot override"):
        generation = _generation()
        state = _shape_state()
        original = mmm.grow_hdiv_mmm_by_superposition
        try:
            mmm.grow_hdiv_mmm_by_superposition = lambda **kwargs: generation
            mmm.optimize_mmm_topology(
                lego_options={},
                build_shape_state=lambda observed: state,
                linearize_shape_step=object(),
                deformation_factory=object(),
                rebuild_shape_model=object(),
                evaluate_shape_model=object(),
                move_limit=0.1,
                shape_options={"move_limit": 0.2},
            )
        finally:
            mmm.grow_hdiv_mmm_by_superposition = original
