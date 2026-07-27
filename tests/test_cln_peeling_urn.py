import numpy as np
import pytest

pytest.importorskip("torch")

from radia.urn import CLNPeelingConfig, train_cln_peeling_urn
from radia.urn.cln_peeling_urn import (
    _compose_stage,
    _ensure_nonempty_mask,
    _fit_composite_seed,
    _fit_paired_stage,
    _peel_tail,
    _tail_trust_weight,
    _weighted_s_domain_rmse,
)


def test_cln_exact_peel_inverts_stage_composition():
    series = np.array([2.0 + 0.2j, 3.0 + 0.4j, 4.0 + 0.8j])
    shunt = np.array([8.0 + 0.5j, 9.0 + 0.7j, 10.0 + 0.9j])
    tail = np.array([20.0 + 1.0j, 25.0 + 1.5j, 30.0 + 2.0j])

    driving_point = _compose_stage(series, shunt, tail)
    peeled, peeled_y = _peel_tail(driving_point, series, shunt)

    np.testing.assert_allclose(peeled, tail, rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(peeled_y, 1.0 / tail, rtol=1.0e-12, atol=1.0e-12)


@pytest.mark.parametrize(
    ("initial", "expected_true", "expected_false"),
    [
        (np.zeros(4, dtype=bool), 1, 3),
        (np.ones(4, dtype=bool), 3, 1),
    ],
)
def test_cln_split_mask_keeps_both_branches(initial, expected_true, expected_false):
    logits = np.array([-2.0, -1.0, 1.0, 2.0])

    mask = _ensure_nonempty_mask(initial, logits)

    assert np.count_nonzero(mask) >= expected_true
    assert np.count_nonzero(~mask) >= expected_false


def test_cln_peeling_smoke_freezes_even_and_odd_branches():
    freqs = np.logspace(2, 5, 14)
    omega = 2.0 * np.pi * freqs
    z_series = 2.0 + 1j * omega * 1.0e-5
    z_shunt = 1.0 / (0.02 / (1.0 + 1j * omega * 8.0e-5))
    z_tail = 30.0 + 0.0j
    z = _compose_stage(z_series, z_shunt, z_tail)
    cfg = CLNPeelingConfig(
        n_stages=1,
        branch_epochs=20,
        branch_lr=8.0e-3,
        branch_restarts=1,
        pair_epochs=30,
        pair_polish_epochs=10,
        pair_lr=8.0e-3,
        pair_restarts=1,
        branch_sparsity_weight=0.0,
        positive_real_weight=0.1,
        min_tail_sensitivity=0.0,
        residual_real_tolerance=100.0,
        max_stage_relative_degradation=10.0,
        seed=37,
    )

    model = train_cln_peeling_urn(freqs, z, cfg, verbose=False)

    assert len(model.stages) == 1
    stage = model.stages[0]
    assert stage.accepted
    assert np.all(stage.split_fraction > 0.0)
    assert np.all(stage.split_fraction < 1.0)
    assert stage.series_impedance.active_count() > 0
    assert stage.shunt_impedance.active_count() > 0

    pred_stored = model.predict(freqs)
    pred_lookahead = model.predict_terminated(freqs, termination="lookahead")
    np.testing.assert_allclose(pred_stored, z, rtol=1.0e-9, atol=1.0e-9)
    assert np.all(np.isfinite(pred_lookahead))
    assert np.isfinite(model.s_domain_rmse_terminated(z, termination="lookahead"))


@pytest.fixture(scope="module")
def smoke_ladder_model():
    freqs = np.logspace(2, 5, 14)
    omega = 2.0 * np.pi * freqs
    z_series = 2.0 + 1j * omega * 1.0e-5
    z_shunt = 1.0 / (0.02 / (1.0 + 1j * omega * 8.0e-5))
    z = _compose_stage(z_series, z_shunt, np.full(freqs.shape, 30.0 + 0.0j))
    cfg = CLNPeelingConfig(
        n_stages=1,
        branch_epochs=20,
        branch_lr=8.0e-3,
        branch_restarts=1,
        pair_epochs=30,
        pair_polish_epochs=10,
        pair_lr=8.0e-3,
        pair_restarts=1,
        branch_sparsity_weight=0.0,
        positive_real_weight=0.1,
        min_tail_sensitivity=0.0,
        residual_real_tolerance=100.0,
        max_stage_relative_degradation=10.0,
        min_trusted_fraction=0.0,
        seed=37,
    )
    model = train_cln_peeling_urn(freqs, z, cfg, verbose=False)
    assert len(model.stages) == 1
    return freqs, z, model


def test_cln_tail_trust_weight_flags_cancellation_and_inactive_tail():
    config = CLNPeelingConfig()
    z_ref = 50.0
    target = np.array([100.0 + 10.0j, 30000.0 + 20000.0j, 80.0 + 5.0j])
    series = np.array([10.0 + 1.0j, 1500.0 + 300.0j, 79.9 + 4.99j])
    tail_y = np.array([1.0e-2 + 1.0e-3j, -1.7e-5 + 1.4e-7j, 2.0e-2 + 0.0j])

    weight = _tail_trust_weight(target, series, tail_y, z_ref, config)

    assert weight[0] > 0.9
    # parallel-resonance band: |Y_tail|*z_ref = 8.5e-4 << margin 5e-3
    assert weight[1] < 0.05
    # series branch nearly cancels the target
    assert weight[2] < 0.5


def test_cln_weighted_s_domain_rmse_ignores_zero_weight_points():
    target = np.array([50.0 + 0.0j, 60.0 + 0.0j, 70.0 + 0.0j])
    fit = target.copy()
    fit[1] = -1000.0 + 0.0j
    weight = np.array([1.0, 0.0, 1.0])

    assert _weighted_s_domain_rmse(fit, target, weight) < 1.0e-12
    assert _weighted_s_domain_rmse(fit, target, None) > 0.1


def test_cln_untrusted_spike_points_do_not_block_acceptance():
    freqs = np.logspace(2, 5, 14)
    omega = 2.0 * np.pi * freqs
    z_series = 2.0 + 1j * omega * 1.0e-5
    z_shunt = 1.0 / (0.02 / (1.0 + 1j * omega * 8.0e-5))
    z = _compose_stage(z_series, z_shunt, np.full(freqs.shape, 30.0 + 0.0j))
    spike = 7
    z[spike] = -3000.0 + 0.5j  # non-passive spike no passive dictionary can fit
    weight = np.ones(freqs.shape)
    weight[spike] = 0.0
    cfg = CLNPeelingConfig(
        n_stages=1,
        branch_epochs=25,
        branch_lr=8.0e-3,
        branch_restarts=1,
        pair_epochs=30,
        pair_polish_epochs=10,
        pair_lr=8.0e-3,
        pair_restarts=1,
        branch_sparsity_weight=0.0,
        positive_real_weight=0.1,
        min_tail_sensitivity=0.0,
        residual_real_tolerance=2.0,
        max_stage_relative_degradation=10.0,
        min_trusted_fraction=0.0,
        seed=53,
    )

    composite = _fit_composite_seed(
        freqs, z, cfg, seed_offset=0, sample_weight=weight
    )
    trusted_stage = _fit_paired_stage(
        freqs, z, composite, cfg, stage_index=0, sample_weight=weight
    )
    assert (
        trusted_stage.metrics["min_parallel_real_normalized"]
        < -cfg.residual_real_tolerance
    )
    assert (
        trusted_stage.metrics["min_parallel_real_trusted"]
        >= -cfg.residual_real_tolerance
    )
    assert trusted_stage.accepted

    unweighted_stage = _fit_paired_stage(
        freqs, z, composite, cfg, stage_index=0, sample_weight=None
    )
    assert not unweighted_stage.accepted


def test_cln_stage_records_trust_weights(smoke_ladder_model):
    _freqs, _z, model = smoke_ladder_model
    stage = model.stages[0]

    np.testing.assert_allclose(stage.sample_weight, 1.0)
    assert np.all(stage.tail_trust_weight > 0.0)
    assert np.all(stage.tail_trust_weight <= 1.0)
    assert 0.0 <= stage.metrics["trusted_fraction"] <= 1.0


def test_cln_composite_warm_start_reproduces_previous_lookahead():
    from dataclasses import replace

    freqs = np.logspace(2, 5, 16)
    omega = 2.0 * np.pi * freqs
    z_series = 2.0 + 1j * omega * 1.0e-5
    z_shunt = 1.0 / (0.02 / (1.0 + 1j * omega * 8.0e-5))
    z = _compose_stage(z_series, z_shunt, np.full(freqs.shape, 30.0 + 0.0j))
    cfg = CLNPeelingConfig(
        branch_epochs=10,
        branch_restarts=1,
        branch_sparsity_weight=0.0,
        seed=61,
    )

    donor = _fit_composite_seed(freqs, z, cfg, seed_offset=0)
    # Different target scale -> different y_ref; the gate rescaling in the
    # warm start must still reproduce the donor branch exactly at epoch 0.
    other_target = 2.5 * z
    frozen_cfg = replace(cfg, branch_epochs=0)
    seeded = _fit_composite_seed(
        freqs,
        other_target,
        frozen_cfg,
        seed_offset=1,
        warm_start_model=donor.model,
    )

    np.testing.assert_allclose(
        seeded.response(freqs), donor.response(freqs), rtol=1.0e-9
    )


def test_cln_two_stage_training_with_global_objective_runs():
    freqs = np.logspace(2, 5, 14)
    omega = 2.0 * np.pi * freqs
    inner = _compose_stage(
        0.5 + 1j * omega * 3.0e-6,
        1.0 / (0.05 / (1.0 + 1j * omega * 2.0e-5)),
        np.full(freqs.shape, 12.0 + 0.0j),
    )
    z = _compose_stage(
        2.0 + 1j * omega * 1.0e-5,
        1.0 / (0.02 / (1.0 + 1j * omega * 8.0e-5)),
        inner,
    )
    cfg = CLNPeelingConfig(
        n_stages=2,
        branch_epochs=15,
        branch_restarts=1,
        pair_epochs=20,
        pair_polish_epochs=5,
        pair_restarts=1,
        branch_sparsity_weight=0.0,
        positive_real_weight=0.1,
        min_tail_sensitivity=0.0,
        residual_real_tolerance=100.0,
        max_stage_relative_degradation=10.0,
        min_trusted_fraction=0.0,
        seed=71,
    )
    assert cfg.global_objective and cfg.composite_warm_start  # defaults

    model = train_cln_peeling_urn(freqs, z, cfg, verbose=False)

    # Stage 0 must freeze; stage 1 exercises the through-ladder objective
    # (frozen outer chain composition) whether or not it survives the gates.
    assert len(model.stages) >= 1
    assert "global_lookahead_s_rmse" in model.stages[0].metrics
    prediction = model.predict_terminated(freqs, termination="lookahead")
    assert np.all(np.isfinite(prediction.real))
    assert np.all(np.isfinite(prediction.imag))


def test_cln_frozen_stage_records_global_lookahead_rmse(smoke_ladder_model):
    freqs, z, model = smoke_ladder_model
    stage = model.stages[0]

    recorded = stage.metrics["global_lookahead_s_rmse"]
    recomputed = model.s_domain_rmse_terminated(z, termination="lookahead")
    assert recorded == pytest.approx(recomputed, rel=1.0e-9)
    # A frozen stage is never the one rejected by the global stopping rule.
    assert "rejected_by_global_degradation" not in stage.metrics


def test_cln_predict_terminated_supports_dense_grid(smoke_ladder_model):
    freqs, _z, model = smoke_ladder_model
    dense = np.logspace(1.5, 5.5, 41)

    for termination in ("lookahead", "constant", "open", "short"):
        prediction = model.predict_terminated(dense, termination=termination)
        assert prediction.shape == dense.shape
        assert np.all(np.isfinite(prediction.real))
        assert np.all(np.isfinite(prediction.imag))

    with pytest.raises(ValueError, match="training grid"):
        model.predict_terminated(dense, termination="stored")
    on_grid = model.predict_terminated(freqs, termination="stored")
    assert np.all(np.isfinite(on_grid.real))


def test_cln_audit_passivity_dense_grid(smoke_ladder_model):
    _freqs, _z, model = smoke_ladder_model

    report = model.audit_passivity(points_per_decade=16, extrapolation_decades=1.0)

    assert report["passive"] is True
    assert report["n_points"] >= 16 * 5
    assert len(report["branches"]) == 3 * len(model.stages)
    for branch in report["branches"]:
        assert branch["finite"]
        assert branch["min_re_z"] >= 0.0
        assert branch["min_re_y"] >= 0.0
    assert report["ladder_lookahead"]["finite"]
    assert report["ladder_lookahead"]["min_re_z"] >= 0.0


def test_cln_rejected_split_does_not_freeze_stage():
    freqs = np.logspace(2, 4, 10)
    z = np.full(freqs.shape, 5.0 + 0.0j)
    cfg = CLNPeelingConfig(
        n_stages=1,
        branch_epochs=4,
        pair_epochs=4,
        pair_polish_epochs=2,
        pair_restarts=1,
        min_tail_sensitivity=2.0,
        seed=41,
    )

    model = train_cln_peeling_urn(freqs, z, cfg, verbose=False)

    assert model.stages == []
    assert len(model.training_log) == 1
    assert model.training_log[0]["accepted"] is False
    np.testing.assert_allclose(model.final_tail_impedance, z)
