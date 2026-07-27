import numpy as np
import pytest

torch = pytest.importorskip("torch")

from radia.urn import (
    CauerLadderURN,
    CauerLadderURNConfig,
    fit_rational_pole_zero,
    train_cauer_ladder_alternating,
    train_cauer_ladder_progressive,
    train_cauer_ladder_tail_then_polish,
)


def _one_section_target(freqs, *, z_ref=50.0):
    omega = 2.0 * np.pi * freqs
    omega_ref = np.sqrt(np.min(omega) * np.max(omega))
    x = 1j * omega / omega_ref
    zbar = 0.08 + 0.20 * x + 1.0 / (0.04 + 0.30 * x)
    return z_ref * zbar, omega_ref


def test_cauer_ladder_is_autodifferentiable():
    freqs = np.logspace(2, 5, 12)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(n_sections=2, omega_ref=omega_ref, z_ref=50.0)
    model = CauerLadderURN(freqs, cfg, z_data=z)
    omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
    z_target = torch.tensor(z, dtype=torch.complex128)

    total, z_loss, y_loss, _parts = model.zy_losses(omega, z_target)
    total.backward()

    assert total.detach().item() > 0.0
    assert z_loss.detach().item() > 0.0
    assert y_loss.detach().item() > 0.0
    for param in model.parameters():
        assert param.grad is not None
        assert torch.all(torch.isfinite(param.grad))


def test_cauer_alternating_training_smoke():
    freqs = np.logspace(2, 5, 16)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(
        n_sections=1,
        omega_ref=omega_ref,
        z_ref=50.0,
        lr=1.0e-2,
        n_restarts=1,
        alternating_cycles=6,
        alternating_block_epochs=4,
        regularization_weight=0.0,
        log_smoothness_weight=0.0,
    )

    model = train_cauer_ladder_alternating(freqs, z, cfg, verbose=False)

    assert np.all(np.isfinite(model.predict(freqs)))
    assert model.training_history[0]["loss"] >= model.training_history[-1]["loss"]


def test_cauer_peeling_initialization_is_finite():
    freqs = np.logspace(2, 5, 16)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(
        n_sections=2,
        omega_ref=omega_ref,
        z_ref=50.0,
        use_peeling_initialization=True,
    )
    model = CauerLadderURN(freqs, cfg, z_data=z)

    model.initialize_from_peeling(z)

    assert np.all(np.isfinite(model.predict(freqs)))
    for item in model.parameter_summary():
        assert item["R_ohm"] >= 0.0
        assert item["L_h"] >= 0.0
        assert item["G_siemens"] >= 0.0
        assert item["C_f"] >= 0.0


def test_rational_pole_zero_fit_predicts_finite_response():
    freqs = np.logspace(2, 5, 16)
    z, omega_ref = _one_section_target(freqs)

    fit = fit_rational_pole_zero(freqs, z, order=1, omega_ref=omega_ref, z_ref=50.0)
    pred = fit.predict(freqs)

    assert np.all(np.isfinite(pred))
    assert fit.poles.size == 1
    assert fit.zeros.size == 1


def test_cauer_rational_initialization_reduces_one_section_error():
    freqs = np.logspace(2, 5, 24)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(
        n_sections=1,
        omega_ref=omega_ref,
        z_ref=50.0,
        rational_order=2,
        rational_cauer_max_nfev=300,
    )
    model = CauerLadderURN(freqs, cfg, z_data=z)
    default_error = np.linalg.norm(model.predict(freqs) - z)

    model.initialize_from_rational_fit(z)
    initialized_error = np.linalg.norm(model.predict(freqs) - z)

    assert np.all(np.isfinite(model.predict(freqs)))
    assert initialized_error < default_error


def test_cauer_least_squares_polish_reduces_one_section_error():
    freqs = np.logspace(2, 5, 24)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(
        n_sections=1,
        omega_ref=omega_ref,
        z_ref=50.0,
        least_squares_max_nfev=200,
    )
    model = CauerLadderURN(freqs, cfg, z_data=z)
    default_error = np.linalg.norm(model.predict(freqs) - z)

    model.fit_to_response_least_squares(z)
    polished_error = np.linalg.norm(model.predict(freqs) - z)

    assert polished_error < default_error


def test_cauer_tail_then_polish_training_smoke():
    freqs = np.logspace(2, 5, 16)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(
        n_sections=2,
        omega_ref=omega_ref,
        z_ref=50.0,
        lr=8.0e-3,
        n_restarts=1,
        alternating_cycles=2,
        alternating_block_epochs=2,
        frozen_outer_sections=1,
        tail_train_cycles=2,
        polish_train_cycles=2,
        regularization_weight=0.0,
        log_smoothness_weight=0.0,
    )

    model = train_cauer_ladder_tail_then_polish(freqs, z, cfg, verbose=False)

    assert model.config.n_sections == 2
    assert np.all(np.isfinite(model.predict(freqs)))
    assert {item["stage"] for item in model.training_history} == {"tail", "polish"}


def test_cauer_progressive_grows_to_requested_depth():
    freqs = np.logspace(2, 5, 12)
    z, omega_ref = _one_section_target(freqs)
    cfg = CauerLadderURNConfig(
        n_sections=2,
        omega_ref=omega_ref,
        z_ref=50.0,
        lr=8.0e-3,
        n_restarts=1,
        progressive_start_sections=1,
        progressive_stage_epochs=1,
        alternating_cycles=2,
        alternating_block_epochs=2,
        regularization_weight=0.0,
        log_smoothness_weight=0.0,
    )

    model = train_cauer_ladder_progressive(freqs, z, cfg, verbose=False)

    assert model.config.n_sections == 2
    assert model.config.total_parameters == 8
    assert np.all(np.isfinite(model.predict(freqs)))
