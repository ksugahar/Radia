import numpy as np
import pytest

torch = pytest.importorskip("torch")

from radia.urn import (
    YAdmittanceURN,
    YAdmittanceURNActiveBasis,
    YAdmittanceURNConfig,
    active_basis_refit_config,
    refit_y_admittance_active_bases,
    s_domain_rmse,
    train_y_admittance_urn,
)
from radia.urn.y_admittance_urn import complex_smooth_l1


def test_y_admittance_default_dictionary_matches_research_meeting_model():
    freqs = np.logspace(2, 6, 8)
    z = 10.0 + 1j * 2.0 * np.pi * freqs * 1.0e-4
    cfg = YAdmittanceURNConfig()

    model = YAdmittanceURN(freqs, cfg, z_data=z)
    omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)

    assert cfg.total_basis_functions == 22
    assert len(model.basis_labels) == 22
    assert model.basis_matrix(omega).shape == (freqs.size, 22)

    gates = model.gates().detach().numpy()
    assert gates.shape == (22,)
    assert np.all(gates > 0.0)


def test_paper_22_basis_uses_research_meeting_training_defaults():
    cfg = YAdmittanceURNConfig.paper_22_basis(n_epochs=30)

    assert cfg.total_basis_functions == 22
    assert cfg.n_epochs == 30
    assert cfg.lr == pytest.approx(5.0e-3)
    assert cfg.sparsity_weight == pytest.approx(1.0e-4)
    assert cfg.gate_init == pytest.approx(1.0)
    assert cfg.active_threshold == pytest.approx(1.0e-3)


def test_s_domain_rmse_matches_manuscript_metric():
    z = np.array([1.0 + 1.0j, 2.0 + 0.5j, 5.0 - 0.25j])
    z_fit = z * np.array([1.0, 1.01, 0.99])

    assert s_domain_rmse(z, z) == pytest.approx(0.0)
    assert s_domain_rmse(z_fit, z) > 0.0


def test_complex_smooth_l1_weighted_matches_unweighted_for_unit_weights():
    residual = torch.tensor(
        [0.5 + 0.2j, -0.01 + 0.03j, 2.0 - 1.0j], dtype=torch.complex128
    )

    unweighted = complex_smooth_l1(residual, 1.0e-2)
    unit = complex_smooth_l1(residual, 1.0e-2, weight=torch.ones(3, dtype=torch.float64))
    zero_mid = complex_smooth_l1(
        residual, 1.0e-2, weight=torch.tensor([1.0, 0.0, 1.0], dtype=torch.float64)
    )

    assert float(unit) == pytest.approx(float(unweighted), rel=1.0e-12)
    assert float(zero_mid) > float(unweighted)  # dropping the smallest residual


def test_active_basis_refit_uses_selected_bases_without_attention():
    freqs = np.logspace(2, 5, 12)
    omega = 2.0 * np.pi * freqs
    tau = 2.0e-4
    y_true = 0.05 / (1.0 + 1j * omega * tau)
    z_true = 1.0 / y_true
    active = [
        YAdmittanceURNActiveBasis(
            basis_index=0,
            basis_type="debye",
            local_index=0,
            importance=1.0,
            gate=0.5,
            parameters={"tau": tau},
        )
    ]

    cfg = active_basis_refit_config(
        active,
        YAdmittanceURNConfig.paper_22_basis(n_epochs=0),
    )
    model = refit_y_admittance_active_bases(freqs, z_true, active, cfg, verbose=False)

    assert cfg.total_basis_functions == 1
    assert model.basis_labels == [("debye", 0)]
    summary = model.parameter_summary()[0]
    assert summary["tau"] == pytest.approx(tau)
    assert summary["gate"] == pytest.approx(0.5)


def test_initialize_from_model_embeds_smaller_dictionary_response():
    freqs = np.logspace(2, 5, 16)
    omega = 2.0 * np.pi * freqs
    y_true = 0.05 / (1.0 + 1j * omega * 2.0e-4)
    z_true = 1.0 / y_true
    small_cfg = YAdmittanceURNConfig(
        n_debye=1,
        n_magnetic_debye=0,
        n_cole_cole=0,
        n_magnetic_cole_cole=0,
        n_inductive_cpe=0,
        n_capacitive_cpe=0,
        n_series_rlc=0,
        gate_init=0.5,
    )
    large_cfg = YAdmittanceURNConfig(
        n_debye=2,
        n_magnetic_debye=0,
        n_cole_cole=0,
        n_magnetic_cole_cole=0,
        n_inductive_cpe=0,
        n_capacitive_cpe=0,
        n_series_rlc=0,
        gate_init=0.5,
    )
    small = YAdmittanceURN(freqs, small_cfg, z_data=z_true)
    large = YAdmittanceURN(freqs, large_cfg, z_data=z_true)

    large.initialize_from_model(small, extra_gate=1.0e-12)

    np.testing.assert_allclose(large.predict(freqs), small.predict(freqs), rtol=1.0e-8)


def test_y_admittance_training_smoke_single_debye_basis():
    freqs = np.logspace(2, 5, 24)
    omega = 2.0 * np.pi * freqs
    y_true = 0.05 / (1.0 + 1j * omega * 2.0e-4)
    z_true = 1.0 / y_true
    cfg = YAdmittanceURNConfig(
        n_debye=1,
        n_magnetic_debye=0,
        n_cole_cole=0,
        n_magnetic_cole_cole=0,
        n_inductive_cpe=0,
        n_capacitive_cpe=0,
        n_series_rlc=0,
        sparsity_weight=0.0,
        lr=2.0e-2,
        n_epochs=20,
        n_restarts=1,
    )

    model = train_y_admittance_urn(freqs, z_true, cfg, verbose=False)

    assert np.all(np.isfinite(model.predict(freqs)))
    assert model.training_history[0]["loss"] > model.training_history[-1]["loss"]
    active = model.active_bases(freqs)
    assert len(active) == 1
    assert active[0].basis_type == "debye"
