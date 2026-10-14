import numpy as np

from radia_mcp.radia_ngsolve.cq_urn import (
    cq_convolve,
    cq_weights_from_laplace,
    cq_time_grid_contract_gate,
    fit_nonnegative_debye,
    make_cq_urn_bridge_artifact,
    periodic_ifft_response,
    relaxation_response,
)
from radia_mcp.radia_ngsolve.knowledge.urn import get_urn_documentation


def test_nonnegative_debye_fit_recovers_sparse_passive_ladder():
    freq_hz = np.logspace(-1.0, 2.0, 80)
    tau_grid = np.array([0.005, 0.01, 0.02, 0.05, 0.1, 0.35, 0.7])
    target = relaxation_response(
        1j * 2.0 * np.pi * freq_hz,
        weights=[0.62, 0.30],
        taus=[0.02, 0.35],
        feedthrough=0.08,
    )

    fit = fit_nonnegative_debye(freq_hz, target, tau_grid, active_threshold=1.0e-6)

    assert fit.relative_error < 1.0e-10
    assert fit.feedthrough >= 0.0
    assert all(w >= -1.0e-14 for w in fit.weights)
    assert fit.active_count == 3


def test_cq_bridge_is_causal_where_periodic_ifft_wraps():
    n_steps = 100
    dt = 0.01
    hit_index = 10

    def laplace(s):
        return relaxation_response(s, weights=[0.62, 0.30], taus=[0.02, 0.35], feedthrough=0.08)

    weights = cq_weights_from_laplace(laplace, dt, n_steps)
    signal = np.zeros(n_steps)
    signal[hit_index:] = 1.0

    y_cq = cq_convolve(weights, signal).real
    y_ifft = periodic_ifft_response(laplace, dt, signal).real

    assert np.max(np.abs(y_cq[:hit_index])) < 1.0e-12
    assert np.max(np.abs(y_ifft[:hit_index])) > 1.0e-3


def test_cq_time_grid_contract_accepts_padding_choices_and_rejects_undersampling():
    matlab_style = cq_time_grid_contract_gate(
        1.0e-3, 100, method="bdf2", contour_samples=100
    )
    padded = cq_time_grid_contract_gate(
        1.0e-3, 100, method="bdf2", contour_samples=256
    )
    assert matlab_style["status"] == "ok"
    assert padded["status"] == "ok"
    assert matlab_style["time_end"] == 0.099
    assert matlab_style["min_real_laplace_node"] > 0.0

    undersampled = cq_time_grid_contract_gate(
        1.0e-3, 100, method="bdf2", contour_samples=64, radius=0.9
    )
    assert undersampled["status"] == "needs_attention"
    assert undersampled["checks"]["contour_covers_time_grid"] is False


def test_teaching_artifact_has_feedback_and_output_contract():
    artifact = make_cq_urn_bridge_artifact(n_steps=80, dt=0.01, hit_index=8)

    assert artifact["pass"] is True
    assert artifact["schema"] == "radia.cq_urn_bridge.v1"
    assert artifact["result_output_schema_id"] == "radia.cq_urn_bridge.timeseries.v1"
    assert artifact["mcp_feedback"]["learning_lanes"]["public"] == "encoded"
    assert artifact["checks"]["cq_causal_before_hit"] is True
    assert artifact["checks"]["ifft_periodic_has_prehit_wraparound"] is True


def test_urn_knowledge_exposes_cq_topic():
    text = get_urn_documentation("cq")

    assert "convolution quadrature" in text
    assert "H(s)" in text
    assert "cq_urn_bridge.ipynb" in text
