"""Synthetic contracts for the coefficient-driven envelope Stop core."""

import numpy as np
import pytest

from radia.vim import EnvelopeStopVector3D
from radia.vim import _envelope_stop as es


def _theta():
    rng = np.random.default_rng(20260828)
    history = rng.uniform(0.0, 1.0, es.Pqm)
    history *= rng.uniform(size=es.Pqm) < 0.02
    history *= 40.0 / max(float(history.sum()), 1e-30)
    static = np.array([60.0, 40.0, 25.0, 15.0, 10.0, 6.0, 4.0,
                       3.0, 2.0, 1.5, 1.0, 0.5, 0.25])
    return np.concatenate([history, static])


THETA = _theta()


def _model(dim=1):
    return EnvelopeStopVector3D(
        THETA, b_max_T=1.25, dim=dim,
        allow_unvalidated_rotation=dim > 1)


def _drive(model, flux):
    state = model.state0()[None, :]
    field = np.zeros_like(flux)
    for index, value in enumerate(flux):
        field[index] = model.forward(value[None, :], state)[0]
        state = model.commit(value[None, :], state)
    return field


def test_requires_explicit_coefficients_range_and_rotation_opt_in():
    with pytest.raises(TypeError):
        EnvelopeStopVector3D()
    with pytest.raises(ValueError, match="b_max_T"):
        EnvelopeStopVector3D(THETA, b_max_T=float("nan"))
    with pytest.raises(ValueError, match="not validated"):
        EnvelopeStopVector3D(THETA, b_max_T=1.25, dim=3)
    model = _model(dim=3)
    assert model.rotational_validation == "unvalidated-explicit-opt-in"
    assert model.b_max() == pytest.approx(1.25)


@pytest.mark.parametrize("bad", [np.zeros(10), np.full(es.DPAR, np.nan), -np.ones(es.DPAR)])
def test_coefficient_schema_fails_loudly(bad):
    with pytest.raises(ValueError):
        EnvelopeStopVector3D(bad, b_max_T=1.0)


def test_uniaxial_oddness_restart_and_scalar_design_agree():
    model = _model()
    flux = np.concatenate([
        np.linspace(0.0, 1.0, 50), np.linspace(1.0, -1.0, 100),
        np.linspace(-1.0, 0.7, 80)])[:, None]
    field = _drive(model, flux)[:, 0]
    np.testing.assert_allclose(_drive(model, -flux)[:, 0], -field, atol=1e-9)
    np.testing.assert_allclose(field, es.envelope_stop_scalar_H(flux[:, 0], THETA), atol=1e-9)

    cut = 77
    state = model.state0()[None, :]
    for value in flux[:cut]:
        state = model.commit(value[None, :], state)
    tail = _drive_from_state(model, flux[cut:], state)
    np.testing.assert_allclose(tail[:, 0], field[cut:], atol=1e-12)


def _drive_from_state(model, flux, state):
    field = np.zeros_like(flux)
    for index, value in enumerate(flux):
        field[index] = model.forward(value[None, :], state)[0]
        state = model.commit(value[None, :], state)
    return field


def test_analytic_tangent_matches_independent_finite_difference():
    model = _model(dim=3)
    rng = np.random.default_rng(31)
    flux = rng.normal(size=(5, 3))
    flux *= rng.uniform(0.1, 1.0, (5, 1)) / np.linalg.norm(flux, axis=1)[:, None]
    states = np.tile(model.state0(), (5, 1))
    analytic = model.tangent_analytic(flux, states)
    finite_difference = model.tangent(flux, states)
    relative = (np.linalg.norm(analytic - finite_difference, axis=(1, 2)) /
                np.maximum(np.linalg.norm(finite_difference, axis=(1, 2)), 1e-30))
    assert float(relative.max()) < 1e-5


def test_first_law_defect_shrinks_with_step_refinement():
    model = _model()
    minima = []
    for samples in (120, 480):
        flux = (0.8 * np.sin(np.linspace(0.0, 4.0 * np.pi, samples)))[:, None]
        state = model.state0()[None, :]
        previous_h = previous_energy = None
        minimum = np.inf
        for index, value in enumerate(flux):
            field = model.forward(value[None, :], state)[0]
            energy = model.energy(value[None, :], state)[0]
            if previous_h is not None:
                work = 0.5 * (field + previous_h) @ (value - flux[index - 1])
                minimum = min(minimum, float(work - (energy - previous_energy)))
            previous_h, previous_energy = field, energy
            state = model.commit(value[None, :], state)
        minima.append(minimum)
    assert minima[1] >= minima[0] - 1e-12
    assert abs(minima[1]) < 0.35 * abs(minima[0])


def test_static_only_inverse_round_trip_and_range_gate():
    theta = np.concatenate([np.zeros(es.Pqm), THETA[es.Pqm:]])
    model = EnvelopeStopVector3D(theta, b_max_T=0.9)
    state = model.state0()[None, :]
    target_b = np.array([[0.44]])
    target_h = model.forward(target_b, state)
    recovered, permeability = model.inverse(target_h, state, B0=0.8 * target_b, tol=1e-10)
    np.testing.assert_allclose(recovered, target_b, atol=1e-8)
    assert np.all(np.isfinite(permeability))

    outside_h = model.forward(np.array([[1.0]]), state)
    with pytest.raises(RuntimeError, match="validated"):
        model.inverse(outside_h, state, B0=np.array([[0.8]]), tol=1e-10)
