import numpy as np
import pytest

from radia.harmonic_balance import (
    hysteresis_cycle_metrics,
    periodic_phase,
    project_odd_sine_harmonics,
    solve_odd_harmonic_balance,
    synthesize_odd_sine_series,
)
from radia.planar_materials import PlayHysteresis


def test_odd_sine_projection_is_exact_on_an_integer_period_window():
    harmonics = (1, 3, 5, 7)
    expected = np.array([1.2, -0.35, 0.08, 0.015])
    phase = periodic_phase(2048, period_count=3)
    waveform = synthesize_odd_sine_series(expected, harmonics, phase)

    actual = project_odd_sine_harmonics(
        waveform,
        harmonics,
        samples_per_period=2048,
        period_count=3,
    )

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-15)


@pytest.mark.parametrize("harmonics", [(0, 1), (1, 2), (3, 1), (1, 1)])
def test_projection_rejects_ambiguous_or_non_odd_harmonics(harmonics):
    with pytest.raises(ValueError):
        project_odd_sine_harmonics(
            np.zeros(32), harmonics, samples_per_period=32
        )


def test_projection_rejects_a_duplicated_period_endpoint():
    with pytest.raises(ValueError, match="do not duplicate"):
        project_odd_sine_harmonics(
            np.zeros(33), (1, 3), samples_per_period=32
        )


def test_nonlinear_material_requires_odd_harmonics_beyond_the_fundamental():
    harmonics = (1, 3, 5, 7)
    h_fundamental = 7.25e-4

    result = solve_odd_harmonic_balance(
        lambda h: 2.0 * np.tanh(h / 1.0e-3),
        lambda _b: np.array([h_fundamental, 0.0, 0.0, 0.0]),
        np.array([h_fundamental, 0.0, 0.0, 0.0]),
        harmonics=harmonics,
        damping=1.0,
        relative_tolerance=1.0e-13,
        samples_per_period=4096,
    )

    assert result.converged
    assert result.iterations == 1
    assert abs(result.b_coefficients[2]) > 1.0e-3
    assert abs(result.b_coefficients[3]) > 1.0e-4
    fundamental_only = synthesize_odd_sine_series(
        np.array([result.b_coefficients[0]]), (1,), periodic_phase(4096)
    )
    full_waveform = 2.0 * np.tanh(h_fundamental * np.sin(periodic_phase(4096)) / 1.0e-3)
    assert np.linalg.norm(full_waveform - fundamental_only) / np.linalg.norm(full_waveform) > 0.03


def test_harmonic_balance_reports_nonconvergence_instead_of_nominal_mode_count():
    with pytest.raises(RuntimeError, match="did not converge"):
        solve_odd_harmonic_balance(
            lambda h: h,
            lambda b: b + 1.0,
            np.zeros(2),
            harmonics=(1, 3),
            max_iterations=3,
        )


def test_play_state_is_committed_only_after_acceptance_and_cycle_loss_is_positive():
    eta = np.arange(10, dtype=float) * 300.0
    play = PlayHysteresis(eta=eta, w=np.full(10, 140.0))
    committed = play.fresh_state(1)
    trial_field = np.array([1800.0])

    first_trial = play.M(trial_field, committed)
    second_trial = play.M(trial_field, committed)
    np.testing.assert_array_equal(committed, play.fresh_state(1))
    np.testing.assert_allclose(second_trial, first_trial, rtol=0.0, atol=0.0)

    peak = 2400.0
    field = np.concatenate(
        [
            np.linspace(peak, -peak, 129, endpoint=False),
            np.linspace(-peak, peak, 129),
        ]
    )
    flux = []
    for value in field:
        h = np.array([value])
        flux.append(4.0e-7 * np.pi * play.M(h, committed)[0])
        committed = play.advance(h, committed)
    metrics = hysteresis_cycle_metrics(field, np.asarray(flux))

    assert metrics.passive_orientation
    assert metrics.signed_energy_density_j_per_m3 > 1000.0
    assert metrics.closure_relative < 0.02


def test_cycle_metrics_preserve_orientation_instead_of_taking_absolute_value():
    field = np.array([1.0, 0.0, -1.0, 0.0, 1.0])
    flux = np.array([0.0, 1.0, 0.0, -1.0, 0.0])
    forward = hysteresis_cycle_metrics(field, flux)
    reverse = hysteresis_cycle_metrics(field[::-1], flux[::-1])

    assert forward.signed_energy_density_j_per_m3 == pytest.approx(
        -reverse.signed_energy_density_j_per_m3
    )
    assert forward.passive_orientation is not reverse.passive_orientation
