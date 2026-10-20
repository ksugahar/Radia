"""Curve diagnostics only: a tangent jump is not a solver failure verdict."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'bh_audit', ROOT / 'validation_test/feec/bh_saturation_audit.py')
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
MU0 = 4e-7 * np.pi


def test_shared_payload_preserves_cpu_field_energy_and_barrier():
    from scipy.interpolate import PPoly
    h = np.array([0., 10., 100., 500., 2000., 10000.])
    b = MU0*(h+1e6*np.tanh(h/1000.))
    helpers = AUDIT.load_curve_helpers()
    fields, energy, _ = helpers['_bh_inverse_funcs'](h, b)
    payload = helpers['_bh_inverse_payload'](h, b)
    poly = PPoly(payload['coefficients'].reshape(-1, 4).T, payload['knots'])
    amplitudes = np.array([0., 1., .2*payload['cap'], payload['cap'], 1.01*payload['cap']])
    clipped = np.minimum(amplitudes, payload['cap'])
    over = np.maximum(amplitudes-payload['cap'], 0.)
    expected_h = poly(clipped)+payload['barrier_slope']*over
    expected_energy = poly.antiderivative()(clipped)-poly.antiderivative()(0.)
    expected_energy += payload['last_field']*over+.5*payload['barrier_slope']*over**2
    secant, _ = fields(amplitudes)
    np.testing.assert_allclose(secant*amplitudes, expected_h, rtol=2e-14, atol=1e-12)
    np.testing.assert_allclose(energy(amplitudes), expected_energy, rtol=2e-14, atol=1e-12)
    assert payload['initial_reluctivity'] == poly.derivative()(0.)
    assert not np.shares_memory(payload['coefficients'], poly.c)


def test_synthetic_table_cap_below_interpolated_peak():
    report = AUDIT.saturation_report(AUDIT.regression_table())
    assert report['interpolated_M_nondecreasing'] is False
    assert report['first_decreasing_node_interval'] == [30000, 100000]
    assert .60 < report['barrier_to_peak_ratio'] < .61
    assert 29000 < report['peak_H_Am'] < 30000
    assert report['minimum_interpolated_chi'] < -.8
    assert report['tangent_jump_ratio'] > 1e6
    assert report['inverse_error'] is None


def test_shipped_table_is_monotone_including_interpolation():
    table = np.loadtxt(ROOT / 'src/radia/panels/samples/em_sample_bh.txt')
    report = AUDIT.saturation_report(table)
    assert report['rows'] == 100
    assert report['interpolated_M_nondecreasing'] is True
    assert report['endpoint_M_Am'] == report['peak_M_Am']
    assert report['first_decreasing_node_interval'] is None


def test_linear_positive_susceptibility():
    h = np.array([0., 1., 10., 100.])
    report = AUDIT.saturation_report(np.column_stack((h, MU0 * 101 * h)))
    assert report['interpolated_M_nondecreasing'] is True
    assert report['minimum_interpolated_chi'] == pytest.approx(100.)
    assert report['peak_M_Am'] == pytest.approx(10000.)


def test_constant_zero_m_is_not_negative_susceptibility():
    h = np.array([0., 1., 10., 100.])
    report = AUDIT.saturation_report(np.column_stack((h, MU0 * h)))
    assert report['interpolated_M_nondecreasing'] is True
    # A constant M table is not invertible; this is explicit, not a solver fallback.
    assert report['inverse_error'] is not None


def test_monotone_node_m_gives_monotone_interpolated_m():
    """PCHIP of B(H) let M fall between these rising nodes (chi < -0.1); the law now interpolates M."""
    h = np.array([0., 1., 2., 3.])
    b = MU0 * np.array([0., 2., 3., 4.])
    assert np.all(np.diff(b / MU0 - h) >= -1e-14)
    report = AUDIT.saturation_report(np.column_stack((h, b)))
    assert report['first_decreasing_node_interval'] is None
    assert report['interpolated_M_nondecreasing'] is True
    assert report['minimum_interpolated_chi'] >= -1e-12
    assert report['peak_M_Am'] == pytest.approx(report['endpoint_M_Am'])


@pytest.mark.parametrize('table', [[[0, 0], [0, 1]], [[0, 0], [1, float('nan')]]])
def test_invalid_table_fails_loudly(table):
    with pytest.raises(ValueError):
        AUDIT.saturation_report(table)


def test_recorded_native_repeat_reaches_the_true_residual():
    """Replay evidence, not a substitute for running the native regression module."""
    path = ROOT / 'validation_test/feec/results/bh_saturation_20260914/qualified/result.json'
    report = json.loads(path.read_text(encoding='utf-8'))
    assert report['completed'] is True and report['passed'] is True
    assert report['exit_code'] == 0
    assert report['pip_check'] == 'No broken requirements found.'
    assert len(report['tests']) == 9
    assert all(item['outcome'] == 'passed' for item in report['tests'])
    deep_jacobi = []
    for solve in report['solves']:
        assert 'exception' not in solve
        stats = solve['statistics']
        assert stats['nonlinear_converged_final_stage'] is True
        if stats.get('nonlinear_solver') == 'picard-mass-riesz':
            continue  # This existing test exercises a different convergence contract.
        residual = stats['nonlinear_final_relative_residual']
        assert np.isfinite(residual) and 0 <= residual <= 1e-6
        if 'deep_saturation' in solve['test'] and solve['preconditioner'] == 'jacobi':
            deep_jacobi.append(solve)
    assert len(deep_jacobi) == 1
