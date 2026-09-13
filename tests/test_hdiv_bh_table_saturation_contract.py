"""Curve diagnostics only: a tangent jump is not a solver failure verdict."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'bh_audit', ROOT / 'validation_test/feec/bh_saturation_audit.py')
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
MU0 = 4e-7 * np.pi


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


def test_monotone_node_m_does_not_guarantee_monotone_pchip_b():
    h = np.array([0., 1., 2., 3.])
    b = MU0 * np.array([0., 2., 3., 4.])
    assert np.all(np.diff(b / MU0 - h) >= -1e-14)
    report = AUDIT.saturation_report(np.column_stack((h, b)))
    assert report['first_decreasing_node_interval'] is None
    assert report['interpolated_M_nondecreasing'] is False
    assert report['minimum_interpolated_chi'] < -.1
    assert report['peak_M_Am'] > report['endpoint_M_Am']


@pytest.mark.parametrize('table', [[[0, 0], [0, 1]], [[0, 0], [1, float('nan')]]])
def test_invalid_table_fails_loudly(table):
    with pytest.raises(ValueError):
        AUDIT.saturation_report(table)
