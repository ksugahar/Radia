"""Guards: explicit finite-volume circuit inputs and measurement identity.

Assumes: CPU native from this source; no physical oracle or stored baseline.
These are public fail-loud contracts, not fin accuracy or surrogate claims.
"""
from dataclasses import replace
import numpy as np
import pytest
from radia.volume_peec import build_volume_peec
from radia.volume_fin import solve_volume_fin, measure_volume_fin
from radia.fin_observables import FinMeasurementRegions


def topology(sigma=5.8e7):
    return build_volume_peec(([0., .001],)*3, np.ones((1, 1, 1), dtype=bool), sigma,
        source_faces=[(0, 0, 0, 0, -1)], sink_faces=[(0, 0, 0, 0, 1)])


@pytest.mark.parametrize('frequency,drive', [(-1., 1.), (float('nan'), 1.), (1., 0.), (1., complex(0, float('inf')))])
def test_invalid_circuit_request_raises(frequency, drive):
    with pytest.raises(ValueError, match='Finite'):
        solve_volume_fin(topology(), frequency_hz=frequency,
                         drive_current_peak_A=drive, quadrature_order=32)


def test_measurement_rejects_wrong_material_identity_and_rotated_roi():
    model = topology()
    solution = solve_volume_fin(model, frequency_hz=5000.,
                                drive_current_peak_A=1., quadrature_order=32)
    roi = FinMeasurementRegions((0., 0., 0.), (0., 1., 0.), (1., 0., 0.),
                                 .0002, .0008, (0., .001))
    assert not solution.branch_current_peak_A.flags.writeable
    assert not solution.inductance_H.flags.writeable
    with pytest.raises(ValueError, match='identity mismatch'):
        measure_volume_fin(topology(1e7), solution, roi, quadrature_order=32, convention='explicit')
    rotated = replace(roi, fin_axis=(2**-.5, 2**-.5, 0.), span_axis=(0., 0., 1.))
    with pytest.raises(ValueError, match='Cartesian'):
        measure_volume_fin(model, solution, rotated, quadrature_order=32, convention='explicit')
    with pytest.raises(ValueError, match='convention'):
        measure_volume_fin(model, solution, roi, quadrature_order=32, convention='')


@pytest.mark.parametrize('nodes', [np.array([[-1, 0], [0, 2]]), np.array([[3, 0], [0, 2]])])
def test_invalid_connectivity_rejected_before_native(nodes):
    model = replace(topology(), branches=nodes)
    with pytest.raises(ValueError, match='branch node index'):
        solve_volume_fin(model, frequency_hz=5000., drive_current_peak_A=1., quadrature_order=32)
