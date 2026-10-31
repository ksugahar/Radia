"""Guards: circuit inputs, measurement identity, material/ROI/phase bookkeeping.

Assumes: CPU native from this source; explicit boxes and peak complex drive.
The mixed-material measurement contract uses q64 and the existing linear
RELATIVE_LIMIT only as an algebra-consistency budget. The companion independent
Green diagnostics observed q64 differences about 6e-12 (not an error bound).
Heat has no quadrature: its budget is an operation-count roundoff allowance.
Neither budget is a new physical criterion, surrogate tolerance or baseline.
These contracts do not establish fin formulation/field accuracy.
"""
from dataclasses import replace
import itertools
import json
import numpy as np
import pytest
from radia._residual_gate import RELATIVE_LIMIT
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


def _input_cell_roi_heat(edges, sigma, supports, axes, currents, slabs):
    """Independent box arithmetic; no production partition/intersection/cell map.

    J_i=I_i/A_i is uniform on the declared support, and perpendicular Cartesian
    components have no Joule cross term. Integrate each component over EACH
    input material cell and explicit ROI slab using its input conductivity.
    """
    power = dict.fromkeys(slabs, 0.)
    for i, support in enumerate(supports):
        area = 1.
        for axis in range(3):
            if axis != axes[i]:
                area *= support[1, axis]-support[0, axis]
        coefficient = .5*(currents[i].real**2+currents[i].imag**2)/area**2
        for index in itertools.product(range(2), repeat=3):
            cell_lo = [edges[a][index[a]] for a in range(3)]
            cell_hi = [edges[a][index[a]+1] for a in range(3)]
            for name, boxes in slabs.items():
                for roi_lo, roi_hi in boxes:
                    volume = 1.
                    for axis in range(3):
                        left = max(support[0, axis], cell_lo[axis], roi_lo[axis])
                        right = min(support[1, axis], cell_hi[axis], roi_hi[axis])
                        volume *= max(0., right-left)
                    power[name] += coefficient*volume/sigma[index]
    return power


def test_mixed_volume_region_measurement_contract(record_testsuite_property):
    # Different widths in every coordinate expose width/height interchange.
    edges = tuple(np.array(e)*1e-3 for e in (
        (-1.1, -.15, .9), (-.8, .2, 1.35), (-.4, .05, .65)))
    sigma = np.empty((2, 2, 2))
    for index in itertools.product(range(2), repeat=3):
        sigma[index] = 5.8e7 if sum(index) % 2 == 0 else 1e7
    model = build_volume_peec(edges, np.ones((2, 2, 2), dtype=bool), sigma,
        source_faces=[(0, 0, 0, 0, -1)], sink_faces=[(1, 1, 1, 2, 1)])
    drive = .7+1.1j  # Removing conjugation cannot hide behind a real drive.
    solved = solve_volume_fin(model, frequency_hz=12000.,
                              drive_current_peak_A=drive, quadrature_order=64)
    origin = (.00007, .00003, 0.)
    root_width, tip = .00027, .00083
    span = (-.00066, .00057)
    roi = FinMeasurementRegions(origin, (0., 1., 0.), (1., 0., 0.),
                                 root_width, tip, span)
    measured = measure_volume_fin(model, solved, roi, quadrature_order=64,
        convention='Declared Cartesian whole conductor; body/outside mutual retained; peak drive')
    # Define the input regions independently of ROI.intervals/point_region.
    xlo, xhi = origin[0]+span[0], origin[0]+span[1]
    y0, y1, y2 = origin[1]-root_width, origin[1]+root_width, origin[1]+tip
    slabs = {}
    for name, lower, upper in (('body', -np.inf, y0), ('root', y0, y1),
                               ('fin_shaft', y1, y2), ('tip', y2, np.inf)):
        slabs[name] = [((xlo, lower, -np.inf), (xhi, upper, np.inf))]
    slabs['outside_window'] = [((-np.inf,)*3, (xlo, np.inf, np.inf)),
                               ((xhi, -np.inf, -np.inf), (np.inf,)*3)]
    expected = _input_cell_roi_heat(edges, sigma, model.support_boxes_m, model.axes,
                                     solved.branch_current_peak_A, slabs)
    names = tuple(slabs)
    total_heat = sum(expected.values())
    # Conservative gamma_k bound for this finite, non-sliver box fixture:
    # 64 elementary operations per branch/input-cell intersection path.
    k = 64*len(model.branches)*sigma.size
    gamma_k = k*np.finfo(float).eps/(1-k*np.finfo(float).eps)
    axis_norms = [float(np.linalg.norm(solved.branch_current_peak_A[model.axes == a]))
                  for a in range(3)]
    for name, value in {
            'sigma_S_per_m': sigma.tolist(), 'drive_peak_A': [drive.real, drive.imag],
            'axis_current_norms_A': axis_norms, 'ROI': roi.as_dict(),
            'actual_heat_W': measured['mean_power_W'], 'independent_heat_W': expected,
            'heat_roundoff_gamma_k': gamma_k,
            'green_partition_difference_H': measured['green_partition_difference_H'],
            'original_partial_L_H': solved.partial_inductance_H,
            'heat_partition_difference_W': measured['heat_partition_difference_W']}.items():
        record_testsuite_property(name, json.dumps(value, allow_nan=False))
    assert all(value > 0 for value in axis_norms), axis_norms
    assert all(value > 0 for value in expected.values()), expected
    np.testing.assert_allclose(
        [measured['mean_power_W'][name] for name in names], [expected[name] for name in names],
        rtol=gamma_k, atol=gamma_k*total_heat,
        err_msg=f'Regional per-cell material/intersection heat contract: {names}')
    assert abs(measured['heat_partition_difference_W']) <= gamma_k*total_heat, (
        measured['heat_partition_difference_W'], gamma_k*total_heat)
    green_budget = RELATIVE_LIMIT*abs(solved.partial_inductance_H)
    assert abs(measured['green_partition_difference_H']) <= green_budget, (
        'Original restricted Green partition contract',
        measured['green_partition_difference_H'], green_budget)
