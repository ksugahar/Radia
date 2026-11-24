"""Guards: common-grid back-EMF comparison and explicit lane diagnostics.

Assumes: uniformly sampled synthetic periodic three-phase linkage signals.
"""

import math

import numpy as np
import pytest

from radia_mcp.motor.back_emf_cross_lane_gate import back_emf_cross_lane_gate


LIMITS = {'grid_deg_e': 1e-7, 'amplitude_relative': 0.01, 'phase_deg_e': 0.2,
          'rms_relative': 0.01, 'thd_absolute': 0.01}


def _lane():
    theta = np.arange(48) * 2 * np.pi / 48
    phase = theta[:, None] + np.radians([0, -120, 120])
    return {
        'mechanical_angles_deg': np.degrees(theta / 3).tolist(),
        'phase_flux_wb': (0.08 * np.cos(phase)).tolist(),
        'pole_pairs': 3, 'electrical_offset_deg': 0, 'speed_rpm': 900,
        'stack_basis': 'per_metre', 'stack_length_m': 0.17,
        'linkage_polarity': 'positive_winding', 'model_dimension': '2d',
        'emf_method': 'central_difference_of_linkage',
    }


def test_equivalent_time_grid_and_stack_basis():
    """Guards: equivalence across time/angle and total/per-metre inputs.

    Assumes: constant positive mechanical speed and one electrical period.
    """
    first, second = _lane(), _lane()
    second['times_s'] = [angle / (900 * 6) for angle in second.pop('mechanical_angles_deg')]
    second.update(stack_basis='total', already_scaled=True)
    second['phase_flux_wb'] = (np.asarray(second['phase_flux_wb']) * second['stack_length_m']).tolist()
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'pass'
    assert result['metrics_per_metre']['linkage'][0]['amplitude'] == pytest.approx([0.08] * 3)
    assert result['metrics_per_metre']['emf'][0]['phase_deg_e'] == pytest.approx([-90, 150, 30])


def test_harmonic_rms_and_thd_match_analytic_sampled_signal():
    """Guards: THD excludes DC and uses RMS harmonic energy by phase.

    Assumes: sampled cosine fundamental and third harmonic on a complete grid.
    """
    lane = _lane()
    theta = np.radians(3 * np.asarray(lane['mechanical_angles_deg']))[:, None]
    phases = theta + np.radians([0, -120, 120])
    lane['phase_flux_wb'] = (0.01 + 0.08 * np.cos(phases) + 0.02 * np.cos(3 * phases)).tolist()
    result = back_emf_cross_lane_gate(lane, lane, LIMITS)
    assert result['status'] == 'pass'
    flux = result['metrics_per_metre']['linkage'][0]
    emf = result['metrics_per_metre']['emf'][0]
    assert flux['rms'] == pytest.approx([math.sqrt(0.01**2 + (0.08**2 + 0.02**2) / 2)] * 3)
    assert flux['thd'] == pytest.approx([0.25] * 3)
    step = math.tau / len(theta)
    assert emf['thd'] == pytest.approx([0.25 * math.sin(3 * step) / math.sin(step)] * 3)


def test_duplicate_endpoint_requires_explicit_drop():
    """Guards: deliberate removal of a duplicate electrical endpoint.

    Assumes: the closing sample repeats the first linkage row.
    """
    first, second = _lane(), _lane()
    second['mechanical_angles_deg'].append(120)
    second['phase_flux_wb'].append(second['phase_flux_wb'][0])
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'invalid_input'
    assert result['checks']['periodic_grid'] is False
    second['drop_duplicate_endpoint'] = True
    assert back_emf_cross_lane_gate(first, second, LIMITS)['lanes'][1]['dropped_duplicate_endpoint']
    assert back_emf_cross_lane_gate(first, second, LIMITS)['status'] == 'pass'


def test_pole_count_ambiguity_and_double_scaling():
    """Guards: declared pole pairs and unique stack conversion.

    Assumes: a factor-of-two pole-pair discrepancy is ambiguous.
    """
    first, second = _lane(), _lane()
    second['pole_pairs'] = 6
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'invalid_input'
    assert 'pole_count_vs_pole_pairs_ambiguity' in result['findings']
    second['pole_pairs'] = 3
    for lane in (first, second):
        lane.update(stack_basis='total', already_scaled=True)
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'invalid_input'
    assert result['checks']['stack_normalization'] is False


def test_backward_difference_voltage_applies_half_step():
    """Guards: native backward-difference voltage is compared at its midpoint.

    Assumes: native terminal voltage uses the positive linkage derivative.
    """
    first, second = _lane(), _lane()
    flux = np.asarray(second['phase_flux_wb'])
    dt = 60 / (900 * 3 * len(flux))
    second.update(emf_method='backward_difference_native', emf_polarity='terminal_drop',
                  phase_emf_v=((flux - np.roll(flux, 1, axis=0)) / dt).tolist())
    limits = {**LIMITS, 'amplitude_relative': 0.001}
    result = back_emf_cross_lane_gate(first, second, limits)
    lane = result['lanes'][1]
    assert result['status'] == 'pass'
    assert not lane['input_emf_phase_before_pass']
    assert result['checks']['input_emf_phase_2']
    assert lane['emf_sample_shift_s'] == pytest.approx(-dt / 2)
    assert lane['emf_sample_shift_deg_e'] == pytest.approx(-3.75)
    assert not lane['input_emf_amplitude_before_pass']
    assert result['checks']['input_emf_amplitude_2']
    assert lane['input_emf_amplitude_correction_factor'] == pytest.approx(math.cos(math.pi / len(flux)))
    second['emf_method'] = 'reported_emf'
    assert back_emf_cross_lane_gate(first, second, LIMITS)['status'] == 'mismatch'


@pytest.mark.parametrize(('swap', 'sign', 'finding'), [
    (True, 1, 'phase_order_abc_vs_acb'),
    (False, -1, 'phase_linkage_sign_inversion'),
])
def test_phase_order_and_sign_findings(swap, sign, finding):
    """Guards: swapped phases and inverted signs are visible mismatches.

    Assumes: the first lane establishes the expected phase order and sign.
    """
    first, second = _lane(), _lane()
    flux = np.asarray(second['phase_flux_wb'])
    second['phase_flux_wb'] = (sign * flux[:, [0, 2, 1] if swap else [0, 1, 2]]).tolist()
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'mismatch'
    assert finding in result['findings']


def test_dimensionality_declares_model_difference_and_grid_offset_is_checked():
    """Guards: dimensionality findings are advisory, grid offsets are not.

    Assumes: the two synthetic waveforms are identical.
    """
    first, second = _lane(), _lane()
    second['model_dimension'] = '3d'
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'pass'
    assert 'model_difference_declared_required' in result['findings']
    second['electrical_offset_deg'] = math.degrees(0.1)
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'invalid_input'
    assert result['checks']['common_electrical_grid'] is False


def test_reported_voltage_factor_and_double_amplitude():
    """Guards: accurate-phase reported voltage can still fail amplitude.

    Assumes: reported voltage is the analytic derivative of sinusoidal linkage.
    """
    first, second = _lane(), _lane()
    theta = np.radians(3 * np.asarray(second['mechanical_angles_deg']))[:, None]
    phase = theta + np.radians([0, -120, 120])
    voltage = 0.08 * math.tau * 45 * np.sin(phase)
    second.update(emf_method='reported_emf', phase_emf_v=voltage.tolist())
    limits = {**LIMITS, 'amplitude_relative': 0.001}
    result = back_emf_cross_lane_gate(first, second, limits)
    assert result['status'] == 'pass'
    assert not result['lanes'][1]['input_emf_amplitude_before_pass']
    second['phase_emf_v'] = (2 * voltage).tolist()
    result = back_emf_cross_lane_gate(first, second, limits)
    assert result['status'] == 'mismatch'
    assert result['checks']['input_emf_phase_2']
    assert not result['checks']['input_emf_amplitude_2']


def test_speed_mismatch_and_stack_scaling_diagnosis():
    """Guards: mismatched speeds and suspected stack mistakes are named.

    Assumes: per-metre values are compared at a common constant speed.
    """
    first, second = _lane(), _lane()
    second['speed_rpm'] *= 2
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert result['status'] == 'invalid_input'
    assert not result['checks']['speed_rpm']
    assert 'speed_differs' in result['findings']
    second['speed_rpm'] = first['speed_rpm']
    second['phase_flux_wb'] = (np.asarray(second['phase_flux_wb']) * second['stack_length_m']).tolist()
    assert 'stack_scaling_suspected' in back_emf_cross_lane_gate(first, second, LIMITS)['findings']


def test_label_and_data_inversion_are_both_reported():
    """Guards: declared polarity changes do not hide inverted data.

    Assumes: the lane signals are exact negatives over a common grid.
    """
    first, second = _lane(), _lane()
    second['linkage_polarity'] = 'negative_winding'
    second['phase_flux_wb'] = (-np.asarray(second['phase_flux_wb'])).tolist()
    result = back_emf_cross_lane_gate(first, second, LIMITS)
    assert {'phase_linkage_sign_inversion', 'linkage_polarity_convention_differs'} <= set(result['findings'])
