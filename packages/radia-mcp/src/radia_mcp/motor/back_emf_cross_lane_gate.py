"""Compare two periodic three-phase linkage and back-EMF exports."""

import math


_METHODS = {'central_difference_of_linkage', 'backward_difference_native', 'reported_emf'}
_TOLERANCES = ('grid_deg_e', 'amplitude_relative', 'phase_deg_e', 'rms_relative', 'thd_absolute')


def _wrap(degrees):
    return (degrees + 180.0) % 360.0 - 180.0


def _positive(value, name):
    number = float(value)
    if isinstance(value, bool) or not math.isfinite(number) or number <= 0:
        raise ValueError(f'{name} must be finite and positive')
    return number


def _array(value, shape, name):
    array = [tuple(float(component) for component in row) for row in value]
    if len(array) != shape[0] or any(len(row) != shape[1] or not all(map(math.isfinite, row)) for row in array):
        raise ValueError(f'{name} must be finite with shape {shape}')
    return array


def _flag(lane, name):
    value = lane.get(name, False)
    if not isinstance(value, bool):
        raise ValueError(f'{name} must be a boolean')
    return value


def _prepare(lane, poles, tolerance):
    omega = _positive(lane['speed_rpm'], 'speed_rpm') * (math.tau / 60)
    has_angles = 'mechanical_angles_deg' in lane
    if has_angles == ('times_s' in lane):
        raise ValueError('provide mechanical_angles_deg or times_s, exactly once')
    coordinate = [float(value) for value in lane['mechanical_angles_deg' if has_angles else 'times_s']]
    if len(coordinate) < 8 or not all(map(math.isfinite, coordinate)):
        raise ValueError('at least eight finite angle or time samples are required')
    offset = float(lane['electrical_offset_deg'])
    if not math.isfinite(offset):
        raise ValueError('electrical_offset_deg must be finite')
    angle = coordinate if has_angles else [math.degrees(omega * value) for value in coordinate]
    theta = [poles * value + offset for value in angle]
    n = len(theta)
    flux = _array(lane['phase_flux_wb'], (n, 3), 'phase_flux_wb')
    method = lane['emf_method']
    if method not in _METHODS:
        raise ValueError('unknown emf_method')
    voltage = None
    if method != 'central_difference_of_linkage' or 'phase_emf_v' in lane:
        voltage = _array(lane['phase_emf_v'], (n, 3), 'phase_emf_v')
    if not all(map(math.isfinite, theta)) or any(b <= a for a, b in zip(theta, theta[1:])):
        raise ValueError('electrical angles must increase strictly')
    dropped = abs(theta[-1] - theta[0] - 360) <= tolerance
    if dropped:
        if not _flag(lane, 'drop_duplicate_endpoint'):
            raise ValueError('duplicate periodic endpoint')
        theta, flux = theta[:-1], flux[:-1]
        if voltage is not None:
            voltage = voltage[:-1]
    n = len(theta)
    step_deg = 360 / n
    if tolerance >= step_deg / 2:
        raise ValueError('grid tolerance must be smaller than half an electrical sample interval')
    if n < 8 or any(abs(value - theta[0] - i * step_deg) > tolerance for i, value in enumerate(theta)):
        raise ValueError('grid must span one uniformly sampled electrical period')
    basis = lane['stack_basis']
    if basis not in ('per_metre', 'total'):
        raise ValueError('stack_basis must be per_metre or total')
    length = _positive(lane['stack_length_m'], 'stack_length_m')
    scaled = _flag(lane, 'already_scaled')
    if scaled and basis != 'total':
        raise ValueError('already_scaled requires total stack basis')
    factor = 1 / length if basis == 'total' else 1.0
    polarity = lane['linkage_polarity']
    if polarity not in ('positive_winding', 'negative_winding'):
        raise ValueError('linkage_polarity must name the signed winding convention')
    if lane['model_dimension'] not in ('2d', '3d'):
        raise ValueError('model_dimension must be 2d or 3d')
    voltage_sign = lane.get('emf_polarity', 'generator')
    if voltage_sign not in ('generator', 'terminal_drop'):
        raise ValueError('emf_polarity must be generator or terminal_drop')
    flux = [tuple(value * factor for value in row) for row in flux]
    step_s = math.radians(step_deg) / (poles * omega)
    emf = [tuple((previous - following) / (2 * step_s)
                 for previous, following in zip(flux[index - 1], flux[(index + 1) % n]))
           for index in range(n)]
    shift_s = -step_s / 2 if method == 'backward_difference_native' else 0.0
    if voltage is not None:
        voltage_factor = -factor if voltage_sign == 'terminal_drop' else factor
        voltage = [tuple(value * voltage_factor for value in row) for row in voltage]
    if not all(math.isfinite(value) for row in [*flux, *emf, *(voltage or [])] for value in row):
        raise ValueError('stack and speed conversion produced nonfinite samples')
    return theta, flux, emf, voltage, shift_s, step_s, dropped, factor


def _metrics(values, start_deg):
    n = len(values)
    cosine = [math.cos(math.tau * i / n) for i in range(n)]
    sine = [math.sin(math.tau * i / n) for i in range(n)]
    result = {'amplitude': [], 'phase_deg_e': [], 'rms': [], 'thd': []}
    for phase in zip(*values):
        scale = max(map(abs, phase))
        if not scale:
            raise ValueError('each phase must have a resolved nonzero fundamental')
        unit = [value / scale for value in phase]
        real = 2 * math.fsum(value * basis for value, basis in zip(unit, cosine)) / n
        imag = -2 * math.fsum(value * basis for value, basis in zip(unit, sine)) / n
        amplitude = math.hypot(real, imag)
        if amplitude <= math.ulp(1.0):
            raise ValueError('each phase must have a resolved nonzero fundamental')
        mean = math.fsum(unit) / n
        residual = math.fsum((value - mean - real * c + imag * s)**2
                             for value, c, s in zip(unit, cosine, sine)) / n
        peak = amplitude * scale
        if not math.isfinite(peak):
            raise ValueError('fundamental amplitude exceeds the supported numeric range')
        result['amplitude'].append(peak)
        result['phase_deg_e'].append(_wrap(math.degrees(math.atan2(imag, real)) - start_deg))
        result['rms'].append(math.sqrt(math.fsum(value**2 for value in unit) / n) * scale)
        result['thd'].append(math.sqrt(2 * residual) / amplitude)
    return result


def _compare(first, second, limits):
    differences = {}
    tolerances = {'amplitude': limits['amplitude_relative'], 'phase_deg_e': limits['phase_deg_e'],
                  'rms': limits['rms_relative'], 'thd': limits['thd_absolute']}
    for key in tolerances:
        pairs = tuple(zip(first[key], second[key]))
        if key == 'phase_deg_e':
            delta = [_wrap(b - a) for a, b in pairs]
        elif key in ('amplitude', 'rms'):
            delta = [(b - a) / max(a, b) for a, b in pairs]
        else:
            delta = [b - a for a, b in pairs]
        differences[key] = {'difference': delta, 'pass': all(abs(value) <= tolerances[key] for value in delta)}
    return differences


def _diagnose_phase(first, second, tolerance):
    def error(order, offset):
        return max(abs(_wrap(second[j] - a - offset)) for a, j in zip(first, order))

    direct_error = error((0, 1, 2), 0)
    inverted_error = error((0, 1, 2), 180)
    swapped_error = error((0, 2, 1), 0)
    swapped_inverted_error = error((0, 2, 1), 180)
    acb = min(swapped_error, swapped_inverted_error) < min(direct_error, inverted_error)
    if acb and min(swapped_error, swapped_inverted_error) <= tolerance:
        return True, bool(swapped_inverted_error <= tolerance)
    return False, bool(inverted_error < direct_error and inverted_error <= tolerance)


def back_emf_cross_lane_gate(first_lane: dict, second_lane: dict, tolerances: dict) -> dict:
    """Gate per-metre abc linkage and generator EMF on a shared electrical grid.

    Supply electrical_offset_deg, speed_rpm, stack_length_m, stack_basis,
    linkage_polarity, model_dimension and emf_method for each lane. For times,
    mechanical zero is at time zero and the electrical offset sets the reference.
    Native/reported methods require phase_emf_v on the input sample grid.
    The main comparison uses periodic central differences for both lanes.
    """
    checks = {}
    findings = []
    result = {'schema': 'radia-motor-back-emf-cross-lane/v1', 'checks': checks, 'findings': findings}
    stage = 'tolerances'
    try:
        limits = {name: _positive(tolerances[name], name) for name in _TOLERANCES}
        checks[stage] = True
        stage = 'pole_pairs'
        poles = [lane['pole_pairs'] for lane in (first_lane, second_lane)]
        if any(type(p) is not int or p < 1 for p in poles):
            raise ValueError('pole_pairs must be positive integers')
        if poles[0] != poles[1]:
            if max(poles) == 2 * min(poles):
                findings.append('pole_count_vs_pole_pairs_ambiguity')
            raise ValueError('pole_pairs differ')
        checks[stage] = True
        stage = 'stack_normalization'
        scaled = [_flag(lane, 'already_scaled') for lane in (first_lane, second_lane)]
        if all(scaled):
            raise ValueError('both lanes claim already-scaled linkage')
        for lane in (first_lane, second_lane):
            _positive(lane['stack_length_m'], 'stack_length_m')
        checks[stage] = True
        stage = 'periodic_grid'
        prepared = [_prepare(lane, poles[0], limits['grid_deg_e']) for lane in (first_lane, second_lane)]
        checks[stage] = True
        stage = 'common_electrical_grid'
        a, b = (item[0] for item in prepared)
        if len(a) != len(b) or any(abs(_wrap(y - x)) > limits['grid_deg_e'] for x, y in zip(a, b)):
            raise ValueError('electrical grids differ after angle mapping')
        checks[stage] = True
        stage = 'fundamentals_resolved'
        flux, emf = ([_metrics(item[column], item[0][0]) for item in prepared] for column in (1, 2))
        checks[stage] = True
        detail = {name: _compare(*values, limits) for name, values in (('linkage', flux), ('emf', emf))}
        result['comparisons'] = detail
        result['metrics_per_metre'] = {'linkage': flux, 'emf': emf}
        acb, inverted = _diagnose_phase(flux[0]['phase_deg_e'], flux[1]['phase_deg_e'], limits['phase_deg_e'])
        polarity_matches = first_lane['linkage_polarity'] == second_lane['linkage_polarity']
        if acb:
            findings.append('phase_order_abc_vs_acb')
        if inverted or not polarity_matches:
            findings.append('phase_linkage_sign_inversion' if inverted else 'linkage_polarity_convention_differs')
        checks['phase_order'] = not acb
        checks['linkage_polarity'] = polarity_matches and not inverted
        for name, comparison in detail.items():
            checks[f'{name}_metrics'] = all(item['pass'] for item in comparison.values())
        stage = 'input_emf_phase'
        result['lanes'] = []
        for index, (angles, _, _, voltage, shift, dt, dropped, factor) in enumerate(prepared):
            offset = shift / dt * 360 / len(angles)
            lane_report = {'emf_sample_shift_s': shift, 'emf_sample_shift_deg_e': offset,
                           'dropped_duplicate_endpoint': dropped, 'per_metre_factor': factor}
            if voltage is not None:
                native = _metrics(voltage, angles[0])['phase_deg_e']
                before = [_wrap(x - y) for x, y in zip(native, emf[index]['phase_deg_e'])]
                after = [_wrap(value - offset) for value in before]
                passed = all(abs(value) <= limits['phase_deg_e'] for value in after)
                lane_report.update(input_emf_phase_before_deg_e=before,
                                   input_emf_phase_after_deg_e=after,
                                   input_emf_phase_before_pass=all(abs(value) <= limits['phase_deg_e'] for value in before))
                checks[f'{stage}_{index + 1}'] = passed
            result['lanes'].append(lane_report)
        if first_lane['model_dimension'] != second_lane['model_dimension']:
            findings.append('model_difference_declared_required')
            result['model_difference_examples'] = ['end leakage', 'axial fringing']
        result['status'] = 'pass' if all(checks.values()) else 'mismatch'
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        checks[stage] = False
        result.update(status='invalid_input', error=str(exc))
    return result
