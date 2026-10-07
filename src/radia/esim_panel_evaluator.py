"""Per-run ESIM acceleration; final certification always uses direct cells.

Linear interpolation in log H preserves endpoint passivity. Midpoint probes
control empirical interpolation error, not a rigorous global error bound.
Consequently interpolation is never used for final constitutive certification.
"""
from __future__ import annotations

import hashlib
import json
import math
import time

import numpy as np

from radia.esim_cell_problem import require_esim_converged

def _solve_cell(solver, field):
    started = time.perf_counter()
    result = require_esim_converged(
        solver.solve(float(field)),
        f'ESIM direct cell at H={field:g}')
    value = complex(result['Z'])
    if not np.isfinite(value) or value.real < 0:
        raise RuntimeError('ESIM cell returned nonfinite or nonpassive impedance')
    return value, time.perf_counter() - started


def finite_cell_configuration(solver):
    """Fingerprint inputs/state of the fixed finite-slab/cylinder cell law."""
    names = (
        'half_thickness', 'sigma', 'frequency', 'omega', 'rho', 'n_nodes',
        'geometry', 'linear_mu_r', 'use_complex_mu', 'mu_initial',
        'delta', 'xi', 'mesh_points', 'n_elements')
    missing = [name for name in (*names, 'bh_interp', 'mu_interp') if not hasattr(solver, name)]
    if missing:
        raise ValueError('Unsupported ESIM finite-cell configuration; missing attributes: ' + ', '.join(missing))
    state = {name: getattr(solver, name) for name in names}
    state['cell_type'] = f'{type(solver).__module__}.{type(solver).__qualname__}'
    if solver.bh_interp is not None:
        state['bh'] = {name: getattr(solver.bh_interp, name) for name in
                       ('H_data', 'B_data', 'mu_initial', '_saturation_magnetization')}
    if solver.mu_interp is not None:
        state['complex_mu'] = {name: getattr(solver.mu_interp, name) for name in (
            'frequency', 'mode', 'mu_initial', 'mu_prime_r', 'mu_double_prime_r',
            'H_data', 'mu_prime_r_data',
            'mu_double_prime_r_data', 'f_data', 'mu_prime_r_f_data',
            'mu_double_prime_r_f_data') if hasattr(solver.mu_interp, name)}
    return state


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, complex):
        return {'real': value.real, 'imag': value.imag}
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f'Unsupported cell identity value {type(value).__name__}')


class PanelESIMEvaluator:
    """Owned fixed-law evaluator, with no persistent/cross-run table reuse.

    ``configuration`` may supply a complete identity for an explicit test law.
    Production uses the same finite-cell fingerprint in every mode.
    A mutation invalidates the evaluator loudly instead of reusing old data.
    Certification deliberately bypasses even the exact-value cache.
    """
    def __init__(self, solver, *, mode='direct', interpolation_tol=1e-4,
                 configuration=None, max_table_cells=4096):
        if mode not in ('direct', 'table'):
            raise ValueError('ESIM evaluator must be direct or table')
        if not np.isfinite(interpolation_tol) or not 0 < interpolation_tol < 1:
            raise ValueError('ESIM interpolation tolerance must be finite in (0,1)')
        if isinstance(max_table_cells, bool) or int(max_table_cells) != max_table_cells or max_table_cells < 3:
            raise ValueError('ESIM table cell budget must be an integer >=3')
        self.solver = solver
        self.mode = mode
        self.tolerance = float(interpolation_tol)
        self.max_table_cells = int(max_table_cells)
        self.configuration = configuration or (lambda: finite_cell_configuration(solver))
        self.identity = self._identity()
        self.cache = {}
        self.knots = []
        self.counts = dict(iteration_cell_calls=0, table_cell_calls=0,
                           certification_cell_calls=0, cache_hits=0,
                           interpolation_hits=0, interval_probes=0)
        self.cell_seconds = 0.0
        self.evaluation_seconds = 0.0
        self.certification_seconds = 0.0
        self.max_accepted_probe_error = 0.0

    def _identity(self):
        payload = dict(cell=self.configuration(),
                       solve_controls={'tol': 1e-6, 'max_iter': 50, 'relaxation': .5},
                       field_floor_A_per_m=1e-3)
        raw = json.dumps(payload, sort_keys=True, default=_json_value,
                         allow_nan=False, separators=(',', ':')).encode()
        return hashlib.sha256(raw).hexdigest()

    def _fields(self, fields):
        if self._identity() != self.identity:
            raise RuntimeError('ESIM cell configuration changed; create a new evaluator')
        values = np.asarray(fields, dtype=float)
        if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError('ESIM panel fields must be a finite nonnegative 1D array')
        return np.maximum(values, 1e-3)

    def _direct(self, field, category):
        self.counts[category] += 1
        value, elapsed = _solve_cell(self.solver, field)
        self.cell_seconds += elapsed
        return value

    def _direct_many(self, fields, category):
        return [self._direct(field, category) for field in fields]

    def _cached_pair(self, fields):
        missing = list(dict.fromkeys(float(h) for h in fields if float(h) not in self.cache))
        if self.counts['table_cell_calls'] + len(missing) > self.max_table_cells:
            raise RuntimeError('ESIM adaptive table cell budget exhausted')
        if missing:
            values = self._direct_many(missing, 'table_cell_calls')
            self.cache.update(zip(missing, values))
        self.counts['cache_hits'] += len(fields) - len(missing)
        return [self.cache[float(h)] for h in fields]

    def _cached(self, field, category):
        field = float(field)
        if field in self.cache:
            self.counts['cache_hits'] += 1
            return self.cache[field]
        if category == 'table_cell_calls' and self.counts[category] >= self.max_table_cells:
            raise RuntimeError('ESIM adaptive table cell budget exhausted')
        value = self._direct(field, category)
        self.cache[field] = value
        return value

    def _refine(self, low, high, depth=0):
        zl, zh = self._cached_pair((low, high))
        middle = math.exp((math.log(low) + math.log(high)) / 2)
        if middle <= low or middle >= high or depth >= 40:
            raise RuntimeError('ESIM table refinement exhausted floating-point/depth resolution')
        exact = self._cached(middle, 'table_cell_calls')
        self.counts['interval_probes'] += 1
        difference = exact - (zl + zh) / 2
        error = max(abs(difference.real), abs(difference.imag)) / max(abs(exact), 1e-30)
        if error <= self.tolerance:
            self.max_accepted_probe_error = max(self.max_accepted_probe_error, error)
            return [low, high]
        left = self._refine(low, middle, depth + 1)
        right = self._refine(middle, high, depth + 1)
        return left[:-1] + right

    def _extend(self, low, high):
        if self.knots and low >= self.knots[0] and high <= self.knots[-1]:
            return
        if self.knots:
            low, high = min(low, self.knots[0]), max(high, self.knots[-1])
        if low == high:
            self._cached(low, 'table_cell_calls')
            self.knots = [low]
            return
        # Seed at no more than one decade spacing; probe/adapt each interval.
        count = max(1, math.ceil(math.log10(high) - math.log10(low)))
        if count > self.max_table_cells:
            raise RuntimeError('ESIM field range exceeds table cell budget')
        seeds = np.geomspace(low, high, count + 1)
        seeds[0], seeds[-1] = low, high
        knots = []
        for a, b in zip(seeds[:-1], seeds[1:]):
            knots.extend(self._refine(float(a), float(b))[:-1])
        self.knots = knots + [high]

    def evaluate(self, fields):
        values = self._fields(fields)
        started = time.perf_counter()
        try:
            if not len(values):
                return np.empty(0, dtype=complex)
            if self.mode == 'direct':
                return np.array(self._direct_many(values, 'iteration_cell_calls'))
            self._extend(float(values.min()), float(values.max()))
            logs = np.log(self.knots)
            impedances = np.array([self.cache[h] for h in self.knots])
            output = []
            for field in values:
                if float(field) in self.cache:
                    self.counts['cache_hits'] += 1
                    output.append(self.cache[float(field)])
                else:
                    self.counts['interpolation_hits'] += 1
                    output.append(complex(np.interp(math.log(field), logs, impedances.real),
                                          np.interp(math.log(field), logs, impedances.imag)))
            return np.array(output)
        finally:
            self.evaluation_seconds += time.perf_counter() - started

    def certify(self, fields):
        values = self._fields(fields)
        started = time.perf_counter()
        try:
            return np.array(self._direct_many(values, 'certification_cell_calls'))
        finally:
            self.certification_seconds += time.perf_counter() - started

    def diagnostics(self):
        return dict(mode=self.mode, configuration_sha256=self.identity,
                    interpolation_tolerance=self.tolerance,
                    interpolation_error_control='empirical-midpoint-probes',
                    final_certification='direct-all-panels', field_floor_A_per_m=1e-3,
                    table_knots=len(self.knots), max_table_cells=self.max_table_cells,
                    max_accepted_probe_error=self.max_accepted_probe_error,
                    **self.counts, direct_cell_calls=sum(self.counts[k] for k in
                    ('iteration_cell_calls', 'table_cell_calls', 'certification_cell_calls')),
                    cell_seconds=self.cell_seconds, evaluation_seconds=self.evaluation_seconds,
                    certification_seconds=self.certification_seconds, workers=1,
                    cell_time_kind='sum-of-direct-cell-wall-times')
