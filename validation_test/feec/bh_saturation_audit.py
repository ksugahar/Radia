"""Audit the actual inverse-BH helpers without importing NGSolve or native Radia.

This is curve diagnostics, not a nonlinear-solver acceptance verdict. It does
not alter input tables, validation policy, or the production constitutive law.
"""
import argparse
import ast
import hashlib
import importlib.metadata
import json
from math import pi
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'src/radia/vim/_nonlinear.py'


def load_curve_helpers(source=SOURCE):
    """Execute selected production scalar helpers, not a second implementation."""
    tree = ast.parse(Path(source).read_text(encoding='utf-8'))
    functions = {'_validate_bh_table', '_bh_table_funcs', '_bh_inverse_funcs'}
    constants = {'_MU0', '_CHI_DIFF_FLOOR'}
    nodes = [node for node in tree.body if
             (isinstance(node, ast.FunctionDef) and node.name in functions)
             or (isinstance(node, ast.Assign) and any(
                 isinstance(target, ast.Name) and target.id in constants for target in node.targets))]
    namespace = {'np': np, 'pi': pi}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), namespace)
    for name in functions | constants:
        if name not in namespace:
            raise RuntimeError('Missing production helper: ' + name)
    return namespace


def regression_table():
    path = ROOT / 'validation_test/feec/test_hdiv_vim_energy_newton.py'
    tree = ast.parse(path.read_text(encoding='utf-8'))
    nodes = [node for node in tree.body if isinstance(node, ast.Assign) and any(
        isinstance(target, ast.Name) and target.id in {'_H', '_B'} for target in node.targets)]
    namespace = {'np': np}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)
    return np.column_stack((namespace['_H'], namespace['_B']))


def saturation_report(table):
    helpers = load_curve_helpers()
    helpers['_validate_bh_table'](table, context='saturation audit')
    h, b = np.asarray(table, dtype=float).T
    mu0 = helpers['_MU0']
    _, curve, derivative, hmax, endpoint = helpers['_bh_table_funcs'](h, b)
    # PCHIP is piecewise cubic. Its derivative extrema and M stationary points
    # can be evaluated within every polynomial interval, without a sampling grid.
    stationary = derivative.solve(mu0, extrapolate=False)
    stationary = stationary[np.isfinite(stationary) & (stationary >= h[0]) & (stationary <= hmax)]
    candidates = np.unique(np.concatenate((h, stationary)))
    magnetization = curve(candidates) / mu0 - candidates
    peak_index = int(np.argmax(magnetization))
    extrema = derivative.derivative().roots(extrapolate=False)
    extrema = extrema[np.isfinite(extrema) & (extrema >= h[0]) & (extrema <= hmax)]
    probes = np.unique(np.concatenate((h, extrema)))
    chi = derivative(probes) / mu0 - 1
    tolerances = 128 * np.finfo(float).eps * np.maximum(1., np.abs(chi + 1))
    node_m = b / mu0 - h
    drops = np.flatnonzero(np.diff(node_m) < -128 * np.finfo(float).eps * max(1., float(np.max(abs(node_m)))))
    result = dict(schema='radia.bh-saturation-audit.v1', rows=len(h),
                  source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  endpoint_M_Am=endpoint, peak_M_Am=float(magnetization[peak_index]),
                  peak_H_Am=float(candidates[peak_index]),
                  minimum_interpolated_chi=float(np.min(chi)),
                  chi_roundoff_tolerance=float(tolerances[np.argmin(chi)]),
                  interpolated_M_nondecreasing=bool(np.all(chi >= -tolerances)),
                  first_decreasing_node_interval=(h[drops[0]:drops[0]+2].tolist() if len(drops) else None),
                  inverse_error=None)
    try:
        fields, _, maximum = helpers['_bh_inverse_funcs'](h, b)
        cap = maximum * (1 - 1e-9)
        _, tangents = fields(np.array([np.nextafter(cap, -np.inf), np.nextafter(cap, np.inf)]))
        result.update(barrier_M_Am=cap, barrier_to_peak_ratio=cap/result['peak_M_Am'],
                      tangent_below=float(tangents[0]), tangent_above=float(tangents[1]),
                      tangent_jump_ratio=float(tangents[1]/tangents[0]))
    except ValueError as exc:
        result['inverse_error'] = str(exc)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sample = ROOT / 'src/radia/panels/samples/em_sample_bh.txt'
    fixture = ROOT / 'validation_test/feec/test_hdiv_vim_energy_newton.py'
    report = {'regression_table': saturation_report(regression_table()),
              'shipped_sample': saturation_report(np.loadtxt(sample)),
              'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in (fixture, sample)},
              'scipy_version': importlib.metadata.version('scipy'),
              'scope': 'Constitutive curve properties only; not a solver convergence verdict.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
