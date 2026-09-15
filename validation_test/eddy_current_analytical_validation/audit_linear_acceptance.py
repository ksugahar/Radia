"""Audit saved linear measurements and their recovered source identities.

This is NOT a new solver run or a certification of an installed native binary.
The recovered original must first match the recorded raw SHA256. Only CRLF/LF
conversion is then allowed when comparing it with the current checkout.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIMITS = dict(loss_relative_error=.02, reaction_relative_error=.02,
              complex_tangential_field_relative_l2=.02,
              bem_loss_refinement_relative_change=.015,
              bem_power_balance_relative_error=.01,
              heat_projection_relative_error=1e-10,
              fem_exterior_relative_change=.005,
              fem_power_balance_relative_error=1e-8)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def audit(archive):
    evidence = HERE / 'ih_sibc_closure_20260915.json'
    saved = json.loads(evidence.read_text(encoding='utf-8'))
    if {c['bore'] for c in saved['cases']} != {0., .0075} or len(saved['cases']) != 2:
        raise ValueError('Require both solid and bored cases')
    sources = {}
    for name in ('axisym_sibc_reference.py', 'validate_sibc_reaction_power.py',
                 'bem_sibc_solver.py', 'bem_loop_extension.py', 'workpiece_surface.py'):
        is_reference = name in ('axisym_sibc_reference.py', 'validate_sibc_reaction_power.py')
        original = archive / ('validation' if is_reference else 'src/radia') / name
        current = (HERE if is_reference else ROOT / 'src/radia') / name
        old, new = original.read_bytes(), current.read_bytes()
        if sha(old) != saved['source_sha256'][name]:
            raise ValueError(f'Recovered original hash mismatch: {name}')
        if old.replace(b'\r\n', b'\n') != new.replace(b'\r\n', b'\n'):
            raise ValueError(f'Numerical source changed: {name}; rerun solvers')
        sources[name] = dict(recorded_raw_sha256=sha(old), current_raw_sha256=sha(new),
                             lf_sha256=sha(new.replace(b'\r\n', b'\n')))
    cases = []
    for case in saved['cases']:
        fem, coarse, fine = case['fem'][-1], case['bem'][-2], case['bem'][-1]
        expected = (np.interp(fine['z'], fem['z'], fem['H_z_real'])
                    + 1j*np.interp(fine['z'], fem['z'], fem['H_z_imag']))
        actual = np.array(fine['H_z_real']) + 1j*np.array(fine['H_z_imag'])
        metrics = dict(
            loss_relative_error=abs(fine['P_surface']/fem['P_surface']-1),
            reaction_relative_error=abs(fine['P_reaction_complete']/fem['P_reaction']-1),
            complex_tangential_field_relative_l2=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected)),
            bem_loss_refinement_relative_change=abs(fine['P_surface']/coarse['P_surface']-1),
            bem_power_balance_relative_error=fine['complete_power_relative_error'],
            heat_projection_relative_error=fine['heat_projection_relative_error'],
            fem_exterior_relative_change=abs(fem['P_surface']/case['fem'][0]['P_surface']-1),
            fem_power_balance_relative_error=max(r['power_balance_relative_error'] for r in case['fem']))
        checks = {k: math.isfinite(v) and 0 <= v <= LIMITS[k] for k, v in metrics.items()}
        cases.append(dict(bore=case['bore'], metrics=metrics, checks=checks, passed=all(checks.values())))
    return dict(scope=saved['scope'], evidence_sha256=sha(evidence.read_bytes()),
                audit_source_sha256=sha(Path(__file__).read_bytes()),
                audit_kind='saved-measurement and recovered-source audit; no solver rerun',
                source_identity=sources, limits=LIMITS, cases=cases,
                recorded_native_sha256=saved['source_sha256']['_radia_pybind.pyd'],
                installed_native_certified=False,
                excluded=['nonlinear steel', 'temperature-dependent BH',
                          'general linear-material/frequency/geometry envelope',
                          'SIBC approximation error versus volume-resolved conductor FEM',
                          'installed-package or real-file Simulink delivery acceptance'],
                passed=all(c['passed'] for c in cases))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.archive)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))
    if not result['passed']:
        raise SystemExit('Linear acceptance audit failed')
