"""Repository-owned linear SIBC acceptance; no external user-case inputs.

Scope: axisymmetric solid/bored copper cylinders, uniform scalar SIBC,
prescribed circular current source, P1 dense BEM, total loop field.
FEM is the independent axisymmetric A-form reference. BEM is a 3D surface
discretization of the same axisymmetric geometry, not a 3D volume FEM solve.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import ngsolve as ng
import radia.bem_sibc_solver as sibc
import radia.bem_loop_extension as loop
import radia.workpiece_surface as surface
import radia._radia_pybind as native
from axisym_sibc_reference import solve as fem_solve
from validate_sibc_reaction_power import solve as bem_solve


def metrics(fem, coarse, fine):
    reference = (np.interp(fine['z'], fem['z'], fem['H_z_real'])
                 + 1j*np.interp(fine['z'], fem['z'], fem['H_z_imag']))
    actual = np.array(fine['H_z_real']) + 1j*np.array(fine['H_z_imag'])
    return {
        'loss_relative_error': abs(fine['P_surface']/fem['P_surface']-1),
        'reaction_relative_error': abs(fine['P_reaction_complete']/fem['P_reaction']-1),
        'complex_tangential_field_relative_l2': float(np.linalg.norm(actual-reference)/np.linalg.norm(reference)),
        'bem_loss_refinement_relative_change': abs(fine['P_surface']/coarse['P_surface']-1),
        'bem_power_balance_relative_error': fine['complete_power_relative_error'],
        'heat_projection_relative_error': fine['heat_projection_relative_error'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--recheck', action='store_true',
                        help='Recheck saved measurements with tighter gates; do not rerun solvers')
    args = parser.parse_args()
    ng.SetNumThreads(4)
    t0 = time.perf_counter()
    report = dict(scope=__doc__, host=platform.node(), python=sys.version,
                  ngsolve=ng.__version__, threads=4, cases=[])
    limits = dict(loss_relative_error=.02, reaction_relative_error=.02,
                  complex_tangential_field_relative_l2=.02,
                  bem_loss_refinement_relative_change=.015,
                  bem_power_balance_relative_error=.01,
                  heat_projection_relative_error=1e-10,
                  fem_exterior_relative_change=.005,
                  fem_power_balance_relative_error=1e-8)
    report['limits'] = limits
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.recheck:
        report = json.loads(args.output.read_text(encoding='utf-8'))
        for key, limit in limits.items():
            if limit > report['limits'][key]:
                raise ValueError('Recheck cannot relax an existing acceptance limit')
        for module in (sibc, loop, surface, native):
            source = Path(module.__file__)
            if hashlib.sha256(source.read_bytes()).hexdigest() != report['source_sha256'][source.name]:
                raise ValueError('Numerical source changed; rerun the solver acceptance')
        report['limits'] = limits
        if len(report['cases']) != 2 or {case['bore'] for case in report['cases']} != {0., .0075}:
            raise ValueError('Both solid and bored acceptance cases are required')
        for case in report['cases']:
            case['metrics'] = metrics(case['fem'][-1], case['bem'][-2], case['bem'][-1])
            case['metrics']['fem_exterior_relative_change'] = abs(case['fem'][-1]['P_surface']/case['fem'][0]['P_surface']-1)
            case['metrics']['fem_power_balance_relative_error'] = max(row['power_balance_relative_error'] for row in case['fem'])
            case['checks'] = {key: value <= limits[key] for key,value in case['metrics'].items()}
            case['passed'] = all(case['checks'].values())
        report['passed'] = all(case['passed'] for case in report['cases'])
        report['gate_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
        if not report['passed']:
            raise RuntimeError('Tightened IH acceptance failed')
        print('Saved solver measurements pass the tightened gates; numerical source hashes match',flush=True)
        return
    with ng.TaskManager():
        for bore in (0., .0075):
            fem = [fem_solve(bore=bore, outer=radius) for radius in (.15, .3)]
            bem = []
            for maxh in (.003, .002, .0015):
                row = bem_solve(maxh, bore=bore, coil_a=.0005)
                bem.append(row)
                print(json.dumps(dict(bore=bore, maxh=maxh, P=row['P_surface'],
                      balance=row['complete_power_relative_error'])), flush=True)
            measured = metrics(fem[-1], bem[-2], bem[-1])
            measured['fem_exterior_relative_change'] = abs(fem[-1]['P_surface']/fem[0]['P_surface']-1)
            measured['fem_power_balance_relative_error'] = max(r['power_balance_relative_error'] for r in fem)
            checks = {key: float(value) <= limits[key] for key,value in measured.items()}
            case = dict(bore=bore, fem=fem, bem=bem, metrics=measured, checks=checks,
                        passed=all(checks.values()))
            report['cases'].append(case)
            print(json.dumps(dict(bore=bore, metrics=measured, checks=checks)), flush=True)
            args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    paths = [Path(__file__), Path(__file__).with_name('axisym_sibc_reference.py'),
             Path(__file__).with_name('validate_sibc_reaction_power.py'),
             Path(sibc.__file__), Path(loop.__file__), Path(surface.__file__), Path(native.__file__)]
    report['source_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    report['elapsed_s'] = time.perf_counter()-t0
    report['passed'] = all(case['passed'] for case in report['cases'])
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    if not report['passed']:
        raise RuntimeError('IH SIBC acceptance failed; inspect checks in JSON')


if __name__ == '__main__':
    main()
