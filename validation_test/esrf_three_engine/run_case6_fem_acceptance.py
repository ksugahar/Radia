"""Re-solve both FEM routes with the installed wheel of an accepted HDiv run."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time
import traceback

import ngsolve as ng
import numpy as np
import radia as rad

from case6_acceptance_contract import require_accepted_hdiv, verify_relocated_identity
from esrf_coil_yoke import (average_observation_field, build_radia_coil_source,
                           core_selector, get_case, observation_volume_quadrature)
from run_coil_yoke_three_engine import (_comparison_gate, _is_converged_result,
                                      _load_shared_engines, _pairwise, _require_mesh_contract,
                                      _sha256, _validated_field)
from radia.esrf_examples import get_esrf_bh_table
from radia.kelvin_identify_ngsolve import detect_kelvin_offset, has_kelvin_identification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hdiv-result', type=Path, required=True)
    parser.add_argument('--iron-mesh', type=Path, required=True,
                        help='Restaged iron input matching the accepted HDiv SHA-256')
    parser.add_argument('--wheel', type=Path, required=True,
                        help='Original candidate wheel matching the accepted HDiv SHA-256')
    parser.add_argument('--fem-mesh', type=Path, required=True)
    parser.add_argument('--fem-report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--engine', choices=('reduced_a', 'mixed_total_reduced_omega', 'both'), default='both')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    hdiv = json.loads(args.hdiv_result.read_text(encoding='utf-8'))
    require_accepted_hdiv(hdiv)
    direct = json.loads(importlib.metadata.distribution('radia').read_text('direct_url.json') or '{}')
    if direct.get('dir_info', {}).get('editable'):
        raise RuntimeError('Use the non-editable repair wheel')
    package = Path(rad.__file__).resolve().parent
    runtime_identity = verify_relocated_identity(
        hdiv['implementation'], package, args.iron_mesh, args.wheel)
    fem_sha = _sha256(args.fem_mesh)
    if fem_sha != 'dfc12b84f80db1fa17fb5012b6f87c072e8aa1fb3c085b61269f047a879a54ac':
        raise RuntimeError('FEM mesh is not the audited ESRF6 mesh')
    mesh_contract = _require_mesh_contract(args.fem_report)
    ng.SetNumThreads(args.threads)
    with ng.TaskManager():
        mesh = ng.Mesh(str(args.fem_mesh))
    if not has_kelvin_identification(mesh):
        raise RuntimeError('Missing Kelvin identification')
    centre = tuple(map(float, detect_kelvin_offset(mesh)))
    case = get_case(6)
    points, field_points = observation_volume_quadrature(6, half_width_m=hdiv['observation_half_width_m'])
    if not np.array_equal(points, hdiv['observation_points_m']):
        raise RuntimeError('Observation points changed')
    rad.UtiDelAll()
    coil, source = build_radia_coil_source(6)
    bh = np.asarray(get_esrf_bh_table(6), dtype=float)
    if source != hdiv['source'] or not np.array_equal(bh, hdiv['bh_table']):
        raise RuntimeError('Coil source or BH table differs from HDiv')
    engines = _load_shared_engines()
    # Importing a checkout adapter must not redirect the selected installed package.
    if Path(rad.__file__).resolve().parent != package or any(Path(p).resolve() != package for p in rad.__path__):
        raise RuntimeError('Adapter redirected the installed package')
    common = dict(nonlinear=True, order=2, nonlinear_tolerance=2e-5,
                  nonlinear_maximum_iterations=200, nonlinear_verbose=True,
                  kelvin_center=centre, kelvin_radius=case.kelvin_radius_m,
                  points=field_points, observation_points=field_points)
    settings = {
        'reduced_a': dict(linear_solver='bddc', relax=0.1, anderson_depth=0),
        'mixed_total_reduced_omega': dict(source_trace_tolerance=0.05, relaxation=0.3,
                                        anderson_depth=2, source_projection_order=2,
                                        bonus_intorder=4, exact_exterior_source=False),
    }
    fields = {'hdiv_mmm': _validated_field(hdiv['B_T'], len(points)).tolist()}
    report = dict(schema='radia.validation.esrf6-repaired-three-engine.v1',
                  machine=platform.node(), started_utc=datetime.now(timezone.utc).isoformat(),
                  argv=sys.argv, implementation=hdiv['implementation'],
                  runtime_identity=runtime_identity,
                  hdiv_result_sha256=_sha256(args.hdiv_result), fem_mesh_sha256=fem_sha,
                  fem_report_sha256=_sha256(args.fem_report), fem_contract=mesh_contract,
                  runner_sha256=_sha256(Path(__file__)),
                  contract_sha256=_sha256(Path(require_accepted_hdiv.__code__.co_filename)),
                  adapter_sha256=_sha256(Path(engines.__file__)),
                  source=source, bh_table=bh.tolist(), observation_points_m=points.tolist(),
                  observation_half_width_m=hdiv['observation_half_width_m'],
                  threads=args.threads, fem_order=2, settings=settings,
                  nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=200,
                  fields_T=fields, engines={}, completed=False, passed=False,
                  acceptance_scope='ESRF6 BDM1, nominal excitation, core gap vector-field agreement <=3%',
                  limitations=['Not an absolute-error certificate or deep-saturation validation.',
                               'The historical wheel uses default Hodge quadrature and mixed bonus=4; '
                               'source quadrature convergence is not certified by algebraic residuals.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')

    save()
    start = time.perf_counter()
    try:
        for name, solve in [('reduced_a', engines.solve_reduced_a),
                            ('mixed_total_reduced_omega', engines.solve_omega)]:
            if args.engine not in ('both', name):
                continue
            print('START '+name, flush=True)
            samples, diagnostic = solve(mesh, coil, bh, **common, **settings[name])
            fields[name] = average_observation_field(_validated_field(samples, len(field_points)), len(points)).tolist()
            report['engines'][name] = diagnostic
            report['elapsed_s'] = time.perf_counter()-start
            save()
            converged = _is_converged_result(diagnostic)
            print('DONE '+name+' '+json.dumps(dict(converged=converged, runtime_s=diagnostic['runtime_s'])), flush=True)
            if not converged:
                raise RuntimeError(name+' nonlinear solve did not converge')
        if args.engine == 'both':
            values = {k: np.asarray(v) for k, v in fields.items()}
            pairs, maximum, passed = _comparison_gate(values, core_selector(6, points), 0.03)
            report.update(pairwise_core=pairs, maximum_core_pairwise_relative_rms=maximum,
                          pairwise_raw=_pairwise(values, np.ones(len(points), dtype=bool)),
                          passed=passed, completed=True)
        else:
            report['completed'] = True
        print('COMPLETE '+str(args.output), flush=True)
    except Exception as exc:
        report['exception'] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        report['elapsed_s'] = time.perf_counter()-start
        save()
    return 0 if report['passed'] or args.engine != 'both' else 2


if __name__ == '__main__':
    raise SystemExit(main())
