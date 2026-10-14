"""Isolate ESRF6 nonlinear convergence in an identified installed repair wheel."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import traceback

import ngsolve as ng
import numpy as np
import radia as rad
from radia import vim
from radia.esrf_examples import get_esrf_bh_table

from esrf_coil_yoke import (
    average_observation_field, build_radia_coil_source,
    observation_volume_quadrature,
)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iron', type=Path, required=True)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--order', type=int, choices=(1, 2), default=1)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--mode', choices=('mass-riesz', 'picard-warmstart', 'picard'),
                        default='mass-riesz')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    dist = importlib.metadata.distribution('radia')
    direct = json.loads(dist.read_text('direct_url.json') or '{}')
    if direct.get('dir_info', {}).get('editable'):
        raise RuntimeError('An installed non-editable candidate wheel is required')
    package = Path(rad.__file__).resolve().parent
    native = package / '_radia_pybind.pyd'
    identities = {
        'iron': (args.iron, 'fdd1387281e550e27eab435cfac3ec54e0b9d1339438699fdc50c852dfcafc74'),
        'wheel': (args.wheel, 'e9e6b9f09d105552330ac9b2b6d870187d825d11bf7a78f407553c8c7a866e25'),
        'native': (native, '932e49626c05c1e3d459f20649632d4bfa0adf9b37f325dace4093d74aeeb7b5'),
    }
    verified = {}
    for name, (path, expected) in identities.items():
        actual = sha(path)
        if actual != expected:
            raise RuntimeError(f'{name} identity mismatch: {actual}')
        verified[name] = {'path': str(path.resolve()), 'sha256': actual}
    python_sources = {str(p.relative_to(package)): sha(p)
                      for p in sorted(package.rglob('*.py'))}
    options = dict(order=args.order, gram_eps=1e-12, leaf=64, eta=2.0,
                   curve_order=None, curve_gauss=8, image=None,
                   tol=1e-8, maxit=12000, nl_tol=2e-5, nl_maxit=80,
                   preconditioner='mass-riesz', nonlinear_solver='energy-newton',
                   newton_warmstart='linear')
    if args.mode == 'picard-warmstart':
        options['newton_warmstart'] = 'picard'
    if args.mode == 'picard':
        options['nonlinear_solver'] = 'picard-mass-riesz'
    ng.SetNumThreads(args.threads)
    rad.UtiDelAll()
    coil, source = build_radia_coil_source(6)
    bh = np.asarray(get_esrf_bh_table(6), dtype=float)
    points, evaluation_points = observation_volume_quadrature(6, half_width_m=2e-5)
    report = dict(schema='radia.validation.esrf6-repaired-hdiv.v1',
                  source_commit='4d85e72cce70e02417abcd7972b509dfc2bac632',
                  case=6, mode=args.mode, options=options, threads=args.threads,
                  argv=sys.argv, implementation=dict(version=rad.__version__,
                  radia_path=str(package), ngsolve=ng.__version__,
                  python=sys.version, identities=verified, python_sources=python_sources,
                  driver_sha256=sha(__file__),
                  source_adapter_sha256=sha(sys.modules['esrf_coil_yoke'].__file__)),
                  source=source, bh_table=bh.tolist(), observation_points_m=points.tolist(),
                  observation_half_width_m=2e-5, completed=False, accepted=False,
                  three_method_acceptance='HOLD: HDiv-only diagnostic; mixed quadrature unresolved')
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')

    save()
    start = time.perf_counter()
    print(f'START case6 {args.mode} BDM{args.order}, threads={args.threads}', flush=True)
    try:
        mesh = ng.Mesh(str(args.iron))
        with ng.TaskManager():
            result = vim.Solve(mesh, H_ext=rad.RadiaField(coil, 'h'), bh_table=bh, **options)
            report['solve_returned_s'] = time.perf_counter()-start
            stats = dict(result['nonlinear_solve_stats'])
            report['nonlinear_stats'] = stats
            report['ndof'] = int(result['ndof'])
            report['M_avg_Am'] = np.asarray(result['M_avg']).tolist()
            report['M_elements_Am'] = np.asarray(result['M']).tolist()
            report['timings'] = {k: float(v) for k, v in result.items()
                                 if k.endswith('_wall_s') and np.isscalar(v)}
            save()
            print('SOLVE '+json.dumps(stats), flush=True)
            demag = np.asarray(vim.FieldFromSolution(result, evaluation_points, algorithm='direct'))
        coil_h = np.asarray(rad.Fld(coil, 'h', evaluation_points))
        report['B_T'] = average_observation_field(4e-7*np.pi*(demag+coil_h), len(points)).tolist()
        residual = stats.get('nonlinear_final_relative_residual')
        report['accepted'] = bool(residual is not None and np.isfinite(residual)
                                  and residual <= options['nl_tol']
                                  and stats.get('nonlinear_converged_final_stage') is True)
        report['residual_contract'] = ('true assembled energy-Newton residual' if residual is not None
                                       else 'not reported by this path; diagnostic only, not accepted')
        report['completed'] = True
    except Exception as exc:
        from radia.vim import _solve
        report['nonlinear_stats'] = dict(_solve._LAST_NONLINEAR_SOLVE_STATS)
        report['exception'] = dict(type=type(exc).__name__, message=str(exc),
                                   traceback=traceback.format_exc())
        raise
    finally:
        report['elapsed_s'] = time.perf_counter()-start
        save()
        print('RESULT '+str(args.output), flush=True)


if __name__ == '__main__':
    main()
