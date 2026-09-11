"""Re-solve mixed Omega with controlled assembly quadrature (p <= 3).

A factory module supplies create_case(mesh_path), returning mesh, H_s,
Kelvin H_s, radius, center and mu_r. Sources are constructed once per run.
This diagnostic lane never grants three-engine field acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time

import ngsolve as ng
import numpy as np

from diagnostics import audit_energy, constraint_violation


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def identity(factory, solver):
    import radia
    root = Path(radia.__file__).resolve().parent
    files = {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob('*.py'))}
    natives = {p.name: digest(p) for p in root.glob('_radia_pybind.*')
               if p.suffix in ('.pyd', '.so')}
    if not natives:
        raise RuntimeError('No Radia native extension identified')
    return {'version': getattr(radia, '__version__', None), 'module': str(root),
            'native': natives, 'python_files': files, 'ngsolve': ng.__version__,
            'solver_file': str(Path(solver.__file__).resolve()),
            'solver_sha256': digest(solver.__file__), 'factory_sha256': digest(factory),
            'runner_sha256': digest(__file__),
            'diagnostics_sha256': digest(Path(__file__).with_name('diagnostics.py'))}


def check_runtime(mode):
    import radia
    dist = importlib.metadata.distribution('radia')
    direct = json.loads(dist.read_text('direct_url.json') or '{}')
    if mode == 'wheel':
        if direct.get('dir_info', {}).get('editable', False):
            raise RuntimeError('wheel mode rejects editable distributions')
        actual = Path(radia.__file__).resolve()
        expected = Path(dist.locate_file('radia/__init__.py')).resolve()
        if not actual.samefile(expected):
            raise RuntimeError('Imported Radia is not the installed wheel')
    return {'mode': mode, 'direct_url': direct}


def gates(result):
    residual = result['linear_residual']
    relative = residual['free_dofs']['relative']
    audit = result['assembled_energy']
    ident = audit['fixed_rule_identity']
    scale = max(abs(ident['W']), abs(ident['J']), abs(ident['source_offset']), 1e-30)
    checks = {
        'free_residual': relative is not None and np.isfinite(relative) and relative <= 1e-8,
        'same_rule_energy': abs(ident['W_minus_J_minus_offset']) <= 1e-10 * scale,
        'assembly_reconstruction': abs(audit['default_reconstruction_difference'])
        <= 1e-10 * max(abs(audit['energy']), 1e-30),
    }
    for name, row in residual['blocks'].items():
        value = row['relative'] if row is not None else None
        checks['residual_' + name] = (value is not None and np.isfinite(value) and value <= 1e-8)
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--factory', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('research', 'wheel'), required=True)
    parser.add_argument('--research-solver', type=Path)
    parser.add_argument('--orders', nargs='+', type=int, choices=(1, 2, 3), default=[1, 2, 3])
    parser.add_argument('--bonuses', nargs='+', type=int, default=[4, 8, 12])
    parser.add_argument('--evaluation-order', type=int, default=16)
    parser.add_argument('--source-order', type=int, default=3)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if min(args.bonuses) < 0 or args.evaluation_order < 1:
        parser.error('quadrature orders must be nonnegative (evaluation >= 1)')
    if args.mode == 'wheel' and args.research_solver:
        parser.error('wheel mode forbids Python solver overrides')
    runtime = check_runtime(args.mode)
    if args.research_solver:
        solver = load_module(args.research_solver, 'radia._quadrature_research_solver')
    else:
        from radia import kelvin_solver as solver
    before = identity(args.factory, solver)
    ng.SetNumThreads(args.threads)
    case = load_module(args.factory, '_omega_case').create_case(args.mesh)
    mesh, h_s, h_ext = case['mesh'], case['H_s'], case['H_ext']
    rows, nesting = [], []
    payload = {'schema': 'radia.validation.omega-quadrature.v1', 'host': platform.node(),
               'runtime': runtime, 'implementation': before, 'controls': vars(args).copy(),
               'mesh_sha256': digest(args.mesh), 'mesh_elements': mesh.ne,
               'case': case['controls'], 'rows': rows, 'nesting': nesting,
               'embedding_quadrature_order': 16,
               'assembly_scope': 'bonus applies to volume and interface terms; not load-only',
               'acceptance': 'HOLD: diagnostics are not three-engine field acceptance'}
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        payload['source_unchanged'] = identity(args.factory, solver) == before
        payload['mesh_unchanged'] = digest(args.mesh) == payload['mesh_sha256']
        dependency = case['controls'].get('coil_source_file')
        if dependency:
            payload['source_unchanged'] &= digest(dependency) == case['controls']['coil_source_sha256']
        args.output.write_text(json.dumps(payload, indent=2, default=str, allow_nan=False), encoding='utf-8')
        if not payload['source_unchanged'] or not payload['mesh_unchanged']:
            raise RuntimeError('Source or mesh changed during validation')

    with ng.TaskManager():
        source = solver.project_source_total_hodge(mesh, h_s, ('iron',), order=args.source_order)
        payload['source_harmonic_norm'] = source['relative_harmonic_norm']
        for bonus in sorted(set(args.bonuses)):
            results = {}
            for order in sorted(set(args.orders)):
                start = time.perf_counter()
                result = solver.solve_magnetostatic_mixed_total_reduced_omega_kelvin(
                    mesh, h_s, source['potential'], case['radius'], case['center'],
                    mu_r_by_material={'iron': case['mu_r']}, reduced_materials=('air',),
                    total_materials=('iron', 'kelvin'), interface_boundary='iron_air_interface',
                    order=order, bonus_intorder=bonus, dirichlet_bbbnd='GND', inverse='pardiso',
                    kelvin_mats=('kelvin',), kelvin_interface_boundary='kelvin_int',
                    kelvin_source_h=h_ext, total_source_h=source['harmonic_field'],
                    total_source_materials=('iron',), return_system=True)
                solve_s = time.perf_counter() - start
                audit_energy(result, mesh, h_s, h_ext, source['harmonic_field'], bonus, args.evaluation_order)
                rows.append({'order': order, 'bonus': bonus, 'solve_s': solve_s,
                             'total_s': time.perf_counter()-start, 'ndof': result['fes'].ndof,
                             'linear_residual': result['linear_residual'],
                             'energy': result['assembled_energy'], 'gates': gates(result)})
                results[order] = result
                save()
                print('solved', bonus, order, rows[-1]['gates'], flush=True)
            for low in results:
                for high in results:
                    if low <= high:
                        check = constraint_violation(results[low], results[high], mesh)
                        nesting.append({'bonus': bonus, 'low': low, 'high': high, **check})
                        save()
            del results
    payload['completed'] = True
    save()
    return 0 if all(all(row['gates'].values()) for row in rows) else 2


if __name__ == '__main__':
    raise SystemExit(main())
