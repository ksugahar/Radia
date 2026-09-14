"""Run the energy-Newton regression from an installed candidate and retain telemetry."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import subprocess
import time

import ngsolve as ng
import numpy as np
import pytest
import radia


ROOT = Path(__file__).resolve().parent


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    dependency_check = subprocess.run([sys.executable, '-m', 'pip', 'check'],
                                      capture_output=True, text=True, check=True)
    package = Path(radia.__file__).resolve().parent
    wheel = ROOT / 'radia-4.95.91-cp312-cp312-win_amd64.whl'
    native = package / '_radia_pybind.pyd'
    if sha(wheel) != 'e9e6b9f09d105552330ac9b2b6d870187d825d11bf7a78f407553c8c7a866e25':
        raise RuntimeError('Unexpected wheel')
    if sha(native) != '932e49626c05c1e3d459f20649632d4bfa0adf9b37f325dace4093d74aeeb7b5':
        raise RuntimeError('Unexpected native binary')
    if not package.is_relative_to(ROOT / 'venv'):
        raise RuntimeError('Candidate is not installed in the dedicated venv')
    import zipfile
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            if name.startswith('radia/') and name.endswith('.py'):
                if archive.read(name) != (package / name[6:]).read_bytes():
                    raise RuntimeError('Installed Python differs from wheel: ' + name)
    test = ROOT / 'test_hdiv_vim_energy_newton.py'
    report = dict(schema='radia.energy-newton-native-audit.v1', host=platform.node(),
                  python=sys.version, executable=sys.executable, package=str(package),
                  versions={name: importlib.metadata.version(name) for name in
                            ('radia', 'ngsolve', 'numpy', 'scipy', 'pytest', 'threadpoolctl')},
                  pip_check=dependency_check.stdout.strip(),
                  wheel_sha256=sha(wheel), native_sha256=sha(native),
                  test_sha256=sha(test), runner_sha256=sha(__file__),
                  threads=2, solves=[], tests=[], completed=False)
    ng.SetNumThreads(2)

    def save():
        (ROOT / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')

    audit = report
    class Recorder:
        def pytest_runtest_setup(self, item):
            original = item.module.Solve
            if getattr(original, '_audit_wrapper', False):
                return

            def solve(*args, **kwargs):
                start = time.perf_counter()
                entry = {'test': item.nodeid, 'preconditioner': kwargs.get('preconditioner', 'auto'),
                         'solver': kwargs.get('nonlinear_solver', 'default')}
                try:
                    result = original(*args, **kwargs)
                    stats = result.get('nonlinear_solve_stats', {})
                    entry.update(M_avg=np.asarray(result['M_avg']).tolist(),
                                 statistics={k: v.item() if isinstance(v, np.generic) else v
                                             for k, v in stats.items()
                                             if isinstance(v, (str, bool, int, float, np.generic))})
                    return result
                except Exception as exc:
                    entry['exception'] = {'type': type(exc).__name__, 'message': str(exc)}
                    raise
                finally:
                    entry['elapsed_s'] = time.perf_counter() - start
                    report['solves'].append(entry)
                    save()

            solve._audit_wrapper = True
            # Each setup refreshes the current test name without nesting wrappers.
            solve._original = original
            item.module.Solve = solve

        def pytest_runtest_teardown(self, item):
            wrapped = item.module.Solve
            if getattr(wrapped, '_audit_wrapper', False):
                item.module.Solve = wrapped._original

        def pytest_runtest_logreport(self, report):
            if report.when == 'call' or report.failed or report.skipped:
                audit['tests'].append(dict(nodeid=report.nodeid, phase=report.when,
                                           outcome=report.outcome, seconds=report.duration))
                save()

    save()
    start = time.perf_counter()
    code = pytest.main([str(test), '-vv', '--capture=tee-sys',
                        '--junitxml=' + str(ROOT / 'junit.xml')], plugins=[Recorder()])
    report.update(completed=True, exit_code=int(code), elapsed_s=time.perf_counter() - start,
                  passed=code == 0 and len(report['tests']) == 9
                  and all(test['outcome'] == 'passed' for test in report['tests']))
    save()
    return int(code) if code else (0 if report['passed'] else 2)


if __name__ == '__main__':
    raise SystemExit(main())
