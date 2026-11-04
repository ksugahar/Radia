"""Run the HDiv GPU device tests against an installed exact-SHA radia wheel.

Compute-lane entry: nothing is compiled and no checkout binary is loaded. The
checkout supplies only the test files; every module under test comes from the
wheel. A failed, skipped or missing case fails the lane, and the result JSON is
written either way. Parity of one discretization per test docstring; this is not
independent accuracy evidence.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
GATE = 'RADIA_HDIV_NATIVE_TEST_DIR'
SCHEMA = 'radia.hdiv-installed-lane.v1'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lane_tests(root=ROOT):
    """Every test gated on the native directory belongs to this lane."""
    return sorted(p.relative_to(root).as_posix() for p in (root / 'tests').glob('test_*.py')
                  if GATE in p.read_text(encoding='utf-8'))


def read_cases(junit):
    cases = []
    for case in ET.parse(junit).getroot().iter('testcase'):
        entry = dict(module=case.get('classname'), name=case.get('name'),
                     time_s=float(case.get('time', 0)))
        for tag in ('failure', 'error', 'skipped'):
            node = case.find(tag)
            if node is not None:
                entry[tag] = (node.get('message') or '')[:2000]
        cases.append(entry)
    return cases


def lane_failures(returncode, tests, cases):
    """A green pytest exit is not enough: a skipped lane must not read as a pass."""
    failures = [] if returncode == 0 else [f'pytest returned {returncode}']
    for test in tests:
        module = Path(test).with_suffix('').as_posix().replace('/', '.')
        if not any((c['module'] or '').startswith(module) for c in cases):
            failures.append(f'no case ran from {test}')
    for case in cases:
        for tag in ('failure', 'error', 'skipped'):
            if tag in case:
                failures.append(f"{tag}: {case['module']}::{case['name']}")
    return failures


def check_installed(source_commit):
    """Bind the imported package to the venv and both native receipts to the source."""
    import radia
    import radia._radia_pybind as core
    package = Path(radia.__file__).resolve().parent
    if not package.is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(f'radia imported outside the run venv: {package} not under {sys.prefix}')
    binaries = {}
    for binary in (Path(core.__file__), package / ('_hdiv_cuda' + sysconfig.get_config_var('EXT_SUFFIX'))):
        receipt = json.loads(Path(str(binary) + '.build.json').read_text(encoding='utf-8'))
        if receipt.get('source_commit') != source_commit:
            raise RuntimeError(f"{binary.name} source {receipt.get('source_commit')} != {source_commit}")
        binaries[binary.name] = digest(binary)
    return package, binaries


def check_netgen(overlay):
    """The netgen of this venv must hold the patched library; test_hdiv_pyramid_material assumes it."""
    from importlib.util import find_spec
    receipt = json.loads((overlay / 'netgen-overlay.json').read_text(encoding='utf-8'))
    package = Path(find_spec('netgen').submodule_search_locations[0]).resolve()
    if not package.is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(f'netgen imported outside the run venv: {package} not under {sys.prefix}')
    actual = digest(package / receipt['library'])
    if actual != receipt['library_sha256']:
        raise RuntimeError(f"{receipt['library']} is {actual}, overlay {receipt['library_sha256']}")
    return {key: receipt[key] for key in ('upstream_sha', 'patch_sha256', 'library_sha256')}


def query_gpu(required):
    names = subprocess.check_output(['nvidia-smi', '--query-gpu=name,driver_version,memory.total',
                                     '--format=csv,noheader'], text=True).strip().splitlines()
    if not any(required in name for name in names):
        raise RuntimeError(f'Required GPU {required!r} not present: {names}')
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--wheel', required=True, type=Path)
    parser.add_argument('--require-gpu', required=True)
    parser.add_argument('--netgen-overlay', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = dict(schema=SCHEMA, source_commit=args.source_commit, wheel=args.wheel.name,
                  wheel_sha256=digest(args.wheel), passed=False, cases=[],
                  threads={k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'RADIA_NUM_THREADS', 'MKL_NUM_THREADS')},
                  logical_cpus=os.cpu_count())
    try:
        result['gpu'] = query_gpu(args.require_gpu)
        package, result['binaries'] = check_installed(args.source_commit)
        result['netgen'] = check_netgen(args.netgen_overlay)
        import ngsolve
        import pytest
        result['versions'] = dict(python=sys.version.split()[0], ngsolve=ngsolve.__version__, pytest=pytest.__version__)
        tests = lane_tests()
        result['tests'] = {test: digest(ROOT / test) for test in tests}
        os.environ[GATE] = str(package)
        junit = args.out / 'junit.xml'
        os.chdir(ROOT)
        # tests/conftest.py puts the checkout src first; this lane must import the wheel.
        result['pytest_returncode'] = int(pytest.main(
            tests + ['-q', '-rs', '--noconftest', '-p', 'no:cacheprovider', '-o', 'addopts=', '-o', 'junit_family=xunit1',
                     '--basetemp', str(args.out / 'pytest-tmp'), '--junitxml', str(junit)]))
        result['cases'] = read_cases(junit)
        result['failures'] = lane_failures(result['pytest_returncode'], tests, result['cases'])
        result['passed'] = not result['failures']
    except Exception as error:
        result['error'] = f'{type(error).__name__}: {error}'
    result['elapsed_s'] = time.monotonic() - start
    (args.out / 'result.json').write_text(json.dumps(result, indent=1) + '\n', encoding='utf-8')
    print(json.dumps({k: result.get(k) for k in ('passed', 'failures', 'error', 'elapsed_s')}))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
