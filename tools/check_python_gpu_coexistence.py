"""Check same-source Python CUDA module coexistence before native wheel packaging.

This imports real payloads in two fresh processes, one per complete module order.
It never constructs a resident object or requests a solve. Exact-wheel numerical
acceptance remains separate. Generate source metadata with setup.py egg_info.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import sysconfig

ROOT = Path(__file__).resolve().parents[1]
MODULES = ('_gpu_csr', '_gpu_resident_amg', '_gpu_ams', '_hdiv_cuda')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resident_types(modules):
    """Require distinct native Python types while preserving both public names."""
    amg, ams = (getattr(modules[name], 'Resident', None)
                for name in ('_gpu_resident_amg', '_gpu_ams'))
    if (not isinstance(amg, type) or not isinstance(ams, type)
            or amg is ams or amg.__name__ != 'Resident' or ams.__name__ != 'Resident'):
        raise RuntimeError('AMG and AMS must expose distinct Resident types')
    return {'amg': amg.__module__ + '.' + amg.__name__,
            'ams': ams.__module__ + '.' + ams.__name__, 'distinct': True}


def source_identity(source_commit):
    if not re.fullmatch(r'[0-9a-f]{40}', source_commit):
        raise ValueError('An exact source commit is required')
    git = ['git', '-c', 'safe.directory=' + ROOT.as_posix(), '-C', str(ROOT)]
    head = subprocess.check_output(git + ['rev-parse', 'HEAD'], text=True).strip()
    dirty = subprocess.check_output(git + ['status', '--porcelain', '--untracked-files=all'], text=True)
    if head != source_commit or dirty.strip():
        raise RuntimeError('Coexistence requires the clean exact source checkout')


def payload_receipt(module, package, source_commit):
    binary = Path(module.__file__).resolve()
    if binary.parent != package:
        raise RuntimeError('Native module loaded outside the selected source payload')
    receipt = Path(str(binary) + '.build.json')
    record = json.loads(receipt.read_text(encoding='utf-8'))
    if (record.get('schema') != 'radia.native-build-provenance.v1'
            or record.get('source_commit') != source_commit
            or record.get('source_dirty') is not False
            or record.get('binary_name') != binary.name
            or record.get('binary_sha256') != sha(binary)):
        raise RuntimeError('Native CPU/GPU payload provenance mismatch: ' + binary.name)
    return {'binary_name': binary.name, 'binary_sha256': sha(binary),
            'receipt_sha256': sha(receipt), 'source_commit': source_commit}


def check_order(order, source_commit):
    source_identity(source_commit)
    package = (ROOT / 'src/radia').resolve()
    if 'radia' in sys.modules:
        raise RuntimeError('Coexistence must own the initial Radia import')
    spec = importlib.util.spec_from_file_location('coexistence_payload', ROOT / 'tools/python_gpu_payload.py')
    payload = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(payload)
    for name in MODULES:
        binary = package / (name + sysconfig.get_config_var('EXT_SUFFIX'))
        record = json.loads(Path(str(binary) + '.build.json').read_text(encoding='utf-8'))
        payload.verify_payload(binary.name, binary.read_bytes(), record, source_commit,
                               architectures=('75', '80') if sys.platform == 'win32' else ('100',))
        for relative, digest in record['source_hashes'].items():
            if sha(ROOT / relative) != digest:
                raise RuntimeError('Compiled GPU source closure mismatch: ' + relative)
    # Explicit source-payload contract, not an installed-wheel acceptance route.
    sys.path.insert(0, str(ROOT / 'src'))
    import radia
    if Path(radia.__file__).resolve() != package / '__init__.py':
        raise RuntimeError('Unexpected Radia source import')
    dist = metadata.distribution('radia')
    if (Path(dist.locate_file('radia/__init__.py')).resolve() != package / '__init__.py'
            or dist.version != radia.__version__):
        raise RuntimeError('Generate matching source metadata with setup.py egg_info')
    from radia.gpu import _native
    from radia._cuda_runtime import prepare_cuda_runtime
    loaded = {}
    binaries = {}
    for name in order:
        module = _native(name)  # Retains the production CPU/GPU same-source gate.
        if _native(name) is not module:
            raise RuntimeError('Repeated native import changed module identity: ' + name)
        loaded[name] = module
        binaries[name] = payload_receipt(module, package, source_commit)
    types = resident_types(loaded)
    import radia._radia_pybind as core
    binaries['_radia_pybind'] = payload_receipt(core, package, source_commit)
    runtime = prepare_cuda_runtime()
    libraries = {name: {key: value for key, value in item.items() if key != 'path'}
                 for name, item in runtime['libraries'].items()}
    source_identity(source_commit)
    return dict(passed=True, order=list(order), binaries=binaries, resident_types=types,
                runtime_pins=runtime['pins'], runtime_libraries=libraries)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--order', choices=('forward', 'reverse'), help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.report.exists():
        raise FileExistsError('Coexistence report must be fresh')
    if ROOT == args.report.resolve() or ROOT in args.report.resolve().parents:
        raise ValueError('Keep coexistence evidence outside the source checkout')
    result = dict(schema='radia.python-gpu-coexistence.v1', source_commit=args.source_commit,
                  passed=False, numerical_qualified=False, solve_requested=False,
                  resident_construction_requested=False, device_context_checked=False)
    try:
        if args.order:
            order = MODULES if args.order == 'forward' else MODULES[::-1]
            result.update(check_order(order, args.source_commit))
        else:
            source_identity(args.source_commit)
            orders = []
            for order in ('forward', 'reverse'):
                report = args.report.with_name(args.report.stem + '.' + order + '.json')
                command = [sys.executable, '-I', '-B', str(Path(__file__).resolve()),
                           '--source-commit', args.source_commit, '--report', str(report), '--order', order]
                process = subprocess.run(command, capture_output=True, text=True, timeout=240)
                if not report.is_file():
                    raise RuntimeError('Coexistence worker produced no report: ' + order)
                record = json.loads(report.read_text(encoding='utf-8'))
                orders.append(record)
                result['orders'] = orders
                if process.returncode or record.get('passed') is not True:
                    raise RuntimeError('Coexistence failed in ' + order + ': ' + str(record.get('error')))
            if orders[0]['binaries'] != orders[1]['binaries']:
                raise RuntimeError('Payload identity changed between import orders')
            result['passed'] = True
    except Exception as error:
        result['error'] = (type(error).__name__ + ': ' + str(error)).replace(str(ROOT), '<source>').replace(str(Path(sys.prefix)), '<runtime>')
    finally:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
