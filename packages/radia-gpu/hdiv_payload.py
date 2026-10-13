"""Stage explicitly selected, exact-commit bounded HDiv CUDA artifacts."""
from pathlib import Path
import hashlib
import json
import shutil

WRAPPERS = ('solveHDivElements.m', 'affineHDivQuadrature.m', 'curvedHDivQuadrature.m')


def stage(root, package, commit, hmatrix, projection):
    root, package = Path(root), Path(package)
    receipts = {}
    expected = set()
    for name, folder in [('hmatrix', hmatrix), ('projection', projection)]:
        folder = Path(folder)
        receipt = json.loads((folder / 'build.json').read_text(encoding='utf-8'))
        if (not receipt.get('accepted') or receipt.get('source_dirty', True)
                or receipt.get('source_commit') != commit
                or receipt.get('taskmanager_started', True)):
            raise RuntimeError('HDiv requires clean exact-commit native provenance: ' + name)
        for relative, digest in receipt['source_sha256'].items():
            if Path(relative).is_absolute() or '..' in Path(relative).parts:
                raise RuntimeError('Invalid HDiv source path')
            source = root / relative
            if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != digest:
                raise RuntimeError('HDiv source mismatch: ' + relative)
        binary = folder / ('radia_cuda_' + name + '.mexw64')
        digest = hashlib.sha256(binary.read_bytes()).hexdigest()
        if digest != receipt['binary_sha256']:
            raise RuntimeError('HDiv binary mismatch: ' + name)
        target = package / 'matlab' / binary.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(binary, target)
        expected.add(target.relative_to(package).as_posix())
        receipts[name] = dict(source_commit=commit, source_dirty=False,
                              binary_sha256=digest, source_sha256=receipt['source_sha256'],
                              matlab=receipt['matlab'], taskmanager_started=False)
    for name in WRAPPERS:
        relative = Path('matlab/+radia/+gpu') / name
        target = package / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / 'packages/radia-gpu' / relative, target)
        expected.add(relative.as_posix())
    record = dict(schema='radia.gpu.hdiv-payload.v1', artifacts=receipts,
                  scope='Experimental bounded affine/P2 tetrahedra, HDiv orders 1/2, linear isotropic material',
                  limitations=['Dense near-rule Gram; general adaptive high-order HACApK unsupported.',
                               'Images, internal interfaces, constraints and nonlinear material unaccepted.',
                               'User-installed compatible NGSolve native libraries and Python DLL required; no Python solver adapter.'])
    (package / 'hdiv-build.json').write_text(json.dumps(record, indent=2) + '\n')
    expected.add('hdiv-build.json')
    return expected
