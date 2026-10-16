"""Owned fresh-Engine transaction/schema validation; requires separate resource admission."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test-artifact', type=Path, required=True)
    parser.add_argument('--production-artifact', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reviewed-block-units', type=Path, required=True,
                        help='IH-reviewed JSON list in complete assembler block order')
    parser.add_argument('--binary-source-sha', required=True,
                        help='Exact reviewed native build commit; independently bound to harness HEAD')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    head = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', '-C', str(root), 'status', '--porcelain'], text=True).strip():
        raise RuntimeError('Clean reviewed source required')
    binary_source = subprocess.check_output(
        ['git', '-C', str(root), 'rev-parse', args.binary_source_sha + '^{commit}'], text=True).strip()
    if binary_source != args.binary_source_sha:
        raise RuntimeError('Full exact binary source SHA required')
    allowed_diagnostic_changes = {
        'validation_test/ih_cuda/checkTransactionSnapshot.m',
        'validation_test/ih_cuda/run_transaction_snapshot.py',
    }
    changed = set(subprocess.check_output(
        ['git', '-C', str(root), 'diff', '--name-only', binary_source, head], text=True).splitlines())
    if not changed <= allowed_diagnostic_changes:
        raise RuntimeError('Unreviewed binary/harness source differences: ' + ', '.join(sorted(changed)))
    receipts = []
    for path, test in [(args.test_artifact, True), (args.production_artifact, False)]:
        receipt = json.loads((path / 'build.json').read_text())
        if (receipt.get('accepted') is not True or receipt.get('source_dirty') is not False
                or receipt.get('source_commit') != binary_source
                or receipt.get('transaction_test_build') is not test
                or receipt.get('distribution_eligible') is not (not test)):
            raise RuntimeError('Receipt source/classification mismatch')
        if sha(path / 'radia_cuda_ih.mexw64') != receipt['binary_sha256']:
            raise RuntimeError('MEX hash mismatch')
        for name, digest in receipt['source_sha256'].items():
            if sha(root / name) != digest:
                raise RuntimeError('Build source inventory mismatch: ' + name)
        receipts.append(receipt)
    units = json.loads(args.reviewed_block_units.read_text())
    if units != ['A', 'V', 'A', 'A', 'A', 'A', 'A/m', 'V', 'V']:
        raise RuntimeError('Reviewed complete-block unknown units mismatch')
    fixture = root / 'validation_test/ih_cuda'
    wrappers = root / 'packages/radia-gpu/matlab'
    identity = dict(source_commit=head, validation_harness_sha=head, binary_source_sha=binary_source,
                    allowed_diagnostic_changes=sorted(changed),
                    route='official Engine; fresh operation-owned session',
                    receipts=receipts, block_units_sha256=sha(args.reviewed_block_units),
                    row_units=['V', 'A', 'A', 'A', 'A m^2', 'A m^2', 'A m', 'V', 'V'],
                    gauge_interpretation='Scaled Lagrange multipliers; not independently measured physical potentials',
                    source_sha256={str(p.relative_to(root)): sha(p) for p in [
                        Path(__file__), fixture / 'checkTransactionSnapshot.m',
                        fixture / 'makeTemperatureBHConfig.m',
                        wrappers / '+radia/+simulink/makeCompletePEECBEMConfig.m']})
    if args.out.exists():
        raise RuntimeError('New owned output directory required')
    import matlab.engine
    eng = matlab.engine.start_matlab('-nodesktop -nosplash')
    try:
        eng.addpath(str(wrappers), str(fixture), nargout=0)
        identity['matlab'] = eng.version()
        identity['owned_matlab_pid'] = int(eng.feature('getpid'))
        eng.workspace['reviewed_units'] = units
        eng.workspace['radia_transaction_owned_session'] = True
        eng.workspace['test_artifact'] = str(args.test_artifact.resolve())
        eng.workspace['production_artifact'] = str(args.production_artifact.resolve())
        eng.workspace['output_directory'] = str(args.out.resolve())
        eng.eval('transaction_result=checkTransactionSnapshot(test_artifact,production_artifact,output_directory,reviewed_units);', nargout=0)
        result = json.loads((args.out / 'result.json').read_text(encoding='utf-8'))
        if result.get('accepted') is not True:
            raise RuntimeError('Native transaction/schema result rejected; see owned result.json')
    finally:
        eng.quit()
        if args.out.is_dir():
            (args.out / 'launch-identity.json').write_text(json.dumps(identity, indent=2) + '\n')


if __name__ == '__main__':
    main()
