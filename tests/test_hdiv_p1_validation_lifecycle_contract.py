"""Guards: invalid pins/admission fail before Radia import; no premature adoption.

Assumes: synthetic source-only control inputs, not numerical or GPU acceptance.
"""
import builtins
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_driver():
    path = ROOT/'validation_test/hmatrix_cuda/validate_p1_charge_diagonal.py'
    spec = importlib.util.spec_from_file_location('p1_validation_contract', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('failure', ['missing_cpu_pins', 'cpu_digest', 'oracle_digest', 'expired_gpu'])
def test_invalid_frame_never_imports_radia(tmp_path, monkeypatch, failure):
    driver = load_driver()
    source_sha = '1'*40
    def git_output(argv, **kwargs):
        if 'rev-parse' in argv:
            return source_sha+'\n'
        if 'status' in argv:
            return ''
        raise AssertionError(f'unexpected external command: {argv}')
    monkeypatch.setattr(driver.subprocess, 'check_output', git_output)
    imported = []
    original = builtins.__import__
    def guarded_import(name, *args, **kwargs):
        if name == 'radia' or name.startswith('radia.'):
            imported.append(name)
            raise AssertionError('Radia import before invalid-frame rejection')
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', guarded_import)
    output = tmp_path/'result.json'
    args = ['--backend', 'cpu' if failure.startswith(('missing', 'cpu')) else 'gpu',
            '--git', str(tmp_path/'git'), '--mesh', str(tmp_path/'mesh.vol'),
            '--contract', str(tmp_path/'labels.json'), '--receipt', str(tmp_path/'receipt.json'),
            '--receipt-sha256', '0'*64, '--receipt-provenance', str(tmp_path/'provenance.json'),
            '--out', str(output)]
    if failure == 'cpu_digest':
        binary = tmp_path/'cpu.pyd'
        binary.write_bytes(b'synthetic')
        args += ['--source-root', str(tmp_path), '--cpu-source-sha', source_sha,
                 '--cpu-native', str(binary), '--cpu-native-sha256', '0'*64]
    elif failure in ('oracle_digest', 'expired_gpu'):
        oracle = tmp_path/'oracle.json'
        oracle.write_text(json.dumps(dict(schema=driver.SCHEMA, backend='cpu', execution_complete=True,
                                         cases=[dict(case=name) for name in driver.CASES])), encoding='utf-8')
        digest = hashlib.sha256(oracle.read_bytes()).hexdigest()
        admission = tmp_path/'admission.json'
        admission.write_text(json.dumps(dict(gpu_execution_approved=True, source_sha=source_sha,
                                             oracle_sha256=digest, deadline_utc='2000-01-01T00:00:00+00:00')), encoding='utf-8')
        args += ['--oracle', str(oracle), '--oracle-sha256', '0'*64 if failure=='oracle_digest' else digest,
                 '--admission', str(admission), '--admission-sha256', hashlib.sha256(admission.read_bytes()).hexdigest(),
                 '--nvidia-smi', str(tmp_path/'nvidia-smi')]
    assert driver.main(args) == 1, f'{failure}: actual completed; expected failure'
    result = json.loads(output.read_bytes())
    assert imported == [], f'actual imports={imported}; expected no Radia import'
    assert result['execution_complete'] is False
    assert result['numerically_accepted'] is False
    assert result['criterion_status'] == 'pending-user-approval'
    assert 'adoption_result' not in result


def test_existing_evidence_is_never_overwritten(tmp_path):
    driver = load_driver()
    output = tmp_path/'result.json'
    output.write_bytes(b'preserved raw evidence')
    args = ['--backend', 'cpu', '--git', str(tmp_path/'git'), '--mesh', str(tmp_path/'mesh'),
            '--contract', str(tmp_path/'contract'), '--receipt', str(tmp_path/'receipt'),
            '--receipt-sha256', '0'*64, '--receipt-provenance', str(tmp_path/'provenance'), '--out', str(output)]
    with pytest.raises(FileExistsError, match='Preserve'):
        driver.main(args)
    assert output.read_bytes() == b'preserved raw evidence'
