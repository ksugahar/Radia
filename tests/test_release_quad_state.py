import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[1] / 'tools' / 'release_quad.py'
SPEC = importlib.util.spec_from_file_location('release_quad_state', TOOL)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def state(**targets):
    return {'schema': 'test', 'commit': 'a', 'package_sha256': 'hash', 'targets': targets}


def test_stale_snapshot_does_not_erase_another_host(tmp_path):
    path = tmp_path / 'state.json'
    module._write_simulink_state(path, state(lab={'status': 'passed'}), 'lab')
    stale = state(lab={'status': 'passed'}, mdx1={'status': 'passed'})
    module._write_simulink_state(path, state(lab={'status': 'failed'}), 'lab')
    module._write_simulink_state(path, stale, 'mdx1')
    assert json.loads(path.read_text())['targets'] == {
        'lab': {'status': 'failed'}, 'mdx1': {'status': 'passed'}}


def test_identity_mismatch_preserves_evidence(tmp_path):
    path = tmp_path / 'state.json'
    original = state(lab={'status': 'passed'})
    module._write_simulink_state(path, original, 'lab')
    changed = {**original, 'commit': 'b'}
    with pytest.raises(ValueError, match='identity mismatch'):
        module._write_simulink_state(path, changed, 'lab')
    assert json.loads(path.read_text()) == original


def test_replace_failure_keeps_valid_json_and_releases_lock(tmp_path, monkeypatch):
    path = tmp_path / 'state.json'
    original = state(lab={'status': 'passed'})
    module._write_simulink_state(path, original, 'lab')
    with monkeypatch.context() as patch:
        def fail(*args):
            raise OSError('replace failed')
        patch.setattr(module.os, 'replace', fail)
        with pytest.raises(OSError, match='replace failed'):
            module._write_simulink_state(path, state(lab={'status': 'failed'}), 'lab')
    assert json.loads(path.read_text()) == original
    assert not list(tmp_path.glob('*.tmp'))
    module._write_simulink_state(path, state(lab={'status': 'failed'}), 'lab')
    assert json.loads(path.read_text())['targets']['lab']['status'] == 'failed'


def test_independent_processes_preserve_both_targets(tmp_path):
    code = '''
import importlib.util, pathlib, sys
spec = importlib.util.spec_from_file_location('candidate', sys.argv[1])
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
for index in range(10):
    m._write_simulink_state(pathlib.Path(sys.argv[2]),
        {'schema': 'test', 'commit': 'a', 'package_sha256': 'hash',
         'targets': {sys.argv[3]: {'iteration': index}}}, sys.argv[3])
'''
    path = tmp_path / 'state.json'
    workers = []
    try:
        for host in ('mdx1', 'mdx2'):
            workers.append(subprocess.Popen([sys.executable, '-c', code, str(TOOL), str(path), host]))
        for worker in workers:
            assert worker.wait(timeout=60) == 0
    finally:
        for worker in workers:
            if worker.poll() is None:
                worker.kill()
                worker.wait()
    assert json.loads(path.read_text())['targets'] == {
        'mdx1': {'iteration': 9}, 'mdx2': {'iteration': 9}}
