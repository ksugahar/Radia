"""Guards: readiness remains manual, bounded and fail-loud; owned Engine closes.

Assumes: existing GPU runner context is verified on dispatch; fake Engine here
checks lifecycle/failure contracts only and cannot certify MATLAB or hardware.
"""
import importlib.util
import json
from pathlib import Path
import sys
import types

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def test_readiness_is_manual_only_with_fixed_existing_targets():
    data = yaml.load((ROOT / '.github/workflows/gpu-runner-readiness.yml').read_text(), Loader=yaml.BaseLoader)
    assert set(data['on']) == {'workflow_dispatch'}
    assert data['on']['workflow_dispatch']['inputs']['target']['options'] == ['gpu-a100', 'gpu-rtx']
    assert (ROOT / 'packages/radia-simulink/tools/gpu_readiness/probe_matlab_readiness.m').is_file()
    assert data['permissions'] == {'contents': 'read'}
    for job, target in [('a100', 'gpu-a100'), ('rtx', 'gpu-rtx')]:
        spec = data['jobs'][job]
        assert target in spec['runs-on']
        assert "github.event_name == 'workflow_dispatch'" in spec['if']
        assert spec['timeout-minutes'] == '10'
        probe = next(s for s in spec['steps'] if s.get('name', '').startswith('Probe'))
        assert probe['timeout-minutes'] == '5'
        assert 'if ($LASTEXITCODE -ne 0)' in probe['run']
        upload = next(s for s in spec['steps'] if s.get('uses') == 'actions/upload-artifact@v4')
        assert upload['if'] == 'always()'
        assert upload['with']['path'].splitlines() == ['${{ env.PROBE_EVIDENCE }}/dispatch.json', '${{ env.PROBE_EVIDENCE }}/summary.json']


@pytest.mark.parametrize('failure', [None, 'start', 'probe', 'shutdown'])
def test_engine_owns_session_and_never_reports_failed_lifecycle_as_pass(monkeypatch, tmp_path, failure):
    spec = importlib.util.spec_from_file_location('readiness_engine', ROOT / 'packages/radia-simulink/tools/gpu_readiness/probe_engine.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls = []

    class Owned:
        def addpath(self, *args, **kwargs):
            calls.append('addpath')

        def probe_matlab_readiness(self, *args, **kwargs):
            calls.append('probe')
            if failure == 'probe':
                raise RuntimeError('missing license')
            return json.dumps({'passed': True, 'gpu': {'uuid': 'test-device'}})

        def quit(self):
            calls.append('quit')
            if failure == 'shutdown':
                raise RuntimeError('shutdown failed')

    def start(options):
        calls.append('start')
        if failure == 'start':
            raise RuntimeError('startup failed')
        return Owned()

    matlab = types.ModuleType('matlab')
    engine = types.ModuleType('matlab.engine')
    engine.__file__ = 'fake-engine'
    engine.start_matlab = start
    matlab.engine = engine
    monkeypatch.setitem(sys.modules, 'matlab', matlab)
    monkeypatch.setitem(sys.modules, 'matlab.engine', engine)
    monkeypatch.setattr(module.importlib.metadata, 'version', lambda name: '26.1')
    out = tmp_path / 'engine.json'
    monkeypatch.setattr(sys, 'argv', ['probe', '--matlab-root', str(tmp_path), '--gpu', 'test', '--probe-dir', str(tmp_path), '--result', str(out)])
    assert module.main() == (0 if failure is None else 1)
    record = json.loads(out.read_text())
    assert record['passed'] is (failure is None)
    assert calls.count('start') == 1
    assert calls.count('quit') == (0 if failure == 'start' else 1)
    if failure is None:
        assert record['shutdown'] == 'complete'
