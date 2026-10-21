"""Artifact provenance rejects stale, missing, modified and extra wheels."""
import importlib.util
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('gpu_artifact', ROOT / 'tools/ci_gpu_artifact.py')
gpu = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gpu)


def test_exact_wheel_hash_and_source(tmp_path):
    wheel = tmp_path / 'package.whl'
    wheel.write_bytes(b'wheel')
    sha = 'a' * 40
    gpu.seal(tmp_path, sha)
    gpu.verify(tmp_path, sha)
    with pytest.raises(ValueError, match='source mismatch'):
        gpu.verify(tmp_path, 'b' * 40)
    wheel.write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='hash mismatch'):
        gpu.verify(tmp_path, sha)
    wheel.unlink()
    with pytest.raises(ValueError, match='inventory mismatch'):
        gpu.verify(tmp_path, sha)


def test_absent_artifact_fails(tmp_path, monkeypatch):
    import json
    monkeypatch.setenv('GITHUB_REPOSITORY', 'test/repo')
    monkeypatch.setenv('GH_TOKEN', 'fixture')
    monkeypatch.setattr(gpu, 'request', lambda route: json.dumps({'workflow_runs': []}).encode())
    with pytest.raises(ValueError, match='Missing successful exact-SHA artifact'):
        gpu.fetch(tmp_path, 'a' * 40)


def test_redirect_strips_credential_without_network():
    req = gpu.urllib.request.Request('https://api.github.com/repos/test/repo/actions/artifacts/1/zip',
        headers={'Authorization': 'Bearer secret', 'Accept': 'application/zip'})
    redirected = gpu.SafeRedirect().redirect_request(req, None, 302, 'Found', {},
        'https://storage.example/signed')
    assert redirected.get_header('Authorization') is None
    assert redirected.full_url == 'https://storage.example/signed'


def test_schedule_chooses_newest_complete_successful_main_sha(monkeypatch):
    import json
    old, newest, incomplete = 'a' * 40, 'b' * 40, 'c' * 40
    artifacts = [dict(name=name, expired=False, workflow_run=dict(head_branch='main', head_sha=sha))
        for sha in (old, newest) for name in ('gpu-simulink-wheel', 'gpu-optuna-wheel')]
    artifacts.append(dict(name='gpu-simulink-wheel', expired=False,
        workflow_run=dict(head_branch='main', head_sha=incomplete)))
    routes = []
    def request(route):
        routes.append(route)
        if route.startswith('actions/artifacts?'):
            return json.dumps({'artifacts': artifacts if route.endswith('&page=1') else []}).encode()
        return json.dumps([{'sha': sha} for sha in (incomplete, newest, old)]).encode()
    monkeypatch.setattr(gpu, 'request', request)
    monkeypatch.setattr(gpu, 'resolve', lambda sha: [] if sha == newest else pytest.fail('wrong source'))
    assert gpu.newest_main_sha() == newest
    assert any('page=2' in route for route in routes)


def test_schedule_fails_without_complete_pair(monkeypatch):
    monkeypatch.setattr(gpu, 'request', lambda route: b'{"artifacts": []}')
    with pytest.raises(ValueError, match='No main SHA'):
        gpu.newest_main_sha()
