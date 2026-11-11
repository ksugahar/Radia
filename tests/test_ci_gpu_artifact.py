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


def test_netgen_library_is_sealed_like_a_wheel(tmp_path):
    library = tmp_path / 'nglib.dll'
    library.write_bytes(b'library')
    (tmp_path / 'netgen-overlay.json').write_text('{}')
    sha = 'a' * 40
    gpu.seal(tmp_path, sha)
    gpu.verify(tmp_path, sha)
    library.write_bytes(b'stock')
    with pytest.raises(ValueError, match='hash mismatch'):
        gpu.verify(tmp_path, sha)


def test_absent_artifact_fails(tmp_path, monkeypatch):
    import json
    monkeypatch.setenv('GITHUB_REPOSITORY', 'test/repo')
    monkeypatch.setenv('GH_TOKEN', 'fixture')
    monkeypatch.setattr(gpu, 'request', lambda route: json.dumps({'workflow_runs': []}).encode())
    with pytest.raises(ValueError, match='Missing successful exact-SHA artifact'):
        gpu.fetch(tmp_path, 'a' * 40)


def test_radia_wheel_comes_from_native_release_at_the_exact_sha(monkeypatch):
    import json
    sha = 'a' * 40
    routes = []
    def request(route):
        routes.append(route)
        if route.startswith('actions/workflows/'):
            return json.dumps({'workflow_runs': [dict(id=8, head_sha='b' * 40), dict(id=7, head_sha=sha)]}).encode()
        return json.dumps({'artifacts': [dict(id=9, name='gpu-radia-wheel', expired=False),
                                         dict(id=10, name='gpu-netgen-overlay', expired=False)]}).encode()
    monkeypatch.setattr(gpu, 'request', request)
    selected = gpu.resolve(sha, gpu.RADIA_WHEELS)
    assert [(run['id'], found['id'], artifact) for run, found, artifact in selected] == [
        (7, 9, 'gpu-radia-wheel'), (7, 10, 'gpu-netgen-overlay')]
    assert routes[0] == f'actions/workflows/build-test.yml/runs?head_sha={sha}&status=success&per_page=100'
    routes.clear()
    monkeypatch.setattr(gpu, 'request', lambda route: routes.append(route) or json.dumps({'workflow_runs': []}).encode())
    with pytest.raises(ValueError, match='Missing successful exact-SHA artifact: gpu-simulink-wheel'):
        gpu.resolve(sha)
    assert routes[0] == f'actions/workflows/radia-simulink.yml/runs?head_sha={sha}&event=push&status=success&per_page=100'


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
        for sha in (old, newest) for name in ('gpu-simulink-wheel',)]
    artifacts.append(dict(name='unrelated-wheel', expired=False,
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


def test_schedule_fails_without_required_wheel(monkeypatch):
    monkeypatch.setattr(gpu, 'request', lambda route: b'{"artifacts": []}')
    with pytest.raises(ValueError, match='No main SHA'):
        gpu.newest_main_sha()


def test_compute_bundle_round_trip_and_refusals(tmp_path):
    """Guards: the build-host bundle restores the exact trees on the compute host and never escapes or overwrites."""
    root, temp, sha = tmp_path / 'build', tmp_path / 'temp', 'a' * 40
    (root / 'build-x/sub').mkdir(parents=True)
    (temp / 'projection-1').mkdir(parents=True)
    (root / 'build-x/sub/gateway.mexw64').write_bytes(b'mex')
    (temp / 'projection-1/build.json').write_text('{}')
    sealed = root / 'bundle'
    gpu.bundle(sealed, root, temp, ['checkout:build-x', 'temp:projection-1'])
    gpu.seal(sealed, sha)
    for index, entry in enumerate(['checkout:../x', 'checkout:.', 'temp:missing', 'other:x']):
        with pytest.raises(ValueError):
            gpu.bundle(tmp_path / f'bad{index}', root, temp, [entry])
    compute, scratch = tmp_path / 'compute', tmp_path / 'scratch'
    compute.mkdir()
    scratch.mkdir()
    with pytest.raises(ValueError, match='source mismatch'):
        gpu.unbundle(sealed, 'b' * 40, compute, scratch)
    with pytest.raises(ValueError, match='would overwrite'):
        gpu.unbundle(sealed, sha, root, temp)
    gpu.unbundle(sealed, sha, compute, scratch)
    assert (compute / 'build-x/sub/gateway.mexw64').read_bytes() == b'mex'
    assert (scratch / 'projection-1/build.json').read_text() == '{}'
