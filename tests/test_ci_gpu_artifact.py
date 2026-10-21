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
    class Response:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def read(self): return json.dumps({'workflow_runs': []}).encode()
    monkeypatch.setattr(gpu.urllib.request, 'urlopen', lambda *a, **kw: Response())
    with pytest.raises(ValueError, match='Missing successful exact-SHA artifact'):
        gpu.fetch(tmp_path, 'a' * 40)
