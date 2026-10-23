"""Runtime contract failures do not import an incompatible native library."""
import json

import pytest

from radia_mcp.gmsh import runtime


@pytest.mark.parametrize('version', ['4.15.2', '6.0.0.dev1', '5.0.0'])
def test_unsupported_distribution_fails_before_native_import(monkeypatch, version):
    monkeypatch.setattr(runtime.importlib.metadata, 'version', lambda _: version)
    with pytest.raises(RuntimeError, match='Gmsh 5 development'):
        runtime.qualify_runtime()


def test_missing_provenance_is_explicit(monkeypatch):
    monkeypatch.setattr(runtime.importlib.metadata, 'version', lambda _: '5.0.0.dev20261018')
    monkeypatch.delenv('RADIA_GMSH_PROVENANCE', raising=False)
    with pytest.raises(RuntimeError, match='RADIA_GMSH_PROVENANCE'):
        runtime.qualify_runtime()


def test_wrong_version_in_provenance_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(runtime.importlib.metadata, 'version', lambda _: '5.0.0.dev20261018')
    record = tmp_path / 'provenance.json'
    record.write_text(json.dumps({'distribution_version': '5.0.0.dev20261017', 'wheel_sha256': 'a'*64}))
    monkeypatch.setenv('RADIA_GMSH_PROVENANCE', str(record))
    with pytest.raises(RuntimeError, match='provenance/version mismatch'):
        runtime.qualify_runtime()
