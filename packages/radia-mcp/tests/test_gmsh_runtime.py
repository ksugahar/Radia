"""Runtime contract failures do not import an incompatible native library."""
import json

import pytest

from radia_mcp.gmsh import runtime


@pytest.mark.parametrize('version', ['4.15.2', '6.0.0.dev1'])
def test_unsupported_distribution_fails_before_native_import(monkeypatch, version):
    monkeypatch.setattr(runtime.importlib.metadata, 'version', lambda _: version)
    with pytest.raises(RuntimeError, match='Gmsh 5 required'):
        runtime.qualify_runtime()


@pytest.mark.parametrize('version', ['5.0.0.dev20261018', '5.0.0'])
def test_missing_provenance_is_explicit(monkeypatch, version):
    monkeypatch.setattr(runtime.importlib.metadata, 'version', lambda _: version)
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


@pytest.mark.parametrize('kind', ['NodeData', 'ElementData', 'ElementNodeData'])
@pytest.mark.parametrize('tag', [2**28, 2147483665])
def test_native_field_tag_limit_before_access(tmp_path, kind, tag):
    from types import SimpleNamespace
    path = tmp_path / 'field.msh'
    block = '15 1 1' if kind == 'ElementNodeData' else '1'
    path.write_text(f'$MeshFormat\n4.2 0 8\n$EndMeshFormat\n${kind}\n'
                    f'1\n"field"\n1\n0\n3\n0\n1\n1\n1 1 {tag} {tag}\n'
                    f'{block}\n{tag}\n7\n$End{kind}\n')
    calls = []
    native = SimpleNamespace(open=lambda p: calls.append(p), merge=lambda p: None,
                             write=lambda p: None, view=SimpleNamespace(write=lambda *a: None))
    runtime.install_io_guards(native)
    if tag > 2**31 - 1:
        with pytest.raises(ValueError, match='2147483665.*exceeds 2147483647'):
            native.open(path)
        assert calls == []
    else:
        native.open(path)
        assert calls == [path]
