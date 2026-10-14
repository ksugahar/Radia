"""Projection headers/libs follow the selected interpreter, not system Netgen."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest


def test_projection_dependency_selection_and_identity(tmp_path,monkeypatch):
    path=Path(__file__).resolve().parents[1]/'tools/build_projection_cuda.py'
    spec=importlib.util.spec_from_file_location('projection_build',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    base=tmp_path/'python';net=tmp_path/'isolated/site-packages/netgen'
    for relative in ('include/Python.h','libs/python312.lib'):
        p=base/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'python')
    for relative in ('__init__.py','include/nginterface.h','lib/libngsolve.lib','lib/ngcore.lib','lib/nglib.lib'):
        p=net/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'isolated')
    monkeypatch.setattr(module.sys,'base_prefix',str(base))
    monkeypatch.setattr(module.importlib.util,'find_spec',lambda name:SimpleNamespace(origin=str(net/'__init__.py')))
    monkeypatch.setattr(module.importlib.metadata,'version',lambda name:'test-version')
    prefix,selected,identity=module.projection_dependencies()
    assert prefix==base and selected==net
    assert str(net/'lib/libngsolve.lib') in identity['native_dependency_sha256']
    assert str(base/'include/Python.h') in identity['native_dependency_sha256']
    (net/'lib/libngsolve.lib').write_bytes(b'changed')
    assert module.projection_dependencies()[2]!=identity
    (net/'lib/nglib.lib').unlink()
    with pytest.raises(FileNotFoundError):module.projection_dependencies()
