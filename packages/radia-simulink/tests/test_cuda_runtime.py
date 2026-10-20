"""GPU dependency admission and path ownership without loading CUDA."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest


@pytest.fixture
def runtime(monkeypatch):
    path=Path(__file__).parents[1]/'src/radia_simulink/_cuda_runtime.py'
    spec=importlib.util.spec_from_file_location('runtime_contract',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    monkeypatch.setattr(module.sys,'platform','win32')
    monkeypatch.setattr(module,'_loaded_path',lambda _:None)
    return module


def distributions(runtime,monkeypatch,tmp_path):
    result={}
    for name,version in runtime._VERSIONS.items():
        root=tmp_path/name;root.mkdir();(root/'runtime.dll').write_bytes(b'fixture')
        result[name]=SimpleNamespace(version=version,files=[Path('runtime.dll')],locate_file=lambda p,root=root:root/p)
    monkeypatch.setattr(runtime.metadata,'distribution',result.__getitem__)
    return result


def test_exact_runtime_registers_once_without_loading_cuda(runtime,monkeypatch,tmp_path):
    distributions(runtime,monkeypatch,tmp_path);calls=[]
    monkeypatch.setattr(runtime.os,'add_dll_directory',lambda p:calls.append(p) or SimpleNamespace(close=lambda:None))
    assert calls==[]
    runtime.prepare_cuda_runtime();runtime.prepare_cuda_runtime()
    assert len(calls)==len(runtime._VERSIONS)
    assert all(Path(p).is_relative_to(tmp_path) for p in calls)


def test_different_runtime_fails_before_registering_paths(runtime,monkeypatch,tmp_path):
    installed=distributions(runtime,monkeypatch,tmp_path)
    installed['nvidia-cublas-cu12'].version='12.9.0'
    monkeypatch.setattr(runtime.os,'add_dll_directory',lambda _:pytest.fail('unqualified runtime'))
    with pytest.raises(ImportError,match='requires'):runtime.prepare_cuda_runtime()
    assert runtime._handles==[]


def test_partial_directory_failure_releases_its_handles(runtime,monkeypatch,tmp_path):
    distributions(runtime,monkeypatch,tmp_path);closed=[];calls=[]
    def add(path):
        calls.append(path)
        if len(calls)==2:raise OSError('injected directory failure')
        return SimpleNamespace(close=lambda:closed.append(path))
    monkeypatch.setattr(runtime.os,'add_dll_directory',add)
    with pytest.raises(OSError,match='injected'):runtime.prepare_cuda_runtime()
    assert closed==calls[:1] and runtime._handles==[]


def test_existing_foreign_dll_fails_before_path_registration(runtime,monkeypatch,tmp_path):
    distributions(runtime,monkeypatch,tmp_path)
    monkeypatch.setattr(runtime,'_loaded_path',lambda _:tmp_path/'other-provider/runtime.dll')
    monkeypatch.setattr(runtime.os,'add_dll_directory',lambda _:pytest.fail('foreign provider admitted'))
    with pytest.raises(ImportError,match='already loaded'):runtime.prepare_cuda_runtime()
    assert runtime._handles==[]
