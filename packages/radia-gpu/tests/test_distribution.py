"""Guard the installed resource contract without requiring a MATLAB license."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def load_api():
    spec=importlib.util.spec_from_file_location("gpu_contract_api",ROOT/"src/radia_gpu/__init__.py")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_installed_manifest_rejects_tampering_and_path_escape(tmp_path,monkeypatch):
    api=load_api()
    monkeypatch.setattr(api,"files",lambda package:tmp_path)
    (tmp_path/"matlab").mkdir()
    payload=tmp_path/"matlab/kernel.m"
    payload.write_text("checked payload")
    manifest=dict(schema="radia.gpu.distribution.v1",version=api.__version__,source_commit="a"*40,
                  capabilities=["test payload"],runtime_requirements=[],
                  files={"matlab/kernel.m":hashlib.sha256(payload.read_bytes()).hexdigest()})
    (tmp_path/"manifest.json").write_text(json.dumps(manifest))
    assert api.verify()["passed"]
    payload.write_text("changed payload")
    with pytest.raises(RuntimeError,match="hash mismatch"):
        api.verify()
    manifest["files"]={"../escape":"0"*64}
    (tmp_path/"manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError,match="invalid payload path"):
        api.verify()


def test_package_exposes_only_gpu_solver_and_no_cpu_runtime_dependency():
    import tomllib
    metadata=tomllib.loads((ROOT/"pyproject.toml").read_text())
    assert metadata["project"]["dependencies"]==[]
    api=load_api()
    assert "ngsolve" not in api.__dict__ and "matlab" not in api.__dict__
    assert "checkout" not in api.matlab_path.__doc__.lower()


@pytest.mark.parametrize('member',['../escape','radia_gpu/matlab/radia_mex.mexw64',
                                  'radia_gpu/matlab/external.dll','radia/__init__.py'])
def test_wheel_rejects_escape_cpu_solver_and_external_payload(tmp_path,member):
    import zipfile
    spec=importlib.util.spec_from_file_location('gpu_wheel_contract',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel=tmp_path/'radia_gpu-0.1.0a1-py3-none-win_amd64.whl'
    with zipfile.ZipFile(wheel,'w') as archive:
        archive.writestr(member,b'forbidden')
    with pytest.raises(ValueError):
        verifier.verify(wheel)
