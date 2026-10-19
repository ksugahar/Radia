"""Guard the installed resource contract without requiring a MATLAB license."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def load_api():
    spec=importlib.util.spec_from_file_location("gpu_contract_api",ROOT/"src/radia_simulink/__init__.py")
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
    checked = api.verify()
    assert checked["passed"]
    assert checked['manifest_sha256'] == hashlib.sha256((tmp_path/'manifest.json').read_bytes()).hexdigest()
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
    extras=metadata["project"]["optional-dependencies"]
    assert set(extras["exchange"]).issubset(extras["test"])
    api=load_api()
    assert "ngsolve" not in api.__dict__ and "matlab" not in api.__dict__
    assert "checkout" not in api.matlab_path.__doc__.lower()


@pytest.mark.parametrize('member',['../escape','radia_simulink/matlab/radia_mex.mexw64',
                                  'radia_simulink/matlab/external.dll','radia/__init__.py'])
def test_wheel_rejects_escape_cpu_solver_and_external_payload(tmp_path,member):
    import zipfile
    spec=importlib.util.spec_from_file_location('gpu_wheel_contract',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel=tmp_path/'radia_simulink-0.1.0a1-py3-none-win_amd64.whl'
    with zipfile.ZipFile(wheel,'w') as archive:
        archive.writestr(member,b'forbidden')
    with pytest.raises(ValueError):
        verifier.verify(wheel)


def test_installed_acceptance_rejects_a_different_build_of_the_same_source():
    spec=importlib.util.spec_from_file_location('gpu_installed_identity',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel = dict(source_commit='a'*40, version='0.1.0', manifest_sha256='b'*64)
    verifier.verify_installed_identity(dict(wheel), wheel)
    installed = dict(wheel, manifest_sha256='c'*64)
    with pytest.raises(ValueError, match='differs from release wheel'):
        verifier.verify_installed_identity(installed, wheel)
def test_runtime_and_package_versions_agree():
    import ast
    import tomllib
    metadata = tomllib.loads((ROOT/'pyproject.toml').read_text())
    module = ast.parse((ROOT/'src/radia_simulink/__init__.py').read_text())
    assignment = next(node for node in module.body if isinstance(node, ast.Assign)
                      and any(isinstance(target, ast.Name) and target.id=='__version__' for target in node.targets))
    assert ast.literal_eval(assignment.value) == metadata['project']['version']

