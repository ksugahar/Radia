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


def test_package_keeps_solver_runtime_dependencies_optional_and_lazy():
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


@pytest.mark.parametrize('component,field', [('hmatrix','fault_injection'),('projection','fault_injection'),('projection','lifecycle_diagnostic')])
@pytest.mark.parametrize('mode', [False, True, None, 0, 'False'])
def test_hdiv_mex_wheel_requires_explicit_false_diagnostic_mode(component, field, mode):
    spec=importlib.util.spec_from_file_location('hdiv_mex_receipt_verifier',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    manifest=dict(source_commit='a'*40,files={
        'matlab/radia_cuda_hmatrix.mexw64':'b'*64,
        'matlab/radia_cuda_projection.mexw64':'c'*64})
    artifacts={name:dict(source_dirty=False,source_commit='a'*40,
        taskmanager_started=False,fault_injection=False,
        binary_sha256=manifest['files']['matlab/radia_cuda_'+name+'.mexw64'])
        for name in ('hmatrix','projection')}
    artifacts['projection']['lifecycle_diagnostic']=False
    if mode is None:
        artifacts[component].pop(field)
    else:
        artifacts[component][field]=mode
    receipt=dict(schema='radia.gpu.hdiv-payload.v1',artifacts=artifacts)
    if mode is False:
        verifier.verify_hdiv_receipt(manifest,receipt)
    else:
        with pytest.raises(ValueError,match='HDiv CUDA binary receipt mismatch'):
            verifier.verify_hdiv_receipt(manifest,receipt)





@pytest.mark.parametrize('wheel_tag,native_tag,accepted',[
    ('cp312-cp312','cp312',True),('py3-none','cp312',False),
    ('cp313-cp313','cp312',False),('cp312-cp311','cp312',False)])
def test_python_cuda_wheel_requires_exact_native_abi(tmp_path,wheel_tag,native_tag,accepted):
    import zipfile
    spec=importlib.util.spec_from_file_location('wheel_native_contract',ROOT/'verify_wheel.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    name=f'_csr_cuda.{native_tag}-win_amd64.pyd'
    wheel=tmp_path/f'radia_simulink-0.1.5rc1-{wheel_tag}-win_amd64.whl'
    digest=hashlib.sha256(b'packaging-only synthetic binary').hexdigest()
    manifest=dict(source_commit='a'*40,files={name:digest})
    receipt=dict(schema='radia.csr-python-build.v1',source_commit='a'*40,binary=name,
                 binary_sha256=digest,matlab_engine_user_dependency=False)
    with zipfile.ZipFile(wheel,'w') as z:
        z.writestr('radia_simulink/'+name,b'packaging-only synthetic binary')
        z.writestr('radia_simulink/python-build.json',json.dumps(receipt))
    with zipfile.ZipFile(wheel) as z:
        if accepted:module.verify_python_binding(wheel,z.namelist(),z,manifest)
        else:
            with pytest.raises(ValueError,match='ABI'):
                module.verify_python_binding(wheel,z.namelist(),z,manifest)


def test_native_payload_uses_platform_library_root(monkeypatch):
    import runpy
    import setuptools
    captured = {}
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: captured.update(kwargs))
    runpy.run_path(str(ROOT / "setup.py"))
    distribution = captured["distclass"]({"packages": []})
    command = distribution.get_command_obj("build")
    command.ensure_finalized()
    assert command.build_lib == command.build_platlib
    assert command.build_lib != command.build_purelib

@pytest.mark.parametrize('mutation', ['none', 'abi', 'runtime', 'redistribution', 'source', 'diagnostic', 'device-diagnostic', 'device-mode-missing'])
def test_csr_and_hdiv_wheel_share_abi_and_require_hdiv_qualification(tmp_path, mutation):
    import zipfile
    spec=importlib.util.spec_from_file_location('wheel_dual_native_contract', ROOT/'verify_wheel.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    wheel=tmp_path/'radia_simulink-0.1.5rc1-cp312-cp312-win_amd64.whl'
    manifest=dict(source_commit='a'*40, files={})
    with zipfile.ZipFile(wheel, 'w') as archive:
        for component in ('csr', 'hdiv'):
            abi='cp313' if component=='hdiv' and mutation=='abi' else 'cp312'
            name=f'_{component}_cuda.{abi}-win_amd64.pyd'
            payload=b'packaging-only synthetic '+component.encode()
            digest=hashlib.sha256(payload).hexdigest(); manifest['files'][name]=digest
            receipt=dict(schema=f'radia.{component}-python-build.v1', source_commit='a'*40,
                         binary=name, binary_sha256=digest, matlab_engine_user_dependency=False)
            if component=='hdiv':
                receipt.update(fault_injection=mutation=='diagnostic',device_fault_injection=mutation=='device-diagnostic',runtime_import_checked=mutation!='runtime',
                               redistribution_accepted=mutation!='redistribution',
                               redistribution_reference_sha256='b'*64)
                if mutation=='source': receipt['source_commit']='c'*40
                if mutation=='device-mode-missing': receipt.pop('device_fault_injection')
            archive.writestr('radia_simulink/'+name, payload)
            record='python-build.json' if component=='csr' else 'hdiv-python-build.json'
            archive.writestr('radia_simulink/'+record, json.dumps(receipt))
    with zipfile.ZipFile(wheel) as archive:
        if mutation=='none': module.verify_python_binding(wheel,archive.namelist(),archive,manifest)
        else:
            with pytest.raises(ValueError):
                module.verify_python_binding(wheel,archive.namelist(),archive,manifest)
