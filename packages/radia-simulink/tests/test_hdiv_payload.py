"""Release staging must preserve exact native identity and reject stale builds."""
from pathlib import Path
import hashlib
import importlib.util
import json
import pytest

spec = importlib.util.spec_from_file_location('hdiv_payload', Path(__file__).parents[1] / 'hdiv_payload.py')
payload = importlib.util.module_from_spec(spec)
spec.loader.exec_module(payload)


def test_hdiv_only_python_wheel_uses_the_cpython_abi(tmp_path, monkeypatch):
    import ast
    import os
    import types
    source = Path(__file__).parents[1] / 'setup.py'
    tree = ast.parse(source.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'bdist_wheel')
    namespace = dict(_bdist_wheel=object, ROOT=tmp_path, json=json, os=os,
                     sys=types.SimpleNamespace(platform='win32', version_info=(3,12)),
                     payload_module=lambda:payload)
    exec(compile(ast.Module(body=[cls],type_ignores=[]),str(source),'exec'),namespace)
    monkeypatch.setenv('RADIA_HDIV_PYTHON_RUNTIME',str(tmp_path/'runtime.json'))
    instance = namespace['bdist_wheel']()
    assert instance.get_tag() == ('cp312','cp312','win_amd64')
    monkeypatch.delenv('RADIA_HDIV_PYTHON_RUNTIME')
    assert instance.get_tag() == ('py3','none','win_amd64')


def test_runtime_requirements_are_interface_specific():
    common, interfaces = payload.runtime_contract(True)
    assert common == ['Windows x64']
    assert not any('MATLAB' in item or 'Toolbox' in item for item in interfaces['python_gpu'])
    assert 'MATLAB R2026a' in interfaces['matlab_gpu']
    assert 'Simulink' in interfaces['simulink_gpu']
    assert 'python_gpu' not in payload.runtime_contract(False)[1]


def test_native_notices_are_preserved_in_staged_resources(tmp_path):
    source = Path(__file__).parents[1]
    target = tmp_path/'package'
    assert payload.stage_documentation(source,target) == set(payload.DOCUMENTATION_FILES)
    for name in payload.DOCUMENTATION_FILES:
        assert (source/name).read_bytes() == (target/name).read_bytes()


@pytest.mark.parametrize('fault', [None, 'commit', 'binary', 'taskmanager'])
def test_exact_fmm_artifact_staging(tmp_path, fault):
    root=tmp_path/'source';root.mkdir()
    (root/'native.cu').write_bytes(b'generic public FMM')
    folder=tmp_path/'build';folder.mkdir()
    binary=folder/'radia_cuda_fmm.mexw64';binary.write_bytes(b'fixture artifact')
    receipt=dict(schema='radia.fmm-cuda-build.v1',accepted=True,source_commit='release',
                 source_dirty=False,taskmanager_used=False,python_adapter_required=False,
                 source_sha256={'native.cu':hashlib.sha256((root/'native.cu').read_bytes()).hexdigest()},
                 binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),matlab='fixture')
    if fault=='commit':receipt['source_commit']='stale'
    if fault=='binary':binary.write_bytes(b'changed')
    if fault=='taskmanager':receipt['taskmanager_used']=True
    (folder/'build.json').write_text(json.dumps(receipt))
    package=tmp_path/'package'
    if fault:
        with pytest.raises(RuntimeError):payload.stage_fmm(root,package,'release',folder)
    else:
        assert payload.stage_fmm(root,package,'release',folder)=={'fmm-build.json','matlab/radia_cuda_fmm.mexw64'}


@pytest.mark.parametrize('fault', [None, 'commit', 'dirty', 'binary', 'source', 'taskmanager', 'diagnostic', 'unknown-fault-mode', 'projection-diagnostic', 'projection-mode-missing', 'projection-fault-missing'])
def test_exact_hdiv_artifact_staging(tmp_path, fault):
    root = tmp_path / 'source'
    root.mkdir()
    source = root / 'native.cu'
    source.write_bytes(b'public generic CUDA source')
    wrappers = root / 'packages/radia-simulink/matlab/+radia/+gpu'
    wrappers.mkdir(parents=True)
    for name in payload.WRAPPERS:
        (wrappers / name).write_text('% wrapper')
    for name in payload.SIMULINK_WRAPPERS:
        path = root / 'packages/radia-simulink' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('% selected shared GPU workflow')
    folders = []
    for name in ('hmatrix', 'projection'):
        folder = tmp_path / name
        folder.mkdir()
        binary = folder / ('radia_cuda_' + name + '.mexw64')
        binary.write_bytes(name.encode())
        receipt = dict(accepted=True, source_dirty=False, source_commit='a'*40,
                       taskmanager_started=False, fault_injection=False, matlab='test fixture',
                       source_sha256={'native.cu': hashlib.sha256(source.read_bytes()).hexdigest()},
                       binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
        if name == 'hmatrix' and fault == 'unknown-fault-mode':receipt.pop('fault_injection')
        if name == 'projection':
            receipt['lifecycle_diagnostic']=False
            if fault == 'projection-fault-missing':receipt.pop('fault_injection')
            if fault == 'projection-diagnostic':receipt['lifecycle_diagnostic']=True
            if fault == 'projection-mode-missing':receipt.pop('lifecycle_diagnostic')
            if fault == 'commit': receipt['source_commit'] = 'b'*40
            if fault == 'dirty': receipt['source_dirty'] = True
            if fault == 'binary': receipt['binary_sha256'] = '0'*64
            if fault == 'source': receipt['source_sha256']['native.cu'] = '0'*64
            if fault == 'taskmanager': receipt['taskmanager_started'] = True
            if fault == 'diagnostic': receipt['fault_injection'] = True
        (folder / 'build.json').write_text(json.dumps(receipt))
        folders.append(folder)
    target = tmp_path / 'package'
    if fault:
        with pytest.raises(RuntimeError): payload.stage(root, target, 'a'*40, *folders)
    else:
        files = payload.stage(root, target, 'a'*40, *folders)
        assert files == {'hdiv-build.json','matlab/radia_cuda_hmatrix.mexw64',
                         'matlab/radia_cuda_projection.mexw64'} | {
            'matlab/+radia/+gpu/'+name for name in payload.WRAPPERS} | set(payload.SIMULINK_WRAPPERS)
        for name in ('HMatrixOperator.m','HDivBHMaterial.m','PreparedHDivNonlinear.m'):
            assert 'matlab/+radia/+gpu/'+name in files
        for folder, name in zip(folders, ('hmatrix', 'projection')):
            binary = 'radia_cuda_' + name + '.mexw64'
            assert (target / 'matlab' / binary).read_bytes() == (folder / binary).read_bytes()


@pytest.mark.parametrize('fault',[None,'commit','diagnostic','python-diagnostic','python-mode-missing','device-diagnostic','device-mode-missing','binary','abi','source','runtime','runtime-abi','redistribution','engine','numerics','nan','boolean'])
def test_hdiv_python_requires_exact_abi_runtime_and_redistribution(tmp_path,monkeypatch,fault):
    import sysconfig
    monkeypatch.setattr(sysconfig,'get_config_var',lambda key:'.cp312-win_amd64.pyd')
    root=tmp_path/'source';root.mkdir();(root/'native.cu').write_bytes(b'public test source')
    folder=tmp_path/'build';folder.mkdir();name='_hdiv_cuda.cp312-win_amd64.pyd'
    native=folder/name;native.write_bytes(b'packaging-only synthetic native file');digest=hashlib.sha256(native.read_bytes()).hexdigest()
    binding=dict(fault_injection=False,device_fault_injection=False,schema='radia.hdiv-python-build.v1',binary=name,binary_sha256=digest,extension_suffix='.cp312-win_amd64.pyd',ngsolve_version='6.2.2608',pybind11_version='3.1.0',pybind_header_sha256='1'*64,matlab_engine_user_dependency=False)
    build=dict(accepted=True,source_dirty=False,source_commit='a'*40,fault_injection=False,python_binding=binding,source_sha256={'native.cu':hashlib.sha256((root/'native.cu').read_bytes()).hexdigest()})
    runtime=dict(schema='radia.hdiv-python-runtime.v1',accepted=True,source_commit='a'*40,binary_sha256=digest,ngsolve_version='6.2.2608',runtime_import_checked=True,matlab_engine_user_dependency=False,redistribution_accepted=True,redistribution_reference='synthetic fixture only; not a licensing decision',runtime_identity_sha256='2'*64,cases=[dict(backend='gpu',coefficient_error=1e-13,independent_original_residual=1e-14)])
    if fault=='commit':build['source_commit']='b'*40
    if fault=='diagnostic':build['fault_injection']=True
    if fault=='python-diagnostic':binding['fault_injection']=True
    if fault=='python-mode-missing':binding.pop('fault_injection')
    if fault=='device-diagnostic':binding['device_fault_injection']=True
    if fault=='device-mode-missing':binding.pop('device_fault_injection')
    if fault=='binary':native.write_bytes(b'changed')
    if fault=='abi':binding['extension_suffix']='.cp311-win_amd64.pyd'
    if fault=='source':(root/'native.cu').write_bytes(b'changed')
    if fault=='runtime':runtime['accepted']=False
    if fault=='runtime-abi':runtime['ngsolve_version']='stale'
    if fault=='redistribution':runtime['redistribution_accepted']=False
    if fault=='engine':runtime['matlab_engine_user_dependency']=True
    if fault=='numerics':runtime['cases'][0]['coefficient_error']=1e-4
    if fault=='nan':runtime['cases'][0]['coefficient_error']=float('nan')
    if fault=='boolean':runtime['cases'][0]['coefficient_error']=False
    (folder/'build.json').write_text(json.dumps(build));proof=tmp_path/'runtime.json';proof.write_text(json.dumps(runtime));package=tmp_path/'package'
    if fault:
        with pytest.raises(RuntimeError):payload.stage_python(root,package,'a'*40,folder,proof)
        assert not (package/name).exists()
    else:
        assert payload.stage_python(root,package,'a'*40,folder,proof)=={name,'hdiv-python-build.json'}
        receipt=json.loads((package/'hdiv-python-build.json').read_text())
        assert receipt['runtime_import_checked'] and receipt['redistribution_accepted']
        assert 'redistribution_reference' not in receipt
        assert len(receipt['redistribution_reference_sha256'])==64
