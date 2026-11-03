"""Release staging must preserve exact native identity and reject stale builds."""
from pathlib import Path
import hashlib
import importlib.util
import json
import pytest

spec = importlib.util.spec_from_file_location('hdiv_payload', Path(__file__).parents[1] / 'hdiv_payload.py')
payload = importlib.util.module_from_spec(spec)
spec.loader.exec_module(payload)


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


@pytest.mark.parametrize('fault', [None, 'commit', 'dirty', 'binary', 'source', 'mexcuda', 'schema', 'receipt-missing'])
def test_exact_hdiv_device_mex_staging(tmp_path, fault):
    """Guards: only a clean exact-commit plain MEX enters the wheel, with its MATLAB wrappers.

    Assumes: the receipt fields are those build_hdiv_mex.py writes.
    """
    root = tmp_path / 'source'
    wrappers = root / 'packages/radia-simulink/matlab/+radia/+gpu'
    wrappers.mkdir(parents=True)
    for name in payload.DEVICE_WRAPPERS:
        (wrappers / name).write_text('% device wrapper')
    for name in payload.SIMULINK_WRAPPERS:
        path = root / 'packages/radia-simulink' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('% selected shared GPU workflow')
    source = root / 'packages/radia-simulink/native.cpp'
    source.write_bytes(b'public plain MEX gateway')
    folder = tmp_path / 'build'
    folder.mkdir()
    binary = folder / payload.DEVICE_BINARY
    binary.write_bytes(b'fixture artifact')
    receipt = dict(schema='radia.hdiv-mex-build.v1', source_commit='a'*40, source_dirty=False, mexcuda=False,
                   source_sha256={'packages/radia-simulink/native.cpp': hashlib.sha256(source.read_bytes()).hexdigest()},
                   binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                   cuda_object=dict(name='hdiv.obj', sha256='c'*64))
    if fault == 'commit': receipt['source_commit'] = 'b'*40
    if fault == 'dirty': receipt['source_dirty'] = True
    if fault == 'binary': receipt['binary_sha256'] = '0'*64
    if fault == 'source': source.write_bytes(b'changed')
    if fault == 'mexcuda': receipt['mexcuda'] = True
    if fault == 'schema': receipt['schema'] = 'unknown'
    if fault != 'receipt-missing':
        Path(str(binary) + '.build.json').write_text(json.dumps(receipt))
    target = tmp_path / 'package'
    if fault:
        with pytest.raises((RuntimeError, FileNotFoundError)):
            payload.stage_device(root, target, 'a'*40, folder)
        assert not (target / 'hdiv-device-build.json').exists()
        return
    files = payload.stage_device(root, target, 'a'*40, folder)
    assert files == {'matlab/' + payload.DEVICE_BINARY, 'hdiv-device-build.json'} | {
        'matlab/+radia/+gpu/' + name for name in payload.DEVICE_WRAPPERS} | set(payload.SIMULINK_WRAPPERS), files
    assert (target / 'matlab' / payload.DEVICE_BINARY).read_bytes() == binary.read_bytes()
    public = json.loads((target / 'hdiv-device-build.json').read_text())
    assert public['binary_sha256'] == receipt['binary_sha256'] and public['mexcuda'] is False, public


def test_hdiv_python_payload_is_owned_by_radia():
    """Guards: MEX wheels cannot stage the moved Radia Python domain binding.

    Assumes: Radia owns public Python GPU provenance and runtime checks.
    """
    assert not hasattr(payload, "stage_python")
