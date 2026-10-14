"""Release staging must preserve exact native identity and reject stale builds."""
from pathlib import Path
import hashlib
import importlib.util
import json
import pytest

spec = importlib.util.spec_from_file_location('hdiv_payload', Path(__file__).parents[1] / 'hdiv_payload.py')
payload = importlib.util.module_from_spec(spec)
spec.loader.exec_module(payload)


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


@pytest.mark.parametrize('fault', [None, 'commit', 'dirty', 'binary', 'source', 'taskmanager'])
def test_exact_hdiv_artifact_staging(tmp_path, fault):
    root = tmp_path / 'source'
    root.mkdir()
    source = root / 'native.cu'
    source.write_bytes(b'public generic CUDA source')
    wrappers = root / 'packages/radia-gpu/matlab/+radia/+gpu'
    wrappers.mkdir(parents=True)
    for name in payload.WRAPPERS:
        (wrappers / name).write_text('% wrapper')
    folders = []
    for name in ('hmatrix', 'projection'):
        folder = tmp_path / name
        folder.mkdir()
        binary = folder / ('radia_cuda_' + name + '.mexw64')
        binary.write_bytes(name.encode())
        receipt = dict(accepted=True, source_dirty=False, source_commit='a'*40,
                       taskmanager_started=False, matlab='test fixture',
                       source_sha256={'native.cu': hashlib.sha256(source.read_bytes()).hexdigest()},
                       binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
        if name == 'projection':
            if fault == 'commit': receipt['source_commit'] = 'b'*40
            if fault == 'dirty': receipt['source_dirty'] = True
            if fault == 'binary': receipt['binary_sha256'] = '0'*64
            if fault == 'source': receipt['source_sha256']['native.cu'] = '0'*64
            if fault == 'taskmanager': receipt['taskmanager_started'] = True
        (folder / 'build.json').write_text(json.dumps(receipt))
        folders.append(folder)
    target = tmp_path / 'package'
    if fault:
        with pytest.raises(RuntimeError): payload.stage(root, target, 'a'*40, *folders)
    else:
        files = payload.stage(root, target, 'a'*40, *folders)
        assert len(files) == 6
        for folder, name in zip(folders, ('hmatrix', 'projection')):
            binary = 'radia_cuda_' + name + '.mexw64'
            assert (target / 'matlab' / binary).read_bytes() == (folder / binary).read_bytes()
