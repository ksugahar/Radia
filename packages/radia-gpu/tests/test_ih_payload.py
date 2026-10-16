"""The research IH distribution must reject stale or incomplete provenance."""
from pathlib import Path
import hashlib
import importlib.util
import json
import pytest

spec=importlib.util.spec_from_file_location('ih_payload',Path(__file__).parents[1]/'ih_payload.py')
payload=importlib.util.module_from_spec(spec);spec.loader.exec_module(payload)


@pytest.mark.parametrize('fault',[None,'commit','dirty','binary','source',
                                  'python','taskmanager','features','sample','missing-source','test-build','ineligible','missing-eligibility'])
def test_ih_exact_research_staging(tmp_path,fault):
    root=tmp_path/'source';folder=tmp_path/'build';folder.mkdir()
    sources={}
    for relative in payload.REQUIRED_SOURCES:
        p=root/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'generic kernel fixture')
        sources[relative]=hashlib.sha256(p.read_bytes()).hexdigest()
    for name in payload.MATLAB_FILES:
        p=root/'packages/radia-gpu/matlab'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('% fixture')
    for name in payload.SAMPLES+('TEMPERATURE_BH_EXPERIMENT.md',):
        p=root/'validation_test/ih_cuda'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('% fixture')
    doc=root/'packages/radia-gpu/docs/IH_TEMPERATURE_BH.md';doc.parent.mkdir();doc.write_text('Research only')
    binary=folder/'radia_cuda_ih.mexw64';binary.write_bytes(b'fixture artifact')
    receipt=dict(accepted=True,transaction_test_build=False,distribution_eligible=True,source_commit='a'*40,source_dirty=False,
                 taskmanager_regions=False,python_runtime_required=False,cpu_fallback=False,
                 features=sorted(payload.FEATURES),source_sha256=sources,
                 binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                 matlab='fixture',cuda_targets=['fixture'])
    if fault=='test-build':receipt['transaction_test_build']=True
    if fault=='ineligible':receipt['distribution_eligible']=False
    if fault=='missing-eligibility':receipt.pop('distribution_eligible')
    if fault=='commit':receipt['source_commit']='b'*40
    if fault=='dirty':receipt['source_dirty']=True
    if fault=='binary':binary.write_bytes(b'changed')
    if fault=='source':sources[next(iter(sources))]='0'*64
    if fault=='python':receipt['python_runtime_required']=True
    if fault=='taskmanager':receipt['taskmanager_regions']=True
    if fault=='features':receipt['features']=[]
    if fault=='missing-source':sources.pop('packages/radia-gpu/native/ih_esim_cell.cuh')
    if fault=='sample':(root/'validation_test/ih_cuda'/payload.SAMPLES[0]).unlink()
    (folder/'build.json').write_text(json.dumps(receipt))
    out=tmp_path/'package'
    if fault:
        with pytest.raises(RuntimeError):payload.stage(root,out,'a'*40,folder)
        assert not out.exists()
    else:
        files=payload.stage(root,out,'a'*40,folder)
        assert files=={p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file()}
        public=json.loads((out/'ih-build.json').read_text())
        assert public['research_only'] and not public['experimental_validation']
        assert public['binary_sha256']==hashlib.sha256((out/'matlab/radia_cuda_ih.mexw64').read_bytes()).hexdigest()
