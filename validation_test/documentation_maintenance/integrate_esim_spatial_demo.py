"""Integrate recovered mdx artifacts and regenerate canonical citations/sidecar."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'packages/radia-mcp/src'))
from radia_mcp.document_meta.tools import document_meta_write_notebook_result_json
from render_notebook_bibliography import render


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def integrate(stage):
    doc=ROOT/'docs/ih_esim_benchmark'
    folder=doc/'spatial_demo_data'
    folder.mkdir(exist_ok=True)
    manifest=json.loads((stage/'recovery_hashes.json').read_bytes())
    for name, expected in manifest.items():
        assert digest(stage/name)==expected, name
    names=['compute_record.json','vol_check.json','spatial_metrics.json']
    for f in (10000,50000,100000):
        names += [f'f{f}.json',f'f{f}_qsurf.sol',f'f{f}.log']
    for name in names: shutil.copyfile(stage/name,folder/name)
    nbpath=doc/'esim_spatial_demo.ipynb'
    shutil.copyfile(stage/nbpath.name,nbpath)
    nb=json.loads(nbpath.read_bytes())
    text=''.join(nb['cells'][0]['source'])
    opening,marker,theory=text.partition('## Model and equations')
    if marker:
        nb['cells'][0]['source']=opening.splitlines(keepends=True)
        nb['cells'].append(dict(cell_type='markdown',id='esim-spatial-theory',metadata={},
                               source=(marker+theory).splitlines(keepends=True)))
        nbpath.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
    render(nbpath)
    run=json.loads((folder/'compute_record.json').read_bytes())
    execution=dict(hostname=run['hostname'],solver_versions=run['versions'],
                   runtime=json.loads((stage/'visualization_runtime.json').read_bytes()),
                   visual_qa='Saved scene/state audit performed; rendered visual QA pending (browser URL security boundary).',
                   scope='Fresh three-frequency weak-coupled ESIM solves and notebook reconstruction on mdx2; no thermal or volume-FEM validation.')
    receipt=document_meta_write_notebook_result_json(str(nbpath),overwrite=True)
    assert 'error' not in receipt,receipt
    sidepath=nbpath.with_name(nbpath.stem+'_result.json')
    side=json.loads(sidepath.read_bytes())
    side['sidecar_writer_versions']=side.pop('versions',{})
    side['versions']=run['versions']
    side['execution']=execution
    side['artifacts_sha256']={name:digest(folder/name) for name in names}
    sidepath.write_text(json.dumps(side,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    audit=ROOT/'validation_test/documentation_maintenance/esim_spatial_20260915'
    audit.mkdir(exist_ok=True)
    for name in ('recovery_hashes.json','run_cases.py','visualization_runtime.json','notebook_execution.log','missing_ocp.log'):
        if (stage/name).exists(): shutil.copyfile(stage/name,audit/name)
    receipt=dict(generated_at_utc=datetime.now(timezone.utc).isoformat(),
        recovered_sha256=manifest, notebook_sha256=digest(nbpath), execution=execution,
        helper_source=dict(path='docs/ih_esim_benchmark/esim_spatial_demo.py',
                           sha256=digest(doc/'esim_spatial_demo.py'),
                           source=(doc/'esim_spatial_demo.py').read_text(encoding='utf-8')),
        remote_cleanup='Pending evidence commit; owner: documentation maintenance; trigger: commit and hash verification.')
    (audit/'receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    # Add a discovery link without touching any historical code/output cell.
    oldpath=doc/'esim_showcase.ipynb'
    old=json.loads(oldpath.read_bytes())
    cell_id='esim-spatial-demo-discovery'
    if not any(c.get('id')==cell_id for c in old['cells']):
        old['cells'].insert(1,dict(cell_type='markdown',id=cell_id,metadata={},source=[
            '## Inspect the coil, field and heating in 3D\n\n',
            '[Open the ESIM spatial demo](esim_spatial_demo.ipynb) for fresh mdx2 results at 10, 50 and 100 kHz, ',
            'saved WebGUI scenes, and the distinction between local ESIM heating and the normalized heat-transfer artifact.\n']))
        oldpath.write_text(json.dumps(old,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
        oldside=doc/'esim_showcase_result.json'
        d=json.loads(oldside.read_bytes())
        d['notebook_sha256']=digest(oldpath)
        d['discovery_link_update']={'date':'2026-09-15','code_and_saved_outputs_unchanged':True}
        oldside.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('recovered',type=Path)
    integrate(p.parse_args().recovered)
