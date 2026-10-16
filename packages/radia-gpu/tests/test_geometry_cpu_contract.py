"""Reject mismatched numerical evidence before loading independent CPU solvers."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT=Path(__file__).resolve().parents[3]/'validation_test/ams_cuda/geometry_cpu_validate.py'

@pytest.mark.parametrize('fault,expected',[
    ('binary','loaded MEX/build identity'),('mesh','mesh differs'),
    ('source','Full build source SHA'),('path','actually loaded module'),
    ('test_receipt','distribution eligibility'),('missing_eligibility','distribution eligibility'),
    ('labels','Dirichlet/boundary labels'),('dirichlet','Dirichlet/boundary labels')])
def test_geometry_evidence_identity_failure_is_saved(tmp_path,fault,expected):
    vol=tmp_path/'fixture.vol';vol.write_bytes(b'Input contract only; never a solver mesh')
    mat=tmp_path/'fixture.mat';module=tmp_path/'module.json';receipt=tmp_path/'receipt.json';out=tmp_path/'result.json'
    report=dict(geometry=dict(mesh_sha256=hashlib.sha256(vol.read_bytes()).hexdigest()),
                internal_parity_accepted=True,registry_restored=True,loaded_binary='C:/candidate/radia_cuda_ams.mexw64',
                dirichlet='.*',boundary_labels=['surface'])
    loaded=dict(path=report['loaded_binary'],sha256='b'*64)
    build=dict(source_commit='a'*40,source_dirty=False,sha256={'build-ams-cuda/radia_cuda_ams.mexw64':'b'*64},
               distribution_eligible=True,geometry_lifecycle_test_build=False)
    if fault=='binary':loaded['sha256']='c'*64
    if fault=='mesh':report['geometry']['mesh_sha256']='c'*64
    if fault=='source':build['source_commit']='short'
    if fault=='path':loaded['path']='C:/other/radia_cuda_ams.mexw64'
    if fault=='test_receipt':build.update(distribution_eligible=False,geometry_lifecycle_test_build=True)
    if fault=='missing_eligibility':build.pop('distribution_eligible')
    if fault=='labels':report['boundary_labels']=['surface','surface']
    if fault=='dirichlet':report.pop('dirichlet')
    mat.with_suffix('.mat.json').write_text(json.dumps(report));module.write_text(json.dumps(loaded));receipt.write_text(json.dumps(build))
    run=subprocess.run([sys.executable,str(SCRIPT),'--vol',str(vol),'--mat',str(mat),
        '--module',str(module),'--receipt',str(receipt),'--out',str(out)],capture_output=True,text=True)
    assert run.returncode!=0
    result=json.loads(out.read_text())
    assert not result['accepted'] and not result['timing_accepted'] and not result['accepted_for_adoption']
    assert expected in result['error']
