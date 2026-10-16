"""Privacy and completeness are mandatory for IH evidence export."""
import importlib.util
from pathlib import Path
import pytest
spec=importlib.util.spec_from_file_location('ih_public',Path(__file__).parents[1]/'ih_public_evidence.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

@pytest.mark.parametrize('fault',[None,'unknown','nested','nan','pathversion','incomplete','power'])
def test_ih_public_evidence(fault):
    r=dict(schema='radia.gpu.ih-installed-qualification.v1',accepted=True,source_commit='a'*40,toolkit_commit='b'*40,wheel_sha256='c'*64,manifest_sha256='d'*64,binary_sha256='e'*64,matlab_version='26.1.0.3346908 (R2026a) Update 5',toolkit_accepted=True,cpu=dict(accepted=True,relative_errors=dict.fromkeys(module.ERRORS,1e-12),true_relative_residual=1e-12,scope='private free text'),lifecycle=dict.fromkeys(module.LIFECYCLE,True),loaded_mex='C:/private/example.mexw64',hostname='private-host')
    if fault=='unknown':r['unrecognized_path']='C:/private'
    if fault=='nested':r['cpu']['host_path']='C:/private'
    if fault=='nan':r['cpu']['relative_errors']['heat']=float('nan')
    if fault=='pathversion':r['matlab_version']='C:/private'
    if fault=='incomplete':r['lifecycle'].pop('reset')
    if fault=='power':r['power']={'accepted':True}
    if fault:
        with pytest.raises(ValueError):module.ih_public(r)
    else:
        public=module.ih_public(r)
        assert 'loaded_mex' not in public and 'hostname' not in public
        assert 'scope' not in public['cpu']
        assert public['release_qualified'] is False
