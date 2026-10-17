"""Trust and source-change negative contracts; no heavy validation or network."""
import copy
from datetime import datetime,timezone,timedelta
import importlib.util
import io
import json
from pathlib import Path
import zipfile
import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('validation_gate',ROOT/'tools/validation_source_gate.py')
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
NOW=datetime(2026,10,17,tzinfo=timezone.utc)
XML=b'<testsuite><testcase name="physics"/></testsuite>'


def proof():
    receipt={'schema':gate.SCHEMA,'profile':gate.PROFILE,'source_sha':'a'*40,'run_id':12,'run_attempt':1,
             'fingerprint':'fingerprint','accepted':True,'reusable':True,'case_count':1,'junit_sha256':gate.digest(XML)}
    run={'id':12,'run_attempt':1,'path':gate.WORKFLOW,'head_branch':'main','status':'completed','conclusion':'success',
         'event':'push','head_sha':'a'*40,'updated_at':NOW.isoformat()}
    return receipt,run


def archive(receipt,xml=XML):
    data=io.BytesIO()
    with zipfile.ZipFile(data,'w') as z:
        z.writestr('receipt.json',json.dumps(receipt,indent=2));z.writestr('solver-heavy.xml',xml)
    return data.getvalue()


def artifact(data):
    return {'name':gate.ARTIFACT,'expired':False,'digest':'sha256:'+gate.digest(data),'workflow_run':{'id':12},
            'expires_at':(NOW+timedelta(days=1)).isoformat()}


def test_trusted_receipt_and_verbatim_evidence(monkeypatch,tmp_path):
    receipt,run=proof();data=archive(receipt)
    calls=[];monkeypatch.setattr(gate.subprocess,'run',lambda command,**kw:calls.append(command))
    actual,xml,raw=gate.verified_receipt(tmp_path,artifact(data),run,data,'fingerprint',NOW)
    assert actual==receipt and xml==XML and raw==json.dumps(receipt,indent=2).encode()
    assert '--is-ancestor' in calls[0] and calls[0][-1]=='origin/main'


@pytest.mark.parametrize('field,value',[('head_branch','feature'),('conclusion','failure'),('event','pull_request'),
                                         ('path','.github/workflows/untrusted.yml'),('run_attempt',2)])
def test_untrusted_run_rejected(monkeypatch,tmp_path,field,value):
    receipt,run=proof();run[field]=value;data=archive(receipt)
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:None)
    with pytest.raises(ValueError):gate.verified_receipt(tmp_path,artifact(data),run,data,'fingerprint',NOW)


@pytest.mark.parametrize('change',['archive','expired','expirytime','old','wrong-input','wrong-sha','failed-evidence','missing-digest','boolean-count','unknown-expiry','nonboolean-expiry'])
def test_bad_artifact_rejected(monkeypatch,tmp_path,change):
    receipt,run=proof();xml=XML
    if change=='wrong-sha':receipt['source_sha']='b'*40
    if change=='boolean-count':receipt['case_count']=True
    if change=='failed-evidence':xml=b'<testsuite><testcase><skipped/></testcase></testsuite>';receipt['junit_sha256']=gate.digest(xml)
    data=archive(receipt,xml);meta=artifact(data)
    if change=='archive':data+=b'corrupt'
    if change=='expired':meta['expired']=True
    if change=='unknown-expiry':meta.pop('expired')
    if change=='nonboolean-expiry':meta['expired']=0
    if change=='expirytime':meta['expires_at']=NOW.isoformat()
    if change=='old':run['updated_at']=(NOW-timedelta(days=31)).isoformat()
    if change=='missing-digest':meta.pop('digest')
    fingerprint='different' if change=='wrong-input' else 'fingerprint'
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:None)
    with pytest.raises(ValueError):gate.verified_receipt(tmp_path,meta,run,data,fingerprint,NOW)


def test_main_ancestry_required(monkeypatch,tmp_path):
    receipt,run=proof();data=archive(receipt)
    def fail(*a,**k):raise gate.subprocess.CalledProcessError(1,['git'])
    monkeypatch.setattr(gate.subprocess,'run',fail)
    with pytest.raises(gate.subprocess.CalledProcessError):gate.verified_receipt(tmp_path,artifact(data),run,data,'fingerprint',NOW)


def test_source_scope_includes_physics_inputs_but_not_workflow_routing(tmp_path,monkeypatch):
    files={'src/solver.py':'physics','validation_test/case/input.vol':'mesh','tests/conftest.py':'fixture',
           'packages/demo/tests/fixture.h5':'fixture-bytes','Build.ps1':'build','.github/workflows/sparsesolv.yml':'routing','tools/validation_source_gate.py':'gate'}
    for path,text in files.items():
        p=tmp_path/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
    manifest=tmp_path/'tests/test_tier_manifest.json'
    manifest.write_text(json.dumps({'profiles':{'solver-heavy':{'paths':['validation_test/case/test.py']}}}))
    test=tmp_path/'validation_test/case/test.py';test.write_text('assert True');files[str(test.relative_to(tmp_path))]='assert True'
    monkeypatch.setattr(gate,'git',lambda root,*args:'\0'.join(files) if args[0]=='ls-files' and '--others' not in args else '')
    initial,_=gate.source_fingerprint(tmp_path,{'deps':'1'})
    (tmp_path/'.github/workflows/sparsesolv.yml').write_text('new routing')
    assert gate.source_fingerprint(tmp_path,{'deps':'1'})[0]==initial
    for path in ['src/solver.py','validation_test/case/input.vol','Build.ps1','tests/conftest.py','packages/demo/tests/fixture.h5']:
        p=tmp_path/path;p.write_text(files[path]+'changed')
        assert gate.source_fingerprint(tmp_path,{'deps':'1'})[0]!=initial
        p.write_text(files[path])
    assert gate.source_fingerprint(tmp_path,{'deps':'2'})[0]!=initial


@pytest.mark.parametrize('lookup', ['miss','error','force','unknown'])
def test_miss_unknown_force_execute_validation(tmp_path,monkeypatch,lookup):
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    def fingerprint(*a):
        if lookup=='unknown':raise ValueError('unknown')
        return 'fingerprint',1
    monkeypatch.setattr(gate,'source_fingerprint',fingerprint)
    class Client:
        def find(self,*a):
            if lookup=='force':pytest.fail('force must not look up cache')
            if lookup=='error':raise OSError('offline')
            return None
    calls=[]
    def execute(command,**kw):
        calls.append(command);(tmp_path/'out/solver-heavy.xml').write_bytes(XML)
        return type('Result',(),{'returncode':0})()
    monkeypatch.setattr(gate.subprocess,'run',execute)
    assert gate.run(tmp_path,tmp_path/'out',lookup=='force',Client())==0
    report=json.loads((tmp_path/'out/verification.json').read_text())
    assert calls and report['mode']=='execute' and not report['current_artifact_certified']


def test_failed_execution_never_passes(tmp_path,monkeypatch):
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    monkeypatch.setattr(gate,'source_fingerprint',lambda *a:('fingerprint',1))
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:type('Result',(),{'returncode':1})())
    assert gate.run(tmp_path,tmp_path/'out')==1
    assert not json.loads((tmp_path/'out/verification.json').read_text())['accepted']


def test_workflow_keeps_ordinary_and_release_gates_separate():
    workflow=(ROOT/'.github/workflows/sparsesolv.yml').read_text()
    assert '--profile solver-numerics' in workflow and '-RadiaOnly -BuildJobs 2' in workflow
    assert 'tools/validation_source_gate.py' in workflow and 'force_validation:' in workflow
    assert "steps.heavy-validation.outputs.mode == 'execute'" in workflow
    assert 'actions: read' in workflow

def test_redirect_does_not_forward_token_to_artifact_storage():
    request=gate.urllib.request.Request('https://api.github.com/repos/example/repo/actions/artifacts/1/zip',headers={'Authorization':'Bearer secret'})
    redirect=gate._SafeRedirect().redirect_request(request,None,302,'redirect',{},'https://storage.example/signed-archive')
    assert not redirect.has_header('Authorization')
    with pytest.raises(ValueError):gate._SafeRedirect().redirect_request(request,None,302,'redirect',{},'http://storage.example/archive')

def test_verified_reuse_keeps_original_identity_and_bytes(tmp_path,monkeypatch):
    receipt,prior=proof();data=archive(receipt);raw=json.dumps(receipt,indent=2).encode()
    monkeypatch.setattr(gate,'git',lambda *a:'b'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    monkeypatch.setattr(gate,'source_fingerprint',lambda *a:('fingerprint',1))
    class Client:
        def find(self,*a):return {'id':99},receipt,XML,raw
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:pytest.fail('reused validation executed'))
    assert gate.run(tmp_path,tmp_path/'out',client=Client())==0
    report=json.loads((tmp_path/'out/verification.json').read_text())
    assert report['current']['source_sha']=='b'*40 and report['original']['source_sha']=='a'*40
    assert not report['current_artifact_certified'] and report['mode']=='verified-reuse'
    assert (tmp_path/'out/receipt.json').read_bytes()==raw

def test_old_xml_cannot_turn_empty_execution_into_pass(tmp_path,monkeypatch):
    out=tmp_path/'out';out.mkdir();(out/'solver-heavy.xml').write_bytes(XML)
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    monkeypatch.setattr(gate,'source_fingerprint',lambda *a:('fingerprint',1))
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:type('Result',(),{'returncode':0})())
    assert gate.run(tmp_path,out)==1
    assert not json.loads((out/'verification.json').read_text())['accepted']


@pytest.mark.parametrize('input_path', ['src/extra.py', 'validation_test/case/generated.h5',
                                       'packages/demo/tests/new_fixture.vol'])
def test_unknown_untracked_inputs_prevent_reuse(tmp_path, monkeypatch, input_path):
    manifest=tmp_path/'tests/test_tier_manifest.json'
    manifest.parent.mkdir()
    manifest.write_text(json.dumps({'profiles':{'solver-heavy':{'paths':[]}}}))
    def files(root, *args):
        return input_path if '--others' in args else 'src/solver.py'
    monkeypatch.setattr(gate, 'git', files)
    with pytest.raises(ValueError, match='untracked validation input'):
        gate.source_fingerprint(tmp_path, {})


@pytest.mark.parametrize('failure', [None, 'schema', 'source_commit', 'source_dirty',
                                     'source_change_fingerprint_sha256', 'binary_name',
                                     'binary_bytes', 'binary_sha256'])
def test_current_native_sidecar_is_hash_and_source_bound(tmp_path, monkeypatch, failure):
    binary=tmp_path/'src/radia/_radia_pybind.pyd'
    binary.parent.mkdir(parents=True)
    data=b'current test binary';binary.write_bytes(data)
    record=dict(schema='radia.native-build-provenance.v1',source_commit='a'*40,
                source_dirty=False,source_change_fingerprint_sha256=gate.digest(b'\0'),
                binary_name=binary.name,binary_bytes=len(data),binary_sha256=gate.digest(data))
    if failure:record[failure]='wrong'
    Path(str(binary)+'.build.json').write_text(json.dumps(record))
    monkeypatch.setattr(gate, 'git', lambda *a:'a'*40)
    if failure:
        with pytest.raises(ValueError, match='native provenance mismatch'):
            gate.validate_native_sidecar(tmp_path)
    else:
        gate.validate_native_sidecar(tmp_path)
