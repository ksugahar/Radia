"""Guards: explicit role bindings and Radia-only forward editable operations.

Assumes: remote transport and installers are observed mocks; binding validation
and intent persistence have independent executable tests in test_editable_intent.
"""
from __future__ import annotations
import importlib.util
import inspect
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

_TOOL=Path(__file__).resolve().parents[1]/'tools/release_quad.py'
_SPEC=importlib.util.spec_from_file_location('radia_release_quad_intent',_TOOL)
release_quad=importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)
intent=release_quad.editable_intent

@pytest.fixture(autouse=True)
def explicit_roles(monkeypatch,tmp_path):
    monkeypatch.setenv(intent.INTENT_FILE_ENV,str(tmp_path/'intent.json'))
    for role in ('100','VALIDATION'):
        host='build-test' if role=='100' else 'validation-test'
        values={'SSH':host,'HOSTNAME':host,'PYTHON':f'C:/Runtime/{host}/python.exe',
                'QA_PYTHON':f'C:/QA/{host}/Scripts/python.exe',
                'EDITABLE_REPO':f'C:/Sources/{host}','EDITABLE_COMMIT':'a'*40,
                'EDITABLE_PATH_ALIASES':'[]'}
        for key,value in values.items():monkeypatch.setenv(f'RADIA_RELEASE_{key}_{role}',value)


def response(counts,rc=0):
    return rc,dict(schema=intent.SCHEMA,rows=[],counts=counts,record_file='test.json',interpreter='C:/Synthetic/python.exe',exit_code=rc),''


def test_remote_verifier_is_the_shipped_tool_without_private_default_aliases():
    assert release_quad.REMOTE_EDITABLE_VERIFY==_TOOL.with_name('editable_intent.py').read_text(encoding='utf-8')
    assert release_quad._PATH_ALIASES==()


def test_verify_editable_visits_both_normal_roles(monkeypatch):
    calls=[]
    monkeypatch.setattr(release_quad,'_verify_role_editable',lambda role:calls.append(role) or 0)
    assert release_quad.cmd_verify_editable(None)==0
    assert set(calls)=={'100','validation'} and len(calls)==2


@pytest.mark.parametrize('role',['100','validation'])
def test_normal_verification_limits_the_remote_cli_to_radia(monkeypatch,role):
    calls=[]
    monkeypatch.setattr(release_quad,'_verify_live_editable',lambda target:0)
    monkeypatch.setattr(release_quad,'_remote_editable_intent',lambda target,argv:calls.append((target,argv)) or response(dict(match=1,drift=0,unverified=0)))
    assert release_quad._verify_role_editable(role)==0
    target,argv=calls[0]
    assert target['role']==role
    assert argv==['verify','--package','radia','--expect','radia='+target['source']]


@pytest.mark.parametrize('counts,rc,expected',[(dict(match=1,drift=0,unverified=0),0,0),(dict(match=0,drift=1,unverified=0),1,1),(dict(match=0,drift=0,unverified=1),1,1),(dict(match=1,drift=0,unverified=0),1,1)])
def test_remote_drift_unverified_or_failed_exit_are_not_success(monkeypatch,counts,rc,expected):
    monkeypatch.setattr(release_quad,'_verify_live_editable',lambda target:0)
    monkeypatch.setattr(release_quad,'_remote_editable_intent',lambda *args:response(counts,rc))
    report={}
    assert release_quad._verify_role_editable('validation',report)==expected
    assert report['drift']==counts['drift'] and report['unverified']==counts['unverified']


def test_failed_live_check_stops_before_intent_transport(monkeypatch):
    monkeypatch.setattr(release_quad,'_verify_live_editable',lambda target:4)
    monkeypatch.setattr(release_quad,'_remote_editable_intent',lambda *args:pytest.fail('No transport after failed payload check'))
    assert release_quad._verify_role_editable('100')==1


@pytest.mark.parametrize('role',['100','validation'])
def test_record_current_binds_one_role_and_one_package(monkeypatch,role):
    seen=[]
    monkeypatch.setattr(release_quad,'_remote_editable_intent',lambda target,argv:seen.append((target,argv)) or (0,{'results':[],'record_file':'r','exit_code':0},''))
    assert release_quad._record_release_intent_remote(role,'reviewed update')==0
    target,argv=seen[0]
    assert target['role']==role
    assert argv[:4]==['repoint','--record-current','--package','radia']
    assert '--source' not in argv and '--rollback' not in argv


def args_for(role='validation',**updates):
    values=dict(host=role,package=['radia'],source=[f'C:/Sources/{"build-test" if role=="100" else "validation-test"}'],reason='reviewed update',record_current=False,rollback=False,dry_run=False)
    return SimpleNamespace(**(values|updates))


@pytest.mark.parametrize('role',['100','validation'])
def test_repoint_dry_run_exposes_bound_plan_without_any_action(monkeypatch,capsys,role):
    monkeypatch.setattr(release_quad,'run',lambda *a,**k:pytest.fail('dry-run network action'))
    assert release_quad.cmd_repoint(args_for(role,dry_run=True))==0
    plan=json.loads(capsys.readouterr().out)
    assert plan['role']==role and plan['mode']=='editable'
    assert plan['binding']['role']==('build' if role=='100' else 'validation')
    assert plan['install'][-1]==args_for(role).source[0]+'[gpu]'


@pytest.mark.parametrize('updates',[{'source':[]},{'source':['C:/Wrong']},{'reason':''},{'record_current':True},{'rollback':True},{'package':['radia-mcp']},{'host':'retired-test'}])
def test_repoint_refuses_incomplete_wrong_source_and_foreign_package_requests(monkeypatch,updates):
    monkeypatch.setattr(release_quad,'_deploy_editable_role',lambda *a:pytest.fail('Invalid plan mutated runtime'))
    monkeypatch.setattr(release_quad,'_remote_editable_intent',lambda *a:pytest.fail('Invalid plan recorded intent'))
    assert release_quad.cmd_repoint(args_for(**updates))==2


@pytest.mark.parametrize('role',['100','validation'])
def test_repoint_routes_only_selected_normal_role(monkeypatch,role):
    seen=[]
    monkeypatch.setattr(release_quad,'_deploy_editable_role',lambda actual,reason:seen.append((actual,reason)) or 3)
    assert release_quad.cmd_repoint(args_for(role))==3
    assert seen==[(role,'reviewed update')]


@pytest.mark.parametrize('live_rc',[0,4])
def test_record_current_checks_payload_before_recording(monkeypatch,live_rc):
    seen=[]
    monkeypatch.setattr(release_quad,'_verify_live_editable',lambda target:seen.append('live') or live_rc)
    monkeypatch.setattr(release_quad,'_record_release_intent_remote',lambda *a:seen.append('record') or 0)
    monkeypatch.setattr(release_quad,'_verify_role_editable',lambda *a:seen.append('verify') or 0)
    assert release_quad.cmd_repoint(args_for(record_current=True,source=[]))==(4 if live_rc else 0)
    assert seen==(['live'] if live_rc else ['live','record','verify'])


def test_done_source_never_uses_an_intent_entry(tmp_path):
    record=intent.load_intent()
    intent.set_entry(record,'radia',{'source':'C:/Unrelated','commit':'b'*40})
    intent.save_intent(record)
    assert release_quad._done_release_source(SimpleNamespace())==str(release_quad.REPO)
    selected=tmp_path/'explicit-release'
    assert release_quad._done_release_source(SimpleNamespace(release_source=str(selected)))==str(selected.resolve())


def test_repoint_never_stops_foreign_processes_or_uninstalls_packages():
    for name in ('cmd_repoint','cmd_restore_editable','cmd_verify_editable','_record_release_intent_remote','_verify_role_editable','_remote_editable_intent','_editable_install_program'):
        source=inspect.getsource(getattr(release_quad,name))
        for forbidden in ('pip uninstall','uninstall -y','Stop-Process','taskkill'):
            assert forbidden not in source,(name,forbidden)
