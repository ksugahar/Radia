"""Guards: V-N skips, missing cases and interrupted runs cannot become PASS.

Assumes: synthetic pytest cases in clean temporary Git repositories, not physics
models or release certification. Native identity/resource admission are producer
contracts outside this receipt plugin; elapsed times are observations only.
"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('validation_receipt', ROOT/'tools/validation_execution_receipt.py')
receipt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(receipt)


def repo(tmp_path, body='def test_case():\n    assert True\n'):
    source = tmp_path/'source'
    source.mkdir()
    (source/'validation_test').mkdir()
    (source/'validation_test/test_case.py').write_text(body, encoding='utf-8')
    (source/'pytest.ini').write_text('[pytest]\n', encoding='utf-8')
    def git(*args):
        return subprocess.check_output(['git','-c',f'core.hooksPath={source / ".git/empty-hooks"}','-c','user.name=Synthetic','-c',
                                       'user.email=synthetic@example.invalid',*args], cwd=source).decode().strip()
    git('init','-q');git('add','.');git('commit','-qm','Synthetic source')
    plan = dict(schema='radia.validation-execution-plan.v1',source_sha=git('rev-parse','HEAD'),
                source_tree=git('rev-parse','HEAD^{tree}'),run_id='synthetic',attempt=1,
                profile=dict(name='synthetic-V-N',sha256='1'*64),
                expected_nodeids=['validation_test/test_case.py::test_case'])
    plan_path=tmp_path/'plan.json';plan_path.write_text(json.dumps(plan),encoding='utf-8')
    return source,plan_path,tmp_path/'receipt.json',tmp_path/'outputs'


def execute(paths, extra=()):
    source,plan,out,outputs=paths
    env=dict(os.environ,PYTHONPATH=str(ROOT/'tools'),PYTEST_DISABLE_PLUGIN_AUTOLOAD='1',
             PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',PYTHONUTF8='1')
    command=[sys.executable,'-m','pytest','-p','validation_execution_receipt','--noconftest',
             '--rootdir',str(source),'--validation-plan',str(plan),
             '--validation-source-root',str(source),'--validation-receipt',str(out),
             '--validation-output-root',str(outputs),'validation_test','-q',*extra]
    process=subprocess.run(command,cwd=source,env=env,capture_output=True,text=True,
                           encoding='utf-8',errors='replace',timeout=40)
    data=json.loads(out.read_text(encoding='utf-8')) if out.exists() else None
    return process,data


@pytest.mark.parametrize('body,status',[
 ('def test_case():\n    assert True\n','passed'),
 ('def test_case():\n    assert False\n','failed'),
 ('import pytest\ndef test_case():\n    pytest.skip("visible")\n','incomplete'),
 ('import pytest\n@pytest.mark.xfail(reason="known")\ndef test_case():\n    assert False\n','incomplete'),
 ('import pytest\npytest.skip("module skipped", allow_module_level=True)\ndef test_case(): pass\n','incomplete'),
 ('import pytest\n@pytest.mark.xfail(reason="expected failure")\ndef test_case():\n    assert True\n','incomplete'),
 ('import pytest\n@pytest.fixture\ndef broken():\n    raise RuntimeError("setup failed")\ndef test_case(broken): pass\n','failed'),
 ('import pytest\n@pytest.fixture\ndef broken():\n    yield\n    raise RuntimeError("teardown failed")\ndef test_case(broken): pass\n','failed'),
])
def test_real_pytest_outcomes_preserve_completeness(tmp_path,body,status):
    result,data=execute(repo(tmp_path,body))
    assert data['status']==status,(result.stdout,result.stderr,data)
    assert (result.returncode==0)==(status=='passed')
    assert data['qualification_scope'].startswith('source-nightly-only')
    assert data['state']=='finished'
    assert bool(data['passed_pytest_cases']) == (status=='passed')
    if status=='incomplete': assert data['executed']==[]


def test_missing_expected_case_is_incomplete(tmp_path):
    paths=repo(tmp_path);plan=json.loads(paths[1].read_text())
    plan['expected_nodeids'].append('validation_test/test_missing.py::test_missing')
    paths[1].write_text(json.dumps(plan))
    result,data=execute(paths)
    assert result.returncode!=0 and data['status']=='incomplete'
    assert data['missing_collection']==['validation_test/test_missing.py::test_missing']
    assert data['executed']==['validation_test/test_case.py::test_case']


def test_unexpected_case_is_blocked_before_execution(tmp_path):
    body='def test_case(): pass\ndef test_unreviewed():\n    raise RuntimeError("must not execute")\n'
    result,data=execute(repo(tmp_path,body))
    assert result.returncode!=0 and data['status']=='incomplete'
    assert data['executed']==[]
    assert data['unexpected']==['validation_test/test_case.py::test_unreviewed']


@pytest.mark.parametrize('extra',[('--collect-only',),('-k','never_selected')])
def test_collect_only_and_deselection_never_pass(tmp_path,extra):
    paths=repo(tmp_path)
    result,data=execute(paths,extra)
    assert result.returncode!=0 and data['status']=='incomplete'
    assert data['missing_execution']==['validation_test/test_case.py::test_case']


def test_interrupted_process_leaves_incomplete_checkpoint(tmp_path):
    result,data=execute(repo(tmp_path,'import os\ndef test_case():\n    os._exit(9)\n'))
    assert result.returncode==9
    assert data['status']=='incomplete' and data['state']=='running'
    assert data['missing_execution']==['validation_test/test_case.py::test_case']


def test_workspace_change_is_failed_without_reset(tmp_path):
    body='from pathlib import Path\ndef test_case():\n    Path("pytest.ini").write_text("changed")\n'
    paths=repo(tmp_path,body);result,data=execute(paths)
    assert result.returncode!=0 and data['status']=='failed' and data['workspace_changed']
    assert (paths[0]/'pytest.ini').read_text()=='changed'


def test_checkpoint_and_numerical_file_hashes_are_observations(tmp_path):
    paths=repo(tmp_path);record=receipt.Receipt(paths[1],paths[0],paths[2],paths[3])
    record.collection(['validation_test/test_case.py::test_case'])
    for phase in ('setup','call','teardown'):
        record.phase('validation_test/test_case.py::test_case',phase,'passed',duration_s=0.001)
    assert json.loads(paths[2].read_text())['status']=='incomplete'
    paths[3].mkdir();numeric=paths[3]/'observed.json';numeric.write_text('{"value": 2}')
    assert record.finish(0)==0
    data=json.loads(paths[2].read_text())
    assert data['numerical_outputs'][0]['sha256']==receipt.sha256(numeric)
    assert data['numerical_output_status']=='observed'
    assert data['phase_durations_s']['validation_test/test_case.py::test_case']['call']==0.001
    with pytest.raises(ValueError,match='already exists'):receipt.Receipt(paths[1],paths[0],paths[2],paths[3])


@pytest.mark.parametrize('change',['empty','duplicate','unsafe','source','profile','dirty'])
def test_invalid_or_unbound_plans_fail_loudly(tmp_path,change):
    paths=repo(tmp_path);plan=json.loads(paths[1].read_text())
    if change=='empty':plan['expected_nodeids']=[]
    elif change=='duplicate':plan['expected_nodeids']*=2
    elif change=='unsafe':plan['expected_nodeids']=['validation_test/../escape.py::test_case']
    elif change=='source':plan['source_sha']='0'*40
    elif change=='profile':plan['profile']['sha256']='not-a-hash'
    elif change=='dirty':(paths[0]/'untracked.txt').write_text('dirty')
    paths[1].write_text(json.dumps(plan))
    with pytest.raises(ValueError):receipt.Receipt(paths[1],paths[0],paths[2],paths[3])


def test_evidence_cannot_overwrite_source(tmp_path):
    source,plan,_,outputs=repo(tmp_path)
    with pytest.raises(ValueError,match='outside source'):
        receipt.Receipt(plan,source,source/'receipt.json',outputs)


def test_output_scope_cannot_include_source_or_receipt(tmp_path):
    source,plan,out,_=repo(tmp_path)
    with pytest.raises(ValueError,match='output scope'):
        receipt.Receipt(plan,source,out,tmp_path)


def test_parent_traversal_cannot_redirect_receipt_into_source(tmp_path):
    source,plan,_,outputs=repo(tmp_path)
    with pytest.raises(ValueError,match='outside source'):
        receipt.Receipt(plan,source,tmp_path/'nonexistent/../source/receipt.json',outputs)


def test_missing_phase_reports_never_pass(tmp_path):
    paths=repo(tmp_path);record=receipt.Receipt(paths[1],paths[0],paths[2],paths[3])
    record.collection(['validation_test/test_case.py::test_case'])
    record.phase('validation_test/test_case.py::test_case','call','passed')
    assert record.finish(0)!=0
    assert json.loads(paths[2].read_text())['status']=='incomplete'


def test_source_identity_cannot_change_after_plan_binding(tmp_path,monkeypatch):
    paths=repo(tmp_path);record=receipt.Receipt(paths[1],paths[0],paths[2],paths[3])
    record.collection(['validation_test/test_case.py::test_case'])
    for phase in ('setup','call','teardown'):
        record.phase('validation_test/test_case.py::test_case',phase,'passed')
    original=receipt.git
    def changed(root,*args):
        return '2'*40 if args==('rev-parse','HEAD') else original(root,*args)
    monkeypatch.setattr(receipt,'git',changed)
    assert record.finish(0)!=0
    assert record.data['source_identity_changed'] and record.data['status']=='failed'


def test_duplicate_phase_is_not_a_successful_retry(tmp_path):
    paths=repo(tmp_path);record=receipt.Receipt(paths[1],paths[0],paths[2],paths[3])
    record.collection(['validation_test/test_case.py::test_case'])
    for phase in ('setup','call','call','teardown'):
        record.phase('validation_test/test_case.py::test_case',phase,'passed')
    assert record.finish(0)!=0 and record.data['status']=='incomplete'


def test_duplicate_json_keys_are_refused(tmp_path):
    plan=tmp_path/'duplicate.json';plan.write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError,match='duplicate plan key'):receipt.read_plan(plan)
