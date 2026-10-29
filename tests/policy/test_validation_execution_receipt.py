"""Guards: V-N skips, missing cases and interrupted runs cannot become PASS.

Assumes: synthetic pytest cases in clean temporary Git repositories, not physics
models or release certification. Native identity/resource admission are producer
contracts outside this receipt plugin; elapsed times are observations only.
"""
import importlib.util
import hashlib
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
    assert data['numerical_outputs'][0]['sha256']==hashlib.sha256(b'{"value": 2}').hexdigest()
    assert data['numerical_outputs'][0]['bytes']==len(b'{"value": 2}')
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


def test_modified_copy_rootdir_cannot_certify_bound_checkout(tmp_path):
    paths = repo(tmp_path, 'def test_case():\n    assert False\n')
    copy = tmp_path/'copy'
    copy.mkdir()
    (copy/'validation_test').mkdir()
    (copy/'validation_test/test_case.py').write_text('def test_case(): pass\n')
    result, data = execute(paths, ('--rootdir', str(copy)))
    assert result.returncode != 0 and data is None
    assert 'rootdir must equal' in result.stderr


def test_untracked_or_external_item_is_rejected(tmp_path):
    from types import SimpleNamespace
    paths = repo(tmp_path)
    record = receipt.Receipt(paths[1], paths[0], paths[2], paths[3])
    for path in (paths[0]/'validation_test/untracked.py', tmp_path/'external.py'):
        path.write_text('def test_case(): pass')
        item = SimpleNamespace(path=path, nodeid='validation_test/untracked.py::test_case')
        with pytest.raises(pytest.UsageError, match='tracked source'):
            record.collection([item.nodeid], [item])


def test_zero_exit_interruption_is_rejected_by_consumer(tmp_path):
    paths = repo(tmp_path, 'import os\ndef test_case():\n    os._exit(0)\n')
    result, data = execute(paths)
    assert result.returncode == 0 and data['state'] == 'running'
    with pytest.raises(ValueError, match='finished'):
        receipt.verify_receipt(paths[2], paths[1])


def test_consumer_verifies_finished_identity_and_plan(tmp_path):
    paths = repo(tmp_path)
    with pytest.raises(OSError):
        receipt.verify_receipt(paths[2], paths[1])
    result, data = execute(paths)
    assert result.returncode == 0
    assert receipt.verify_receipt(paths[2], paths[1]) == data
    for key, value in [('effective_exit_code', 1), ('plan_sha256', '0'*64),
                       ('identity', {}), ('state', 'running'), ('executed', [])]:
        changed = dict(data, **{key: value})
        paths[2].write_text(json.dumps(changed))
        with pytest.raises(ValueError):
            receipt.verify_receipt(paths[2], paths[1])


def _publish_probe(tmp_path, monkeypatch, *, windows=True, oversleep=False):
    """Use real checkpoint files and an independent, controllable I/O clock."""
    from types import SimpleNamespace
    paths = repo(tmp_path)
    record = receipt.Receipt(paths[1], paths[0], paths[2], paths[3])
    clock = dict(now=0.0, sleeps=[])
    def sleep(seconds):
        clock['sleeps'].append(seconds)
        clock['now'] += 1.0 if oversleep else seconds
    monkeypatch.setattr(receipt, '_WINDOWS', windows)
    monkeypatch.setattr(receipt, 'time', SimpleNamespace(
        monotonic=lambda: clock['now'], sleep=sleep))
    record.started = 0.0
    return record, paths[2], clock


def _windows_publish_error(code):
    error = OSError('injected atomic publication failure')
    error.winerror = code
    return error


@pytest.mark.parametrize('code', [5, 32])
def test_transient_windows_publish_records_actual_retries(tmp_path, monkeypatch, code):
    record, path, clock = _publish_probe(tmp_path, monkeypatch)
    previous = record.data['write_retries']
    replace = receipt.os.replace
    calls = []
    def transient(source, target):
        calls.append(clock['now'])
        if len(calls) <= 2:
            raise _windows_publish_error(code)
        return replace(source, target)
    monkeypatch.setattr(receipt.os, 'replace', transient)
    record.write()
    data = json.loads(path.read_text())
    assert len(calls) == 3, f'actual replacement attempts={calls}, expected 3'
    assert data['write_retries'] == previous + 2, (data['write_retries'], previous + 2)
    assert record.data['write_retries'] == data['write_retries']
    assert not path.with_name(path.name + '.pending').exists()


def test_permanent_windows_publish_keeps_checkpoint_and_pending(tmp_path, monkeypatch):
    record, path, clock = _publish_probe(tmp_path, monkeypatch)
    before = path.read_bytes()
    previous = record.data['write_retries']
    calls = []
    error = _windows_publish_error(5)
    def permanent(*args):
        calls.append(clock['now'])
        raise error
    monkeypatch.setattr(receipt.os, 'replace', permanent)
    with pytest.raises(OSError) as raised:
        record.write()
    assert raised.value is error
    assert len(calls) == 10, f'actual replacement attempts={calls}, expected 10'
    assert clock['now'] <= 1.0, (clock['now'], 'at most 1 second')
    assert path.read_bytes() == before, 'failed publication changed prior checkpoint'
    pending = json.loads(path.with_name(path.name + '.pending').read_text())
    assert pending['write_retries'] == previous + 9, (pending['write_retries'], previous + 9)


def test_publish_deadline_stops_before_another_attempt(tmp_path, monkeypatch):
    record, path, clock = _publish_probe(tmp_path, monkeypatch, oversleep=True)
    calls = []
    previous = record.data['write_retries']
    error = _windows_publish_error(32)
    def blocked(*args):
        calls.append(clock['now'])
        raise error
    monkeypatch.setattr(receipt.os, 'replace', blocked)
    with pytest.raises(OSError) as raised:
        record.write()
    assert raised.value is error
    assert calls == [0.0], f'actual attempts={calls}, expected only initial attempt'
    assert record.data['write_retries'] == previous
    assert path.with_name(path.name + '.pending').exists()


def test_slow_retry_encoding_cannot_start_replace_after_deadline(tmp_path, monkeypatch):
    record, path, clock = _publish_probe(tmp_path, monkeypatch)
    before = path.read_bytes()
    calls = []
    error = _windows_publish_error(5)
    dump = receipt.json.dump
    def slow_retry(data, stream, **kwargs):
        dump(data, stream, **kwargs)
        if calls:
            clock['now'] += 1.0
    def blocked(*args):
        calls.append(clock['now'])
        raise error
    monkeypatch.setattr(receipt.json, 'dump', slow_retry)
    monkeypatch.setattr(receipt.os, 'replace', blocked)
    with pytest.raises(OSError) as raised:
        record.write()
    assert raised.value is error
    assert calls == [0.0], f'replace started after retry deadline: {calls}'
    assert path.read_bytes() == before
    assert path.with_name(path.name + '.pending').exists()


@pytest.mark.parametrize('windows,code', [(True, 2), (True, None), (False, 5)])
def test_unrelated_publish_errors_raise_immediately(tmp_path, monkeypatch, windows, code):
    record, path, clock = _publish_probe(tmp_path, monkeypatch, windows=windows)
    calls = []
    previous = record.data['write_retries']
    error = _windows_publish_error(code)
    def invalid(*args):
        calls.append(clock['now'])
        raise error
    monkeypatch.setattr(receipt.os, 'replace', invalid)
    with pytest.raises(OSError) as raised:
        record.write()
    assert raised.value is error
    assert calls == [0.0], f'actual attempts={calls}, expected one immediate failure'
    assert clock['sleeps'] == [], f'non-retryable failure slept: {clock["sleeps"]}'
    assert record.data['write_retries'] == previous
    assert path.with_name(path.name + '.pending').exists()
