"""Guards: receipt-before-GUI, complete owned-tree cleanup, streaming eligibility.
Assumes: Windows native ownership tests run on the supported viewer host.
"""
import asyncio
import ctypes as c
from ctypes import wintypes as w
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace as NS

import pytest
from radia_mcp.mesh_viewer import mcp as adapter, post_view
from radia_mcp.mesh_viewer._worker_process import run_worker, transaction_deadline
from radia_mcp.gmsh.runtime import require_msh42


@pytest.mark.parametrize('value', ['nan', 'inf', '29', '901', 'invalid'])
def test_invalid_deadline_does_not_launch(value, monkeypatch):
    monkeypatch.setenv('RADIA_MESH_VIEWER_TIMEOUT_SECONDS', value)
    with pytest.raises(ValueError, match='30..900'): transaction_deadline()


def test_post_receipt_precedes_only_one_msh_gui(tmp_path, monkeypatch):
    source=tmp_path/'input.vol';source.write_bytes(b'source')
    destination=tmp_path/'post';destination.mkdir()
    output=destination/'verified.msh'; events=[]
    digest=post_view.sha256(source)
    receipt={}
    async def conversion(argv, **kwargs):
        events.append('conversion');await asyncio.sleep(.01)
        output.write_bytes(b'mesh')
        receipt.update(schema='radia.mesh-post-view.v1',output=str(output),
            sha256=post_view.sha256(output),source_sha256=digest,group_scope='materials',
            solver_input=False,staging='staging',physical_groups=[[3,'iron']])
        Path(argv[-1]).write_text(json.dumps(receipt),encoding='utf-8')
        events.append('verified')
    async def preflight(argv, **kwargs):
        assert events==['conversion','verified']
        assert '--post-output' not in argv and str(output) in argv
        assert kwargs['env']['PYTHONIOENCODING']=='utf-8' and kwargs['env']['PYTHONUTF8']=='1'
        assert argv[argv.index('--expected-sha256')+1]==post_view.sha256(output)
        events.append('preflight');return NS(stdout=b'{}')
    def gui(argv, **kwargs):events.append('GUI');return NS(pid=123)
    monkeypatch.setattr(adapter,'require_interactive_desktop',lambda:{})
    monkeypatch.setattr(post_view,'run_worker',conversion)
    monkeypatch.setattr(adapter,'run_worker',preflight)
    monkeypatch.setattr(adapter.subprocess,'Popen',gui)
    monkeypatch.setattr(adapter.subprocess,'DETACHED_PROCESS',8,raising=False)
    monkeypatch.setattr(adapter.subprocess,'CREATE_NEW_PROCESS_GROUP',0x200,raising=False)
    async def scenario():
        task=asyncio.create_task(adapter.open_mesh_viewer(str(source),str(output),'check',digest,'materials'))
        await asyncio.sleep(.001);assert not task.done()
        return await task
    result=asyncio.run(scenario())
    assert events==['conversion','verified','preflight','GUI']
    assert result['status']=='launched'
    assert result['post_view']=={**receipt,'worker_receipt':result['post_view']['worker_receipt']}
    assert Path(result['post_view']['worker_receipt']).is_file()
    assert result['post_view']['sha256']==post_view.sha256(output)


@pytest.mark.parametrize('broken', [
    'schema','output','source_sha256','group_scope','solver_input','sha256','missing target'])
def test_real_receipt_binding_rejects_wrong_worker_result(tmp_path,monkeypatch,broken):
    source=tmp_path/'input.vol';source.write_bytes(b'source')
    destination=tmp_path/'post';destination.mkdir()
    output=destination/'view.msh'
    digest=post_view.sha256(source)
    async def worker(argv,**kwargs):
        if broken!='missing target':output.write_bytes(b'mesh')
        receipt={'schema':'radia.mesh-post-view.v1','output':str(output),
            'source_sha256':digest,'group_scope':'materials','solver_input':False,
            'sha256':'a'*64,'staging':'staging','physical_groups':[[3,'iron']]}
        changes={'schema':'other.schema','output':str(destination/'other.msh'),
            'source_sha256':'0'*64,'group_scope':'boundaries','solver_input':True,
            'sha256':'not-a-sha256'}
        if broken in changes:receipt[broken]=changes[broken]
        Path(argv[-1]).write_text(json.dumps(receipt),encoding='utf-8')
    async def forbidden_preflight(*args,**kwargs):pytest.fail('Preflight before receipt binding')
    monkeypatch.setattr(post_view,'run_worker',worker)
    monkeypatch.setattr(adapter,'run_worker',forbidden_preflight)
    monkeypatch.setattr(adapter,'require_interactive_desktop',lambda:{})
    monkeypatch.setattr(adapter.subprocess,'Popen',lambda *a,**k:pytest.fail('GUI before receipt binding'))
    # The real prepare_post_view_async must reject the worker's otherwise-successful receipt.
    with pytest.raises(RuntimeError,match='Unbound conversion receipt'):
        asyncio.run(adapter.open_mesh_viewer(str(source),str(output),'check',digest,'materials'))


def test_launcher_check_rejects_mismatching_expected_sha(tmp_path,monkeypatch,capsys):
    from radia_mcp.mesh_viewer import launcher
    target=tmp_path/'view.msh'
    target.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n',encoding='utf-8')
    monkeypatch.setattr(launcher,'require_interactive_desktop',lambda:{})
    monkeypatch.setattr(launcher,'runtime_identity',lambda path:{})
    monkeypatch.setattr(launcher,'open_mesh',lambda path:pytest.fail('GUI during --check'))
    assert launcher.main([str(target),'--check','--expected-sha256','f'*64])!=0
    assert 'Verified MSH hash changed' in capsys.readouterr().err


@pytest.mark.parametrize('gate',['inventory','labels','physical groups','element count','source hash','receipt binding'])
def test_conversion_gate_failure_never_launches(tmp_path,monkeypatch,gate):
    source=tmp_path/'input.vol';source.write_bytes(b'source')
    async def fail(*args,**kwargs):raise RuntimeError('conversion '+gate)
    monkeypatch.setattr(adapter,'require_interactive_desktop',lambda:{})
    monkeypatch.setattr(adapter,'prepare_post_view_async',fail)
    monkeypatch.setattr(adapter.subprocess,'Popen',lambda *a,**k:pytest.fail('GUI before gate success'))
    with pytest.raises(RuntimeError,match=gate):
        asyncio.run(adapter.open_mesh_viewer(str(source),str(tmp_path/'out.msh'),'check','b'*64,'materials'))


def test_native_eligibility_streams_and_checks_late_fields(tmp_path,monkeypatch):
    path=tmp_path/'field.msh'
    header='$MeshFormat\n4.2 0 8\n$EndMeshFormat\n'
    field='$NodeData\n1\n"field"\n1\n0\n3\n0\n1\n1\n1 1 {tag} {tag}\n1\n{tag}\n7\n$EndNodeData\n'
    path.write_text(header+(' '*1023+'\n')*512+field.format(tag=1))
    monkeypatch.setattr(Path,'read_text',lambda *a,**k:pytest.fail('whole-file read'))
    require_msh42(path)
    with path.open('a') as out:out.write(field.format(tag=2**31+17))
    with pytest.raises(ValueError,match='2147483647'):require_msh42(path)


def test_streaming_duplicate_tags_still_fail(tmp_path):
    path=tmp_path/'duplicates.msh'
    path.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$NodeData\n'
                    '1\n"field"\n1\n0\n3\n0\n1\n2\n1 2 1 1\n2\n1 1\n7 8\n$EndNodeData\n')
    with pytest.raises(ValueError,match='duplicate'):require_msh42(path)


def _alive(pid):
    k=c.WinDLL('kernel32',use_last_error=True)
    k.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];k.OpenProcess.restype=w.HANDLE
    k.WaitForSingleObject.argtypes=[w.HANDLE,w.DWORD];k.CloseHandle.argtypes=[w.HANDLE]
    handle=k.OpenProcess(0x100000|0x1000,False,pid)
    if not handle:return False
    try:return k.WaitForSingleObject(handle,0)==258
    finally:k.CloseHandle(handle)


@pytest.mark.skipif(sys.platform!='win32',reason='Windows process-tree contract')
@pytest.mark.parametrize('mode',['timeout','cancel'])
def test_owned_grandchild_dies_and_unrelated_sentinel_survives(tmp_path,mode):
    pidfile=tmp_path/'pids.json'
    script=("import subprocess,sys,json,time,os; "
            "p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)']); "
            "open(sys.argv[1],'w').write(json.dumps([os.getpid(),p.pid]));time.sleep(60)")
    sentinel=subprocess.Popen([str(Path(sys.base_prefix)/'python.exe'),'-c','import time;time.sleep(60)'],
                              creationflags=subprocess.CREATE_NO_WINDOW)
    async def scenario():
        task=asyncio.create_task(run_worker([sys.executable,'-c',script,str(pidfile)],
            env=dict(os.environ),deadline=time.monotonic()+10,phase='ownership test'))
        end=time.monotonic()+3
        while not pidfile.exists():
            if task.done():await task
            if time.monotonic()>end:raise RuntimeError('worker did not reach readiness')
            await asyncio.sleep(.01)
        if mode=='cancel':
            task.cancel()
            with pytest.raises(asyncio.CancelledError):await task
        else:
            with pytest.raises(RuntimeError,match='timed out'):await task
    try:
        asyncio.run(scenario())
        pids=json.loads(pidfile.read_text())
        assert all(not _alive(pid) for pid in pids),pids
        assert sentinel.poll() is None
    finally:
        sentinel.terminate();sentinel.wait(timeout=5)


@pytest.mark.skipif(sys.platform!='win32',reason='Windows worker encoding contract')
def test_real_outer_japanese_error_is_not_replaced(monkeypatch):
    monkeypatch.delenv('PYTHONIOENCODING',raising=False)
    monkeypatch.delenv('PYTHONUTF8',raising=False)
    text='指定されたモジュールが見つかりません。'
    with pytest.raises(RuntimeError) as error:
        asyncio.run(run_worker([sys.executable,'-c','import sys;print('+repr(text)+',file=sys.stderr);sys.exit(1)'],
            env=dict(os.environ),deadline=time.monotonic()+30,phase='outer encoding'))
    assert text in str(error.value) and '\ufffd' not in str(error.value)


@pytest.mark.skipif(sys.platform!='win32',reason='Windows suspended assignment contract')
def test_failed_assignment_does_not_run_child(tmp_path,monkeypatch):
    from radia_mcp.mesh_viewer import _worker_process as workers
    marker=tmp_path/'must-not-exist'
    def fail(self,pid):raise RuntimeError('assignment failure')
    monkeypatch.setattr(workers._Job,'attach_resume',fail)
    with pytest.raises(RuntimeError,match='assignment failure'):
        asyncio.run(run_worker([sys.executable,'-c','from pathlib import Path;Path('+repr(str(marker))+').touch()'],
                              env=dict(os.environ),deadline=time.monotonic()+30,phase='assignment'))
    assert not marker.exists()
