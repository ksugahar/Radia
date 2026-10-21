"""Observe an owned MATLAB process with a separately staged Frida wheel.

No system service, shared installation, API replacement or CUDA probing is used.
This is a trace collector, not a blanket proof of CUDA API coverage. Opaque
export-table calls, hook errors or missing positive controls stay inconclusive.
"""
import argparse
import ctypes as c
import hashlib
import json
import os
from pathlib import Path
import sys
import threading
import time


def assert_same_job(pid):
    k=c.WinDLL('kernel32',use_last_error=True)
    k.OpenJobObjectW.argtypes=[c.c_uint32,c.c_int,c.c_wchar_p];k.OpenJobObjectW.restype=c.c_void_p
    k.OpenProcess.argtypes=[c.c_uint32,c.c_int,c.c_uint32];k.OpenProcess.restype=c.c_void_p
    k.AssignProcessToJobObject.argtypes=[c.c_void_p,c.c_void_p];k.AssignProcessToJobObject.restype=c.c_int
    k.QueryInformationJobObject.argtypes=[c.c_void_p,c.c_int,c.c_void_p,c.c_uint32,c.c_void_p]
    k.QueryInformationJobObject.restype=c.c_int
    k.CloseHandle.argtypes=[c.c_void_p]
    name=os.environ.get('RADIA_OWNED_JOB_NAME')
    if not name or not name.startswith('Local\\RadiaCsrTrace-'):
        raise RuntimeError('Named owned bounded supervisor job required')
    job=k.OpenJobObjectW(0x0001|0x0004,False,name)
    if not job:raise c.WinError(c.get_last_error())
    process=None
    try:
        data=c.create_string_buffer(8+8*256)
        def members():
            if not k.QueryInformationJobObject(job,3,data,len(data),None):raise c.WinError(c.get_last_error())
            count=c.c_uint32.from_buffer(data,4).value
            return list((c.c_size_t*count).from_buffer(data,8))
        if os.getpid() not in members():raise RuntimeError('Controller does not belong to the declared owned job')
        if pid not in members():
            # Frida's helper can spawn outside the caller's job. The target is
            # still suspended; join the same existing limits before any resume.
            process=k.OpenProcess(0x0100|0x0001,False,pid)
            if not process or not k.AssignProcessToJobObject(job,process):raise c.WinError(c.get_last_error())
        if pid not in members():raise RuntimeError('Trace target is not in the bounded owned job')
    finally:
        if process:k.CloseHandle(process)
        k.CloseHandle(job)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tools',type=Path,required=True)
    parser.add_argument('--matlab-root',type=Path,required=True)
    parser.add_argument('--workload',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--native-binary',type=Path,required=True)
    parser.add_argument('--backend',choices=('cpu','gpu','observer','visibility'),required=True)
    parser.add_argument('--cuda-visible-devices')
    args=parser.parse_args()
    if args.cuda_visible_devices is not None:os.environ['CUDA_VISIBLE_DEVICES']=args.cuda_visible_devices
    if args.output.exists():raise RuntimeError('Never overwrite a previous trace')
    scratch=args.output.parent/'trace-temp';scratch.mkdir(exist_ok=True)
    os.environ['TEMP']=os.environ['TMP']=str(scratch)
    sys.path.insert(0,str(args.tools))
    import frida
    device=frida.get_local_device()
    events=[];errors=[];output=[];finished=threading.Event();pid=None;session=None
    device.on('output',lambda process,fd,data:output.append(dict(pid=process,fd=fd,text=data.decode('utf-8',errors='replace'))))
    result=dict(schema='radia.csr-context-trace.v1',accepted=False,backend=args.backend,
                observer='Frida',observer_version=frida.__version__,dropped_events=0,
                context_gate='inconclusive',source='trace_cuda_context.js',cuda_visible_devices=args.cuda_visible_devices)
    def message(value,data):
        if value.get('type')=='send':
            if len(events)<100000:events.append(value['payload'])
            else:result['dropped_events']+=1
        else:errors.append(value)
    result['native_binary_sha256']=hashlib.sha256(args.native_binary.read_bytes()).hexdigest()
    result['collector_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    started=time.monotonic()
    try:
        executable=args.matlab_root/'bin/win64/MATLAB.exe'
        command=[str(executable),'-nodesktop','-nosplash','-singleCompThread','-batch',
                 "run('"+str(args.workload).replace("'","''")+"')"]
        pid=device.spawn(command,stdio='pipe')
        result['target_pid']=pid
        assert_same_job(pid)
        session=device.attach(pid)
        def detached(reason,crash=None):
            result['detach_reason']=reason
            result['crash_report']=getattr(crash,'report',None)
            finished.set()
        session.on('detached',detached)
        script_path=Path(__file__).with_suffix('.js')
        text=script_path.read_text(encoding='utf-8')
        result['script_sha256']=hashlib.sha256(text.encode()).hexdigest()
        script=session.create_script(text);script.on('message',message);script.load()
        if not any(e.get('kind')=='ready' for e in events):
            # Messages are asynchronous; readiness is required before resume.
            for _ in range(50):
                if any(e.get('kind')=='ready' for e in events):break
                time.sleep(.02)
        if not any(e.get('kind')=='ready' for e in events):raise RuntimeError('Observer not ready before resume')
        device.resume(pid)
        if not finished.wait(180):raise TimeoutError('Owned trace target exceeded internal watchdog')
        context_calls=[e for e in events if e.get('kind')=='context_api']
        opaque=[e for e in events if e.get('name')=='cuGetExportTable' and e.get('kind')=='observed_api']
        coverage_errors=[e for e in events if e.get('kind')=='coverage_error']
        loaded=any('radia_cuda_csr_pcg' in e.get('name','').lower() and Path(e['path']).samefile(args.native_binary) for e in events if e.get('kind')=='module')
        if hashlib.sha256(args.native_binary.read_bytes()).hexdigest()!=result['native_binary_sha256']:
            raise RuntimeError('Native binary changed during tracing')
        expected=args.workload.with_suffix('.result.json')
        numerical=json.loads(expected.read_text()) if expected.is_file() else None
        complete=(loaded or args.backend in ('observer','visibility')) and not errors and not coverage_errors and not result['dropped_events'] and numerical and numerical.get('accepted') is True
        result.update(numerical=numerical,context_calls=len(context_calls),opaque_export_tables=len(opaque),coverage_errors=coverage_errors)
        if complete and args.backend=='visibility' and numerical.get('visible_device_count')==0:
            result.update(accepted=True,context_gate='process-visible-device-count-zero')
        elif complete and args.backend=='gpu' and context_calls:
            result.update(accepted=True,context_gate='positive-control-detected')
        elif complete and args.backend=='observer' and not context_calls and not opaque:
            result.update(accepted=True,context_gate='observer-control-no-observed-context')
        elif complete and args.backend=='cpu' and not context_calls and not opaque:
            result.update(accepted=True,context_gate='no-observed-context-API-in-measured-interval')
        # The result is scoped to exports/dynamically resolved APIs. Pair it with
        # a matching GPU positive control before any overall CPU-context verdict.
    except BaseException as error:
        result['error']=repr(error)
        raise
    finally:
        if pid is not None and not finished.is_set():
            try:device.kill(pid)
            except Exception as error:result['owned_target_cleanup_error']=repr(error)
        if session is not None:
            try:session.detach()
            except Exception:pass
        result.update(events=events,observer_errors=errors,target_output=output,elapsed_seconds=time.monotonic()-started)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
