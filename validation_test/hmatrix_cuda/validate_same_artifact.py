"""Validate an existing 100-built artifact using independent CPU validation inputs.
No build/install and no Python solver snapshot prerequisite for GPU users.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,time,os

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);p.add_argument('--profile',action='store_true');a=p.parse_args();j=a.job
    record=dict(schema='radia.hmatrix-same-artifact.v1',accepted=False,scope='Existing numerical/operator and bounded physical material subgates; not general HDiv-MMM acceptance',true_residual=None,true_residual_reason='Per-case original-equation residuals are reported in material_solver, not a single suite residual',performance_acceptance=False)
    eng=None
    try:
        manifest=json.loads((j/'manifest.json').read_text(encoding='utf-8'))
        for name,h in manifest.items():
            if sha(j/name)!=h:raise RuntimeError('Staged artifact/input hash mismatch: '+name)
        record['sha256']=manifest
        if (j/'runtime').is_dir():os.environ['PATH']=str(j/'runtime')+';'+os.environ['PATH']
        os.environ['MKL_NUM_THREADS']='1';os.environ['MKL_THREADING_LAYER']='SEQUENTIAL'
        gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.free,utilization.gpu','--format=csv,noheader'],text=True).strip()
        record['gpu_before']=gpu
        fields=gpu.split(',')
        if len(fields)!=5 or int(fields[-1].strip().split()[0])>5:raise RuntimeError('GPU is busy or ambiguous')
        import matlab.engine
        started=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');record['engine_start_seconds']=time.perf_counter()-started
        pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.addpath(str(j),str(j/'matlab'),nargout=0);record['matlab']=eng.version()
        started=time.perf_counter();eng.checkResidentHMatrix(str(j/'inputs.mat'),str(j/'matlab.json'),nargout=0);record['suite_seconds']=time.perf_counter()-started
        result=json.loads((j/'matlab.json').read_text(encoding='utf-8'));assert result['accepted'];record['result']=result
        for name,h in manifest.items():assert sha(j/name)==h
        if a.profile:
            import ctypes,statistics
            def idle():
                def cpu_ticks():
                    idle=ctypes.c_ulonglong();kernel=ctypes.c_ulonglong();user=ctypes.c_ulonglong()
                    if not ctypes.windll.kernel32.GetSystemTimes(ctypes.byref(idle),ctypes.byref(kernel),ctypes.byref(user)):raise OSError('GetSystemTimes failed')
                    return idle.value,kernel.value+user.value
                before=cpu_ticks();time.sleep(2);after=cpu_ticks()
                total=after[1]-before[1]
                if total<=0:raise RuntimeError('CPU sample unavailable')
                cpu=100*(1-(after[0]-before[0])/total)
                class Memory(ctypes.Structure):
                    _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(x,ctypes.c_ulonglong) for x in ['total','available','total_page','available_page','total_virtual','available_virtual','extended']]
                memory=Memory();memory.length=ctypes.sizeof(memory)
                if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory)):raise OSError('GlobalMemoryStatusEx failed')
                gpu=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu','--format=csv,noheader,nounits'],text=True).strip()
                if cpu>=10 or int(gpu)>5 or memory.available<8*2**30:raise RuntimeError('Idle timing budget unavailable')
                return dict(cpu_percent=cpu,gpu_percent=int(gpu),free_host_bytes=memory.available)
            record['idle_before']=idle()
            eng.profileResidentHMatrix(str(j/'inputs.mat'),str(j/'profile.json'),nargout=0)
            record['idle_after']=idle()
            profile=json.loads((j/'profile.json').read_text(encoding='utf-8'));assert profile['accepted']
            record['operator_profile']=profile
            record['median_operator_whole_seconds']={mode:statistics.median(x['whole_seconds'] for x in profile['samples'] if x['mode']==mode) for mode in ('rebuild','resident')}
        record['accepted']=True
        record['timing_limitations']='Suite duration includes assertions and reference comparisons; not transfer-inclusive whole-solver performance'
        record['memory_limitations']='Exact owned device counters only; no global/process peak'
    except BaseException as e:record['error']=f'{type(e).__name__}: {e}';raise
    finally:
        if eng is not None:eng.quit()
        (j/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
