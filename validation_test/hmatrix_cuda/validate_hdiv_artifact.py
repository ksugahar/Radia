"""Same-artifact native HDiv validation; Python only orchestrates independent fixtures."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,statistics,subprocess,sys,time

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def idle():
    def ticks():
        a,b,c=(ctypes.c_ulonglong() for _ in range(3))
        if not ctypes.windll.kernel32.GetSystemTimes(ctypes.byref(a),ctypes.byref(b),ctypes.byref(c)):raise OSError('CPU sample failed')
        return a.value,b.value+c.value
    a=ticks();time.sleep(2);b=ticks();cpu=100*(1-(b[0]-a[0])/(b[1]-a[1]))
    gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,utilization.gpu,memory.free,memory.total','--format=csv,noheader,nounits'],text=True).strip()
    fields=gpu.split(',')
    class Memory(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(x,ctypes.c_ulonglong) for x in ['total','available','total_page','available_page','total_virtual','available_virtual','extended']]
    memory=Memory();memory.length=ctypes.sizeof(memory)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory)):raise OSError('Host memory sample failed')
    if cpu>10 or int(fields[2])>5 or memory.available<8*2**30 or float(fields[3])<max(2048,.2*float(fields[4])):raise RuntimeError(f'Idle GPU timing budget unavailable: CPU={cpu}, GPU={gpu}, host_free_bytes={memory.available}')
    return dict(cpu_percent=cpu,gpu=gpu,free_host_bytes=memory.available,limitation='System samples do not exclude every instantaneous external load')

def main():
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);p.add_argument('--warm',type=int,default=2);a=p.parse_args();j=a.job.resolve();eng=None
    r=dict(accepted=False,release_accepted=False,general_hdiv_validation_complete=False,scope='Bounded public affine tetra and P2 sphere, existing near-rule dense Gram, linear material',cases=[],timing_classification='Exact host diagnostic, not idle-mdx performance acceptance',memory_limitations='Exact owned buffers and available-memory samples; no cuBLAS, host or global allocation peak')
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        manifest=json.loads((j/'manifest.json').read_text());r['sha256']=manifest
        for name,h in manifest.items():
            if sha(j/name)!=h:raise RuntimeError('Staged hash mismatch: '+name)
        os.environ.update(MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL')
        prefix=Path(sys.prefix);site=prefix/'Lib/site-packages';os.environ['PATH']=';'.join(map(str,[j/'runtime',prefix/'Scripts',prefix/'Library/bin',prefix/'bin',site/'netgen',site/'netgen_mesher.libs',site/'ngsolve.libs',site/'ngsolve_openblas',Path(sys.base_prefix)]))+';'+os.environ['PATH']
        r['idle_before']=idle()
        import matlab.engine
        start=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');r['owned_engine_start_seconds']=time.perf_counter()-start;r['matlab']=eng.version()
        eng.addpath(str(j),str(j/'matlab'),nargout=0)
        pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        from scipy.io import loadmat
        for case in json.loads((j/'cases.json').read_text()):
            row=dict(case=case['name'],samples=[]);r['cases'].append(row)
            eng.workspace['mesh_path']=str(j/case['mesh']);eng.workspace['fixture_path']=str(j/case['fixture']);eng.workspace['degree']=float(case['order'])
            eng.eval('cpu=load(fixture_path);',nargout=0)
            eng.workspace['assembly_mode']=case.get('assembly_mode','direct')
            for k in range(a.warm+1):
                idle();output=j/(case['name']+f'-{k}.mat');eng.workspace['output_path']=str(output)
                eng.eval("d=gpuDevice;free_before=d.AvailableMemory;solution=radia.gpu.solveHDivElements(mesh_path,degree,[1 .2 .1],2.5,1e-11,1000,assembly_mode);free_after=d.AvailableMemory;coefficient_error=norm(solution.coefficients-cpu.expected)/norm(cpu.expected);physical_error=norm(solution.mean_magnetization-cpu.mean_operator*cpu.expected)/norm(cpu.mean_operator*cpu.expected);independent_true_residual=norm(cpu.rhs-cpu.system*solution.coefficients)/norm(cpu.rhs);save(output_path,'solution','free_before','free_after','coefficient_error','physical_error','independent_true_residual','-v7');assert(max([coefficient_error physical_error])<1e-8);assert(independent_true_residual<1e-10);assert(radia_cuda_hmatrix('live_count')==0);",nargout=0)
                raw=loadmat(output,simplify_cells=True);solution=raw['solution'];sample={key:float(raw[key]) for key in ['coefficient_error','physical_error','independent_true_residual','free_before','free_after']};sample.update({key:solution[key] for key in ['whole_seconds','iterations','true_relative_residual','operator_owned_device_bytes','solver_owned_device_bytes','setup_peak_owned_device_bytes','combined_peak_owned_device_bytes']});sample['phase']='first_payload' if k==0 else 'warm';row['samples'].append(sample)
            row['warm_median_whole_seconds']=statistics.median(s['whole_seconds'] for s in row['samples'][1:])
        r['idle_after']=idle()
        for name,h in manifest.items():
            if sha(j/name)!=h:raise RuntimeError('Artifact changed during validation')
        r['accepted']=True
    except BaseException as e:r['error']=f'{type(e).__name__}: {e}';raise
    finally:
        if eng is not None:
            try:eng.quit()
            except BaseException as e:r['release_error']=str(e)
        (j/'result.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':main()
