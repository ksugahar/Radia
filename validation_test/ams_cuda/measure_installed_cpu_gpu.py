"""Idle same-host installed AMS CPU1/CPU16/GPU1 nonlinear whole-API ABBA."""
from pathlib import Path
import argparse,contextlib,io,hashlib,importlib.util,json,os,socket,statistics,subprocess,sys,time
START=time.perf_counter()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(a):
    if a.cycles<2:raise ValueError('At least two ABBA cycles are required')
    j=a.job.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);eng=None
    record=dict(schema='radia.ams-installed-cpu-gpu-performance.v1',accepted=False,timing_accepted=False,host=socket.gethostname(),samples=[],cpu_blas_threads=1,gpu_host_setup_threads=1,
        timing_boundary='Fresh nonlinear mesh/space/operator/RHS/preconditioner/solve/result recovery per API call. GPU includes uploads and download; independent checks and output serialization excluded.',
        initialization='CPU reference and first initialized GPU API are recorded separately. Warm ABBA follows both; no OS/device cold-reset or fresh-process cold-performance claim.',
        memory_definition='50ms absolute process RSS/private samples and GPU free endpoints; includes loaded runtimes and validation/serialization around the timed API (GPU module-probe subprocess included), not an API-only peak or transient global GPU peak.')
    try:
        manifest=json.loads((j/'manifest.json').read_text())
        for name,digest in manifest.items():
            p=(j/name).resolve()
            if not p.is_relative_to(j) or sha(p)!=digest:raise RuntimeError('Input manifest mismatch: '+name)
        record['payload_sha256']=manifest
        owner=json.loads((j/'owner.json').read_text());checked=json.loads((j/'vol_check.json').read_text(encoding='utf-8-sig'))
        if not checked['passed'] or checked['input_sha256']!=sha(j/'model.vol') or owner['vol_sha256']!=sha(j/'model.vol') or checked['contract_sha256']!=sha(j/'cube-labels.json'):raise RuntimeError('Versioned mesh check mismatch')
        spec=importlib.util.spec_from_file_location('wheel_verifier',j/'verify_wheel.py');verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
        wheel=verifier.verify(a.wheel)
        installed=json.loads(subprocess.check_output([str(a.installed_python),'-I','-c','import json,radia_simulink;print(json.dumps(radia_simulink.verify()))'],text=True));verifier.verify_installed_identity(installed,wheel)
        payload=Path(installed['matlab_path']);record['wheel']=wheel;record['installed']=installed
        mex_hash=json.loads((payload.parent/'manifest.json').read_text())['files']['matlab/radia_cuda_ams.mexw64']
        os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL',OPENBLAS_NUM_THREADS='1')
        import netgen
        net=Path(netgen.__file__).parent;paths=[j/'runtime',net,net.parent/'ngsolve_openblas',Path(sys.prefix)/'bin',Path(sys.base_prefix)]
        os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH'];handles=[os.add_dll_directory(str(p)) for p in paths if p.is_dir()]
        import numpy as np,ngsolve as ng,matlab.engine
        from scipy.io import loadmat
        from measure_native_cpu_gpu import idle,Memory
        from nonlinear_resident import cpu_run
        from ams_tuning import checked_errors
        sys.path.insert(0,str(j/'cpu'));import _radia_ams_cpu_reference as bridge
        receipt_path=j/'bridge-provenance/build.json';receipt=json.loads(receipt_path.read_text())
        if sha(bridge.__file__)!=receipt['sha256']['build-ams-cuda/_radia_ams_cpu_reference.pyd']:raise RuntimeError('CPU bridge differs from verified build receipt')
        for name,digest in receipt['source_sha256'].items():
            p=(j/'bridge-provenance/source'/name).resolve()
            if not p.is_relative_to(j/'bridge-provenance/source') or sha(p)!=digest:raise RuntimeError('CPU reference source identity mismatch')
        record['cpu_reference_build']=dict(source_commit=receipt['source_commit'],source_dirty=receipt['source_dirty'],receipt_sha256=sha(receipt_path),verified_source_count=len(receipt['source_sha256']))
        record['versions']=dict(python=sys.version,ngsolve=ng.__version__,numpy=np.__version__,logical_cpus=os.cpu_count(),cpu_bridge_sha256=sha(bridge.__file__),mesh_sha256=sha(j/'model.vol'),runner_sha256=sha(__file__))
        record['cpu_hardware']=json.loads(subprocess.check_output(['powershell','-NoProfile','-Command','Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json'],text=True))
        record['gpu_hardware']=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader,nounits'],text=True).strip()
        record['cpu_diagnostics']='Python Newton iteration messages captured in memory inside the API timer; no diagnostic disk I/O in timing.'
        record['process_init_including_input_verification_seconds']=time.perf_counter()-START;record['idle_before']=idle()
        def cpu(threads):
            ng.SetNumThreads(threads);timer=time.perf_counter()
            mesh=ng.Mesh(str(j/'model.vol'));fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT();b2=ng.curl(u)*ng.curl(u)
            form=ng.BilinearForm(fes);form+=ng.Variation((.5*b2+b2*b2+.5*u*u)*ng.dx)
            rhs=ng.LinearForm(fes);rhs+=ng.CoefficientFunction((20,0,0))*v*ng.dx
            with ng.TaskManager():rhs.Assemble()
            free=np.asarray(list(fes.FreeDofs()),dtype=bool);rhs.vec.FV().NumPy()[~free]=0
            gradient,_=fes.CreateGradient();xyz=np.array([[p[0],p[1],p[2]] for p in mesh.ngmesh.Points()])
            with contextlib.redirect_stdout(io.StringIO()):
                value,info=cpu_run(3,fes,form,rhs.vec,gradient,[xyz[:,c] for c in range(3)],bridge,threads)
            seconds=time.perf_counter()-timer
            info.update(ngsolve_threads=threads,ams_setup_coarse_factor_threads=1,whole_api_seconds=seconds)
            return value,info,(fes,form,rhs,free)
        reference,info,validation=cpu(1);record['first_initialized_cpu_api']=info
        fes,form,rhs,free=validation
        def validate(value):
            ng.SetNumThreads(1);probe=ng.GridFunction(fes);probe.vec.FV().NumPy()[:]=np.asarray(value).ravel();applied=probe.vec.CreateVector()
            with ng.TaskManager():form.Apply(probe.vec,applied)
            return checked_errors(value,reference,applied.FV().NumPy(),rhs.vec.FV().NumPy(),free,1e-9)
        record['cpu_reference_check']=validate(reference)
        timer=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');record['matlab_start_seconds']=time.perf_counter()-timer
        eng.restoredefaultpath(nargout=0);eng.addpath(str(j),nargout=0);pid=int(eng.feature('getpid'))
        eng.eval('timer=tic;g=gpuDevice;wait(g);gpu_init_seconds=toc(timer);',nargout=0);record['gpu_init_seconds']=float(eng.workspace['gpu_init_seconds'])
        def gpu(path):
            eng.measure_installed_cpu_gpu(str(j),str(payload),str(path),nargout=0)
            data=loadmat(path,simplify_cells=True);module_path=path.with_suffix('.module.json')
            subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(j/'ams_loaded_module.ps1'),'-ProcessId',str(pid),'-OutputPath',str(module_path)],check=True,capture_output=True)
            module=json.loads(module_path.read_text(encoding='utf-8-sig'));expected=payload/'radia_cuda_ams.mexw64'
            if Path(module['path']).resolve()!=expected.resolve() or module['sha256']!=mex_hash or sha(expected)!=mex_hash:raise RuntimeError('Actual installed AMS module mismatch')
            return data['x'],dict(data['report'],native_module=module)
        idle();value,first=gpu(out/'first-gpu.mat');record['first_initialized_gpu_api']=dict(first,independent=validate(value))
        threads_list=[1]
        if (os.cpu_count() or 1)>=16:threads_list.append(16)
        else:record['cpu16_unavailable']='Fewer than sixteen logical CPUs reported'
        for threads in threads_list:
            for cycle in range(a.cycles):
                for route in ('cpu','gpu','gpu','cpu'):
                    guard=idle();index=len(record['samples']);path=out/f'{threads}-{index}-{route}.mat'
                    with Memory(os.getpid() if route=='cpu' else pid) as memory:
                        if route=='cpu':value,report,_=cpu(threads)
                        else:value,report=gpu(path)
                    seconds=float(report['whole_api_seconds'])
                    if not np.isfinite(seconds) or seconds<=0:raise RuntimeError('Invalid timing')
                    row=dict(cpu_baseline_ngsolve_threads=threads,gpu_host_setup_threads=1,route=route,cycle=cycle,seconds=seconds,report=report,independent=validate(value),host_memory=memory.result(),idle_before=guard)
                    record['samples'].append(row);print(threads,route,seconds,flush=True)
        record['summary']={}
        for threads in threads_list:
            rows=[s for s in record['samples'] if s['cpu_baseline_ngsolve_threads']==threads]
            cpu_median=statistics.median(s['seconds'] for s in rows if s['route']=='cpu');gpu_median=statistics.median(s['seconds'] for s in rows if s['route']=='gpu')
            record['summary'][str(threads)]=dict(cpu_median_seconds=cpu_median,gpu_median_seconds=gpu_median,cpu_over_gpu=cpu_median/gpu_median)
        record['idle_after']=idle()
        def require_finite(value):
            if isinstance(value,dict):
                for v in value.values():require_finite(v)
            elif isinstance(value,(list,tuple)):
                for v in value:require_finite(v)
            elif isinstance(value,np.ndarray):require_finite(value.tolist())
            elif isinstance(value,(int,float,np.number)) and not np.isfinite(value):raise RuntimeError('Nonfinite performance evidence')
        require_finite(record)
        record.update(accepted=True,timing_accepted=socket.gethostname().lower()=='ngsolve')
    except BaseException as e:record['error']=repr(e);raise
    finally:
        if eng is not None:
            try:eng.quit()
            except BaseException as e:record['engine_release_error']=repr(e);record['accepted']=False;record['timing_accepted']=False
        (out/'result.json').write_text(json.dumps(record,indent=2,default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v))+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('job','wheel','installed-python','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--cycles',type=int,default=2);main(p.parse_args())
