"""Same installed wheel: linear HCurl AMS/IC CPU1/16 versus GPU AMS scaling.

GPU users require only MATLAB/mexcuda and the mesh. Python builds independent
CPU operators and references for validation, never GPU input snapshots.
"""
from pathlib import Path
import argparse,gc,hashlib,importlib.util,json,os,socket,statistics,subprocess,sys,time

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def finite(value):
    import numpy as np
    if isinstance(value,dict):
        for v in value.values():finite(v)
    elif isinstance(value,(list,tuple,np.ndarray)):
        for v in value:finite(v)
    elif isinstance(value,(float,int,np.number)) and not np.isfinite(value):raise RuntimeError('Nonfinite evidence')

def main(a):
    if a.cycles<2:raise ValueError('At least two balanced cycles required')
    job=a.job.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);eng=None
    r=dict(schema='radia.ams-ic-installed-scaling.v1',accepted=False,timing_accepted=False,host=socket.gethostname(),cases=[],
        formulation='Affine lowest-order HCurl nograds: curl(u)*curl(v)+u*v, J=(1,0,0) A/m^2, homogeneous tangential boundary, SI',
        comparison='Independent Python CPU AMS and IC at NGSolve 1/16 threads versus installed MATLAB/mexcuda GPU AMS; GPU IC unavailable in this wheel',
        timing_boundary='Fresh mesh load, FE matrix/RHS construction, preconditioner setup and first zero-guess solve. Reuse changes only RHS scale, retaining matrix and preconditioner. GPU setup and solves include transfers and gathered result; validation and serialization excluded.',
        initialization='First initialized call per configuration is recorded separately and excluded from balanced medians; MATLAB/GPU startup separate; no OS/device cold reset claim',
        memory_definition='50ms absolute process RSS/private samples include surrounding validation/serialization and loaded runtimes; GPU native owned bounds and free-memory endpoints, not transient whole-device peak',
        threads_policy=dict(cpu_ngsolve=[1,16],cpu_blas=1,cpu_ams_coarse_factor=1,gpu_host_setup=1),ic_configuration=dict(shift=1.05,auto_shift=False,ordering='native IC default',initial_guess='zero'))
    try:
        if socket.gethostname().lower()!='ngsolve' and not a.diagnostic:raise RuntimeError('Formal timing requires selected A100 host')
        for name,digest in json.loads((job/'manifest.json').read_text()).items():
            p=(job/name).resolve()
            if not p.is_relative_to(job) or sha(p)!=digest:raise RuntimeError('Payload changed '+name)
        os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL')
        import netgen
        net=Path(netgen.__file__).parent
        paths=[job/'runtime',net,net.parent/'ngsolve_openblas',Path(sys.prefix)/'bin',Path(sys.base_prefix)]
        os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH'];dll_handles=[os.add_dll_directory(str(p)) for p in paths if p.is_dir()]
        import numpy as np,ngsolve as ng,matlab.engine
        from scipy.io import loadmat
        from measure_native_cpu_gpu import idle,Memory
        def quiet():
            if a.diagnostic:
                import psutil
                cpu=psutil.cpu_percent(interval=1);raw=subprocess.check_output(['nvidia-smi','--query-gpu=name,utilization.gpu,memory.free,memory.total','--format=csv,noheader,nounits'],text=True).strip()
                name,util,free,total=[v.strip() for v in raw.split(',')];available=psutil.virtual_memory().available
                if float(free)<max(2048,.2*float(total)) or available<8*2**30:raise RuntimeError('Numerical validation memory reserve failed')
                return dict(cpu_percent=cpu,gpu=name,gpu_utilization_percent=float(util),device_free_MiB=float(free),host_available_bytes=available,timing_accepted=False,settling_rejections=[])
            # A completed 16-thread call can leave a short utilization tail.
            # Keep the idle thresholds; record bounded settling, never time load.
            failed=[]
            for attempt in range(5):
                try:return dict(idle(),settling_rejections=failed)
                except RuntimeError as ex:
                    if not str(ex).startswith('Idle measurement guard failed:'):raise
                    failed.append(str(ex))
                    if attempt==4:raise
                    time.sleep(1)
        sys.path.insert(0,str(job/'cpu'));import _radia_ams_cpu_reference as bridge
        import sparsesolv_ngsolve as ic_native
        receipt=json.loads((job/'bridge-provenance/build.json').read_text())
        if sha(bridge.__file__)!=receipt['sha256']['build-ams-cuda/_radia_ams_cpu_reference.pyd']:raise RuntimeError('AMS bridge identity')
        for name,digest in receipt['source_sha256'].items():
            p=(job/'bridge-provenance/source'/name).resolve()
            if not p.is_relative_to(job/'bridge-provenance/source') or sha(p)!=digest:raise RuntimeError('AMS source receipt')
        owner=json.loads((job/'owner.json').read_text())
        if sha(ic_native.__file__)!=owner['cpu_ic_origin']['sha256']:raise RuntimeError('IC native identity')
        spec=importlib.util.spec_from_file_location('wheel_verifier',job/'verify_wheel.py');verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
        wheel=verifier.verify(a.wheel);import radia_gpu
        installed=radia_gpu.verify();verifier.verify_installed_identity(installed,wheel)
        if wheel['version']!='0.1.4' or wheel['wheel_sha256']!=owner['wheel_sha256']:raise RuntimeError('Selected stable wheel identity')
        payload=Path(installed['matlab_path']);mex_hash=json.loads((payload.parent/'manifest.json').read_text())['files']['matlab/radia_cuda_ams.mexw64']
        ic_hash=None
        if a.ic_directory:
            a.ic_directory=a.ic_directory.resolve();ic_build=json.loads((job/'ic-build.json').read_text());ic_hash=ic_build['binary_sha256_by_file']['radia_cuda_ic.mexw64']
            if sha(a.ic_directory/'radia_cuda_ic.mexw64')!=ic_hash:raise RuntimeError('IC candidate build identity mismatch')
            r['ic_candidate']=dict(build_source_commit=ic_build['source_commit'],binary_sha256=ic_hash,source_sha256=ic_build['source_sha256_by_file']['radia_cuda_ic.cu'],distribution_status='Experimental candidate separate from unchanged stable 0.1.4 AMS wheel')
            r['comparison']='Independent CPU AMS/IC at NGSolve1/16 versus stable GPU AMS and experimental GPU IC'
        r.update(wheel=wheel,installed=installed,cpu_ams_build_receipt_sha256=sha(job/'bridge-provenance/build.json'),cpu_ic_origin=owner['cpu_ic_origin'],versions=dict(python=sys.version,ngsolve=ng.__version__,numpy=np.__version__,cpu_native_ngslib_sha256=sha(ng.ngslib.__file__)),protocol_sha256={n:sha(job/n) for n in ('measure_scale_installed.py','measure_scale_installed.m')})
        r['cpu_hardware']=json.loads(subprocess.check_output(['powershell','-NoProfile','-Command','Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json'],text=True))
        if (os.cpu_count() or 0)<16:raise RuntimeError('CPU16 baseline unavailable')
        r['gpu_hardware']=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader,nounits'],text=True).strip()
        r['idle_before']=quiet();timer=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');r['matlab_start_seconds']=time.perf_counter()-timer
        eng.restoredefaultpath(nargout=0);eng.addpath(str(job),str(payload),nargout=0);pid=int(eng.feature('getpid'))
        eng.eval('timer=tic;g=gpuDevice;wait(g);gpu_init_seconds=toc(timer);',nargout=0);r['gpu_init_seconds']=float(eng.workspace['gpu_init_seconds'])
        scales=[1.,2.,.5,1.]
        def build(path,threads):
            ng.SetNumThreads(threads);timer=time.perf_counter();mesh=ng.Mesh(str(path));fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT()
            form=ng.BilinearForm(fes);form+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx
            rhs=ng.LinearForm(fes);rhs+=ng.CoefficientFunction((1,0,0))*v*ng.dx
            preparation=time.perf_counter()-timer;timer=time.perf_counter()
            with ng.TaskManager():form.Assemble();rhs.Assemble()
            assembly=time.perf_counter()-timer;free=np.asarray(list(fes.FreeDofs()),dtype=bool);rhs.vec.FV().NumPy()[~free]=0
            return mesh,fes,form,rhs,free,dict(mesh_space_seconds=preparation,matrix_rhs_assembly_seconds=assembly)
        for case in json.loads((job/'meshes.json').read_text()):
            if a.targets and case['target_dofs'] not in a.targets:continue
            path=job/case['mesh'];check=json.loads((job/case['check']).read_text())
            if not check['passed'] or check['input_sha256']!=sha(path) or check['contract_sha256']!=sha(job/'cube-labels.json'):raise RuntimeError('Mesh/label contract mismatch')
            record=dict(mesh=case,samples=[],first_initialized=[]);r['cases'].append(record)
            # Independent NGSolve sparsecholesky is accuracy truth, never GPU input.
            guard=quiet();timer=time.perf_counter();reference_data=build(path,1);mesh,fes,form,rhs,free,_=reference_data
            timer_factor=time.perf_counter();inverse=form.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky');factor_seconds=time.perf_counter()-timer_factor
            x=rhs.vec.CreateVector()
            with ng.TaskManager():x.data=inverse*rhs.vec
            reference=x.FV().NumPy().copy();del inverse
            record['independent_reference']=dict(backend='NGSolve sparsecholesky',threads=1,factor_seconds=factor_seconds,total_seconds=time.perf_counter()-timer,dofs=fes.ndof,free_dofs=int(free.sum()),idle_before=guard)
            if fes.ndof!=case['actual_dofs']:raise RuntimeError('FE dimension mismatch')
            def validate(values):
                values=np.asarray(values)
                if values.ndim==1:values=values[:,None]
                ng.SetNumThreads(1);results=[]
                for k,scale in enumerate(scales[:values.shape[1]]):
                    value=values[:,k];v=rhs.vec.CreateVector();v.FV().NumPy()[:]=value;applied=v.CreateVector()
                    with ng.TaskManager():applied.data=form.mat*v
                    error=float(np.linalg.norm(value-scale*reference)/np.linalg.norm(scale*reference));residual=float(np.linalg.norm((applied.FV().NumPy()-scale*rhs.vec.FV().NumPy())[free])/np.linalg.norm(scale*rhs.vec.FV().NumPy()[free]))
                    finite(value);finite([error,residual])
                    if error>1e-7 or residual>1e-6:raise RuntimeError(f'Independent parity/true residual failed {error} {residual}')
                    results.append(dict(current_scale=scale,solution_relative=error,true_relative_residual=residual))
                return results
            record['independent_reference']['check']=validate(reference)[0]
            def cpu(kind,threads):
                ng.SetNumThreads(threads);start=time.perf_counter();m,space,aform,b,mask,info=build(path,threads);timer=time.perf_counter()
                if kind=='ams':
                    gradient,_=space.CreateGradient();xyz=np.asarray([[p[0],p[1],p[2]] for p in m.ngmesh.Points()]);pre=bridge.preconditioner(aform.mat,gradient,space.FreeDofs(),*list(xyz.T),1)
                else:
                    with ng.TaskManager():pre=ic_native.ICPreconditioner(aform.mat,freedofs=space.FreeDofs(),shift=1.05)
                info['preconditioner_setup_seconds']=time.perf_counter()-timer;info['setup_seconds']=time.perf_counter()-start
                values=[];reports=[];seconds=[]
                for scale in scales:
                    timer=time.perf_counter();scaled=b.vec.CreateVector();scaled.data=scale*b.vec
                    solver=ng.CGSolver(aform.mat,pre=pre,tol=1e-12,maxiter=20000);solution=b.vec.CreateVector()
                    with ng.TaskManager():solution.data=solver*scaled
                    value=solution.FV().NumPy().copy();seconds.append(time.perf_counter()-timer);values.append(value)
                    reports.append(dict(iterations=int(solver.GetSteps()),ngsolve_threads=threads))
                info.update(fresh_total_seconds=info['setup_seconds']+seconds[0],solve_seconds=seconds,scales=scales,info=reports,preconditioner=kind,matrix_reused=True,preconditioner_reused=True)
                return np.column_stack(values),info
            def gpu(index,route):
                output=out/f"{case['target_dofs']}-{index}.mat";backend='ic' if route=='gpu-ic' else 'ams'
                eng.measure_scale_installed(str(path),str(payload),str(output),backend,str(a.ic_directory or ''),nargout=0)
                data=loadmat(output,simplify_cells=True);module=output.with_suffix('.module.json')
                subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(job/'ams_loaded_module.ps1'),'-ProcessId',str(pid),'-OutputPath',str(module)],check=True,capture_output=True)
                identity=json.loads(module.read_text(encoding='utf-8-sig'))
                if Path(identity['path']).resolve()!=(payload/'radia_cuda_ams.mexw64').resolve() or identity['sha256']!=mex_hash:raise RuntimeError('Actually loaded AMS identity changed')
                report=dict(data['report'],native_module=identity)
                if route=='gpu-ic':
                    ic_identity=output.with_suffix('.ic-module.json')
                    subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(job/'loaded_module.ps1'),'-ProcessId',str(pid),'-ModuleName','radia_cuda_ic.mexw64','-OutputPath',str(ic_identity)],check=True,capture_output=True)
                    actual=json.loads(ic_identity.read_text(encoding='utf-8-sig'))
                    if actual['sha256']!=ic_hash or Path(actual['path']).resolve()!=(a.ic_directory/'radia_cuda_ic.mexw64').resolve():raise RuntimeError('Actually loaded IC candidate identity changed')
                    report['ic_native_module']=actual
                return data['solution'],report
            routes=[('cpu-ams',1),('cpu-ic',1),('gpu-ams',1),('cpu-ams',16),('cpu-ic',16)]
            if a.ic_directory:routes.append(('gpu-ic',1))
            for index,(route,threads) in enumerate(routes):
                gc.collect();guard=quiet()
                with Memory(pid if route.startswith('gpu-') else os.getpid()) as memory:
                    value,info=gpu('first-'+str(index),route) if route.startswith('gpu-') else cpu(route[4:],threads)
                record['first_initialized'].append(dict(route=route,threads=threads,report=info,independent=validate(value),host_memory=memory.result(),idle_before=guard))
                del value;gc.collect()
            for cycle in range(a.cycles):
                for route,threads in routes+list(reversed(routes)):
                    gc.collect();guard=quiet();index=len(record['samples'])
                    with Memory(pid if route.startswith('gpu-') else os.getpid()) as memory:
                        value,info=gpu(index,route) if route.startswith('gpu-') else cpu(route[4:],threads)
                    row=dict(route=route,threads=threads,cycle=cycle,report=info,independent=validate(value),host_memory=memory.result(),idle_before=guard)
                    finite(row);record['samples'].append(row);del value;gc.collect();print(case['target_dofs'],route,threads,info['fresh_total_seconds'],flush=True)
                    (out/'progress.json').write_text(json.dumps(r,indent=2,default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v))+'\n')
            record['summary']={}
            for route,threads in routes:
                rows=[s for s in record['samples'] if s['route']==route and s['threads']==threads]
                record['summary'][route+'-'+str(threads)]=dict(fresh_median_seconds=statistics.median(float(s['report']['fresh_total_seconds']) for s in rows),setup_median_seconds=statistics.median(float(s['report']['setup_seconds']) for s in rows),reuse_median_seconds=statistics.median(float(v) for s in rows for v in np.asarray(s['report']['solve_seconds']).ravel()[1:]),iterations=[s['report']['info'] for s in rows])
            del reference_data,mesh,fes,form,rhs,x,reference;gc.collect()
        r['idle_after']=quiet();finite(r);r.update(accepted=True,timing_accepted=not a.diagnostic)
    except BaseException as ex:r['error']=repr(ex);raise
    finally:
        if eng is not None:
            try:eng.quit()
            except BaseException as ex:r['engine_release_error']=repr(ex);r.update(accepted=False,timing_accepted=False)
        (out/'result.json').write_text(json.dumps(r,indent=2,allow_nan=False,default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v))+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('job','wheel','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--cycles',type=int,default=2);p.add_argument('--targets',nargs='*',type=int);p.add_argument('--diagnostic',action='store_true');p.add_argument('--ic-directory',type=Path);main(p.parse_args())
