"""Independent FE oracle and same-host GPU/GPU AMS tuning acceptance."""
from pathlib import Path
import argparse,ctypes,hashlib,importlib.util,json,os,socket,statistics,subprocess,sys,time

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(a):
    job=a.job.resolve();out=job/'results';out.mkdir(exist_ok=False);eng=None
    record=dict(schema='radia.ams-fem-tuning.v1',accepted=False,timing_accepted=False,host=socket.gethostname(),samples=[],
        memory_definition='Exact added coarse-cache/workspace bytes, device-free endpoints and 50ms MATLAB-process samples; not a transient global-device peak',
        reuse_scope='Fixed affine linear FE matrix, labels, coefficients and hierarchy; explicit source-current scaling. Nonlinear solve uses existing adaptive updates.')
    try:
        manifest=json.loads((job/'manifest.json').read_text())
        for n,h in manifest.items():
            p=(job/n).resolve()
            if not p.is_relative_to(job) or sha(p)!=h:raise RuntimeError('Payload changed: '+n)
        record['payload_sha256']=manifest
        if (job/'owner.json').is_file() and (job/'vol_check.json').is_file():
            owner=json.loads((job/'owner.json').read_text(encoding='utf-8-sig'));volCheck=json.loads((job/'vol_check.json').read_text(encoding='utf-8-sig'))
            if not volCheck.get('passed') or ('vol_sha256' in owner and owner['vol_sha256']!=sha(job/'model.vol')):raise RuntimeError('Checked mesh identity mismatch')
            record['vol_check_passed']=True
        matlabPath=job/'matlab';tunedPath=job/'tuned'
        if bool(a.installed_python)!=bool(a.wheel):raise ValueError('Both --installed-python and --wheel are required')
        if a.installed_python:
            spec=importlib.util.spec_from_file_location('wheel_verifier',job/'verify_wheel.py');verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
            wheel=verifier.verify(a.wheel.resolve())
            installed=json.loads(subprocess.check_output([str(a.installed_python.resolve()),'-I','-c','import json,radia_gpu;print(json.dumps(radia_gpu.verify()))'],text=True))
            verifier.verify_installed_identity(installed,wheel)
            record['installed_wheel']=dict(wheel=wheel,installed=installed)
            matlabPath=Path(installed['matlab_path']).resolve();tunedPath=matlabPath
            wheelManifest=json.loads((matlabPath.parent/'manifest.json').read_text())
            tunedHash=wheelManifest['files']['matlab/radia_cuda_ams.mexw64']
        else:tunedHash=manifest['tuned/radia_cuda_ams.mexw64']
        os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL')
        import netgen
        net=Path(netgen.__file__).parent
        paths=[job/'runtime',net,net.parent/'ngsolve_openblas',Path(sys.prefix)/'bin',Path(sys.base_prefix)]
        os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
        handles=[os.add_dll_directory(str(p)) for p in paths if p.is_dir()]
        import numpy as np,ngsolve as ng,matlab.engine
        from scipy.io import loadmat
        from measure_native_cpu_gpu import idle,Memory
        from nonlinear_resident import cpu_run
        sys.path.insert(0,str(job/'cpu'));import _radia_ams_cpu_reference as bridge
        ng.SetNumThreads(1);record['idle_before']=idle()
        mesh=ng.Mesh(str(job/'model.vol'));fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT()
        free=np.asarray(list(fes.FreeDofs()),dtype=bool);G,_=fes.CreateGradient()
        xyz=np.array([[p[0],p[1],p[2]] for p in mesh.ngmesh.Points()])
        form=ng.BilinearForm(fes);b2=ng.curl(u)*ng.curl(u);form+=ng.Variation((.5*b2+b2*b2+.5*u*u)*ng.dx)
        f=ng.LinearForm(fes);f+=ng.CoefficientFunction((20,0,0))*v*ng.dx
        with ng.TaskManager():f.Assemble()
        f.vec.FV().NumPy()[~free]=0
        nonlinear,cpuInfo=cpu_run(3,fes,form,f.vec,G,[xyz[:,c] for c in range(3)],bridge,1)
        linearForm=ng.BilinearForm(fes);linearForm+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx
        with ng.TaskManager():linearForm.Assemble()
        linearRHS=f.vec.CreateVector();linearRHS.data=.05*f.vec
        inverse=linearForm.mat.Inverse(freedofs=fes.FreeDofs(),inverse='sparsecholesky')
        lin=linearRHS.CreateVector()
        with ng.TaskManager():lin.data=inverse*linearRHS
        linear=lin.FV().NumPy().copy()
        record['cpu_reference']=dict(nonlinear=cpuInfo,linear_backend='NGSolve sparsecholesky',bridge_sha256=sha(bridge.__file__),mesh_sha256=sha(job/'model.vol'),dofs=fes.ndof,elements=mesh.ne)
        record['versions']=dict(python=sys.version,ngsolve=ng.__version__,numpy=np.__version__)
        def validate(value,nonlinearMode,scale=1):
            probe=ng.GridFunction(fes);probe.vec.FV().NumPy()[:]=np.asarray(value).ravel();r=probe.vec.CreateVector()
            with ng.TaskManager():
                if nonlinearMode:form.Apply(probe.vec,r)
                else:r.data=linearForm.mat*probe.vec
            rhs=f.vec.FV().NumPy() if nonlinearMode else linearRHS.FV().NumPy()*scale
            denominator=np.linalg.norm(rhs[free]);true=float(np.linalg.norm((r.FV().NumPy()-rhs)[free])/max(denominator,1e-300))
            expected=nonlinear if nonlinearMode else linear*scale
            error=float(np.linalg.norm(np.asarray(value).ravel()-expected)/max(np.linalg.norm(expected),1e-300))
            if error>1e-7 or true>(1e-9 if nonlinearMode else 1e-10):raise RuntimeError(f'Independent CPU parity/true residual failed: {error}, {true}')
            return dict(solution_relative=error,true_relative_residual=true)
        record['cpu_self_checks']=[validate(nonlinear,True),validate(linear,False)]
        t=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');record['matlab_start_seconds']=time.perf_counter()-t
        eng.restoredefaultpath(nargout=0)
        eng.addpath(str(job),nargout=0);eng.eval('timer=tic;g=gpuDevice;wait(g);gpu_init_seconds=toc(timer);',nargout=0)
        record['gpu_init_seconds']=float(eng.workspace['gpu_init_seconds']);pid=int(eng.feature('getpid'))
        record['cold_boundary']='First API call per variant in one owned initialized MATLAB process is saved separately; not an operating-system/device cold reset. Startup and GPU initialization are separate.'
        for cycle in range(a.cycles):
            for variant in ('baseline','tuned','tuned','baseline'):
                guard=idle();path=out/f'{cycle}-{len(record["samples"])}-{variant}.mat'
                binaryPath=tunedPath if variant=='tuned' else job/'baseline'
                with Memory(pid) as memory:eng.ams_tuning_sample(str(job),str(binaryPath),str(path),str(matlabPath),True,nargout=0)
                modulePath=path.with_suffix('.module.json')
                subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(job/'ams_loaded_module.ps1'),'-ProcessId',str(pid),'-OutputPath',str(modulePath)],check=True,capture_output=True)
                raw=loadmat(path,simplify_cells=True);report=raw['report']
                loaded=Path(report['loaded_binary']).resolve()
                expectedHash=tunedHash if variant=='tuned' else manifest['baseline/radia_cuda_ams.mexw64']
                if loaded!=binaryPath/'radia_cuda_ams.mexw64' or sha(loaded)!=expectedHash:raise RuntimeError('Loaded MEX identity mismatch')
                module=json.loads(modulePath.read_text(encoding='utf-8-sig'));report['native_module']=module
                if Path(module['path']).resolve()!=loaded or module['sha256']!=expectedHash:raise RuntimeError('Actual loaded Windows MEX module mismatch')
                eng.eval('clear radia_cuda_ams;',nargout=0);eng.rmpath(str(binaryPath),nargout=0)
                errors=dict(nonlinear=validate(raw['nonlinear'],True),linear=validate(raw['linear'],False),
                    reused=[validate(raw['reuse'][:,k],False,float(s)) for k,s in enumerate(np.atleast_1d(report['scales']))])
                sample=dict(variant=variant,cycle=cycle,report=report,errors=errors,idle_before=guard,host_memory=memory.result());record['samples'].append(sample)
                print(variant,'nonlinear',report['nonlinear_seconds'],'linear',report['linear_fresh_seconds'],'reuse',np.median(report['reused_seconds']),flush=True)
                (out/'result.json').write_text(json.dumps(record,indent=2,default=encode)+'\n')
        record['summary']={}
        record['first_api_seconds']={v:next({k:float(s['report'][k]) for k in ('nonlinear_seconds','linear_fresh_seconds','prepared_setup_seconds')} for s in record['samples'] if s['variant']==v) for v in ('baseline','tuned')}
        for variant in ('baseline','tuned'):
            rows=[s['report'] for s in record['samples'] if s['variant']==variant]
            record['summary'][variant]={k:statistics.median(float(r[k]) for r in rows) for k in ('nonlinear_seconds','linear_fresh_seconds','prepared_setup_seconds')}
            # Exclude zero RHS and first graph-build solve; retain varying nonzero scales.
            record['summary'][variant]['reused_warm_seconds']=statistics.median(float(np.atleast_1d(r['reused_seconds'])[k]) for r in rows for k in (1,2,3,5))
        record['idle_after']=idle();record.update(accepted=True,timing_accepted=socket.gethostname().lower()=='ngsolve')
    except BaseException as e:record['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2,default=encode)+'\n')
def encode(v):
    if hasattr(v,'tolist'):return v.tolist()
    return str(v)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);p.add_argument('--cycles',type=int,default=2)
    p.add_argument('--installed-python',type=Path);p.add_argument('--wheel',type=Path);main(p.parse_args())
