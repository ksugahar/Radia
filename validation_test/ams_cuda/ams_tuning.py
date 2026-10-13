"""Independent FE oracle and same-host GPU/GPU AMS tuning acceptance."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,statistics,subprocess,sys,time

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
        eng.addpath(str(job),nargout=0);eng.eval('g=gpuDevice;wait(g);',nargout=0);pid=int(eng.feature('getpid'))
        for cycle in range(a.cycles):
            for variant in ('baseline','tuned','tuned','baseline'):
                guard=idle();path=out/f'{cycle}-{len(record["samples"])}-{variant}.mat'
                with Memory(pid) as memory:eng.ams_tuning_sample(str(job),str(job/variant),str(path),nargout=0)
                raw=loadmat(path,simplify_cells=True);report=raw['report']
                loaded=Path(report['loaded_binary']).resolve()
                if loaded!=job/variant/'radia_cuda_ams.mexw64' or sha(loaded)!=manifest[variant+'/radia_cuda_ams.mexw64']:raise RuntimeError('Loaded MEX identity mismatch')
                errors=dict(nonlinear=validate(raw['nonlinear'],True),linear=validate(raw['linear'],False),
                    reused=[validate(raw['reuse'][:,k],False,float(s)) for k,s in enumerate(np.atleast_1d(report['scales']))])
                sample=dict(variant=variant,cycle=cycle,report=report,errors=errors,idle_before=guard,host_memory=memory.result());record['samples'].append(sample)
                print(variant,'nonlinear',report['nonlinear_seconds'],'linear',report['linear_fresh_seconds'],'reuse',np.median(report['reused_seconds']),flush=True)
                (out/'result.json').write_text(json.dumps(record,indent=2,default=encode)+'\n')
        record['summary']={}
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
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);p.add_argument('--cycles',type=int,default=2);main(p.parse_args())
