"""Formal idle A100 same-host CPU/GPU ABBA, original numerical gates unchanged."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,statistics,subprocess,sys,threading,time
START=time.perf_counter()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def idle():
    import psutil
    cpu=psutil.cpu_percent(interval=1)
    raw=subprocess.check_output(['nvidia-smi','--query-gpu=name,utilization.gpu,memory.free,memory.total','--format=csv,noheader,nounits'],text=True).strip()
    name,util,free,total=[x.strip() for x in raw.split(',')]
    available=psutil.virtual_memory().available
    if cpu>10 or float(util)>5 or float(free)<max(2048,.2*float(total)) or available<8*2**30:raise RuntimeError('Idle measurement guard failed: '+raw+f', CPU={cpu}, host_free={available}')
    return dict(cpu_percent=cpu,gpu=name,gpu_utilization_percent=float(util),device_free_MiB=float(free),host_available_bytes=available)

class Memory:
    def __init__(self,pid):
        import psutil
        self.process=psutil.Process(pid);self.samples=[];self.stop=threading.Event()
    def sample(self):
        m=self.process.memory_info();self.samples.append(dict(rss=m.rss,private=getattr(m,'private',m.vms),lifetime_peak_rss=getattr(m,'peak_wset',m.rss)))
    def __enter__(self):
        self.sample()
        def observe():
            while not self.stop.wait(.05):self.sample()
        self.thread=threading.Thread(target=observe,daemon=True);self.thread.start();return self
    def __exit__(self,*args):self.stop.set();self.thread.join();self.sample()
    def result(self):return dict(before=self.samples[0],after=self.samples[-1],sampled_peak_rss=max(s['rss'] for s in self.samples),sampled_peak_private=max(s['private'] for s in self.samples),sample_count=len(self.samples),period_seconds=.05,definition='50ms process samples include baseline/imported runtimes; lifetime peak is not a phase peak; short allocation spikes may be missed')

def main(a):
    job=a.job.resolve();eng=None;r=dict(schema='radia.native-ams-fmm-formal-performance.v1',accepted=False,timing_accepted=False,release_accepted=False,host=socket.gethostname(),date='2026-10-13',samples=[],scope='Bounded independently validated HCurl Newton/AMS, scalar/dipole Laplace FMM and flat P1 boundary operators',classification='Formal same-host NGSolve(A100) CPU/GPU evaluation, explicitly authorized by user',cpu_threads=1,gpu_host_setup_threads=1,cold_definition='First numerical API call per case in one fresh owned process; process imports, MATLAB startup and GPU initialization recorded separately',memory_limits='Exact native owned bounds and device-free endpoint samples plus 50ms host-process samples; no global-device/CUDA runtime transient peak certification')
    try:
        if socket.gethostname().lower()!='ngsolve':raise RuntimeError('Formal GPU host must be NGSolve')
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        manifest=json.loads((job/'manifest.json').read_text())
        for n,h in manifest.items():
            p=(job/n).resolve()
            if not p.is_relative_to(job) or sha(p)!=h:raise RuntimeError('Payload hash mismatch: '+n)
        r['payload_sha256']=manifest
        os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL',OPENBLAS_NUM_THREADS='1')
        import netgen
        net=Path(netgen.__file__).parent
        paths=[job/'runtime',net,net.parent/'ngsolve_openblas',Path(sys.prefix)/'bin',Path(sys.base_prefix)]
        os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
        handles=[os.add_dll_directory(str(p)) for p in paths if p.is_dir()]
        import numpy as np,ngsolve as ng,matlab.engine,psutil
        from scipy.io import loadmat
        from scipy.sparse import coo_matrix
        from ngsolve.bem import LaplaceSL,LaplaceDL
        sys.path.insert(0,str(job/'cpu'));import _radia_ams_cpu_reference as bridge
        from nonlinear_resident import cpu_run
        from fmm_numeric import direct
        ng.SetNumThreads(1)
        r['cpu_import_seconds']=time.perf_counter()-START
        r['versions']=dict(python=sys.version,ngsolve=ng.__version__,numpy=np.__version__)
        r['cpu_reference']=dict(bridge_sha256=sha(bridge.__file__),ngsolve_extension_sha256=sha(ng.ngslib.__file__),nonlinear_cpu_sha256=sha(job/'nonlinear_resident.py'),direct_sum_sha256=sha(job/'fmm_numeric.py'))
        r['idle_before']=idle()
        t=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');r['matlab_start_seconds']=time.perf_counter()-t
        pid=int(eng.feature('getpid'));r['matlab']=eng.version()
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.addpath(str(job),str(job/'matlab'),nargout=0)
        eng.eval('timer=tic;device=gpuDevice;wait(device);gpu_init_seconds=toc(timer);',nargout=0);r['gpu_init_seconds']=float(eng.workspace['gpu_init_seconds'])
        phys=loadmat(job/'physics-input.mat',simplify_cells=True)
        numeric={name:loadmat(job/file,simplify_cells=True) for name,file in [('fmm-small','numeric-input.mat'),('fmm-large','numeric-large.mat')]}
        boundary=loadmat(job/'boundary-reference.mat',simplify_cells=True);ids=np.asarray(boundary['global_dofs'],dtype=int).ravel();density=np.asarray(boundary['x']).ravel()
        def cpu(name):
            if name=='ams':
                mesh=ng.Mesh(str(job/'model.vol'));fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT();b2=ng.curl(u)*ng.curl(u)
                form=ng.BilinearForm(fes);form+=ng.Variation((.5*b2+b2*b2+.5*u*u)*ng.dx)
                rhs=ng.LinearForm(fes);rhs+=ng.CoefficientFunction((20,0,0))*v*ng.dx
                with ng.TaskManager():rhs.Assemble()
                free=np.asarray(list(fes.FreeDofs()),dtype=bool);rhs.vec.FV().NumPy()[~free]=0
                G,_=fes.CreateGradient();xyz=np.array([[p[0],p[1],p[2]] for p in mesh.ngmesh.Points()])
                value,info=cpu_run(3,fes,form,rhs.vec,G,[xyz[:,c] for c in range(3)],bridge,1)
                return [value],info
            if name in numeric:
                s=numeric[name];p,e=direct(s['sources'],s['targets'],s['charge'],s['dipole']);return [p,e],{}
            if name=='physical':
                s=phys;ns=len(s['sphere']);zero=np.zeros((ns,3));weights=s['weights'].reshape(-1,1)
                p,e=direct(s['sphere'],s['targets'],weights,zero)
                dl,_=direct(s['sphere'],s['targets'],np.zeros((ns,1)),s['sphere']*weights)
                ps,es=direct(s['coil'],s['coilTargets'],s['current'],np.zeros((len(s['coil']),9)))
                scale=1.25663706212e-6/(4*np.pi);B=scale*np.column_stack([es[:,5]-es[:,7],es[:,6]-es[:,2],es[:,1]-es[:,3]])
                return [p[:,0]/(4*np.pi),dl[:,0]/(4*np.pi),B,scale*ps,e],{}
            mesh=ng.Mesh(str(job/'boundary.vol'));fes=ng.H1(mesh,order=1);u,v=fes.TnT()
            with ng.TaskManager():
                sl=(LaplaceSL(u.Trace()*ng.ds(bonus_intorder=16),use_fmm=False)*v.Trace()*ng.ds(bonus_intorder=16)).mat
                dl=(LaplaceDL(u.Trace()*ng.ds(bonus_intorder=16),use_fmm=False)*v.Trace()*ng.ds(bonus_intorder=16)).mat
            values=[]
            for m in [sl,dl]:
                i,j,x=m.COO();matrix=coo_matrix((x,(i,j)),shape=(m.height,m.width)).tocsr();values.append(matrix[ids][:,ids]@density)
            return values,{}
        references={};results={}
        def sample(name,route,phase):
            guard=idle()
            if route=='cpu':
                with Memory(os.getpid()) as memory:
                    start=time.perf_counter();values,info=cpu(name);seconds=time.perf_counter()-start
                if name not in references:references[name]=values
            else:
                out=job/f'{name}-{len(r["samples"])}.mat'
                with Memory(pid) as memory:eng.measure_native_sample(str(job),name,str(out),nargout=0)
                raw=loadmat(out,simplify_cells=True);values=raw['values'];values=[values] if name=='ams' else list(values);info=raw['report'];seconds=float(info['whole_seconds'])
                expected=references[name];errors=[float(np.linalg.norm(np.asarray(x).ravel()-np.asarray(y).ravel())/max(np.linalg.norm(y),np.finfo(float).tiny)) for x,y in zip(values,expected)]
                if len(values)!=len(expected):raise RuntimeError('Output count changed')
                # A-axis has analytically zero potential: use absolute error.
                if name=='physical':errors[3]=float(np.linalg.norm(np.asarray(values[3])-np.asarray(expected[3])))
                limit=1e-7 if name=='ams' else 1e-6 if name=='boundary' else 1e-5
                if max(errors)>limit:raise RuntimeError(f'Independent {name} parity failed: {errors}')
                info['independent_relative_errors']=errors
                if name=='ams':
                    mesh=ng.Mesh(str(job/'model.vol'));fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT();b2=ng.curl(u)*ng.curl(u)
                    form=ng.BilinearForm(fes);form+=ng.Variation((.5*b2+b2*b2+.5*u*u)*ng.dx)
                    rhs=ng.LinearForm(fes);rhs+=ng.CoefficientFunction((20,0,0))*v*ng.dx
                    x=ng.GridFunction(fes);x.vec.FV().NumPy()[:]=np.asarray(values[0]).ravel();res=x.vec.CreateVector()
                    with ng.TaskManager():rhs.Assemble();form.Apply(x.vec,res)
                    free=np.asarray(list(fes.FreeDofs()),dtype=bool);true=float(np.linalg.norm((res.FV().NumPy()-rhs.vec.FV().NumPy())[free])/np.linalg.norm(rhs.vec.FV().NumPy()[free]))
                    if true>1e-9:raise RuntimeError('Independent GPU nonlinear true residual failed')
                    info['independent_cpu_true_residual']=true
            row=dict(case=name,backend=route,phase=phase,whole_seconds=seconds,idle_before=guard,host_process_memory=memory.result(),info=info)
            r['samples'].append(row);(job/'result.json').write_text(json.dumps(r,indent=2,default=encode)+'\n');print(name,route,phase,seconds,flush=True)
        for name in ['ams','fmm-small','fmm-large','physical','boundary']:
            sample(name,'cpu','first');sample(name,'gpu','first')
            for cycle in range(a.cycles):
                for route in ['cpu','gpu','gpu','cpu']:sample(name,route,f'ABBA-{cycle}')
            rows=[x for x in r['samples'] if x['case']==name and x['phase']!='first']
            c=statistics.median(x['whole_seconds'] for x in rows if x['backend']=='cpu');g=statistics.median(x['whole_seconds'] for x in rows if x['backend']=='gpu')
            results[name]=dict(cpu_warm_median_seconds=c,gpu_transfer_inclusive_warm_median_seconds=g,cpu_over_gpu_speedup=c/g)
        r['summary']=results;r['idle_after']=idle()
        subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(job/'modules.ps1'),'-ProcessId',str(pid),'-OutputPath',str(job/'modules.json')],check=True,capture_output=True)
        loaded={x['name'].lower():x for x in json.loads((job/'modules.json').read_text(encoding='utf-8-sig'))}
        for name in ['radia_cuda_ams.mexw64','radia_cuda_fmm.mexw64']:
            if loaded[name]['sha256']!=manifest[name]:raise RuntimeError('Actual MATLAB MEX hash mismatch')
        for p in (job/'runtime').glob('*.dll'):
            m=loaded[p.name.lower()]
            if m['sha256']!=sha(p) or Path(m['path']).parent!=job/'runtime':raise RuntimeError('Actual MATLAB native dependency mismatch')
        r['actual_matlab_native_dependencies_verified']=True
        for n,h in manifest.items():
            if sha(job/n)!=h:raise RuntimeError('Payload changed during benchmark')
        r.update(accepted=True,timing_accepted=True)
    except BaseException as e:r['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        (job/'result.json').write_text(json.dumps(r,indent=2,default=encode)+'\n')
def encode(v):
    import numpy as np
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);p.add_argument('--cycles',type=int,default=2);a=p.parse_args()
    if a.cycles<2:p.error('At least two ABBA cycles required')
    main(a)
