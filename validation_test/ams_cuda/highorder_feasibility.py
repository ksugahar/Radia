"""Same-host high-order BDDC+AMS versus generic GPU IC validation fixtures.

This does not implement or certify high-order MATLAB/mexcuda FE assembly.
"""
from pathlib import Path
import argparse,hashlib,json,os,sys,time,statistics,subprocess,socket

def main(args):
    job=args.job.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL')
    import netgen
    net=Path(netgen.__file__).parent
    paths=[job/'runtime',net,net.parent/'ngsolve_openblas',Path(sys.prefix)/'bin',Path(sys.base_prefix)]
    os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
    handles=[os.add_dll_directory(str(p)) for p in paths if p.is_dir()]
    import numpy as np,ngsolve as ng,matlab.engine,psutil
    from scipy.io import savemat,loadmat
    from scipy.sparse import coo_matrix,diags
    sys.path.insert(0,str(job/'cpu'));import sparsesolv_ngsolve as ssn
    sys.path.insert(0,str(job));from measure_native_cpu_gpu import idle,Memory
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    record=dict(schema='radia.ams-highorder-feasibility.v1',host=socket.gethostname(),cases=[],
        scope='CPU BDDC+AMS1/16 vs generic GPU IC on CPU-assembled validation CSR; not GPU FE assembly',
        GPU_FE_assembly_supported=False,unsupported_reason='assembleHCurl supports only affine lowest-order tetrahedra',
        fixed_threads=[1,16],blas_threads=1,ngsolve=ng.__version__,numpy=np.__version__,
        native_cpu_sha256=sha(ssn.__file__),ic_sha256=sha(args.ic_directory/'radia_cuda_ic.mexw64'),
        formal_wheel_accepted=False,all_routes_passed=False)
    engine=None;scales=[1.,2.,.5,1.]
    try:
        if socket.gethostname().lower()!='ngsolve' and not args.cpu_only:raise RuntimeError('Same-host GPU comparison requires selected A100 worker')
        meshfile=job/'10000.vol';contract=json.loads((job/'10000.check.json').read_text())
        if not contract['passed'] or contract['input_sha256']!=sha(meshfile):raise RuntimeError('Checked mesh identity failed')
        for order in (2,3):
            case=dict(order=order,mesh_sha256=sha(meshfile),cpu=[]);record['cases'].append(case)
            ng.SetNumThreads(1);mesh=ng.Mesh(str(meshfile));space=ng.HCurl(mesh,order=order,nograds=True,dirichlet='.*');u,v=space.TnT()
            full=ng.BilinearForm(space,symmetric=True);full+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx
            rhs=ng.LinearForm(space);rhs+=ng.CoefficientFunction((1,0,0))*v*ng.dx
            with ng.TaskManager():full.Assemble();rhs.Assemble()
            free=np.asarray(list(space.FreeDofs()),dtype=bool);rhs.vec.FV().NumPy()[~free]=0
            with ng.TaskManager():inverse=full.mat.Inverse(space.FreeDofs(),inverse='sparsecholesky');ref=rhs.vec.CreateVector();ref.data=inverse*rhs.vec
            reference=ref.FV().NumPy().copy();del inverse
            case.update(total_dofs=space.ndof,free_dofs=int(free.sum()))
            def validate(values):
                values=np.asarray(values);values=values[:,None] if values.ndim==1 else values;checks=[]
                for k in range(values.shape[1]):
                    vec=rhs.vec.CreateVector();vec.FV().NumPy()[:]=values[:,k]
                    applied=vec.CreateVector()
                    with ng.TaskManager():applied.data=full.mat*vec
                    error=float(np.linalg.norm(values[:,k]-scales[k]*reference)/np.linalg.norm(scales[k]*reference))
                    residual=float(np.linalg.norm((applied.FV().NumPy()-scales[k]*rhs.vec.FV().NumPy())[free])/np.linalg.norm(scales[k]*rhs.vec.FV().NumPy()[free]))
                    if not np.isfinite(error+residual) or error>1e-7 or residual>1e-6:raise RuntimeError(f'Independent true residual/parity failed {error} {residual}')
                    checks.append(dict(solution_relative=error,true_relative_residual=residual))
                return checks
            case['reference_check']=validate(reference)
            for threads in (1,16):
                item=dict(threads=threads,preconditioner='NGSolve BDDC with coarsetype=sparsesolv_ams',status='pending',process_rss_before=psutil.Process().memory_info().rss);case['cpu'].append(item)
                try:
                    guard=idle() if not args.cpu_only else {'scope':'Build-host CPU-only feasibility, not GPU performance evidence'}
                    ng.SetNumThreads(threads);timer=time.perf_counter();m=ng.Mesh(str(meshfile));s=ng.HCurl(m,order=order,nograds=True,dirichlet='.*');uu,vv=s.TnT()
                    item['mesh_space_seconds']=time.perf_counter()-timer
                    a=ng.BilinearForm(s,symmetric=True,condense=True);a+=(ng.curl(uu)*ng.curl(vv)+uu*vv)*ng.dx
                    b=ng.LinearForm(s);b+=ng.CoefficientFunction((1,0,0))*vv*ng.dx
                    p=ng.Preconditioner(a,'bddc',coarsetype='sparsesolv_ams',inverse='sparsecholesky',coarseflags={'cycles':1,'lean_coarse':0,'print_level':0})
                    timer=time.perf_counter()
                    with ng.TaskManager():a.Assemble();b.Assemble()
                    item['assembly_and_preconditioner_seconds']=time.perf_counter()-timer
                    item['preconditioner_setup_seconds']=None
                    item['phase_limitation']='BDDC local contributions/coarse setup occur during NGSolve assembly; no fabricated disjoint timer'
                    item['coarse_statistics']=ssn.AMSCoarseStats()
                    values=[];seconds=[];iterations=[]
                    for scale in scales:
                        timer=time.perf_counter();scaled=b.vec.CreateVector();scaled.data=scale*b.vec
                        with ng.TaskManager():
                            condensed=scaled.CreateVector();condensed.data=scaled+a.harmonic_extension_trans*scaled
                            solver=ng.CGSolver(a.mat,pre=p,tol=1e-12,maxiter=20000);x=scaled.CreateVector();x.data=solver*condensed
                            x.data+=a.harmonic_extension*x;x.data+=a.inner_solve*scaled
                        values.append(x.FV().NumPy().copy());seconds.append(time.perf_counter()-timer);iterations.append(int(solver.GetSteps()))
                    item.update(status='passed',solve_seconds=seconds,iterations=iterations,warm_reuse_median_seconds=statistics.median(seconds[1:]),independent=validate(np.column_stack(values)),idle_before=guard)
                except Exception as ex:item.update(status='failed',error=repr(ex));print('CPU feasibility FAILED',order,threads,repr(ex),flush=True)
                finally:item.update(process_rss_after=psutil.Process().memory_info().rss,memory_definition='Process RSS endpoints including retained independent reference; not transient peak')
            if args.cpu_only:
                case['gpu']=dict(status='not_run',reason='CPU-only build-host feasibility');continue
            # Explicit solver-validation fixture. NGSolve owns assembly; this
            # exchange is not a required GPU-user assembly adapter.
            rows,columns,values=full.mat.COO();matrix=coo_matrix((np.asarray(values),(np.asarray(rows),np.asarray(columns))),shape=(space.ndof,space.ndof)).tocsr()
            from scipy.sparse import triu,tril
            if triu(matrix,1).nnz==0 or tril(matrix,-1).nnz==0:
                matrix=matrix+matrix.T-diags(matrix.diagonal());case['coo_storage']='Expanded symmetric triangular storage without numerical approximation'
            projection=diags(free.astype(float));matrix=(projection@matrix@projection+diags((~free).astype(float))).tocsr();matrix.sum_duplicates();matrix.sort_indices()
            probe=np.linspace(.2,1.7,space.ndof);probe[~free]=0;probe_vec=rhs.vec.CreateVector();probe_vec.FV().NumPy()[:]=probe
            applied=probe_vec.CreateVector()
            with ng.TaskManager():applied.data=full.mat*probe_vec
            csr_error=float(np.linalg.norm((matrix@probe-applied.FV().NumPy())[free])/np.linalg.norm(applied.FV().NumPy()[free]))
            if not np.isfinite(csr_error) or csr_error>1e-12:raise RuntimeError('Projected fixture CSR does not reproduce NGSolve operator action')
            case['fixture_csr_action_relative']=csr_error
            if matrix.nnz>=2**31 or matrix.shape[0]>=2**31:raise RuntimeError('GPU CSR int32 capacity')
            arrays=dict(rows=float(space.ndof),rowptr=matrix.indptr.astype(np.int32)[:,None],columns=matrix.indices.astype(np.int32)[:,None],values=matrix.data[:,None],free=free.astype(np.int32)[:,None],rhs=rhs.vec.FV().NumPy()[:,None])
            if sum(getattr(v,'nbytes',8) for v in arrays.values())>=2**31:raise RuntimeError('MAT v7 capacity: use HDF5; never split')
            fixture=out/f'p{order}-fixture.mat';savemat(fixture,arrays);case['fixture_sha256']=sha(fixture)
            if engine is None:
                timer=time.perf_counter();engine=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');record['matlab_start_seconds']=time.perf_counter()-timer
                engine.restoredefaultpath(nargout=0);engine.addpath(str(Path(__file__).parent),nargout=0)
            gpu=dict(status='pending',scope='GPU IC solver-only on independent CPU-assembled validation CSR');case['gpu']=gpu
            try:
                guard=idle();output=out/f'p{order}-gpu.mat'
                with Memory(int(engine.feature('getpid'))) as memory:engine.highorder_ic_fixture(str(fixture),str(args.ic_directory.resolve()),str(output),nargout=0)
                data=loadmat(output,simplify_cells=True);report=dict(data['report']);checks=validate(data['solution'])
                module=out/f'p{order}-module.json'
                subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(job/'loaded_module.ps1'),'-ProcessId',str(int(engine.feature('getpid'))),'-ModuleName','radia_cuda_ic.mexw64','-OutputPath',str(module)],check=True,capture_output=True)
                actual=json.loads(module.read_text(encoding='utf-8-sig'))
                if actual['sha256']!=record['ic_sha256'] or Path(actual['path']).resolve()!=(args.ic_directory/'radia_cuda_ic.mexw64').resolve():raise RuntimeError('Loaded fixed IC identity mismatch')
                gpu.update(status='passed',report=report,independent=checks,idle_before=guard,host_memory=memory.result(),actual_loaded=actual)
            except Exception as ex:gpu.update(status='failed',error=repr(ex));print('GPU feasibility FAILED',order,repr(ex),flush=True)
            (out/'progress.json').write_text(json.dumps(record,indent=2,default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v))+'\n')
        record['evaluated']=True
        record['all_routes_passed']=all(all(x['status']=='passed' for x in c['cpu']) and (args.cpu_only or c['gpu']['status']=='passed') for c in record['cases'])
    except Exception as ex:record['error']=repr(ex);raise
    finally:
        if engine is not None:engine.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2,allow_nan=False,default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v))+'\n')
    if not record['all_routes_passed']:raise RuntimeError('Some high-order routes failed; see explicit per-route evidence, no fallback')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job',required=True,type=Path);p.add_argument('--ic-directory',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--cpu-only',action='store_true');main(p.parse_args())
