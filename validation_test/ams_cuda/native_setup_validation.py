"""Validate native MATLAB AMS setup against independent Python CPU execution.

Python generates validation inputs/oracles only; GPU users pass numeric FEM
operators directly to prepareAMS. All outcomes, including failures, save JSON.
"""
from pathlib import Path
import argparse, json, os, socket, statistics, subprocess, sys, time
import numpy as np
import ngsolve as ng
from netgen.csg import unit_cube
from scipy.sparse import coo_matrix
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'packages/radia-gpu/src'))
from radia_gpu.ams_exchange import write_exchange, file_sha256

def csr(matrix):
    i,j,v=matrix.COO()
    a=coo_matrix((v,(i,j)),shape=(matrix.height,matrix.width)).tocsr()
    a.sort_indices()
    return dict(rows=a.shape[0],cols=a.shape[1],rowptr=a.indptr.astype(np.int32),
                columns=a.indices.astype(np.int32),values=a.data.astype(np.float64))

def main(args):
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.ams-matlab-native-validation.v1',accepted=False,
        host=socket.gethostname(),timing_classification='diagnostic; not formal performance acceptance',
        a100_started=False,python_role='independent CPU oracle and validation orchestration only',
        gpu_runtime_python_adapter_required=False,taskmanager_regions_enabled=False,
        execution_mode='CPU fixture generation only' if args.generate_only else 'native MATLAB validation',
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),cases=[])
    eng=None
    try:
        import netgen
        net=Path(netgen.__file__).parent
        os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
        sys.path.insert(0,str(ROOT/'build-ams-cuda'))
        import _radia_ams_snapshot as bridge
        ng.SetNumThreads(1)
        if not args.generate_only:
            import matlab.engine
            eng=matlab.engine.start_matlab('-nodesktop -nosplash')
            eng.addpath(str(Path(__file__).parent),nargout=0)
        record['build']=json.loads((ROOT/'build-ams-cuda/build.json').read_text())
        for h in args.maxh:
            mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=h))
            print(f'mesh ready: maxh={h}, vertices={mesh.nv}',flush=True)
            fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*')
            u,v=fes.TnT();a=ng.BilinearForm(fes,symmetric=False)
            a+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx;a.Assemble()
            gradient,_=fes.CreateGradient()
            coords=[[mesh.ngmesh.Points()[i+1][j] for i in range(mesh.nv)] for j in range(3)]
            free=np.array(list(fes.FreeDofs()),dtype=bool)
            rng=np.random.default_rng(782);exact=rng.normal(size=fes.ndof)*free
            x=a.mat.CreateRowVector();x.FV().NumPy()[:]=exact
            b=a.mat.CreateColVector();b.data=a.mat*x
            rhs=b.FV().NumPy().copy();rhs[~free]=0;b.FV().NumPy()[:]=rhs
            probes=np.random.default_rng(31415).normal(size=(3,fes.ndof))
            for cycle in args.cycles:
                reference=bridge.snapshot(a.mat,gradient,fes.FreeDofs(),*coords,probes,cycle,1)
                pre=reference.pop('cpu_preconditioner');expected=reference.pop('cpu_probe_results');cpu_setup=reference.pop('setup_diagnostics')
                timer=time.perf_counter();solver=ng.CGSolver(a.mat,pre=pre,tol=1e-9,maxiter=2000)
                solution=b.CreateVector();solution.data=solver*b;cpu_seconds=time.perf_counter()-timer
                cpu=solution.FV().NumPy().copy()
                residual=b.CreateVector();residual.data=a.mat*solution-b
                cpu_res=float(np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(rhs[free]))
                if cpu_res>1e-6:raise RuntimeError('CPU true residual failed')
                model=dict(schema='radia.ams-native-input.v1',A=csr(a.mat),G=csr(gradient),
                    coordinates=np.asarray(coords).T.copy(),free=free.astype(np.int32),cycle=cycle)
                case=dict(maxh=h,dofs=fes.ndof,cycle=cycle,cpu_true_residual=cpu_res,
                          cpu_reference_solve_seconds=cpu_seconds,
                          cpu_reference_setup_seconds=cpu_setup['hierarchy_seconds'],
                          cpu_reference_setup_and_solve_seconds=cpu_setup['hierarchy_seconds']+cpu_seconds,samples=[])
                record['cases'].append(case)
                path=out/f'h{h}-c{cycle}.hdf5'
                write_exchange(path,dict(model=model,snapshot=reference,rhs=rhs,probes=probes,expected=expected))
                case['input_sha256']=file_sha256(path)
                case.update(input_bytes=path.stat().st_size,vertices=mesh.nv,
                            matrix_nonzeros=len(model['A']['values']))
                print(f'CPU/HDF5 ready: dofs={fes.ndof}, cycle={cycle}, bytes={path.stat().st_size}',flush=True)
                if args.generate_only:
                    continue
                warmup=out/f'h{h}-c{cycle}-cold.mat'
                eng.benchmark(str(path),str(warmup),str(ROOT/'build-ams-cuda'),
                              str(ROOT/'packages/radia-gpu/matlab'),1.0,nargout=0)
                cold=loadmat(warmup,simplify_cells=True)
                case['first_call']=dict(native_setup_seconds=float(cold['nativeSetupSeconds']),
                    gpu_solve_seconds=float(cold['coldSeconds']),
                    setup_transfer_solve_download_seconds=float(cold['nativeSetupSeconds']+cold['inclusiveSeconds']))
                for index,threads in enumerate([1,args.threads,args.threads,1]*args.rounds):
                    target=out/f'h{h}-c{cycle}-{index}.mat'
                    eng.benchmark(str(path),str(target),str(ROOT/'build-ams-cuda'),
                                  str(ROOT/'packages/radia-gpu/matlab'),float(threads),nargout=0)
                    result=loadmat(target,simplify_cells=True)
                    actual=np.asarray(result['inclusiveSolution']).reshape(-1)
                    error=float(np.linalg.norm((actual-cpu)[free])/np.linalg.norm(cpu[free]))
                    x.FV().NumPy()[:]=actual;residual.data=a.mat*x-b
                    true_res=float(np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(rhs[free]))
                    probe=float(np.max(result['probeErrors']))
                    sample=dict(setup_threads=threads,native_setup_seconds=float(result['nativeSetupSeconds']),
                        transfer_solve_download_seconds=float(result['inclusiveSeconds']),
                        setup_transfer_solve_download_seconds=float(result['nativeSetupSeconds']+result['inclusiveSeconds']),
                        warmed_solve_median_seconds=float(np.median(result['times'])),
                        cold_solve_seconds=float(result['coldSeconds']),upload_seconds=float(result['uploadSeconds']),
                        hierarchy_relative_difference=float(result['hierarchyMaxRelativeDifference']),
                        hierarchy_bitwise_identical=True,cpu_relative_difference=error,
                        gpu_true_residual=true_res,probe_relative_difference=probe,
                        gpu_live_allocation_delta_bytes=float(result['memoryBytes']),
                        memory_plan=result['memoryPlan'],device=result['device'],matlab=result['matlabVersion'],
                        rejection_checks=result['checks'],native_setup=result['nativeSetupInfo'])
                    case['samples'].append(sample)
                    if error>1e-6 or true_res>1e-6 or probe>1e-8:raise RuntimeError('Native GPU parity failed')
                before=[s for s in case['samples'] if s['setup_threads']==1]
                after=[s for s in case['samples'] if s['setup_threads']==args.threads]
                key='setup_transfer_solve_download_seconds'
                bt=statistics.median(s[key] for s in before);at=statistics.median(s[key] for s in after)
                case['before_after']=dict(before_seconds=bt,after_seconds=at,speedup=bt/at,
                    comparison='same shared native AMS: serial host setup vs OpenMP host setup; same CUDA solve')
        record['accepted']=True
    except BaseException as e:
        record['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        def convert(x):
            if isinstance(x,np.ndarray):return x.tolist()
            if isinstance(x,np.generic):return x.item()
            raise TypeError(type(x).__name__)
        (out/'result.json').write_text(json.dumps(record,indent=2,default=convert,allow_nan=False)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--maxh',type=float,nargs='+',default=[0.12,0.08])
    p.add_argument('--threads',type=int,default=4);p.add_argument('--rounds',type=int,default=1)
    p.add_argument('--generate-only',action='store_true',help='Write independent CPU/HDF5 fixtures without starting MATLAB')
    p.add_argument('--cycles',type=int,nargs='+',choices=[1,7],default=[1,7])
    main(p.parse_args())
