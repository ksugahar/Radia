"""Small AMS CUDA parity acceptance; larger studies belong on idle mdx/A100."""
from __future__ import annotations
import argparse
import importlib.metadata as md
import json
from pathlib import Path
import platform
import socket
import statistics
import subprocess
import sys
import time
import numpy as np
import ngsolve as ng
from netgen.csg import unit_cube
from scipy.io import loadmat
from scipy.sparse import coo_matrix
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'packages/radia-gpu/src'))
from radia_gpu.ams_exchange import write_exchange, read_exchange, file_sha256

def source_identity(root, manifest=None):
    if manifest:
        identity=json.loads(manifest.read_text(encoding='utf-8'))
        commit=identity['source_commit']
        if len(commit)!=40 or any(c not in '0123456789abcdef' for c in commit):
            raise ValueError('Invalid exported Git source identity')
        for relative,digest in identity['source_file_sha256'].items():
            path=(root/relative).resolve()
            if not path.is_relative_to(root.resolve()) or file_sha256(path)!=digest:
                raise RuntimeError('Exported Git source hash mismatch')
        return {'source_commit':commit,'source_dirty':False,
                'provenance_route':'hash-verified exported Git tree'}
    return {'source_commit':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
            'source_dirty':bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip()),
            'provenance_route':'Git checkout'}

def jsonable(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:jsonable(v) for k,v in x.items() if not k.startswith('__')}
    if isinstance(x,(list,tuple)):return [jsonable(v) for v in x]
    if isinstance(x,bytes):return x.decode('utf-8',errors='replace')
    return x

def main(args):
    root=Path(__file__).resolve().parents[2]
    binary=root/'build-ams-cuda'
    gpu_dir=args.gpu_dir.resolve() if args.gpu_dir else binary
    matlab_dir=args.matlab_dir.resolve() if args.matlab_dir else root/'packages/radia-gpu/matlab'
    sys.path.insert(0,str(binary))
    import _radia_ams_snapshot as bridge
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    if args.threads != 1:
        raise ValueError("The mexcuda validation route requires serial NGSolve preprocessing (--threads 1)")
    ng.SetNumThreads(1)
    record={'schema':'radia.ams-cuda-validation.v1','host':socket.gethostname(),'python':sys.executable,
        'platform':platform.platform(),'versions':{p:md.version(p) for p in ('ngsolve','numpy','scipy','h5py')},
        'threads':args.threads,'taskmanager_used':False,'precision':'FP64','cases':[],
        'scope':'experimental assembled-operator AMS; CPU-owned setup, GPU-owned application and coarse solves',
        **source_identity(root,args.source_manifest),
        'build':json.loads((binary/'build.json').read_text()),
        'execution_mode':'cpu-only' if args.cpu_only else 'cpu-gpu-parity','accepted':False}
    if not args.cpu_only:
        record['gpu_binary_sha256']=file_sha256(gpu_dir/'radia_cuda_ams.mexw64')
    result_path=out/'result.json'
    eng=None
    try:
        if not args.cpu_only:
            import matlab.engine
            record['versions']['matlabengine']=md.version('matlabengine')
            timer=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash');record['matlab_start_seconds']=time.perf_counter()-timer
            eng.addpath(str(Path(__file__).parent),nargout=0)
        for h in args.maxh:
            for cycle in (1,7):
                print(f'AMS maxh={h} cycle={cycle}',flush=True)
                timer=time.perf_counter()
                mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=h))
                fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*')
                u,v=fes.TnT();a=ng.BilinearForm(fes,symmetric=False)
                a+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx;a.Assemble()
                gradient,_=fes.CreateGradient()
                assembly=time.perf_counter()-timer
                free=np.array(list(fes.FreeDofs()),dtype=bool)
                coords=[[mesh.ngmesh.Points()[i+1][j] for i in range(mesh.nv)] for j in range(3)]
                probes=np.random.default_rng(31415).normal(size=(3,fes.ndof))
                timer=time.perf_counter()
                snapshot=bridge.snapshot(a.mat,gradient,fes.FreeDofs(),*coords,probes,cycle)
                setup=time.perf_counter()-timer
                expected=snapshot.pop('cpu_probe_results');pre=snapshot.pop('cpu_preconditioner')
                setup_diagnostics=snapshot.pop('setup_diagnostics',{})
                x_exact=np.random.default_rng(782).normal(size=fes.ndof)*free
                native=a.mat.CreateRowVector();native.FV().NumPy()[:]=x_exact
                b=a.mat.CreateColVector();b.data=a.mat*native
                rhs=b.FV().NumPy().copy();rhs[~free]=0;b.FV().NumPy()[:]=rhs
                timer=time.perf_counter();inverse=a.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky');factor=time.perf_counter()-timer
                reference=b.CreateVector();reference.data=inverse*b
                cpu=[];cpu_iterations=[]
                for k in range(6):
                    solver=ng.CGSolver(a.mat,pre=pre,tol=1e-9,maxiter=2000)
                    x=b.CreateVector();timer=time.perf_counter()
                    x.data=solver*b
                    cpu.append(time.perf_counter()-timer);cpu_iterations.append(solver.iterations)
                case=out/f'h{h}-cycle{cycle}';case.mkdir()
                payload=dict(snapshot=snapshot,rhs=rhs,probes=probes,expected=expected)
                src=case/'input.hdf5';exchange=write_exchange(src,payload)
                if args.cpu_only:
                    row,col,val=a.mat.COO();A=coo_matrix((np.asarray(val),(np.asarray(row),np.asarray(col))),shape=(fes.ndof,fes.ndof)).tocsr()
                    actual=x.FV().NumPy().copy();ref=reference.FV().NumPy()
                    residual=float(np.linalg.norm((A@actual-rhs)[free])/np.linalg.norm(rhs[free]))
                    difference=float(np.linalg.norm((actual-ref)[free])/np.linalg.norm(ref[free]))
                    entry={'maxh':h,'cycle':cycle,'dofs':fes.ndof,'free_dofs':int(free.sum()),
                        'assembly_seconds':assembly,'snapshot_setup_and_probes_seconds':setup,'setup_diagnostics':setup_diagnostics,
                        'cpu_direct_factor_seconds':factor,'cpu_ams_seconds':cpu,'cpu_iterations':cpu_iterations,
                        'cpu_warm_median_seconds':statistics.median(cpu[1:]),'cpu_true_residual':residual,
                        'relative_difference_from_direct':difference,'exchange':exchange,'input_sha256':file_sha256(src),
                        'accepted':residual<=1e-6 and difference<=1e-6}
                    record['cases'].append(entry);result_path.write_text(json.dumps(jsonable(record),indent=2)+'\n')
                    if not entry['accepted']:raise RuntimeError('CPU AMS reference check failed')
                    continue
                eng.benchmark(str(src),str(case/'gpu.mat'),str(gpu_dir),str(matlab_dir),nargout=0)
                gpu=loadmat(case/'gpu.mat',simplify_cells=True)
                actual=np.asarray(gpu['solution']).reshape(-1)
                row,col,val=a.mat.COO();A=coo_matrix((np.asarray(val),(np.asarray(row),np.asarray(col))),shape=(fes.ndof,fes.ndof)).tocsr()
                residual=lambda sol:float(np.linalg.norm((A@sol-rhs)[free])/np.linalg.norm(rhs[free]))
                ref=reference.FV().NumPy().copy()
                difference=float(np.linalg.norm((actual-ref)[free])/np.linalg.norm(ref[free]))
                cpu_ams_solution=x.FV().NumPy().copy()
                cpu_ams_difference=float(np.linalg.norm((actual-cpu_ams_solution)[free])/np.linalg.norm(cpu_ams_solution[free]))
                probe=float(np.max(gpu['probeErrors']))
                entry={'maxh':h,'cycle':cycle,'dofs':fes.ndof,'free_dofs':int(free.sum()),'elements':mesh.ne,
                    'assembly_seconds':assembly,'snapshot_setup_and_probes_seconds':setup,'setup_diagnostics':setup_diagnostics,'cpu_direct_factor_seconds':factor,
                    'cpu_ams_seconds':cpu,'cpu_iterations':cpu_iterations,'cpu_warm_median_seconds':statistics.median(cpu[1:]),
                    'gpu':jsonable(gpu),'cpu_true_residual':residual(x.FV().NumPy()),'gpu_true_residual':residual(actual),
                    'relative_difference_from_direct':difference,'max_preconditioner_probe_relative_error':probe,
                    'relative_difference_from_cpu_ams':cpu_ams_difference,
                    'input_sha256':file_sha256(src),'exchange':exchange,
                    'auxiliary_level_rows':[[l['A']['rows'] for l in s['levels']] for s in snapshot['spaces']],
                    'coarse_kinds':[s['levels'][-1]['coarse']['kind'] for s in snapshot['spaces']]}
                entry['gpu'].pop('solution');entry['gpu'].pop('inclusiveSolution')
                entry['gpu_warm_median_seconds']=statistics.median(np.asarray(gpu['times']).reshape(-1))
                entry['speedup_resident_inputs']=entry['cpu_warm_median_seconds']/entry['gpu_warm_median_seconds']
                entry['speedup_including_transfers']=entry['cpu_warm_median_seconds']/float(gpu['inclusiveSeconds'])
                entry['accepted']=bool(entry['gpu_true_residual']<=1e-6 and entry['cpu_true_residual']<=1e-6 and difference<=1e-6 and cpu_ams_difference<=1e-6 and probe<=1e-8)
                record['cases'].append(entry);result_path.write_text(json.dumps(jsonable(record),indent=2)+'\n')
                print(json.dumps({k:entry[k] for k in ('dofs','cycle','gpu_true_residual','max_preconditioner_probe_relative_error','accepted')}),flush=True)
                if not entry['accepted']:raise RuntimeError('AMS CUDA numerical parity failed')
        record['accepted']=True
    except Exception as e:
        record['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        result_path.write_text(json.dumps(jsonable(record),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--maxh',type=float,nargs='+',default=[0.5,0.35]);p.add_argument('--threads',type=int,default=1)
    p.add_argument('--gpu-dir',type=Path);p.add_argument('--matlab-dir',type=Path)
    p.add_argument('--cpu-only',action='store_true',help='Explicit CPU timing route for a host without a GPU')
    p.add_argument('--source-manifest',type=Path,help='Checked exported-Git identity for compute hosts without Git')
    main(p.parse_args())
