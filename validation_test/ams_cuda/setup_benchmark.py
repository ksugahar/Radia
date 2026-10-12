"""Compare identical AMS hierarchies with serial/OpenMP setup; always save JSON.

100 timings are diagnostics. Formal performance acceptance requires idle mdx.
The explicit CPU-only option records transfer-inclusive timing as pending.
"""
from __future__ import annotations
import argparse
import importlib.metadata as md
import json
import os
from pathlib import Path
import socket
import statistics
import sys
import time

import numpy as np
import ngsolve as ng
import psutil
from netgen.csg import unit_cube
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'packages/radia-gpu/src'))
from radia_gpu.ams_exchange import write_exchange, file_sha256
from run import source_identity, jsonable

def same_snapshot(reference, actual):
    """Require identical hierarchy, ordering and FP64 coefficients."""
    if isinstance(reference,np.ndarray):
        if not np.array_equal(reference,actual):
            raise RuntimeError('OpenMP changed AMS hierarchy array')
    elif isinstance(reference,dict):
        if reference.keys()!=actual.keys():raise RuntimeError('Snapshot keys differ')
        for key in reference:same_snapshot(reference[key],actual[key])
    elif isinstance(reference,(list,tuple)):
        if len(reference)!=len(actual):raise RuntimeError('Snapshot lengths differ')
        for left,right in zip(reference,actual):same_snapshot(left,right)
    elif reference!=actual:raise RuntimeError('Snapshot scalar differs')

def main(args):
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(ROOT/'build-ams-cuda'))
    import _radia_ams_snapshot as bridge
    ng.SetNumThreads(1)
    build=json.loads((ROOT/'build-ams-cuda/build.json').read_text(encoding='utf-8'))
    record=dict(schema='radia.ams-setup-comparison.v1',accepted=False,
        host=socket.gethostname(),timing_classification='diagnostic; not formal performance acceptance',
        formal_performance_acceptance='pending: requires idle mdx with a GPU',
        taskmanager_used=False,ngsolve_threads=1,**source_identity(ROOT),
        versions={p:md.version(p) for p in ('ngsolve','numpy','scipy','h5py')},
        build=build,
        memory_measurement=dict(gpu='Live allocation delta after repeated solves; not transient allocator peak',
            host='Process lifetime peak working set; includes previous cases',
            admission='Conservative GPU estimate with reserved free memory'),
        cases=[])
    eng=None
    try:
        for relative,digest in {**build['source_sha256'],**build['sha256']}.items():
            if file_sha256(ROOT/relative)!=digest:raise RuntimeError('Build/source hash mismatch')
        record['parallel_contract_passed']=bool(bridge._check_openmp_contract(args.setup_threads))
        if not args.cpu_only:
            import matlab.engine
            import netgen
            net=Path(netgen.__file__).parent
            os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),
                str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),
                str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
            start=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash')
            record['matlab_start_seconds']=time.perf_counter()-start
            eng.addpath(str(Path(__file__).parent),nargout=0)
        for h in args.maxh:
            mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=h))
            fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*')
            u,v=fes.TnT();a=ng.BilinearForm(fes,symmetric=False)
            a+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx;a.Assemble()
            gradient,_=fes.CreateGradient()
            coords=[[mesh.ngmesh.Points()[i+1][j] for i in range(mesh.nv)] for j in range(3)]
            probes=np.random.default_rng(31415).normal(size=(3,fes.ndof))
            free=np.array(list(fes.FreeDofs()),dtype=bool)
            exact=np.random.default_rng(782).normal(size=fes.ndof)*free
            x=a.mat.CreateRowVector();x.FV().NumPy()[:]=exact
            b=a.mat.CreateColVector();b.data=a.mat*x
            rhs=b.FV().NumPy().copy();rhs[~free]=0;b.FV().NumPy()[:]=rhs
            for cycle in (1,7):
                row=dict(maxh=h,cycle=cycle,dofs=fes.ndof,samples=[])
                record['cases'].append(row)
                reference=None;oracle=None
                # ABBA alternates order to expose warm-up and drift.
                order=[1,args.setup_threads,args.setup_threads,1]*args.rounds
                for index,threads in enumerate(order):
                    start=time.perf_counter()
                    snapshot=bridge.snapshot(a.mat,gradient,fes.FreeDofs(),*coords,probes,cycle,threads)
                    snapshot_wall=time.perf_counter()-start
                    pre=snapshot.pop('cpu_preconditioner')
                    expected=snapshot.pop('cpu_probe_results')
                    diagnostics=snapshot.pop('setup_diagnostics')
                    if reference is None:reference=snapshot;oracle=expected.copy()
                    same_snapshot(reference,snapshot)
                    probe_error=float(np.linalg.norm(expected-oracle)/max(np.linalg.norm(oracle),np.finfo(float).tiny))
                    start=time.perf_counter();solver=ng.CGSolver(a.mat,pre=pre,tol=1e-9,maxiter=2000)
                    solution=b.CreateVector();solution.data=solver*b
                    cpu_reference_solve_seconds=time.perf_counter()-start
                    residual_vector=b.CreateVector();residual_vector.data=a.mat*solution-b
                    residual=residual_vector.FV().NumPy().copy()
                    true_residual=float(np.linalg.norm(residual[free])/np.linalg.norm(rhs[free]))
                    sample=dict(setup_threads=threads,snapshot_wall_seconds=snapshot_wall,
                        setup_seconds=diagnostics['hierarchy_seconds']+diagnostics['export_seconds'],
                        setup=diagnostics,hierarchy_identical=True,probe_relative_difference=probe_error,
                        cpu_true_residual=true_residual,cpu_iterations=solver.iterations,
                        cpu_reference_solve_seconds=cpu_reference_solve_seconds)
                    row['samples'].append(sample)
                    if probe_error>1e-12 or true_residual>1e-6:raise RuntimeError('CPU setup parity failed')
                    directory=out/f'h{h}-c{cycle}-sample{index}';directory.mkdir()
                    path=directory/'input.hdf5'
                    start=time.perf_counter()
                    sample['exchange']=write_exchange(path,dict(snapshot=snapshot,rhs=rhs,probes=probes,expected=expected))
                    sample['hdf5_write_seconds']=time.perf_counter()-start
                    sample['input_sha256']=file_sha256(path)
                    sample['hdf5_bytes']=path.stat().st_size
                    if eng is not None:
                        target=directory/'gpu.mat'
                        eng.benchmark(str(path),str(target),str(ROOT/'build-ams-cuda'),str(ROOT/'packages/radia-gpu/matlab'),nargout=0)
                        gpu=loadmat(target,simplify_cells=True)
                        actual=np.asarray(gpu['inclusiveSolution']).reshape(-1)
                        cpu=solution.FV().NumPy().copy()
                        difference=float(np.linalg.norm((actual-cpu)[free])/np.linalg.norm(cpu[free]))
                        x.FV().NumPy()[:]=actual
                        residual_vector.data=a.mat*x-b
                        r=residual_vector.FV().NumPy().copy()
                        gpu_residual=float(np.linalg.norm(r[free])/np.linalg.norm(rhs[free]))
                        error=float(np.max(gpu['probeErrors']))
                        sample.update(gpu_relative_difference_from_cpu=difference,gpu_true_residual=gpu_residual,
                            gpu_probe_relative_error=error,upload_seconds=float(gpu['uploadSeconds']),
                            gpu_cold_solve_seconds=float(gpu['coldSeconds']),
                            gpu_resident_solve_median_seconds=float(np.median(gpu['times'])),
                            gpu_memory_bytes=float(gpu['memoryBytes']),gpu_memory_plan=jsonable(gpu['memoryPlan']),
                            gpu_device=jsonable(gpu['device']),
                            hdf5_read_seconds=float(gpu['readSeconds']),gpu_inclusive_seconds=float(gpu['inclusiveSeconds']))
                        sample['setup_exchange_transfer_solve_seconds']=(sample['setup_seconds']+sample['hdf5_write_seconds']+
                            sample['hdf5_read_seconds']+sample['gpu_inclusive_seconds'])
                        if difference>1e-6 or gpu_residual>1e-6 or error>1e-8:raise RuntimeError('GPU parity failed')
                    else:sample['setup_exchange_transfer_solve_seconds']=None
                    memory=psutil.Process().memory_info()
                    sample['python_rss_bytes']=memory.rss
                    sample['python_lifetime_peak_rss_bytes']=getattr(memory,'peak_wset',None)
                    (out/'result.json').write_text(json.dumps(jsonable(record),indent=2,allow_nan=False)+'\n',encoding='utf-8')
                row['medians']={str(t):{key:statistics.median(s[key] for s in row['samples'] if s['setup_threads']==t)
                    for key in ('snapshot_wall_seconds','setup_seconds','setup_exchange_transfer_solve_seconds',
                        'cpu_reference_solve_seconds','gpu_resident_solve_median_seconds','gpu_inclusive_seconds',
                        'gpu_memory_bytes')
                    if row['samples'][0].get(key) is not None} for t in (1,args.setup_threads)}
                timing_keys=('snapshot_wall_seconds','setup_seconds','setup_exchange_transfer_solve_seconds')
                before={key:value for key,value in row['medians']['1'].items() if key in timing_keys}
                after={key:value for key,value in row['medians'][str(args.setup_threads)].items() if key in timing_keys}
                row['before_after']=dict(before_setup_threads=1,after_setup_threads=args.setup_threads,
                    before=before,after=after,
                    time_reduction_seconds={key:before[key]-after[key] for key in before},
                    time_reduction_fraction={key:1-after[key]/before[key] for key in before if before[key]>0})
                print(json.dumps(dict(dofs=fes.ndof,cycle=cycle,medians=row['medians'])),flush=True)
        record['accepted']=True
    except Exception as error:
        record['error']=repr(error);raise
    finally:
        if eng is not None:eng.quit()
        (out/'result.json').write_text(json.dumps(jsonable(record),indent=2,allow_nan=False)+'\n',encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--maxh',type=float,nargs='+',default=[0.18])
    p.add_argument('--setup-threads',type=int,choices=range(2,65),default=4)
    p.add_argument('--rounds',type=int,default=2)
    p.add_argument('--cpu-only',action='store_true')
    args=p.parse_args()
    if args.rounds<1:p.error('--rounds must be positive')
    main(args)
