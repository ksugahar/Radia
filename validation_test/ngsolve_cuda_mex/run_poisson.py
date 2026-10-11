"""Compare NGSolve CPU and CUDA MEX on identical assembled Poisson systems."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata as md
import json
import math
from pathlib import Path
import platform
import socket
import statistics
import subprocess
import sys
import time
import matlab.engine
import ngsolve as ng
from netgen.csg import unit_cube
import numpy as np
import psutil
from scipy.io import savemat, loadmat
from scipy.sparse import coo_matrix

def run(args):
    root = Path(__file__).resolve().parents[2]
    out = args.out.resolve(); out.mkdir(parents=True,exist_ok=False)
    ng.SetNumThreads(args.threads)
    cases=[]
    start=time.perf_counter()
    eng=matlab.engine.start_matlab('-nodesktop -nosplash')
    startup=time.perf_counter()-start
    try:
        eng.addpath(str(Path(__file__).parent),nargout=0)
        gateway=eng.check_gateway(str(args.binary_dir.resolve()),nargout=1)
        for index,h in enumerate(args.maxh):
            case=out/f'case-{index}';case.mkdir()
            print(f'Assembling h={h}',flush=True)
            with ng.TaskManager():
                start=time.perf_counter()
                mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=h))
                fes=ng.H1(mesh,order=2,dirichlet='.*')
                u,v=fes.TnT()
                exact=ng.sin(math.pi*ng.x)*ng.sin(math.pi*ng.y)*ng.sin(math.pi*ng.z)
                a=ng.BilinearForm(ng.grad(u)*ng.grad(v)*ng.dx).Assemble()
                f=ng.LinearForm(3*math.pi**2*exact*v*ng.dx).Assemble()
                assembly=time.perf_counter()-start
                free=np.array(list(fes.FreeDofs()),dtype=bool)
                row,col,value=a.mat.COO()
                full=coo_matrix((np.array(value),(np.array(row),np.array(col))),shape=(fes.ndof,fes.ndof)).tocsr()
                matrix=full[free][:,free].tocsr();matrix.sort_indices()
                rhs=f.vec.FV().NumPy()[free].copy()
                start=time.perf_counter()
                inverse=a.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky')
                factor=time.perf_counter()-start
                cpu_direct=[]
                reference=f.vec.CreateVector()
                for _ in range(5):
                    start=time.perf_counter();reference.data=inverse*f.vec
                    cpu_direct.append(time.perf_counter()-start)
                reference_array=reference.FV().NumPy()[free].copy()
                dinv=np.zeros(fes.ndof);dinv[free]=1/matrix.diagonal()
                diagonal=f.vec.CreateVector();diagonal.FV().NumPy()[:]=dinv
                pre=ng.DiagonalMatrix(diagonal)
                cpu_cg=[];iterations=[]
                for _ in range(6):
                    solver=ng.CGSolver(a.mat,pre=pre,tol=1e-10,maxiter=10000)
                    solution=f.vec.CreateVector()
                    start=time.perf_counter();solution.data=solver*f.vec
                    cpu_cg.append(time.perf_counter()-start)
                    iterations.append(solver.iterations)
                cpu_solution=solution.FV().NumPy()[free].copy()
                gfu=ng.GridFunction(fes);gfu.vec.data=reference
                analytic_l2=float(ng.sqrt(ng.Integrate((gfu-exact)**2,mesh,order=8)))
            arrays=dict(rowptr=matrix.indptr.astype(np.int32),columns=matrix.indices.astype(np.int32),
                        values=matrix.data,rhs=rhs,diagonal_inverse=1/matrix.diagonal())
            input_file=case/'system.mat';savemat(input_file,arrays,oned_as='column')
            eng.benchmark_gpu(str(input_file),str(case/'gpu.mat'),str(args.binary_dir.resolve()),nargout=0)
            gpu=loadmat(case/'gpu.mat',simplify_cells=True)
            gpu_solution=np.asarray(gpu['solution']).reshape(-1)
            residual=lambda x:float(np.linalg.norm(matrix@x-rhs)/np.linalg.norm(rhs))
            relative=float(np.linalg.norm(gpu_solution-reference_array)/np.linalg.norm(reference_array))
            cpu_res=residual(cpu_solution);gpu_res=residual(gpu_solution)
            if max(cpu_res,gpu_res)>1e-6 or relative>1e-8:
                raise RuntimeError(f'Poisson parity failed: cpu={cpu_res}, gpu={gpu_res}, delta={relative}')
            row=dict(maxh=h,order=2,ndof=fes.ndof,free_dofs=int(free.sum()),nnz=int(matrix.nnz),
                     assembly_seconds=assembly,cpu_sparsecholesky_factor_seconds=factor,
                     cpu_sparsecholesky_warm_seconds=cpu_direct,
                     cpu_cg_cold_seconds=cpu_cg[0],cpu_cg_warm_seconds=cpu_cg[1:],
                     cpu_cg_warm_median_seconds=statistics.median(cpu_cg[1:]),cpu_cg_iterations=iterations,
                     cpu_true_relative_residual=cpu_res,gpu_true_relative_residual=gpu_res,
                     gpu_vs_cpu_direct_relative_l2=relative,analytic_solution_l2_error=analytic_l2,
                     input_sha256=hashlib.sha256(input_file.read_bytes()).hexdigest(),gpu=gpu['result'])
            row['gpu_resident_speedup_vs_cpu_cg']=row['cpu_cg_warm_median_seconds']/float(row['gpu']['warm_median_seconds'])
            row['gpu_transfer_inclusive_speedup_vs_cpu_cg']=row['cpu_cg_warm_median_seconds']/float(row['gpu']['warm_with_transfers_seconds'])
            cases.append(row)
            (case/'result.json').write_text(json.dumps(row,indent=2,default=convert)+'\n')
            print(json.dumps({k:row[k] for k in ('free_dofs','gpu_true_relative_residual','gpu_vs_cpu_direct_relative_l2','gpu_resident_speedup_vs_cpu_cg')}),flush=True)
        result=dict(schema='radia.ngsolve-cuda-poisson.v1',passed=True,hostname=socket.gethostname(),
                    platform=platform.platform(),python=sys.executable,versions={n:md.version(n) for n in ('ngsolve','netgen-mesher','numpy','scipy','matlabengine')},
                    ngsolve_upstream_cuda_enabled=bool(ng.config.USE_CUDA),cpu_threads=args.threads,
                    matlab_startup_seconds=startup,gateway_checks=dict(gateway),python_peak_rss_bytes=psutil.Process().memory_info().peak_wset,
                    source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
                    build=json.loads((args.binary_dir/'build.json').read_text()),
                    limitations=['Assembled CSR bridge, not an upstream NGSolve matrix-free CUDA backend',
                                 'Jacobi PCG only; IC/AMS/HDiv-MMM not implemented by this experiment',
                                 'CPU/GPU CG stopping norms differ; true residual and direct-reference parity are checked',
                                 'Small installation acceptance on 100, not an A100 scalability study'],cases=cases)
        (out/'result.json').write_text(json.dumps(result,indent=2,default=convert)+'\n')
    finally:eng.quit()

def convert(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--binary-dir',type=Path,required=True)
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--maxh',nargs='+',type=float,default=[0.24,0.16,0.11])
    args=p.parse_args()
    try:
        run(args)
    except Exception as error:
        args.out.mkdir(parents=True,exist_ok=True)
        (args.out/'failed-result.json').write_text(json.dumps(dict(
            schema='radia.ngsolve-cuda-poisson.v1',passed=False,
            hostname=socket.gethostname(),error=str(error),python=sys.executable),indent=2)+'\n')
        raise
