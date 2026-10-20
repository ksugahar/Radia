"""Nontrivial 3D HCurl AMS parity, original residual and reuse-memory validation.

Uses a public in-memory Netgen unit cube. NGSolve owns assembly and the independent
CPU sparsecholesky solve; mexcuda owns the GPU solve. No physical scaling claim.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def main(args):
    import numpy as np
    import ngsolve as ng
    from netgen.csg import unit_cube
    import sys
    source_root=Path(__file__).resolve().parents[2]
    sys.path[:0]=[str(source_root),str(source_root/'packages/radia-simulink/src')]
    from radia_simulink.ams_exchange import write_exchange
    from validation_test.hmatrix_cuda.hdiv_hdf5 import read_hdiv_h5
    import matlab.engine
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    build=args.build.resolve();matlab_dir=args.matlab.resolve()
    receipt=json.loads((build/'build.json').read_text())
    binary=build/'radia_cuda_ams.mexw64'
    digest=hashlib.sha256(binary.read_bytes()).hexdigest()
    if digest!=receipt['sha256']['build-ams-cuda/radia_cuda_ams.mexw64']:
        raise RuntimeError('AMS binary receipt mismatch')
    report=dict(accepted=False,scope='dimensionless linear 3D HCurl unit cube; not nonlinear or scaling acceptance',
                native_source_commit=receipt['source_commit'],binary_sha256=digest,
                helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                matlab_helper_sha256=hashlib.sha256(Path(__file__).with_name('checkSharedAMS.m').read_bytes()).hexdigest(),
                ngsolve=ng.__version__,cpu_threads=args.threads)
    eng=None
    try:
        ng.SetNumThreads(args.threads);start=time.perf_counter()
        with ng.TaskManager():
            mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=args.maxh))
            if set(mesh.GetBoundaries())!={'back','front','left','right','top','bottom'}:
                raise RuntimeError('Unit cube boundary contract changed')
            fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT()
            free=np.asarray(list(fes.FreeDofs()),dtype=bool)
            if not (64<=free.sum() and fes.ndof<20000):
                raise RuntimeError('Case must be nontrivial and below frozen DoF cap')
            form=ng.BilinearForm(fes);form+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx;form.Assemble()
            gradient,_=fes.CreateGradient()
            exact=form.mat.CreateColVector();exact.FV().NumPy()[:]=np.linspace(1,2,fes.ndof)*free
            rhs=exact.CreateVector();rhs.data=form.mat*exact;rhs.FV().NumPy()[~free]=0
            inv=form.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky')
            expected=exact.CreateVector();expected.data=inv*rhs
            def csr(matrix):
                values,columns,rowptr=matrix.CSR()
                return dict(rows=float(matrix.height),cols=float(matrix.width),rowptr=np.array(rowptr,dtype=np.int32,copy=True),columns=np.array(columns,dtype=np.int32,copy=True),values=np.array(values,copy=True))
            coords=np.array([[mesh.ngmesh.Points()[i+1][j] for j in range(3)] for i in range(mesh.nv)])
            model=dict(schema='radia.ams-native-input.v1',A=csr(form.mat),G=csr(gradient),coordinates=coords,free=free.astype(np.int32),cycle=1.)
            write_exchange(out/'input.h5',dict(model=model,rhs=np.array(rhs.FV().NumPy()),expected=np.array(expected.FV().NumPy())))
        report.update(ndof=fes.ndof,free_dofs=int(free.sum()),vertices=mesh.nv,elements=mesh.ne,fe_and_cpu_seconds=time.perf_counter()-start)
        start=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash')
        report['matlab_startup_seconds']=time.perf_counter()-start
        eng.addpath(str(Path(__file__).parent),nargout=0)
        eng.checkSharedAMS(str(out),str(build),str(matlab_dir),float(args.repeats),nargout=0)
        report['gpu']=json.loads((out/'gpu-result.json').read_text())
        solutions=read_hdiv_h5(out/'solutions.h5')[0]['solutions']
        comparisons=[]
        for values in solutions.T:
            probe=exact.CreateVector();probe.FV().NumPy()[:]=values
            residual=exact.CreateVector()
            with ng.TaskManager():residual.data=form.mat*probe
            residual.data-=rhs
            error=float(np.linalg.norm(values-expected.FV().NumPy())/np.linalg.norm(expected.FV().NumPy()))
            rr=float(np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(rhs.FV().NumPy()[free]))
            comparisons.append(dict(cpu_relative_error=error,independent_original_true_residual=rr))
            if error>1e-8 or rr>1e-9:raise RuntimeError('Independent CPU/original equation gate failed')
        report['comparisons']=comparisons;report['accepted']=True
    except BaseException as e:
        report['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True);p.add_argument('--matlab',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--maxh',type=float,default=.12);p.add_argument('--threads',type=int,default=2);p.add_argument('--repeats',type=int,default=100)
    main(p.parse_args())
