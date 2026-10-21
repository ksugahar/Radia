"""CUDA resident AMG lifecycle and repeated-RHS stream-order regression."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--native-directory',type=Path,required=True)
    p.add_argument('--gpu-dll-directory',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();record=dict(passed=False,cases=[])
    dll=os.add_dll_directory(str(args.gpu_dll_directory.resolve()))
    native=None
    try:
        binary=next(args.native_directory.glob('_resident_amg*.pyd'))
        spec=importlib.util.spec_from_file_location('_resident_amg',binary)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        record['binary_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest()
        n=100000
        def csr(n,m,row,col,val):return dict(n=n,m=m,row=np.asarray(row,dtype=np.int32),col=np.asarray(col,dtype=np.int32),val=np.asarray(val,dtype=np.float64))
        a=csr(n,n,np.arange(n+1),np.arange(n),np.full(n,2.))
        prolongation=csr(n,1,np.arange(n+1),np.zeros(n,dtype=int),np.ones(n))
        restriction=csr(1,n,[0,n],np.arange(n),np.ones(n))
        coarse=csr(1,1,[0,1],[0],[2.*n])
        native=module.Resident([dict(a=a,p=prolongation,r=restriction,smoother=np.full(n,.3)),dict(a=coarse,smoother=np.array([.3]))],np.array([1/(2*n)]))
        for i in range(30):
            rhs=np.linspace(1+i,2+i,n)*(1 if i%2 else -1)
            solution,info=native.solve(rhs,1e-12,100)
            residual=float(np.linalg.norm(2*solution-rhs)/np.linalg.norm(rhs))
            record['cases'].append(dict(original_relative_residual=residual,native=info))
            if residual>1e-10:raise RuntimeError('Repeated-RHS original residual failed')
        zero,info=native.solve(np.zeros(n),1e-12,100)
        if np.any(zero) or info['iterations']!=0:raise RuntimeError('Zero RHS failed')
        try:native.solve(np.full(n,np.nan),1e-12,100)
        except ValueError:pass
        else:raise RuntimeError('Nonfinite RHS accepted')
        native.close();native.close()
        try:native.solve(np.ones(n),1e-12,100)
        except RuntimeError:pass
        else:raise AssertionError('Released context accepted')
        # Long rows select the warp kernel; 257 rows exercise its final partial
        # row group. Restriction is also a long row, prolongation a short one.
        from scipy.sparse import diags
        n=257
        operator=diags([np.full(n-abs(k),32. if k==0 else -.1)
                        for k in range(-20,21)],range(-20,21),shape=(n,n)).tocsr()
        def export(a):return csr(*a.shape,a.indptr,a.indices,a.data)
        coarse_value=float(operator.sum())
        native=module.Resident([
            dict(a=export(operator),p=csr(n,1,np.arange(n+1),np.zeros(n),np.ones(n)),
                 r=csr(1,n,[0,n],np.arange(n),np.ones(n)),smoother=np.full(n,1/32.)),
            dict(a=csr(1,1,[0,1],[0],[coarse_value]),smoother=np.array([1/coarse_value]))],
            np.array([1/coarse_value]))
        record['long_row_cases']=[]
        for i in range(3):
            rhs=np.cos(np.arange(n)*.13+i)
            solution,info=native.solve(rhs,1e-12,100)
            reference=np.linalg.solve(operator.toarray(),rhs)
            residual=float(np.linalg.norm(operator@solution-rhs)/np.linalg.norm(rhs))
            error=float(np.linalg.norm(solution-reference)/np.linalg.norm(reference))
            record['long_row_cases'].append(dict(original_relative_residual=residual,relative_error=error,native=info))
            if residual>1e-10 or error>1e-10:raise RuntimeError('Long-row reference mismatch')
        native.close()
        record['ill_conditioned_cases']=[]
        rng=np.random.default_rng(1729)
        from scipy.linalg import solveh_banded
        for n in (1025,2049):
            operator=diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1]).tocsr()
            bands=np.vstack((2*np.ones(n),np.r_[-np.ones(n-1),0.]))
            native=module.Resident([
                dict(a=export(operator),p=csr(n,1,np.arange(n+1),np.zeros(n),np.ones(n)),
                     r=csr(1,n,[0,n],np.arange(n),np.ones(n)),smoother=np.full(n,.3)),
                dict(a=csr(1,1,[0,1],[0],[2.]),smoother=np.array([.3]))],np.array([.5]))
            base=rng.normal(size=n)+.01
            for scale in (1.,1e-60,1e60):
                rhs=np.ascontiguousarray(scale*base)
                solution,info=native.solve(rhs,1e-12,10000)
                assert info['certificate_gate']==1e-11
                reference=solveh_banded(bands,rhs,lower=True)
                residual=float(np.linalg.norm(operator@solution-rhs)/np.linalg.norm(rhs))
                error=float(np.linalg.norm(solution-reference)/np.linalg.norm(reference))
                record['ill_conditioned_cases'].append(dict(n=n,rhs_scale=scale,original_relative_residual=residual,relative_error=error,native=info))
                if residual>1e-11 or error>1e-8:raise RuntimeError('Ill-conditioned original residual/reference failed')
            native.close()
        record.update(passed=True,zero_rhs_passed=True,nonfinite_rejected=True,released_context_rejected=True)
    except Exception as error:record['failure']=repr(error);raise
    finally:
        if native is not None:native.close()
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        dll.close()

if __name__=='__main__':main()
