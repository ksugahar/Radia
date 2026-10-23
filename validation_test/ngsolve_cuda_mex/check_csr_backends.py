"""Bounded same-binary CSR checks; CPU-context tracing remains a separate gate."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


def fixtures():
    from radia.ngsolve_sparse import real_sparse_matrix_type
    import ngsolve as ng
    import numpy as np
    cases=[]
    for n in (3,64):
        rows=[];cols=[];values=[]
        for i in range(n):
            rows.append(i);cols.append(i);values.append(4.)
            if i+1<n:
                rows.extend((i,i+1));cols.extend((i+1,i));values.extend((-1.,-1.))
        a=real_sparse_matrix_type().CreateFromCOO(rows,cols,values,n,n)
        b=a.CreateColVector();b.FV().NumPy()[:]=np.linspace(1.,2.,n)
        with ng.TaskManager():
            inv=a.Inverse(inverse='sparsecholesky');x=b.CreateVector();x.data=inv*b
        av,ci,rp=a.CSR()
        cases.append((a,b,dict(row_pointers=np.array(rp,dtype=np.int32,copy=True),
            columns=np.array(ci,dtype=np.int32,copy=True),values=np.array(av,copy=True),
            rhs=np.array(b.FV().NumPy(),copy=True),inverse_diagonal=np.full(n,.25),
            expected=np.array(x.FV().NumPy(),copy=True))))
    return cases


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=('fixtures','python','matlab-cpu','matlab-gpu'),required=True)
    parser.add_argument('--build-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]
    started=time.perf_counter()
    result=dict(schema='radia.csr-backend-validation.v1',accepted=False,phase=args.phase,
                cuda_context_gate='unverified',cases=[])
    try:
        import numpy as np
        cases=fixtures()
        if args.phase=='fixtures':
            sys.path.insert(0,str(root/'packages/radia-simulink/src'))
            from radia_simulink.ams_exchange import write_exchange
            write_exchange(args.output,dict(cases=[values for _,_,values in cases]))
            return
        receipt=json.loads((args.build_dir/'build.json').read_text())
        result['source_commit']=receipt['source_commit']
        if args.phase=='python':
            # Runtime DLL search is confined to this owned validation process.
            # It is evidence for local execution, not a distribution permission.
            matlab=Path(os.environ['RADIA_VALIDATION_MATLAB_ROOT'])
            handle=os.add_dll_directory(str(matlab/'bin/win64'))
            try:
                sys.path.insert(0,str(root/'packages/radia-simulink/src'))
                import radia_simulink.csr as csr
                module=args.build_dir/receipt['python_binding']['binary']
                assert hashlib.sha256(module.read_bytes()).hexdigest()==receipt['python_binding']['binary_sha256']
                spec=importlib.util.spec_from_file_location('radia_simulink._csr_cuda',module)
                extension=importlib.util.module_from_spec(spec);sys.modules[spec.name]=extension;spec.loader.exec_module(extension)
                result['binary_sha256']=receipt['python_binding']['binary_sha256']
                import ctypes as c
                kernel=c.WinDLL('kernel32',use_last_error=True)
                kernel.GetModuleHandleW.argtypes=[c.c_wchar_p];kernel.GetModuleHandleW.restype=c.c_void_p
                kernel.GetModuleFileNameW.argtypes=[c.c_void_p,c.c_wchar_p,c.c_uint32]
                dll=kernel.GetModuleHandleW('cudart64_12.dll')
                name=c.create_unicode_buffer(4096)
                if not dll or not kernel.GetModuleFileNameW(dll,name,len(name)):
                    raise c.WinError(c.get_last_error())
                runtime=Path(name.value)
                if not runtime.samefile(matlab/'bin/win64/cudart64_12.dll'):
                    raise RuntimeError('Validation loaded a CUDA runtime outside the approved bundle')
                result['runtime_identity']=dict(path=str(runtime),sha256=hashlib.sha256(runtime.read_bytes()).hexdigest())

                for a,b,f in cases:
                    for backend in ('cpu','gpu'):
                        begin=time.perf_counter()
                        # Caller owns NGSolve CPU parallelism, not the adapter.
                        import ngsolve as ng
                        with ng.TaskManager():x,info=csr.solve_spd(a,b,backend=backend,tolerance=1e-12,max_iterations=1000)
                        relative=float(np.linalg.norm(x-f['expected'])/np.linalg.norm(f['expected']))
                        assert relative<=1e-11 and info['true_relative_residual']<=1e-11
                        result['cases'].append(dict(backend=backend,relative_l2=relative,
                            true_relative_residual=info['true_relative_residual'],
                            transfer_inclusive_seconds=time.perf_counter()-begin))
                a,b,f=cases[-1]
                for _ in range(100):csr.solve_spd(a,b,backend='gpu',tolerance=1e-12,max_iterations=1000)
                bad=np.array(f['row_pointers']);bad[-1]-=1
                try:extension.solve(bad,f['columns'],f['values'],f['rhs'],f['inverse_diagonal'],1e-12,1000)
                except RuntimeError:pass
                else:raise AssertionError('Invalid CSR was accepted')
                assert not any(n=='matlab' or n.startswith('matlab.') for n in sys.modules),'Python GPU user process imported MATLAB'
            finally:handle.close()
        else:
            import matlab.engine
            eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
            try:
                eng.addpath(str(args.build_dir),str(root/'packages/radia-simulink/matlab'),str(Path(__file__).parent),nargout=0)
                assert hashlib.sha256(Path(eng.which('radia_cuda_csr_pcg')).read_bytes()).hexdigest()==receipt['binary_sha256']
                mode=args.phase.split('-')[1]
                eng.workspace['fixturePath']=str(args.build_dir/'fixtures.h5');eng.workspace['resultPath']=str(args.output);eng.workspace['backendMode']=mode
                result=json.loads(eng.eval('jsonencode(checkCsrBackends(fixturePath,resultPath,backendMode))'))
                result['binary_sha256']=receipt['binary_sha256'];result['source_commit']=receipt['source_commit']
            finally:eng.quit()
        result['accepted']=True
    except BaseException as error:
        result['error']=repr(error)
        raise
    finally:
        if args.phase!='fixtures':
            result['elapsed_seconds']=time.perf_counter()-started
            args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
