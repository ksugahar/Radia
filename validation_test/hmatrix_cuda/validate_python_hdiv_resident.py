"""Physical HDiv matvec parity through the shared direct Python CUDA component.

The CPU owner exports a bounded deep copy of its normalized HACApK factors.
This validates transfer/packing/B^T G B only; CPU assembly remains explicit and
this is not complete Python HDiv-MMM solver acceptance.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-build',type=Path,required=True)
    p.add_argument('--cpu-root',type=Path,required=True)
    p.add_argument('--cpu-source-commit',required=True)
    p.add_argument('--cpu-native-sha256',required=True)
    p.add_argument('--dependency-dir',type=Path,action='append',required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--topology',nargs='+',choices=('tet','hex'),default=['tet','hex'])
    p.add_argument('--order',nargs='+',type=int,choices=(1,2),default=[1,2])
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    record=dict(accepted=False,scope='actual HDiv matvec with owned normalized HACApK factors; CPU assembly explicit; no complete GPU solver acceptance',cases=[],release_accepted=False)
    handles=[]
    try:
        for directory in args.dependency_dir:
            if not directory.is_dir():raise ValueError('Dependency directory missing')
            handles.append(os.add_dll_directory(str(directory.resolve())))
        if subprocess.check_output(['git','-C',str(args.cpu_root),'rev-parse','HEAD'],text=True).strip()!=args.cpu_source_commit:
            raise ValueError('Independent CPU source mismatch')
        cpu_binary=args.cpu_root/'src/radia/_radia_pybind.pyd'
        if sha(cpu_binary)!=args.cpu_native_sha256:raise ValueError('Independent CPU binary mismatch')
        sys.path.insert(0,str(args.cpu_root/'src'))
        import numpy as np
        import ngsolve as ng
        import radia
        from radia.vim._vim import build_charge_gram
        if Path(radia.__file__).resolve().parent!=cpu_binary.parent.resolve():
            raise ValueError('Unexpected CPU import origin')
        receipt=json.loads((args.native_build/'python-build.json').read_text())
        native=args.native_build/receipt['binary']
        if sha(native)!=receipt['binary_sha256']:raise ValueError('GPU binary mismatch')
        if receipt['ngsolve_version']!=ng.__version__:raise ValueError('NGSolve ABI version mismatch')
        spec=importlib.util.spec_from_file_location('_hdiv_cuda',native)
        gpu=importlib.util.module_from_spec(spec);spec.loader.exec_module(gpu)
        if 'matlab.engine' in sys.modules:raise ValueError('Unexpected MATLAB Engine user dependency')
        root=Path(__file__).resolve().parents[2]
        builder_path=root/'tests/test_hdiv_connected_meshes.py'
        spec=importlib.util.spec_from_file_location('hdiv_connected_fixture',builder_path)
        builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
        ng.SetNumThreads(1)
        record.update(ngsolve=ng.__version__,native_source_commit=json.loads((args.native_build/'build.json').read_text())['source_commit'],
            binary_sha256=sha(native),cpu_source_commit=args.cpu_source_commit,cpu_binary_sha256=sha(cpu_binary),
            fixture_source_sha256=sha(builder_path),validator_sha256=sha(__file__),relative_application_error_limit=1e-11)
        for topology in args.topology:
            for order in args.order:
                start=time.perf_counter()
                with ng.TaskManager():
                    mesh=builder.connected_mesh(topology)
                    mesh_path=args.out/f'{topology}-{order}.vol';mesh.ngmesh.Save(str(mesh_path))
                    contract=args.out/f'{topology}-{order}-labels.json'
                    contract.write_text(json.dumps(dict(schema='radia.vol-label-contract.v1',
                        required=dict(materials=['body'],boundaries=['exterior']),
                        allowed=dict(materials=['body'],boundaries=['exterior']))))
                    checker=args.cpu_root/'packages/cubit-mesh-export/src/cubit_mesh_export/check.py'
                    command=[sys.executable,str(checker),str(mesh_path),'--strict-labels','--contract',str(contract),
                        '--report-json',str(args.out/f'{topology}-{order}-check.json')]
                    if topology=='tet':command.append('--tet-only')
                    subprocess.run(command,check=True)
                    mesh=ng.Mesh(str(mesh_path));fes=ng.HDiv(mesh,order=order)
                    B,G,_=build_charge_gram(fes,eps=1e-10,leafsize=64,eta=2.)
                    n=B.shape[0]
                    if n>256 or fes.ndof>512:raise ValueError('Tiny component validation bound exceeded')
                    snapshot=G._resident_leaf_snapshot(16)
                    if snapshot['dimension']!=n:raise ValueError('Snapshot charge dimension mismatch')
                    obj=gpu._ResidentHDiv(snapshot['leaves'],snapshot['permutation'],snapshot['scale'],fes.ndof,
                        B.indptr.astype(np.uint64),B.indices.astype(np.uint64),np.ascontiguousarray(B.data),backend='gpu')
                    errors=[]
                    try:
                        owned=obj.stats()['owned_device_bytes']
                        for seed in (10,20,30):
                            x=np.random.default_rng(seed).standard_normal(fes.ndof)
                            expected=np.asarray(G.apply_configured_demag(x,True))
                            actual=obj.apply(x)
                            error=float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),np.finfo(float).tiny))
                            errors.append(error)
                            if error>1e-11:
                                record['failed_case']=dict(topology=topology,order=order,seed=seed,error=error)
                                np.savez(args.out/'failed-application.npz',x=x,expected=expected,actual=actual,
                                         sigma=snapshot['scale'],B=B.toarray())
                                raise RuntimeError('Independent CPU physical HDiv matvec mismatch')
                        if obj.stats()['owned_device_bytes']!=owned:raise RuntimeError('Resident allocation changed')
                    finally:obj.close()
                    if not obj.stats()['closed'] or obj.stats()['owned_device_bytes']:
                        raise RuntimeError('Resident was not released')
                record['cases'].append(dict(topology=topology,order=order,dofs=fes.ndof,charges=n,mesh_sha256=sha(mesh_path),
                    relative_errors=errors,owned_device_bytes=owned,closed=True,seconds_diagnostic=time.perf_counter()-start))
        record['accepted']=True
    except BaseException as e:record['error']=repr(e);raise
    finally:
        (args.out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        for handle in handles:handle.close()


if __name__=='__main__':main()
