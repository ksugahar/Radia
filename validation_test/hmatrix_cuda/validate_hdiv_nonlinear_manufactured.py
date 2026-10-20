"""Validate manufactured nonlinear HDiv roots against independent CPU factors.

NGSolve owns mapped FE staging; the analytic cubic material oracle is independent
of the shared GPU Newton core. This is not physical excitation or release acceptance.
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
    record=dict(accepted=False,scope='Manufactured nonlinear connected tet/hex HDiv shared Newton; independent original CPU factor and analytic material oracle; not physical excitation or mixed MMM acceptance',cases=[],release_accepted=False)
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
        root=args.cpu_root
        builder_path=root/'tests/test_hdiv_connected_meshes.py'
        spec=importlib.util.spec_from_file_location('hdiv_connected_fixture',builder_path)
        builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
        ng.SetNumThreads(1)
        record.update(cpu_oracle_arithmetic='independent compensated CPU factors and analytic cubic material',ngsolve=ng.__version__,native_source_commit=json.loads((args.native_build/'build.json').read_text())['source_commit'],
            binary_sha256=sha(native),cpu_source_commit=args.cpu_source_commit,cpu_binary_sha256=sha(cpu_binary),
            fixture_source_sha256=sha(builder_path),validator_sha256=sha(__file__),relative_solution_error_limit=1e-10,independent_original_residual_limit=1e-11)
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
                    from scipy.sparse import coo_matrix
                    u,v=fes.TnT();mass_form=ng.BilinearForm(fes,symmetric=True)
                    mass_form+=ng.InnerProduct(u,v)*ng.dx;mass_form.Assemble()
                    ii,jj,vv=mass_form.mat.COO()
                    mass=coo_matrix((vv,(ii,jj)),shape=(fes.ndof,fes.ndof)).tocsr();mass.sort_indices()
                    obj=gpu._ResidentHDiv(snapshot['leaves'],snapshot['permutation'],snapshot['scale'],fes.ndof,
                        B.indptr.astype(np.uint64),B.indices.astype(np.uint64),np.ascontiguousarray(B.data),backend='gpu')
                    from scipy.sparse import csc_matrix
                    staging=gpu._MaterialStagingContext(str(mesh_path),order)
                    payload=staging.quadrature()
                    E=csc_matrix((payload['evaluation_values'],payload['evaluation_row_indices'],payload['evaluation_column_offsets']),shape=(payload['evaluation_rows'],payload['dofs'])).tocsr();E.sort_indices()
                    weights=np.asarray(payload['weights'])
                    law=dict(knots=np.array([0.,1.,2.]),coefficients=np.array([1.,0.,1.,0.,1.,3.,4.,2.]),integral_base=np.array([0.,.75,6.]),cap=1.5,barrier_slope=7.,last_field=4.875,initial_reluctivity=1.)
                    material=gpu._ResidentBHMaterial(E.indptr.astype(np.uint64),E.indices.astype(np.uint64),np.ascontiguousarray(E.data),weights,fes.ndof,law,backend='gpu')
                    rng=np.random.default_rng(31)
                    target=np.ascontiguousarray(.03*rng.standard_normal(fes.ndof))
                    def gradient(x):
                        field=(E@x).reshape(3,-1).T
                        amplitude=np.linalg.norm(field,axis=1)
                        if np.max(amplitude)>=1.5:raise RuntimeError('Analytic cubic fixture left declared band')
                        return E.T@(weights[:,None]*(1+amplitude**2)[:,None]*field).T.ravel()
                    rhs=np.ascontiguousarray(G._apply_configured_demag_compensated(target,16)+gradient(target))
                    try:
                        for pre in ['reference-diagonal','reference-chebyshev-graph']:
                            solution=obj.solve_nonlinear_material(material,mass.indptr.astype(np.uint64),mass.indices.astype(np.uint64),np.ascontiguousarray(mass.data),rhs,np.zeros(fes.ndof),5e-14,5e-14,40,10000,pre)
                            x=np.asarray(solution['x'])
                            coefficient=float(np.linalg.norm(x-target)/np.linalg.norm(target))
                            residual=float(np.linalg.norm(G._apply_configured_demag_compensated(x,16)+gradient(x)-rhs)/np.linalg.norm(rhs))
                            record['cases'].append(dict(topology=topology,order=order,dofs=fes.ndof,preconditioner=pre,coefficient_error=coefficient,independent_original_residual=residual,reported_residual=solution['relative_residual'],newton_iterations=solution['newton_iterations'],graph_launches=solution['preconditioner_graph_launches'],mesh_sha256=sha(mesh_path)))
                            if coefficient>1e-10 or residual>1e-11:raise RuntimeError('Frozen nonlinear FE accuracy gate failed')
                    finally:
                        material.close();obj.close();staging.close()
                    if obj.stats()['owned_device_bytes'] or not material.stats()['closed']:raise RuntimeError('Resident release failed')
        record['accepted']=True
    except BaseException as e:record['error']=repr(e);raise
    finally:
        (args.out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        for handle in handles:handle.close()


if __name__=='__main__':main()
