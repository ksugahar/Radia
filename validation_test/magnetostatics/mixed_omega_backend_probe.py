"""Compare explicit CSR backends on the production matching-trace operator.

Run with an isolated NGSolve runtime and, only for GPU selection, an installed radia-simulink wheel.
The Python-only formulation is loaded by source path, without borrowing an
incompatible Radia native binary. This is formulation/CSR-binding validation,
not acceptance of the full Radia native distribution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import types


def main():
    from radia.ngsolve_sparse import real_sparse_matrix_type
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--maxh', type=float, default=.65)
    parser.add_argument('--order', type=int, choices=(1,2,3), default=1)
    parser.add_argument('--backend',choices=('cpu','gpu'),default='cpu')
    parser.add_argument('--gpu-dll-directory',type=Path,action='append',default=[])
    args = parser.parse_args()
    started = time.perf_counter()
    import numpy as np
    import ngsolve as ng
    from scipy.sparse import coo_matrix
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, Z
    dlls=[]
    if args.backend=='gpu':
        dlls=[os.add_dll_directory(str(x.resolve())) for x in args.gpu_dll_directory]
        import radia_simulink
        from radia_simulink.csr import solve_spd
    elif args.gpu_dll_directory:raise ValueError('GPU DLL directories require --backend gpu')
    assert ng.__version__ == '6.2.2608', ng.__version__
    root = Path(__file__).resolve().parents[2] / 'src' / 'radia'
    namespace = types.ModuleType('radia')
    namespace.__path__ = [str(root)]
    sys.modules['radia'] = namespace
    from radia.kelvin_solver import solve_magnetostatic_matching_trace_total_reduced_omega as solve
    ng.SetNumThreads(1)
    result = dict(schema='mixed-omega.csr-backend-probe.v1', passed=False,
                  ngsolve=ng.__version__, ngsolve_path=ng.__file__,
                  package_manifest=radia_simulink.verify() if args.backend=='gpu' else None, threads=1, backend=args.backend,
                  formulation_sha256=hashlib.sha256((root/'kelvin_solver.py').read_bytes()).hexdigest(),
                  scope='Python formulation and direct pybind11 CSR; no full Radia native ABI acceptance')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    try:
        total = Box(Pnt(-1,-1,-1),Pnt(1,1,0)); total.mat('total')
        total.faces.name='outer_total'; total.faces.Max(Z).name='interface'
        reduced=Box(Pnt(-1,-1,0),Pnt(1,1,1)); reduced.mat('reduced')
        reduced.faces.name='outer_reduced'; reduced.faces.Min(Z).name='interface'
        mesh=ng.Mesh(OCCGeometry(Glue([total,reduced])).GenerateMesh(maxh=args.maxh))
        fes=ng.H1(mesh,order=args.order)
        source=-ng.x-.4*ng.y-ng.z
        fields=[]
        for cf in (source,-ng.x-.4*ng.y-.01*ng.z+.05*(ng.x**2-ng.y**2),.05*(ng.x**2-ng.y**2)):
            field=ng.GridFunction(fes)
            if args.order == 1:
                for vertex in mesh.vertices:
                    field.vec[fes.GetDofNrs(ng.NodeId(ng.VERTEX,vertex.nr))[0]]=float(cf(mesh(*vertex.point)))
            else:
                field.Set(cf)
            fields.append(field)
        with ng.TaskManager():
            reference=solve(mesh,ng.CF((1,.4,1)),fields[0],
                mu_r_by_material={'total':100.,'reduced':1.},
                reduced_materials=('reduced',),total_materials=('total',),
                interface_boundary='interface',order=args.order,return_system=True,
                solver='cg',cg_preconditioner='local',cg_tolerance=1e-13,
                surface_dirichlet={'total':{'outer_total':fields[1]},'reduced':{'outer_reduced':fields[2]}})
            a=reference['system']['bilinear_form'].mat
            f=reference['system']['linear_form'].vec
            mask=np.array(list(reference['fes'].FreeDofs()),dtype=bool)
            ii,jj,vv=a.COO()
            matrix=coo_matrix((vv,(ii,jj)),shape=(a.height,a.width)).tocsr()
            lift=np.array(reference['solution'].vec.FV().NumPy(),copy=True); lift[mask]=0
            rhs=(np.array(f.FV().NumPy())-matrix@lift)[mask]
            compact=matrix[mask][:,mask].tocoo()
            native=real_sparse_matrix_type().CreateFromCOO(compact.row.tolist(),compact.col.tolist(),compact.data.tolist(),int(mask.sum()),int(mask.sum()))
            b=native.CreateColVector(); b.FV().NumPy()[:]=rhs
            t=time.perf_counter(); inverse=native.Inverse(inverse='sparsecholesky')
            direct=b.CreateVector();direct.data=inverse*b
            cpu=np.array(direct.FV().NumPy(),copy=True)
            cpu_r=float(np.linalg.norm(compact.tocsr()@cpu-rhs)/np.linalg.norm(rhs))
            if not np.isfinite(cpu_r) or cpu_r>1e-6:raise RuntimeError('CPU true residual failed')
            ci=dict(backend='cpu',solver='sparsecholesky',true_relative_residual=cpu_r)
            cpu_seconds=time.perf_counter()-t
        cases=[]
        for repetition in range(3):
            t=time.perf_counter()
            if args.backend=='gpu':gpu,info=solve_spd(native,b,backend='gpu',tolerance=1e-12)
            else:gpu,info=cpu.copy(),ci.copy()
            cases.append(dict(repetition=repetition,transfer_inclusive_seconds=time.perf_counter()-t,info=info,
                cpu_relative_error=float(np.linalg.norm(gpu-cpu)/np.linalg.norm(cpu))))
        reconstructed=ng.GridFunction(reference['fes'])
        reconstructed.vec.FV().NumPy()[:]=lift
        reconstructed.vec.FV().NumPy()[mask]=gpu
        delta=ng.grad(reconstructed)-ng.grad(reference['solution'])
        with ng.TaskManager():
            field_error=float(np.sqrt(ng.Integrate(ng.InnerProduct(delta,delta),mesh)/
                ng.Integrate(ng.InnerProduct(reference['H_cf'],reference['H_cf']),mesh)))
        np.savez(args.out.with_suffix('.npz'),rows=compact.row,columns=compact.col,
                 values=compact.data,rhs=rhs,cpu_solution=cpu,gpu_solution=gpu)
        result.update(elements=mesh.ne,order=args.order,ndof=reference['fes'].ndof,free_dofs=int(mask.sum()),
            cpu_seconds=cpu_seconds,cpu_info=ci,gpu_cases=cases,
            gpu_vs_production_field_relative_rms=field_error,
            production_residual=reference['linear_residual_relative'],
            production_cpu_cg_iterations=reference['iterations'],
            production_cpu_phase_seconds=reference['phase_timings_seconds'],
            elapsed_seconds=time.perf_counter()-started,
            passed=field_error<1e-8 and max(c['cpu_relative_error'] for c in cases)<1e-8)
        assert result['passed'], result
    finally:
        for dll in dlls:dll.close()
        args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
