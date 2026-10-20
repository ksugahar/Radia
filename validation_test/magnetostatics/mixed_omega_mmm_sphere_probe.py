"""Independent HDiv-MMM cross-check of the Kelvin sphere reaction.

Both formulations receive the applied field, not the other's total field.
Passing the iron's total field as H_ext would double count demagnetization.
This is a shared-physics cross-check, not a coupled topology optimization run.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time


def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--runtime-source',type=Path,required=True)
    p.add_argument('--kelvin-result',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--native-dll-directory',type=Path,action='append',default=[])
    args=p.parse_args()
    import numpy as np
    import ngsolve as ng
    from netgen.occ import Sphere,Pnt,OCCGeometry
    sys.path.insert(0,str(args.runtime_source))
    dlls=[os.add_dll_directory(str(x.resolve())) for x in args.native_dll_directory]
    record=dict(schema='mixed-omega.mmm-sphere-crosscheck.v1',passed=False)
    try:
        import radia
        from radia.vim import Solve
        import radia._radia_pybind as native
        ng.SetNumThreads(1)
        assert ng.__version__=='6.2.2608'
        binary_hash=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()
        assert binary_hash=='1e9f8d9b0457b0798f43d3e753301183fe4a79137bdd28413b56eb0eb9ac6071'
        kelvin=json.loads(args.kelvin_result.read_text())
        assert kelvin['passed'] and kelvin['relative_rms']<=.01
        sphere=Sphere(Pnt(0,0,0),kelvin['sphere_radius']);sphere.mat('iron')
        started=time.perf_counter()
        with ng.TaskManager():
            mesh=ng.Mesh(OCCGeometry(sphere).GenerateMesh(maxh=.22))
            result=Solve(mesh,mu_r=100.,H_ext=ng.CF((0,0,1)),order=1,
                linear_solver='cpp-cg',preconditioner='mass-riesz',tol=1e-10,maxit=4000)
            fes=result['gfM'].space
            u,v=fes.TnT();mass=ng.BilinearForm(fes);mass+=ng.InnerProduct(u,v)*ng.dx;mass.Assemble()
            load=ng.LinearForm(fes);load+=ng.CF((0,0,1))*v*ng.dx;load.Assemble()
            coefficients=result['_m_coefficients']
            state=mass.mat.CreateColVector();state.FV().NumPy()[:]=coefficients
            action=mass.mat.CreateColVector();action.data=mass.mat*state
            original_action=np.array(action.FV().NumPy())/99+np.array(
                result['_charge_gram'].apply_configured_demag(np.ascontiguousarray(coefficients),True))
            rhs=np.array(load.vec.FV().NumPy())
            original_residual=float(np.linalg.norm(original_action-rhs)/np.linalg.norm(rhs))
        m=np.array(result['M_avg'])
        exact=np.array([0.,0.,99*3/102])
        kelvin_m=99*np.mean(np.array(kelvin['actual_H'])[:3],axis=0)
        relative=float(np.linalg.norm(m-exact)/np.linalg.norm(exact))
        difference=float(np.linalg.norm(m-kelvin_m)/np.linalg.norm(exact))
        record.update(ngsolve=ng.__version__,source_path=radia.__file__,native_sha256=binary_hash,
            elements=mesh.ne,iterations=int(result['iters']),M_avg=m.tolist(),exact_M=exact.tolist(),
            kelvin_M=kelvin_m.tolist(),mmm_analytic_relative=relative,mmm_kelvin_relative=difference,
            true_relative_residual=original_residual,
            seconds=time.perf_counter()-started,passed=relative<.01 and difference<.01 and original_residual<1e-6,
            scope='Independent same applied-field/material cross-check; not joint or topology solve')
        assert record['passed'],record
    except Exception as error:
        record['failure']=repr(error)
        raise
    finally:
        for dll in dlls:dll.close()
        args.out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
