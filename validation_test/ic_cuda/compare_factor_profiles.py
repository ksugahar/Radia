"""Compare profiled GPU IC actions with the independent Python CPU IC."""
from pathlib import Path
import argparse,json,os,sys,hashlib

def main(a):
    os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    import netgen
    net=Path(netgen.__file__).parent
    handles=[os.add_dll_directory(str(p)) for p in (a.job/'runtime',net,net.parent/'ngsolve_openblas') if p.is_dir()]
    import numpy as np,ngsolve as ng
    from scipy.io import loadmat
    sys.path.insert(0,str(a.job/'cpu'));import sparsesolv_ngsolve as ic
    owner=json.loads((a.job/'owner.json').read_text())
    if hashlib.sha256(Path(ic.__file__).read_bytes()).hexdigest()!=owner['cpu_ic_origin']['sha256']:raise RuntimeError('CPU IC source identity mismatch')
    ng.SetNumThreads(1);records=[]
    for case in json.loads((a.job/'meshes.json').read_text()):
        target=case['target_dofs'];mesh=ng.Mesh(str(a.job/case['mesh']))
        space=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=space.TnT()
        form=ng.BilinearForm(ng.curl(u)*ng.curl(v)*ng.dx+u*v*ng.dx)
        rhs=ng.LinearForm(ng.CoefficientFunction((1,0,0))*v*ng.dx)
        with ng.TaskManager():
            form.Assemble();rhs.Assemble();pre=ic.ICPreconditioner(form.mat,freedofs=space.FreeDofs(),shift=1.05)
            applied=rhs.vec.CreateVector();applied.data=pre*rhs.vec
        expected=applied.FV().NumPy().copy();values=[]
        for folder in (a.baseline,a.optimized):
            value=np.asarray(loadmat(folder/(str(target)+'.mat'),simplify_cells=True)['preconditioned']).ravel()
            if not np.all(np.isfinite(value)):raise RuntimeError('Nonfinite GPU IC factor action')
            values.append(value)
        errors=[float(np.linalg.norm(x-expected)/np.linalg.norm(expected)) for x in values]
        difference=float(np.linalg.norm(values[0]-values[1])/np.linalg.norm(values[0]))
        if not np.all(np.isfinite(errors+[difference])) or max(errors)>1e-8 or difference>1e-10:raise RuntimeError('Independent factor-action parity failed')
        records.append(dict(total_dofs=space.ndof,baseline_cpu_relative=errors[0],optimized_cpu_relative=errors[1],optimized_baseline_relative=difference))
    result=dict(schema='radia.gpu.ic-factor-action-parity.v1',accepted=True,cases=records,relative_limit=1e-8,baseline_difference_limit=1e-10,cpu_threads=1,cpu_origin=owner['cpu_ic_origin'])
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('job','baseline','optimized','out'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
