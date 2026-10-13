"""Independent NGSolve oracle for the Python-free native HCurl/AMS workflow."""
from pathlib import Path
import argparse,hashlib,json,os,socket,subprocess,sys,time
import numpy as np
import ngsolve as ng
from netgen.csg import unit_cube
from scipy.sparse import coo_matrix,csr_matrix
from scipy.io import loadmat

LABELS={'back','front','left','right','top','bottom'}

def scipy_matrix(a):
    i,j,v=a.COO();result=coo_matrix((v,(i,j)),shape=(a.height,a.width)).tocsr();result.sort_indices();return result

def decoded(a):
    return csr_matrix((a['values'],a['columns'],a['rowptr']),shape=(int(a['rows']),int(a['cols'])))

def relative(a,b):
    denominator=np.linalg.norm(b)
    return float(np.linalg.norm(a-b)/max(denominator,np.finfo(float).tiny))

def main(args):
    root=Path(__file__).resolve().parents[2];out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.native-hcurl-ams-validation.v1',accepted=False,host=socket.gethostname(),
        formulation='SI homogeneous nu=1, positive unit gauge mass, J=(1,0,0), tangential homogeneous Dirichlet',
        supported='lowest-order affine tetrahedral HCurl only',python_role='independent numerical oracle only',
        gpu_python_input_required=False,timing_accepted=False,timing_classification='100 diagnostic',
        source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip())
    engine=None
    try:
        import netgen,matlab.engine
        net=Path(netgen.__file__).parent
        os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
        build=json.loads((root/'build-ams-cuda/build.json').read_text())
        for rel,digest in build['source_sha256'].items():
            if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('Source differs from build receipt')
        for rel,digest in build['sha256'].items():
            if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('Binary differs from build receipt')
        record['build']=build
        ng.SetNumThreads(4)
        with ng.TaskManager():generated=ng.Mesh(unit_cube.GenerateMesh(maxh=args.maxh))
        vol=out/'model.vol';generated.ngmesh.Save(str(vol));mesh=ng.Mesh(str(vol))
        if set(mesh.GetBoundaries())!=LABELS:raise RuntimeError('Versioned six-boundary .vol contract failed')
        record['mesh_sha256']=hashlib.sha256(vol.read_bytes()).hexdigest()
        t=time.perf_counter();fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT()
        mass=ng.BilinearForm(fes);mass+=u*v*ng.dx
        a=ng.BilinearForm(fes);a+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx
        current=20 if args.gamma>0 else 1
        f=ng.LinearForm(fes);f+=ng.CoefficientFunction((current,0,0))*v*ng.dx
        with ng.TaskManager():mass.Assemble();a.Assemble();f.Assemble()
        free=np.asarray(list(fes.FreeDofs()),dtype=bool);f.vec.FV().NumPy()[~free]=0
        G,_=fes.CreateGradient();xyz=np.array([[p[0],p[1],p[2]] for p in mesh.ngmesh.Points()])
        sys.path.insert(0,str(root/'build-ams-cuda'));import _radia_ams_cpu_reference as reference
        nonlinear=None
        if args.gamma>0:
            from nonlinear_resident import cpu_run
            nonlinear=ng.BilinearForm(fes);b2=ng.curl(u)*ng.curl(u)
            nonlinear+=ng.Variation((.5*b2+args.gamma*.25*b2*b2+.5*u*u)*ng.dx)
            cpu,cpu_report=cpu_run(3,fes,nonlinear,f.vec,G,[xyz[:,c] for c in range(3)],reference,4)
            record['cpu_solver']=cpu_report
        else:
            pre=reference.preconditioner(a.mat,G,fes.FreeDofs(),*[xyz[:,c] for c in range(3)],1)
            x=f.vec.CreateVector();solver=ng.CGSolver(a.mat,pre=pre,tol=1e-12,maxiter=4000)
            with ng.TaskManager():x.data=solver*f.vec
            cpu=x.FV().NumPy().copy()
        record['cpu_total_seconds']=time.perf_counter()-t
        A=scipy_matrix(a.mat);M=scipy_matrix(mass.mat);grad=scipy_matrix(G);rhs=f.vec.FV().NumPy().copy()
        def true_residual(coefficients):
            if nonlinear is None:return relative((A@coefficients)[free],rhs[free])
            probe=ng.GridFunction(fes);probe.vec.FV().NumPy()[:]=coefficients;r=probe.vec.CreateVector()
            with ng.TaskManager():nonlinear.Apply(probe.vec,r)
            return relative(r.FV().NumPy()[free],rhs[free])
        cpu_res=true_residual(cpu);record['cpu_true_relative_residual']=cpu_res
        record.update(gamma=args.gamma,current_density=[current,0,0])
        if cpu_res>(1e-9 if args.gamma else 1e-10):raise RuntimeError('Independent CPU solve failed true residual')
        startup=time.perf_counter();engine=matlab.engine.start_matlab('-nodesktop -nosplash')
        record['matlab_startup_seconds']=time.perf_counter()-startup
        engine.addpath(str(Path(__file__).parent),nargout=0)
        engine.native_hcurl(str(vol),str(root/'build-ams-cuda'),str(root/'packages/radia-gpu/matlab'),str(out/'native.mat'),float(args.gamma),nargout=0)
        native=loadmat(out/'native.mat',simplify_cells=True);model=native['model']
        actual=decoded(model['A']);actualG=decoded(model['G'])
        if not np.array_equal(actual.indptr,A.indptr) or not np.array_equal(actual.indices,A.indices):raise RuntimeError('Native CSR topology differs from NGSolve')
        metrics=dict(matrix_relative=relative(actual.data,A.data),mass_relative=relative(native['mass_values'],M.data),
            gradient_relative=float(np.linalg.norm((actualG-grad).data)/max(np.linalg.norm(grad.data),np.finfo(float).tiny)),rhs_relative=relative(native['rhs'],rhs),
            solution_relative=relative(native['solution'],cpu),true_relative_residual=true_residual(native['solution']))
        if not np.array_equal(model['free'],free.astype(np.int32)):raise RuntimeError('Native free DOF mask differs')
        record.update(dofs=fes.ndof,elements=mesh.ne,csr_nonzeros=A.nnz,metrics=metrics,gpu=native['report'])
        if max(metrics[k] for k in ('matrix_relative','mass_relative','gradient_relative','rhs_relative'))>1e-12:raise RuntimeError('Independent FE assembly parity failed')
        if metrics['solution_relative']>1e-7 or metrics['true_relative_residual']>(1e-9 if args.gamma else 1e-10):raise RuntimeError('Independent AMS solution/residual parity failed')
        # Compare vector curls with a volume weighted RMS, never field magnitudes.
        maps=native['maps'];dofs=np.asarray(maps['dofs']).reshape(-1,6);curl=np.asarray(maps['curls']).reshape(-1,6,3);volume=np.asarray(maps['volume'])
        diff=np.einsum('eac,ea->ec',curl,(native['solution']-cpu)[dofs]);base=np.einsum('eac,ea->ec',curl,cpu[dofs])
        metrics['B_relative_RMS']=float(np.sqrt(np.sum(volume[:,None]*diff**2)/np.sum(volume[:,None]*base**2)))
        if metrics['B_relative_RMS']>1e-7:raise RuntimeError('Independent vector field parity failed')
        record['accepted']=True
    except BaseException as error:record['error']=repr(error);raise
    finally:
        if engine is not None:engine.quit()
        def encode(value):
            if isinstance(value,np.ndarray):return value.tolist()
            if isinstance(value,np.generic):return value.item()
            raise TypeError(type(value).__name__)
        (out/'result.json').write_text(json.dumps(record,indent=2,default=encode)+'\n',encoding='utf-8')
        print(json.dumps({k:v for k,v in record.items() if k not in ('build','gpu')},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--maxh',type=float,default=.35);p.add_argument('--gamma',type=float,default=0);main(p.parse_args())
