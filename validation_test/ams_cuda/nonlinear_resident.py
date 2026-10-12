"""Small GPU-resident nonlinear magnetostatic correctness study on 100.

Monotone H=(nu0+gamma*|B|^2)B, B=curl(A), with positive gauge mass.
NGSolve owns FE orientation and assembly. Dense quadrature maps intentionally
limit this first correctness case; timings are diagnostic, not scalability claims.
"""
from pathlib import Path
import argparse, json, os, socket, subprocess, sys, time
import numpy as np
import ngsolve as ng
from netgen.csg import unit_cube
from scipy.sparse import coo_matrix
from scipy.io import savemat, loadmat

def csr(a):
    i,j,v=a.COO();m=coo_matrix((v,(i,j)),shape=(a.height,a.width)).tocsr();m.sort_indices()
    return dict(rows=m.shape[0],cols=m.shape[1],rowptr=m.indptr.astype(np.int32),columns=m.indices.astype(np.int32),values=m.data),m

def main(args):
    root=Path(__file__).resolve().parents[2];out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.ams-nonlinear-resident.v1',accepted=False,host=socket.gethostname(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        timing_classification='100 diagnostic; dense quadrature prototype; not performance acceptance',
        python_role='independent FE oracle and validation input generation',a100_used=False,
        material='H=(nu0+gamma*|B|^2)B; nu0=1; gamma=4; positive unit gauge mass',
        gpu_residency='quadrature maps, material, residual, tangent, solution and frozen AMS hierarchy remain GPU arrays; rebuild route downloads tangent values for native host setup',
        limitations=['dense quadrature prototype, not production scaling','AMS workspace currently allocated per linear solve','single cold sample; no speedup acceptance'])
    eng=None
    try:
        import netgen
        net=Path(netgen.__file__).parent
        os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
        sys.path.insert(0,str(root/'build-ams-cuda'));import _radia_ams_snapshot as bridge
        ng.SetNumThreads(1)
        mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=args.maxh));mesh.ngmesh.Save(str(out/'model.vol'))
        fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT();n=fes.ndof
        free=np.asarray(list(fes.FreeDofs()),dtype=bool)
        mass=ng.BilinearForm(fes);mass+=u*v*ng.dx;mass.Assemble()
        linear=ng.BilinearForm(fes);linear+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx;linear.Assemble()
        maps=[];weights=[]
        for el in mesh.Elements():
            fe=fes.GetFE(el);dofs=fes.GetDofNrs(el);trafo=mesh.GetTrafo(el)
            for ip in ng.IntegrationRule(el.type,4):
                mip=trafo(ip);c=np.zeros((3,n));c[:,dofs]=np.asarray(fe.CalcCurlShape(mip)).T
                maps.append(c);weights.append(ip.weight*mip.measure)
        C=np.asarray(maps);w=np.asarray(weights);M=csr(mass.mat)[1].toarray();A,sp=csr(linear.mat)
        reconstructed=M+sum(C[:,k,:].T@(w[:,None]*C[:,k,:]) for k in range(3))
        map_error=np.linalg.norm(reconstructed-sp.toarray())/np.linalg.norm(sp.data)
        record.update(dofs=n,quadrature_points=len(w),map_assembly_relative_error=float(map_error))
        if map_error>1e-12:raise RuntimeError('NGSolve mapped curl reconstruction mismatch')
        form=ng.BilinearForm(fes);B=ng.curl(u);b2=B*B
        form+=ng.Variation((0.5*b2+b2*b2+0.5*u*u)*ng.dx)
        exact=ng.GridFunction(fes);exact.vec.FV().NumPy()[:]=np.random.default_rng(123).normal(size=n)*free*0.07
        rhs=exact.vec.CreateVector();form.Apply(exact.vec,rhs);rhs.FV().NumPy()[~free]=0
        gradient,_=fes.CreateGradient();coords=np.asarray([[mesh.ngmesh.Points()[i+1][k] for i in range(mesh.nv)] for k in range(3)])
        x=ng.GridFunction(fes);history=[];start=time.perf_counter();setup=0
        for outer in range(30):
            residual=x.vec.CreateVector();form.Apply(x.vec,residual);residual.data-=rhs;residual.FV().NumPy()[~free]=0
            rel=np.linalg.norm(residual.FV().NumPy())/np.linalg.norm(rhs.FV().NumPy())
            if rel<=1e-9:break
            form.AssembleLinearization(x.vec)
            t=time.perf_counter();snapshot=bridge.snapshot(form.mat,gradient,fes.FreeDofs(),*coords,np.zeros((1,n)),1,1);setup+=time.perf_counter()-t
            solver=ng.CGSolver(form.mat,pre=snapshot['cpu_preconditioner'],tol=1e-10,maxiter=2000)
            delta=x.vec.CreateVector();delta.data=solver*residual
            check=residual.CreateVector();check.data=form.mat*delta-residual
            linrel=np.linalg.norm(check.FV().NumPy()[free])/np.linalg.norm(residual.FV().NumPy()[free])
            if linrel>1e-6:raise RuntimeError('CPU linear residual')
            alpha=1.;old=form.Energy(x.vec)-ng.InnerProduct(rhs,x.vec);slope=-ng.InnerProduct(residual,delta)
            trial=x.vec.CreateVector()
            while True:
                trial.data=x.vec-alpha*delta
                if form.Energy(trial)-ng.InnerProduct(rhs,trial)<=old+1e-4*alpha*slope:break
                alpha*=.5
                if alpha<2**-20:raise RuntimeError('CPU line search')
            x.vec.data=trial;history.append([outer,float(rel),float(linrel),alpha])
        if rel>1e-9:raise RuntimeError('CPU nonlinear residual')
        cpu=x.vec.FV().NumPy().copy();record['cpu']=dict(total_seconds=time.perf_counter()-start,setup_seconds=setup,history=history,final_residual=float(rel))
        model=dict(schema='radia.ams-native-input.v1',A=A,G=csr(gradient)[0],coordinates=coords.T.copy(),free=free.astype(np.int32),cycle=1)
        rows=np.repeat(np.arange(n),np.diff(sp.indptr));index=(rows+sp.indices*n+1).astype(np.float64)
        savemat(out/'input.mat',dict(model=model,Cx=C[:,0,:],Cy=C[:,1,:],Cz=C[:,2,:],weights=w[:,None],mass=M,rhs=rhs.FV().NumPy()[:,None],nu0=1.,gamma=4.,csrIndex=index[:,None]))
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash');eng.addpath(str(Path(__file__).parent),nargout=0)
        eng.nonlinear_resident(str(out/'input.mat'),str(out/'gpu.mat'),str(root/'build-ams-cuda'),str(root/'packages/radia-gpu/matlab'),nargout=0)
        gpu=loadmat(out/'gpu.mat',simplify_cells=True);record['gpu']=[]
        for solution,report in zip(gpu['solutions'],gpu['reports']):
            actual=np.asarray(solution).reshape(-1);probe=ng.GridFunction(fes);probe.vec.FV().NumPy()[:]=actual
            residual=probe.vec.CreateVector();form.Apply(probe.vec,residual);residual.data-=rhs
            error=float(np.linalg.norm(actual-cpu)/np.linalg.norm(cpu));rr=float(np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(rhs.FV().NumPy()[free]))
            flux=np.einsum('qkn,n->qk',C,actual-cpu);ref=np.einsum('qkn,n->qk',C,cpu)
            fluxerror=float(np.sqrt(np.sum(w[:,None]*flux**2)/np.sum(w[:,None]*ref**2)))
            report['history']=np.asarray(report['history']).tolist();report.update(cpu_solution_relative_error=error,independent_nonlinear_residual=rr,B_relative_RMS_error=fluxerror)
            record['gpu'].append(report)
            if max(error,rr,fluxerror)>1e-7:raise RuntimeError('Independent CPU/GPU parity failed')
        record['accepted']=True
    finally:
        if eng is not None:eng.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--maxh',type=float,default=.35);main(p.parse_args())
