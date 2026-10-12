"""Sparse GPU-resident nonlinear magnetostatics versus independent Python FE.

NGSolve owns orientation and mapped curls. GPU assembly stores six mapped curl
vectors and 36 CSR slots per affine tetrahedron, never a dense global operator.
CPU uses its normal TaskManager-native AMS build, with no hierarchy export.
"""
from pathlib import Path
import argparse, hashlib, importlib.metadata as md, json, os, socket, statistics, subprocess, sys, time
import numpy as np
import ngsolve as ng
from netgen.csg import unit_cube
from scipy.sparse import coo_matrix, csr_matrix
from scipy.io import savemat, loadmat

def csr(a):
    i,j,v=a.COO();m=coo_matrix((v,(i,j)),shape=(a.height,a.width)).tocsr();m.sort_indices()
    return dict(rows=float(m.shape[0]),cols=float(m.shape[1]),rowptr=m.indptr.astype(np.int32),columns=m.indices.astype(np.int32),values=m.data),m

def mapped_tets(mesh,fes,pattern):
    ne=mesh.ne;dofs=np.empty((ne,6),dtype=np.int32);curls=np.empty((ne,6,3));volume=np.empty(ne);slots=np.empty((ne,6,6),dtype=np.int32)
    rule=ng.IntegrationRule(ng.ET.TET,1)
    if len(rule)!=1:raise RuntimeError('Affine tetrahedron integration rule changed')
    for k,el in enumerate(mesh.Elements()):
        if el.type!=ng.ET.TET:raise RuntimeError('Only affine tetrahedra are supported')
        trafo=mesh.GetTrafo(el);local=fes.GetDofNrs(el)
        if trafo.curved or len(local)!=6:raise RuntimeError('Only lowest-order affine HCurl tet6 is supported')
        mip=trafo(rule[0]);dofs[k]=local;curls[k]=np.asarray(fes.GetFE(el).CalcCurlShape(mip));volume[k]=rule[0].weight*mip.measure
        for a,i in enumerate(local):
            start,end=pattern.indptr[i:i+2];row=pattern.indices[start:end];position=np.searchsorted(row,local)
            if np.any(position>=len(row)) or not np.array_equal(row[position],local):raise RuntimeError('Element CSR pattern mismatch')
            slots[k,a]=start+position
    return dofs,curls,volume,slots

def cpu_run(route,fes,form,rhs,gradient,coords,bridge,threads):
    x=ng.GridFunction(fes);free=np.asarray(list(fes.FreeDofs()),dtype=bool);history=[];start=time.perf_counter()
    setup=0;assembly=0;solving=0;builds=0;pre=None;first=0;rebuild_next=False
    for outer in range(30):
        t=time.perf_counter();residual=x.vec.CreateVector()
        with ng.TaskManager():form.Apply(x.vec,residual)
        residual.data-=rhs;residual.FV().NumPy()[~free]=0
        rel=np.linalg.norm(residual.FV().NumPy())/np.linalg.norm(rhs.FV().NumPy())
        print(f'CPU route={route} outer={outer} residual={rel:.3e}',flush=True)
        if rel<=1e-9:break
        with ng.TaskManager():form.AssembleLinearization(x.vec)
        assembly+=time.perf_counter()-t
        if pre is None or route==1 or (route==3 and rebuild_next):
            t=time.perf_counter();pre=bridge.preconditioner(form.mat,gradient,fes.FreeDofs(),*coords,1)
            setup+=time.perf_counter()-t;builds+=1;first=0
        solver=ng.CGSolver(form.mat,pre=pre,tol=1e-12,maxiter=4000)
        delta=x.vec.CreateVector();t=time.perf_counter()
        with ng.TaskManager():delta.data=solver*residual
        solving+=time.perf_counter()-t
        check=residual.CreateVector()
        with ng.TaskManager():check.data=form.mat*delta-residual
        linrel=np.linalg.norm(check.FV().NumPy()[free])/np.linalg.norm(residual.FV().NumPy()[free])
        if linrel>1e-10:raise RuntimeError('CPU linear true residual failed')
        iterations=int(solver.GetSteps())
        if first==0:first=iterations
        rebuild_next=iterations>max(30,2*first)
        with ng.TaskManager():old=form.Energy(x.vec)-ng.InnerProduct(rhs,x.vec)
        alpha=1.;slope=-ng.InnerProduct(residual,delta);trial=x.vec.CreateVector()
        while True:
            trial.data=x.vec-alpha*delta
            with ng.TaskManager():trial_energy=form.Energy(trial)-ng.InnerProduct(rhs,trial)
            if trial_energy<=old+1e-4*alpha*slope:break
            # Near roundoff, an energy subtraction cannot resolve descent.
            # Require actual residual reduction; convergence gates are unchanged.
            if abs(alpha*slope)<=64*np.finfo(float).eps*max(1.,abs(old)):
                trial_r=residual.CreateVector()
                with ng.TaskManager():form.Apply(trial,trial_r)
                trial_r.data-=rhs
                trial_rel=np.linalg.norm(trial_r.FV().NumPy()[free])/np.linalg.norm(rhs.FV().NumPy()[free])
                if trial_rel<.9*rel or trial_rel<=1e-9:break
            alpha*=.5
            if alpha<2**-20:raise RuntimeError('CPU line search failed')
        x.vec.data=trial;history.append([outer,float(rel),iterations,float(linrel),alpha])
    if rel>1e-9:raise RuntimeError('CPU nonlinear residual failed')
    return x.vec.FV().NumPy().copy(),dict(route=route,total_seconds=time.perf_counter()-start,setup_seconds=setup,
        assembly_seconds=assembly,solve_seconds=solving,builds=builds,history=history,final_residual=float(rel),threads=threads)

def main(args):
    root=Path(__file__).resolve().parents[2];out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.ams-nonlinear-sparse-resident.v2',accepted=False,host=socket.gethostname(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        timing_classification='build-host diagnostic; not formal performance acceptance',a100_used=False,
        material='SI: H=(1+4*|B|^2)B; positive unit gauge mass',
        excitation='manufactured FE interpolation of 0.3*(sin(pi*y)*sin(pi*z),sin(pi*x)*sin(pi*z),sin(pi*x)*sin(pi*y))',
        gpu_runtime_python_adapter_required=False,python_role='independent FE oracle and validation input generation',
        routes=args.routes,cpu_threads=args.threads,gpu_host_setup_threads=args.setup_threads,
        assembly_storage='sparse CSR plus 6 mapped curl vectors and 36 CSR slots per tetrahedron',
        gpu_residency='nonlinear CSR buffers and cached AMS/PCG workspace; explicit tangent download only for hierarchy rebuild',
        near_roundoff_line_search='energy decrement <=64*eps*max(1,abs(E)): require 10 percent residual decrease or final nonlinear gate',
        timing_boundary='numeric FE input to final solution; shared FE preparation and MATLAB startup separately; GPU total includes upload',
        nonlinear_tolerance=1e-9,linear_true_tolerance=1e-10,cpu_cg_internal_tolerance=1e-12,independent_parity_limit=1e-7)
    eng=None;preprocessing=time.perf_counter()
    try:
        import netgen
        net=Path(netgen.__file__).parent
        os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
        sys.path.insert(0,str(root/'build-ams-cuda'));import _radia_ams_cpu_reference as bridge
        build=json.loads((root/'build-ams-cuda/build.json').read_text());record['build']=build
        for relative,digest in build['source_sha256'].items():
            if hashlib.sha256((root/relative).read_bytes()).hexdigest()!=digest:raise RuntimeError('Native source changed since build')
        for relative,digest in build['sha256'].items():
            if hashlib.sha256((root/relative).read_bytes()).hexdigest()!=digest:raise RuntimeError('Native binary hash mismatch')
        record['gpu_binary_sha256']=build['sha256']['build-ams-cuda/radia_cuda_ams.mexw64']
        record['versions']={p:md.version(p) for p in ('ngsolve','numpy','scipy','matlabengine')}
        ng.SetNumThreads(args.threads)
        with ng.TaskManager():mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=args.maxh))
        mesh.ngmesh.Save(str(out/'model.vol'));mesh=ng.Mesh(str(out/'model.vol'))
        if set(mesh.GetBoundaries())!={'back','front','left','right','top','bottom'}:raise RuntimeError('Unit-cube boundary contract failed')
        record['mesh_sha256']=hashlib.sha256((out/'model.vol').read_bytes()).hexdigest()
        record['mesh_route']='explicit Netgen alternative -> checked six-boundary .vol -> NGSolve reload'
        fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT();n=fes.ndof
        free=np.asarray(list(fes.FreeDofs()),dtype=bool)
        mass=ng.BilinearForm(fes);mass+=u*v*ng.dx
        linear=ng.BilinearForm(fes);linear+=(ng.curl(u)*ng.curl(v)+u*v)*ng.dx
        with ng.TaskManager():mass.Assemble();linear.Assemble()
        A,pattern=csr(linear.mat);M,mass_pattern=csr(mass.mat)
        if not np.array_equal(pattern.indptr,mass_pattern.indptr) or not np.array_equal(pattern.indices,mass_pattern.indices):raise RuntimeError('Mass/tangent CSR patterns differ')
        dofs,curls,volume,slots=mapped_tets(mesh,fes,pattern)
        local=volume[:,None,None]*np.einsum('eak,ebk->eab',curls,curls)
        reconstructed=M['values'].copy();np.add.at(reconstructed,slots.ravel(),local.ravel())
        error=np.linalg.norm(reconstructed-A['values'])/np.linalg.norm(A['values'])
        if error>1e-12:raise RuntimeError('NGSolve mapped curl reconstruction mismatch')
        del reconstructed,local
        record.update(dofs=n,free_dofs=int(sum(free)),elements=mesh.ne,csr_nonzeros=len(A['values']),map_assembly_relative_error=float(error))
        print(f'Sparse model ready: DOFs={n}, elements={mesh.ne}, nnz={len(A["values"])}',flush=True)
        form=ng.BilinearForm(fes);B=ng.curl(u);b2=B*B
        form+=ng.Variation((0.5*b2+b2*b2+0.5*u*u)*ng.dx)
        exact=ng.GridFunction(fes)
        cf=.3*ng.CoefficientFunction((ng.sin(ng.pi*ng.y)*ng.sin(ng.pi*ng.z),ng.sin(ng.pi*ng.x)*ng.sin(ng.pi*ng.z),ng.sin(ng.pi*ng.x)*ng.sin(ng.pi*ng.y)))
        with ng.TaskManager():exact.Set(cf)
        exact.vec.FV().NumPy()[~free]=0;rhs=exact.vec.CreateVector()
        with ng.TaskManager():form.Apply(exact.vec,rhs)
        rhs.FV().NumPy()[~free]=0;gradient,_=fes.CreateGradient()
        coords=np.asarray([[mesh.ngmesh.Points()[i+1][k] for i in range(mesh.nv)] for k in range(3)])
        record['shared_fe_setup_seconds']=time.perf_counter()-preprocessing
        cpu_samples=[];record["cpu_samples"]=cpu_samples;cpu=None
        for route in args.routes:
            actual,sample=cpu_run(route,fes,form,rhs,gradient,coords,bridge,args.threads);cpu_samples.append(sample)
            if cpu is None:cpu=actual
            elif np.linalg.norm(actual-cpu)/np.linalg.norm(cpu)>1e-7:raise RuntimeError('CPU reuse changes solution')
        record['cpu_samples']=cpu_samples;np.save(out/'cpu_solution.npy',cpu)
        model=dict(schema='radia.ams-native-input.v1',A=A,G=csr(gradient)[0],coordinates=coords.T.copy(),free=free.astype(np.int32),cycle=1.)
        fem=dict(schema='radia.ams-nonlinear-tet6.v1',mass=M,elements=float(mesh.ne),dofs=dofs.reshape(-1,1),slots=slots.reshape(-1,1),
            curls=curls.reshape(-1,1),volume=volume[:,None],free=free.astype(np.int32)[:,None],rhs=rhs.FV().NumPy()[:,None],nu0=1.,gamma=4.)
        savemat(out/'input.mat',dict(model=model,fem=fem,rhs=rhs.FV().NumPy()[:,None],routes=np.asarray(args.routes,dtype=float)))
        record['input_sha256']=hashlib.sha256((out/'input.mat').read_bytes()).hexdigest()
        import matlab.engine
        startup=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash');record['matlab_startup_seconds']=time.perf_counter()-startup
        eng.addpath(str(Path(__file__).parent),nargout=0)
        eng.nonlinear_resident(str(out/'input.mat'),str(out/'gpu.mat'),str(root/'build-ams-cuda'),str(root/'packages/radia-gpu/matlab'),float(args.setup_threads),nargout=0)
        if args.contract:
            eng.resident_contract(str(out/'input.mat'),str(out/'contract.mat'),str(root/'build-ams-cuda'),str(root/'packages/radia-gpu/matlab'),nargout=0)
            record['resident_contract']=loadmat(out/'contract.mat',simplify_cells=True)['report']
        gpu=loadmat(out/'gpu.mat',simplify_cells=True);record['gpu']=[];record['gpu_device']=gpu['gpuInfo']
        for solution,report in zip(gpu['solutions'],gpu['reports']):
            actual=np.asarray(solution).reshape(-1);probe=ng.GridFunction(fes);probe.vec.FV().NumPy()[:]=actual;residual=probe.vec.CreateVector()
            with ng.TaskManager():form.Apply(probe.vec,residual)
            residual.data-=rhs
            solution_error=float(np.linalg.norm(actual-cpu)/np.linalg.norm(cpu))
            rr=float(np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(rhs.FV().NumPy()[free]))
            flux=np.einsum('eak,ea->ek',curls,(actual-cpu)[dofs]);ref=np.einsum('eak,ea->ek',curls,cpu[dofs])
            flux_error=float(np.sqrt(np.sum(volume[:,None]*flux**2)/np.sum(volume[:,None]*ref**2)))
            report['history']=np.asarray(report['history']).tolist();report.update(cpu_solution_relative_error=solution_error,independent_nonlinear_residual=rr,B_relative_RMS_error=flux_error)
            record['gpu'].append(report)
            if max(solution_error,rr,flux_error)>1e-7:raise RuntimeError('Independent CPU/GPU parity failed')
        record['accepted']=True
    except BaseException as error:
        record['error']=repr(error);raise
    finally:
        if eng is not None:eng.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({k:v for k,v in record.items() if k not in ('build','cpu_samples','gpu')},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--maxh',type=float,default=.12)
    p.add_argument('--contract',action='store_true');p.add_argument('--threads',type=int,default=8);p.add_argument('--setup-threads',type=int,default=8)
    p.add_argument('--routes',type=int,nargs='+',choices=[1,2,3],default=[1,2,3,3,2,1]);main(p.parse_args())
