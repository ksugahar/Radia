"""Public P2-sphere HDiv whole-solve parity against independent CPU pybind11."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,sys,shutil

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    for name in ('projection','hmatrix','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--order',type=int,choices=[1,2],default=1)
    p.add_argument('--geometry',choices=['affine','curved'],default='curved');a=p.parse_args()
    if a.geometry=='curved' and a.order==3:p.error('Existing curved charge semantics support orders 1 and 2')
    a.projection=a.projection.resolve();a.hmatrix=a.hmatrix.resolve();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
    result=dict(accepted=False,scope='Bounded public tetrahedral HDiv, dense near-rule Gram, linear isotropic material',geometry=a.geometry,order=a.order,a100_used=False,release_accepted=False,general_hdiv_validation_complete=False,timing_accepted=False)
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        import numpy as np
        import ngsolve as ng
        from netgen.csg import CSGeometry,Sphere,Pnt
        from scipy.io import savemat,loadmat
        import radia.vim._vim as cpu
        ng.SetNumThreads(1)
        root=Path(__file__).resolve().parents[2]
        if a.geometry=='curved':
            geometry=CSGeometry();geometry.Add(Sphere(Pnt(0,0,0),1))
            mesh=ng.Mesh(geometry.GenerateMesh(maxh=3,curvaturesafety=.5,segmentsperedge=1));mesh.Curve(2)
        else:mesh=ng.Mesh(str(root/'tests/fixtures/beam/affine_field_tetra.vol'))
        mesh_path=a.out/'sphere.vol';mesh.ngmesh.Save(str(mesh_path));mesh=ng.Mesh(str(mesh_path))
        with ng.TaskManager():
            fes=ng.HDiv(mesh,order=a.order)
            B,G,mass=cpu.build_charge_gram(fes,intorder=4,curve_order=2 if a.geometry=='curved' else 0,curve_gauss=8,ho_far_factor=1e30,gram_backend='exact-dense',exact_dense_memory_mb=64)
            load=ng.LinearForm(fes);load+=ng.CoefficientFunction((1.,.2,.1))*fes.TestFunction()*ng.dx;load.Assemble();rhs=np.array(load.vec.FV().NumPy())
            reference=G.solve_configured_linear_material_mass_riesz(2.5,rhs,1e-11,1000)
            demag=np.column_stack([B.T@np.asarray(G.matvec_sym(B[:,k].toarray().ravel())) for k in range(fes.ndof)])
            system=demag+2.5*mass.toarray();expected=np.asarray(reference['m'])
            means=[]
            for c in range(3):
                form=ng.LinearForm(fes);form+=fes.TestFunction()[c]*ng.dx;form.Assemble();means.append(np.array(form.vec.FV().NumPy())/ng.Integrate(1,mesh))
        means=np.asarray(means);cpu_residual=float(np.linalg.norm(rhs-system@expected)/np.linalg.norm(rhs))
        if cpu_residual>1e-10:raise RuntimeError('Independent CPU true residual failed')
        fixture=a.out/'cpu.mat';savemat(fixture,dict(B=B,mass=mass,rhs=rhs[:,None],system=system,expected=expected[:,None],mean_operator=means))
        result.update(cells=mesh.ne,dofs=fes.ndof,charges=B.shape[0],cpu_true_residual=cpu_residual,cpu_iterations=int(reference['iters']),cpu_vim_sha256=sha(cpu.__file__),input_sha256=sha(fixture),mesh_sha256=sha(mesh_path),source_sha256=sha(__file__),projection_sha256=sha(a.projection/'radia_cuda_projection.mexw64'),hmatrix_sha256=sha(a.hmatrix/'radia_cuda_hmatrix.mexw64'))
        prefix=Path(sys.base_prefix);dep=Path('C:/temp/hacapk-native-deps-20261013/Library/bin')
        os.environ['PATH']=os.pathsep.join(map(str,[dep,prefix,prefix/'bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs']))+os.pathsep+os.environ['PATH']
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');result['matlab']=eng.version();result['ngsolve']=ng.__version__
        matlab_source=a.out/'matlab';shutil.copytree(root/'packages/radia-gpu/matlab',matlab_source)
        result['matlab_source_sha256']={f.relative_to(matlab_source).as_posix():sha(f) for f in matlab_source.rglob('*.m')}
        eng.addpath(str(a.projection),str(a.hmatrix),str(matlab_source),nargout=0)
        eng.workspace['mesh_path']=str(mesh_path);eng.workspace['degree']=float(a.order);eng.workspace['fixture_path']=str(fixture);eng.workspace['output_path']=str(a.out/'native.mat')
        eng.eval("cpu=load(fixture_path);setup=radia_cuda_projection('assemble_hdiv_elements',mesh_path,degree,[1 .2 .1]);mass_error=norm(setup.mass-full(cpu.mass),'fro')/norm(full(cpu.mass),'fro');charge_error=norm(setup.charge_map-full(cpu.B),'fro')/norm(full(cpu.B),'fro');rhs_error=norm(setup.rhs-cpu.rhs)/norm(cpu.rhs);assert(max([mass_error charge_error rhs_error])<1e-11);solution=radia.gpu.solveHDivElements(mesh_path,degree,[1 .2 .1],2.5,1e-11,1000);coefficient_error=norm(solution.coefficients-cpu.expected)/norm(cpu.expected);physical_error=norm(solution.mean_magnetization-cpu.mean_operator*cpu.expected)/norm(cpu.mean_operator*cpu.expected);independent_true_residual=norm(cpu.rhs-cpu.system*solution.coefficients)/norm(cpu.rhs);assert(radia_cuda_hmatrix('live_count')==0);save(output_path,'solution','mass_error','charge_error','rhs_error','coefficient_error','physical_error','independent_true_residual','-v7');assert(max([coefficient_error physical_error])<1e-8);assert(independent_true_residual<1e-10);",nargout=0)
        native=loadmat(a.out/'native.mat',simplify_cells=True)
        result['native']={key:float(native[key]) for key in ('mass_error','charge_error','rhs_error','coefficient_error','physical_error','independent_true_residual')};result['solution']=native['solution'];result['accepted']=True
    except BaseException as exc:result['error']=f'{type(exc).__name__}: {exc}';raise
    finally:
        if (a.out/'native.mat').exists() and 'native' not in result:
            from scipy.io import loadmat
            native=loadmat(a.out/'native.mat',simplify_cells=True)
            result['native']={key:float(native[key]) for key in ('mass_error','charge_error','rhs_error','coefficient_error','physical_error','independent_true_residual')}
            result['solution']=native['solution']
        if eng is not None:
            try:eng.quit()
            except BaseException as exc:result['release_error']=str(exc)
        (a.out/'result.json').write_text(json.dumps(result,indent=2,default=lambda value:value.tolist())+'\n',encoding='utf-8')

if __name__=='__main__':main()
