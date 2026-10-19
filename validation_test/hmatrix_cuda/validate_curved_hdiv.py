"""Public P2-sphere HDiv whole-solve parity against independent CPU pybind11."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,sys,shutil
from hdiv_validation_inputs import checked_mesh, dependency_directories, thread_record, cpu_receipt, single_thread_blas, observed_blas, sorted_oracle_csr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    for name in ('projection','hmatrix','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--order',type=int,choices=[1,2],default=1)
    p.add_argument('--assembly-mode',choices=['direct','host-cache'],default='direct')
    p.add_argument('--geometry',choices=['affine','curved'],default='curved')
    for name in ('mesh','contract'):
        p.add_argument('--'+name,type=Path,required=True)
        p.add_argument('--'+name+'-sha256',required=True)
    p.add_argument('--dependency-dir',type=Path,action='append',required=True)
    p.add_argument('--cpu-threads',type=int,default=1)
    p.add_argument('--cpu-build-receipt',type=Path,required=True)
    a=p.parse_args()
    if a.geometry=='curved' and a.order==3:p.error('Existing curved charge semantics support orders 1 and 2')
    a.projection=a.projection.resolve();a.hmatrix=a.hmatrix.resolve();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
    result=dict(accepted=False,scope='Bounded public tetrahedral HDiv, dense near-rule Gram, linear isotropic material',geometry=a.geometry,order=a.order,assembly_mode=a.assembly_mode,a100_used=None,release_accepted=False,general_hdiv_validation_complete=False,timing_accepted=False)
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL')
        directories,dll_handles=dependency_directories(a.dependency_dir)
        result['dependency_directories']=list(map(str,directories))
        result['threads']=thread_record(a.cpu_threads)
        result['mesh_check']=checked_mesh(a.mesh,a.contract,a.mesh_sha256,a.contract_sha256,a.out/'vol_check.json')
        import numpy as np
        import ngsolve as ng
        from hdiv_hdf5 import write_hdiv_h5, read_hdiv_h5
        import radia.vim._vim as cpu
        import radia._radia_pybind as cpu_native
        import radia.vim._capabilities as cpu_capabilities
        import matlab.engine
        result['cpu_reference_provenance'] = dict(
            interpreter=sys.executable, python=sys.version,
            ngsolve_version=ng.__version__,
            native_path=str(Path(cpu_native.__file__).resolve()),
            native_sha256=sha(cpu_native.__file__),
            frontend_sha256=sha(cpu.__file__),
            capabilities_sha256=sha(cpu_capabilities.__file__),
            execution='Python/pybind11, caller-owned NGSolve TaskManager',
            scope='Loaded reference files; hashes alone do not certify native build source provenance')
        result['cpu_build_receipt']=cpu_receipt(a.cpu_build_receipt,cpu_native.__file__)
        result['cpu_build_receipt_sha256']=sha(a.cpu_build_receipt)
        ng.SetNumThreads(a.cpu_threads)
        blas_limit,result['observed_blas_after_ngsolve_configuration']=single_thread_blas()
        root=Path(__file__).resolve().parents[2]
        mesh_path=a.out/'checked.vol';shutil.copyfile(a.mesh,mesh_path)
        if sha(mesh_path)!=a.mesh_sha256:raise RuntimeError('Mesh changed during staging')
        mesh=ng.Mesh(str(mesh_path))
        if bool(mesh.GetCurveOrder()>1)!=(a.geometry=='curved'):raise RuntimeError('Declared geometry/mesh curve order mismatch')
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
        result['observed_blas_after_cpu_solve']=observed_blas()
        means=np.asarray(means);cpu_residual=float(np.linalg.norm(rhs-system@expected)/np.linalg.norm(rhs))
        if cpu_residual>1e-10:raise RuntimeError('Independent CPU true residual failed')
        result['cpu_true_residual']=cpu_residual
        result['oracle_sparse_storage']={name:dict(shape=list(value.shape),nnz=int(value.nnz),canonical=bool(value.has_canonical_format),sorted_indices=bool(value.has_sorted_indices)) for name,value in (('B',B),('mass',mass))}
        output_B,output_mass=sorted_oracle_csr(B),sorted_oracle_csr(mass)
        fixture=a.out/'cpu.h5';write_hdiv_h5(fixture,dict(B=output_B,mass=output_mass,rhs=rhs,system=system,expected=expected,mean_operator=means),dict(units='SI',labels=list(mesh.GetBoundaries())))
        result.update(cells=mesh.ne,dofs=fes.ndof,charges=B.shape[0],cpu_true_residual=cpu_residual,cpu_iterations=int(reference['iters']),cpu_vim_sha256=sha(cpu.__file__),input_sha256=sha(fixture),mesh_sha256=sha(mesh_path),source_sha256=sha(__file__),projection_sha256=sha(a.projection/'radia_cuda_projection.mexw64'),hmatrix_sha256=sha(a.hmatrix/'radia_cuda_hmatrix.mexw64'))
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');result['matlab']=eng.version();result['ngsolve']=ng.__version__
        result['owned_matlab_pid']=int(eng.feature('getpid'))
        eng.eval("observedDevice=gpuDevice;observedGpu=struct('name',observedDevice.Name,'compute_capability',observedDevice.ComputeCapability,'selected_device',observedDevice.Index,'backend','MATLAB mexcuda');",nargout=0)
        result['gpu']=json.loads(eng.eval('jsonencode(observedGpu)'))
        result['a100_used']='A100' in result['gpu']['name']
        matlab_source=a.out/'matlab';shutil.copytree(root/'packages/radia-simulink/matlab',matlab_source)
        result['matlab_source_sha256']={f.relative_to(matlab_source).as_posix():sha(f) for f in matlab_source.rglob('*.m')}
        eng.addpath(str(a.projection),str(a.hmatrix),str(matlab_source),nargout=0)
        eng.workspace['mesh_path']=str(mesh_path);eng.workspace['degree']=float(a.order);eng.workspace['fixture_path']=str(fixture);eng.workspace['output_path']=str(a.out/'native.h5')
        eng.workspace['assembly_mode']=a.assembly_mode
        eng.eval("[cpu,fixture_metadata]=radia.gpu.readHDivH5(fixture_path);cpu.rhs=cpu.rhs(:);cpu.expected=cpu.expected(:);assert(isequal(size(cpu.rhs),size(cpu.expected)) && isequal(size(cpu.system),[numel(cpu.rhs) numel(cpu.rhs)]),'radia:HDivShape','CPU vector/system shape mismatch');setup=radia_cuda_projection('assemble_hdiv_elements',mesh_path,degree,[1 .2 .1]);assert(isequal(size(setup.rhs),size(cpu.rhs)),'radia:HDivShape','Native RHS must match CPU column shape');mass_error=norm(setup.mass-full(cpu.mass),'fro')/norm(full(cpu.mass),'fro');charge_error=norm(setup.charge_map-full(cpu.B),'fro')/norm(full(cpu.B),'fro');rhs_error=norm(setup.rhs-cpu.rhs)/norm(cpu.rhs);assert(max([mass_error charge_error rhs_error])<1e-11);solution=radia.gpu.solveHDivElements(mesh_path,degree,[1 .2 .1],2.5,1e-11,1000,assembly_mode);solution.coefficients=solution.coefficients(:);assert(isequal(size(solution.coefficients),size(cpu.expected)),'radia:HDivShape','Coefficient column shape mismatch');coefficient_error=norm(solution.coefficients-cpu.expected)/norm(cpu.expected);mean_reference=cpu.mean_operator*cpu.expected;assert(isequal(size(solution.mean_magnetization),size(mean_reference)),'radia:HDivShape','Mean magnetization/reference shape mismatch');physical_error=norm(solution.mean_magnetization-mean_reference)/norm(mean_reference);independent_true_residual=norm(cpu.rhs-cpu.system*solution.coefficients)/norm(cpu.rhs);assert(radia_cuda_hmatrix('live_count')==0);output=struct('mass_error',mass_error,'charge_error',charge_error,'rhs_error',rhs_error,'coefficient_error',coefficient_error,'physical_error',physical_error,'independent_true_residual',independent_true_residual);output_metadata=fixture_metadata;output_metadata.solution_text=struct;for fn=fieldnames(solution)',value=solution.(fn{1});if isnumeric(value)||islogical(value),output.(['solution_' fn{1}])=value;elseif ischar(value)||isstring(value),output_metadata.solution_text.(fn{1})=value;else,error('radia:HDF5','Unsupported solution field');end;end;radia.gpu.writeHDivH5(output_path,output,output_metadata);assert(max([coefficient_error physical_error])<1e-8);assert(independent_true_residual<1e-10);",nargout=0)
        native,output_metadata=read_hdiv_h5(a.out/'native.h5'); native['solution']={key.removeprefix('solution_'):(value.reshape(-1)[0].item() if value.size==1 else value.tolist()) for key,value in native.items() if key.startswith('solution_')}; native['solution'].update(output_metadata.get('solution_text',{}))
        result['native']={key:float(native[key].reshape(-1)[0]) for key in ('mass_error','charge_error','rhs_error','coefficient_error','physical_error','independent_true_residual')};result['solution']=native['solution'];result['accepted']=True
    except BaseException as exc:result['error']=f'{type(exc).__name__}: {exc}';raise
    finally:
        if (a.out/'native.h5').exists() and 'native' not in result:
            try:
                from hdiv_hdf5 import read_hdiv_h5
                native,output_metadata=read_hdiv_h5(a.out/'native.h5'); native['solution']={key.removeprefix('solution_'):(value.reshape(-1)[0].item() if value.size==1 else value.tolist()) for key,value in native.items() if key.startswith('solution_')}; native['solution'].update(output_metadata.get('solution_text',{}))
                result['native']={key:float(native[key].reshape(-1)[0]) for key in ('mass_error','charge_error','rhs_error','coefficient_error','physical_error','independent_true_residual')}
                result['solution']=native['solution']
            except BaseException as diagnostic_error:
                result['native_diagnostic_read_error']=str(diagnostic_error)
        if eng is not None:
            try:eng.quit();result['owned_engine_quit']=True
            except BaseException as exc:result['release_error']=str(exc);result['accepted']=False
        (a.out/'result.json').write_text(json.dumps(result,indent=2,default=lambda value:value.tolist())+'\n',encoding='utf-8')

if __name__=='__main__':main()
