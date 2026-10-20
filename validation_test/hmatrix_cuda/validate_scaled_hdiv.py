"""Scale parity without constructing a dense DoF-by-DoF CPU oracle.

Python owns only independent CPU reference/certification and Engine orchestration.
GPU setup receives the checked mesh and physical parameters, not Python operators.
"""
from pathlib import Path
import argparse,hashlib,json,os,time,sys,statistics
from validate_hdiv_artifact import idle
from hdiv_validation_inputs import checked_mesh,dependency_directories,cpu_receipt,single_thread_blas,sorted_oracle_csr,process_affinity,thread_record

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('mesh','contract','projection','hmatrix','matlab','out','cpu-build-receipt'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--mesh-sha256',required=True);p.add_argument('--contract-sha256',required=True)
    p.add_argument('--dependency-dir',type=Path,action='append',required=True);p.add_argument('--cpu-threads',type=int,default=16)
    p.add_argument('--preconditioner',choices=['mass-diagonal','mass-cholesky','mass-polynomial'],default='mass-cholesky');p.add_argument('--staging-mode',choices=['dense','sparse-mass','sparse'],default='dense');p.add_argument('--abba-cycles',type=int,default=0);p.add_argument('--assembly-mode',choices=['host-cache','host-cache-moments','hacapk','hacapk-cached','hacapk-device-dense'],default='host-cache')
    p.add_argument('--gpu-pair-cycles',type=int,choices=(0,2,3,4),default=0,help='Same-module HACApK/dense GPU ABBA, independent CPU certification outside timing')
    p.add_argument('--gpu-pair-control',choices=['host-cache-moments','hacapk-cached'],default='host-cache-moments')
    p.add_argument('--hacapk-epsilon',type=float);p.add_argument('--hacapk-leaf-size',type=int);p.add_argument('--hacapk-eta',type=float);p.add_argument('--hacapk-max-rank',type=int);p.add_argument('--hacapk-cache-mib',type=int);p.add_argument('--operator-probes',type=int,default=0,choices=range(0,9),help='Seeded independent CPU B-transpose G B action checks; excluded from timings');a=p.parse_args()
    if a.gpu_pair_cycles and (a.assembly_mode not in ('hacapk-cached','hacapk-device-dense') or a.assembly_mode==a.gpu_pair_control):p.error('GPU pair requires distinct explicit cached HACApK primary/control modes')
    options=dict(epsilon=a.hacapk_epsilon,leaf_size=a.hacapk_leaf_size,eta=a.hacapk_eta,max_rank=a.hacapk_max_rank)
    if a.assembly_mode.startswith('hacapk'):
        if any(v is None for v in options.values()):p.error('HACApK requires explicit epsilon, leaf size, eta and max rank')
        import math
        if not (0<a.hacapk_epsilon<1 and math.isfinite(a.hacapk_eta) and a.hacapk_eta>0 and 1<=a.hacapk_leaf_size<=4096 and 1<=a.hacapk_max_rank<=4096):p.error('Invalid bounded HACApK options')
    elif any(v is not None for v in options.values()):p.error('HACApK options require the hacapk assembly mode')
    else:options={}
    if a.hacapk_cache_mib is not None:
        if a.assembly_mode not in ('hacapk-cached','hacapk-device-dense') or not 1<=a.hacapk_cache_mib<=2048:p.error('Cache MiB requires cached HACApK and 1..2048')
        options['cache_bytes']=a.hacapk_cache_mib*2**20
    if a.abba_cycles not in (0,2,3,4):p.error('ABBA cycles must be zero or 2..4');a.out.mkdir(parents=True,exist_ok=False)
    result=dict(accepted=False,release_accepted=False,formal_reference_accepted=False,scope='Bounded public P2 sphere linear HDiv scaled whole-solve parity',source_sha256=sha(__file__),preconditioner=a.preconditioner,staging_mode=a.staging_mode,assembly_mode=a.assembly_mode,hacapk_options=options)
    from hdiv_process_memory import ProcessMemory
    e=None;cpu_memory=None;gpu_memory=None
    try:
        result['idle_before']=idle()
        dependency_dirs,dll_handles=dependency_directories(a.dependency_dir);checked_mesh(a.mesh,a.contract,a.mesh_sha256,a.contract_sha256,a.out/'mesh-check.json')
        import numpy as np,ngsolve as ng,matlab.engine
        import radia.vim._vim as cpu
        import radia._radia_pybind as native
        from hdiv_hdf5 import write_hdiv_h5,read_hdiv_h5
        result['threads']=thread_record(a.cpu_threads)
        result['memory_limitations']='50ms sampled process RSS/private includes baseline/runtime and may miss spikes; CPU sampler is in controller process, GPU sampler observes separate MATLAB process, unquantified timing bias; exact GPU owned excludes CUDA/runtime/global peaks'
        result['reference']=cpu_receipt(a.cpu_build_receipt,native.__file__);ng.SetNumThreads(a.cpu_threads);limit,result['blas']=single_thread_blas()
        result['cpu_threads']=a.cpu_threads;result['native_sha256']=sha(native.__file__);result['frontend_sha256']=sha(cpu.__file__)
        cpu_memory=ProcessMemory(os.getpid());cpu_memory.__enter__()
        start=time.perf_counter();mesh=ng.Mesh(str(a.mesh))
        with ng.TaskManager():
            fes=ng.HDiv(mesh,order=2);nd=fes.ndof
            B,G,mass=cpu.build_charge_gram(fes,intorder=4,curve_order=2,curve_gauss=8,ho_far_factor=1e30,gram_backend='exact-dense',exact_dense_memory_mb=1024)
            load=ng.LinearForm(fes);load+=ng.CoefficientFunction((1.,.2,.1))*fes.TestFunction()*ng.dx;load.Assemble();rhs=np.array(load.vec.FV().NumPy())
            answer=G.solve_configured_linear_material_mass_riesz(2.5,rhs,1e-11,1000);expected=np.array(answer['m']);means=[]
            for c in range(3):
                form=ng.LinearForm(fes);form+=fes.TestFunction()[c]*ng.dx;form.Assemble();means.append(np.array(form.vec.FV().NumPy())/ng.Integrate(1,mesh))
            means=np.asarray(means)
            def residual(x):return float(np.linalg.norm(rhs-(B.T@np.asarray(G.matvec_sym(np.asarray(B@x).reshape(-1)))+2.5*(mass@x)))/np.linalg.norm(rhs))
            result['cpu_true_residual']=residual(expected)
        result['cpu_setup_assembly_solve_recovery_seconds']=time.perf_counter()-start
        cpu_memory.__exit__(None,None,None);result['cpu_host_process_memory']=cpu_memory.result();cpu_memory=None
        result['cpu_time_limitation']='Retains independent CPU operator for certification; destruction excluded. No CPU/GPU speedup ratio asserted.'
        result.update(dofs=nd,charges=B.shape[0],cells=mesh.ne,cpu_iterations=int(answer['iters']))
        assert result['cpu_true_residual']<1e-10
        oracle=a.out/'cpu.h5';write_hdiv_h5(oracle,dict(B=sorted_oracle_csr(B),mass=sorted_oracle_csr(mass),rhs=rhs,expected=expected,mean_operator=means),dict(units='SI',labels=list(mesh.GetBoundaries())))
        if a.operator_probes:
            if a.assembly_mode not in ('hacapk-cached','hacapk-device-dense'):raise ValueError('Operator probes require explicit cached HACApK')
            rng=np.random.default_rng(20261018);probe_x=rng.standard_normal((nd,a.operator_probes))
            with ng.TaskManager():probe_y=np.column_stack([B.T@np.asarray(G.matvec_sym(np.asarray(B@probe_x[:,k]).reshape(-1))) for k in range(a.operator_probes)])
            probe_path=a.out/'operator-probes.h5';write_hdiv_h5(probe_path,dict(x=probe_x,expected_action=probe_y),dict(units='SI',labels=list(mesh.GetBoundaries()),seed=20261018))
            result['operator_probe_input_sha256']=sha(probe_path)
        result['oracle_sha256']=sha(oracle);result['mesh_sha256']=sha(a.mesh);result['contract_sha256']=sha(a.contract)
        result['idle_before_gpu']=idle()
        t=time.perf_counter();e=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');result['engine_start_seconds']=time.perf_counter()-t
        result['pid']=int(e.feature('getpid'));result['matlab']=e.version();result['matlab_threads']=int(e.maxNumCompThreads());result['matlab_affinity']=process_affinity(result['pid'])
        e.addpath(str(a.projection),str(a.hmatrix),str(a.matlab),str(Path(__file__).resolve().parent),nargout=0)
        result['mex']={}
        for name,directory in (('radia_cuda_projection',a.projection),('radia_cuda_hmatrix',a.hmatrix)):
            loaded=Path(e.which(name)).resolve();expected_path=(directory/(name+'.mexw64')).resolve()
            assert loaded==expected_path
            receipt=json.loads((directory/'build.json').read_text());assert receipt['accepted'] and receipt['binary_sha256']==sha(loaded)
            result['mex'][name]=dict(sha256=sha(loaded),source_commit=receipt['source_commit'])
        e.workspace['mesh_path']=str(a.mesh);e.workspace['oracle_path']=str(oracle);e.workspace['output_path']=str(a.out/'gpu.h5');e.workspace['preconditioner']=a.preconditioner;e.workspace['staging_mode']=a.staging_mode;e.workspace['assembly_mode']=a.assembly_mode
        e.workspace['hacapk_options_json']=json.dumps(options);e.eval('hacapk_options=jsondecode(hacapk_options_json);',nargout=0)
        e.workspace['setup_command']='assemble_hdiv_elements_sparse' if a.staging_mode=='sparse' else ('assemble_hdiv_elements_sparse_mass' if a.staging_mode=='sparse-mass' else 'assemble_hdiv_elements')
        e.eval("[cpu,metadata]=radia.gpu.readHDivH5(oracle_path);setup=radia_cuda_projection(setup_command,mesh_path,2,[1 .2 .1]);mass_error=norm(setup.mass-cpu.mass,'fro')/norm(cpu.mass,'fro');charge_error=norm(setup.charge_map-cpu.B,'fro')/norm(cpu.B,'fro');rhs_error=norm(setup.rhs-cpu.rhs(:))/norm(cpu.rhs(:));assert(max([mass_error charge_error rhs_error])<1e-11);clear setup cpu;",nargout=0)
        gpu_start=time.perf_counter();e.eval('observedDevice=gpuDevice;',nargout=0);result['gpu_initialization_seconds']=time.perf_counter()-gpu_start;result['cold_limitations']='First solve follows Engine/device initialization and native FE parity staging; interpreter and process cold startup are not measured';result['available_gpu_memory_before']=float(e.eval('observedDevice.AvailableMemory'))
        gpu_memory=ProcessMemory(result['pid']);gpu_memory.__enter__()
        for phase in ('first_payload','warm'):
            output_path=a.out/('gpu-'+phase+'.h5');e.workspace['output_path']=str(output_path)
            e.eval("s=radia.gpu.solveHDivElements(mesh_path,2,[1 .2 .1],2.5,1e-11,1000,assembly_mode,preconditioner,staging_mode,false,hacapk_options);assert(radia_cuda_hmatrix('live_count')==0);o=struct('coefficients',s.coefficients,'mean',s.mean_magnetization,'whole_seconds',s.whole_seconds,'assembly_upload_seconds',s.assembly_upload_seconds,'solve_transfer_seconds',s.solve_transfer_seconds,'iterations',s.iterations,'owned_peak',s.combined_peak_owned_device_bytes);radia.gpu.writeHDivH5(output_path,o,metadata);",nargout=0)
            gpu,_=read_hdiv_h5(output_path);x=np.asarray(gpu['coefficients']).reshape(-1);mean=np.asarray(gpu['mean']).reshape(-1);refmean=means@expected
            with ng.TaskManager():true=residual(x)
            errors=dict(coefficient_error=float(np.linalg.norm(x-expected)/np.linalg.norm(expected)),mean_error=float(np.linalg.norm(mean-refmean)/np.linalg.norm(refmean)),true_residual=true)
            sample={k:float(np.asarray(gpu[k]).reshape(-1)[0]) for k in ('whole_seconds','assembly_upload_seconds','solve_transfer_seconds','iterations','owned_peak')};sample.update(errors)
            result[phase]=sample
            result[phase]['operator_stats']=json.loads(e.eval('jsonencode(s.assembly_profile)'))
            assert max(errors['coefficient_error'],errors['mean_error'])<1e-8 and true<1e-10
        gpu_memory.__exit__(None,None,None);result['gpu_host_process_memory']=gpu_memory.result();gpu_memory=None
        if a.operator_probes:
            e.workspace['probe_path']=str(probe_path)
            result['operator_probes']=json.loads(e.eval('jsonencode(checkHDivCPUActions(mesh_path,probe_path,hacapk_options,assembly_mode))'))
            assert result['operator_probes']['accepted'] and result['operator_probes']['count']==a.operator_probes
        if a.gpu_pair_cycles:
            result['gpu_pair_samples']=[]
            primary,control=a.assembly_mode,a.gpu_pair_control
            pairs=[('first_payload',mode) for mode in (primary,control)]+[('warm',mode) for _ in range(a.gpu_pair_cycles) for mode in (primary,control,control,primary)]
            for index,(phase,mode) in enumerate(pairs):
                before=idle();e.workspace['pair_mode']=mode;e.workspace['pair_options_json']=json.dumps(options if mode.startswith('hacapk') else {})
                e.eval('pair_options=jsondecode(pair_options_json);',nargout=0)
                output_path=a.out/('gpu-pair-'+str(index)+'.h5');e.workspace['output_path']=str(output_path)
                with ProcessMemory(result['pid']) as memory:
                    controller=time.perf_counter()
                    e.eval("s=radia.gpu.solveHDivElements(mesh_path,2,[1 .2 .1],2.5,1e-11,1000,pair_mode,preconditioner,staging_mode,false,pair_options);assert(radia_cuda_hmatrix('live_count')==0);o=struct('coefficients',s.coefficients,'mean',s.mean_magnetization,'whole_seconds',s.whole_seconds,'assembly_upload_seconds',s.assembly_upload_seconds,'solve_transfer_seconds',s.solve_transfer_seconds,'iterations',s.iterations,'owned_peak',s.combined_peak_owned_device_bytes);radia.gpu.writeHDivH5(output_path,o,metadata);",nargout=0)
                    output,_=read_hdiv_h5(output_path);controller=time.perf_counter()-controller
                x2=np.asarray(output['coefficients']).reshape(-1);mean2=np.asarray(output['mean']).reshape(-1)
                with ng.TaskManager():true=residual(x2)
                errors=dict(coefficient_error=float(np.linalg.norm(x2-expected)/np.linalg.norm(expected)),mean_error=float(np.linalg.norm(mean2-means@expected)/np.linalg.norm(means@expected)),true_residual=true)
                sample=dict(phase=phase,assembly_mode=mode,outer_controller_seconds=controller,idle_before=before,process_memory=memory.result(),operator_stats=json.loads(e.eval('jsonencode(s.assembly_profile)')))
                sample.update({k:float(np.asarray(output[k]).reshape(-1)[0]) for k in ('whole_seconds','assembly_upload_seconds','solve_transfer_seconds','iterations','owned_peak')});sample.update(errors);result['gpu_pair_samples'].append(sample)
                assert max(errors['coefficient_error'],errors['mean_error'])<1e-8 and true<1e-10
            result['gpu_pair_warm_api_medians']={mode:statistics.median(row['whole_seconds'] for row in result['gpu_pair_samples'] if row['phase']=='warm' and row['assembly_mode']==mode) for mode in (primary,control)}
            result['gpu_pair_control_over_primary']=result['gpu_pair_warm_api_medians'][control]/result['gpu_pair_warm_api_medians'][primary]
            if primary=='hacapk-cached' and control=='host-cache-moments':result['gpu_pair_dense_over_hacapk']=result['gpu_pair_control_over_primary']
            result['gpu_pair_scope']='Same module/input/material/preconditioner. Each sample includes FE setup, assembly, upload, solve, recovery and destruction; excludes validation, Engine, HDF5 and controller transport. First dense call is not process cold; different assembly algorithms and cache policies remain explicit.'
        e.eval('observedDevice=gpuDevice;',nargout=0);result['available_gpu_memory_after']=float(e.eval('observedDevice.AvailableMemory'))
        if a.abba_cycles:
            result['benchmark_samples']=[]
            result['benchmark_first_limitations']='The first benchmark iteration follows independent CPU reference and GPU validation warm-ups; it is not process cold start. Earlier first GPU solve is recorded separately.'
            result['ratio_scope']='Fresh CPU and GPU complete linear workflows, same public mesh/physics, CPU mass-Riesz sparsecholesky versus GPU explicit chosen preconditioner; excludes interpreter/Engine startup, validation, HDF5 serialization and result transport from API totals'
            result['benchmark_memory_scope']='Sampled process baselines include the retained independent CPU certification operator; not standalone application memory peaks'
            def cpu_workflow():
                start=time.perf_counter();mesh2=ng.Mesh(str(a.mesh))
                with ng.TaskManager():
                    fes2=ng.HDiv(mesh2,order=2)
                    B2,G2,M2=cpu.build_charge_gram(fes2,intorder=4,curve_order=2,curve_gauss=8,ho_far_factor=1e30,gram_backend='exact-dense',exact_dense_memory_mb=1024)
                    load2=ng.LinearForm(fes2);load2+=ng.CoefficientFunction((1.,.2,.1))*fes2.TestFunction()*ng.dx;load2.Assemble()
                    answer2=G2.solve_configured_linear_material_mass_riesz(2.5,np.array(load2.vec.FV().NumPy()),1e-11,1000)
                    x2=np.asarray(answer2['m']).copy();iterations=int(answer2['iters']);mean2=[]
                    for c in range(3):
                        form2=ng.LinearForm(fes2);form2+=fes2.TestFunction()[c]*ng.dx;form2.Assemble()
                        mean2.append(float(np.dot(form2.vec.FV().NumPy(),x2)/ng.Integrate(1,mesh2)))
                    del form2,load2,answer2,B2,G2,M2,fes2,mesh2
                return x2,np.asarray(mean2),iterations,time.perf_counter()-start
            order=[('first_payload','cpu'),('first_payload','gpu')]+[('warm',backend) for _ in range(a.abba_cycles) for backend in ('cpu','gpu','gpu','cpu')]
            for index,(phase,backend) in enumerate(order):
                before=idle()
                with ProcessMemory(os.getpid() if backend=='cpu' else result['pid']) as memory:
                    controller=time.perf_counter()
                    if backend=='cpu':x2,mean2,iterations,elapsed=cpu_workflow();owned=None
                    else:
                        output_path=a.out/('benchmark-'+str(index)+'.h5');e.workspace['output_path']=str(output_path)
                        e.eval("s=radia.gpu.solveHDivElements(mesh_path,2,[1 .2 .1],2.5,1e-11,1000,assembly_mode,preconditioner,staging_mode,false,hacapk_options);assert(radia_cuda_hmatrix('live_count')==0);o=struct('coefficients',s.coefficients,'mean',s.mean_magnetization,'whole_seconds',s.whole_seconds,'iterations',s.iterations,'owned_peak',s.combined_peak_owned_device_bytes);radia.gpu.writeHDivH5(output_path,o,metadata);",nargout=0)
                        output,_=read_hdiv_h5(output_path);x2=np.asarray(output['coefficients']).reshape(-1);mean2=np.asarray(output['mean']).reshape(-1)
                        elapsed=float(np.asarray(output['whole_seconds']).reshape(-1)[0]);iterations=int(np.asarray(output['iterations']).reshape(-1)[0]);owned=int(np.asarray(output['owned_peak']).reshape(-1)[0])
                    controller=time.perf_counter()-controller
                with ng.TaskManager():true=residual(x2)
                coefficient_error=float(np.linalg.norm(x2-expected)/np.linalg.norm(expected));mean_error=float(np.linalg.norm(mean2-means@expected)/np.linalg.norm(means@expected))
                sample=dict(phase=phase,backend=backend,user_api_total_seconds=elapsed,outer_controller_seconds=controller,iterations=iterations,coefficient_error=coefficient_error,mean_error=mean_error,true_residual=true,owned_peak_bytes=owned,process_memory=memory.result(),idle_before=before)
                result['benchmark_samples'].append(sample)
                assert coefficient_error<1e-8 and mean_error<1e-8 and true<1e-10
            medians={backend:statistics.median(row['user_api_total_seconds'] for row in result['benchmark_samples'] if row['phase']=='warm' and row['backend']==backend) for backend in ('cpu','gpu')}
            result['benchmark_warm_api_medians']=medians;result['benchmark_cpu_over_gpu']=medians['cpu']/medians['gpu']
            first={backend:next(row['user_api_total_seconds'] for row in result['benchmark_samples'] if row['phase']=='first_payload' and row['backend']==backend) for backend in ('cpu','gpu')}
            result['benchmark_first_api_seconds']=first;result['benchmark_first_cpu_over_gpu']=first['cpu']/first['gpu']
            result['benchmark_warm_outer_medians']={backend:statistics.median(row['outer_controller_seconds'] for row in result['benchmark_samples'] if row['phase']=='warm' and row['backend']==backend) for backend in ('cpu','gpu')}
        result['idle_after']=idle()
        result['native_setup_errors']={k:float(e.eval(k)) for k in ('mass_error','charge_error','rhs_error')};result['accepted']=True
    except BaseException as ex:result['error']=repr(ex);raise
    finally:
        for sampler in (cpu_memory,gpu_memory):
            if sampler is not None:sampler.__exit__(None,None,None)
        if e is not None:e.quit();result['owned_engine_quit']=True
        (a.out/'result.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
