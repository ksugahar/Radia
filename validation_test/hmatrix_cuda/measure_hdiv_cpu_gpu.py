"""Idle same-host independent CPU versus transfer-inclusive MATLAB HDiv timing.

Consumes a hash-verified job payload and an independently built CPU reference.
Python performs CPU validation and orchestration, never GPU solver work.
"""
from pathlib import Path
import argparse, ctypes, hashlib, json, os, statistics, sys, time, math
from validate_hdiv_artifact import idle
from hdiv_process_memory import ProcessMemory
from hdiv_validation_inputs import checked_mesh, dependency_directories, thread_record, oracle_vectors, cpu_receipt, single_thread_blas, observed_blas, process_affinity


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', type=Path, required=True)
    parser.add_argument('--cpu-root', type=Path, required=True)
    parser.add_argument('--cpu-manifest', type=Path, required=True,
                        help='JSON mapping paths relative to cpu-root to SHA256')
    parser.add_argument('--cycles', type=int, default=2)
    parser.add_argument('--wheel', type=Path,
                        help='Require verified installed-wheel identity and use only its MATLAB payload')
    parser.add_argument('--cpu-threads', type=int, default=1, help='NGSolve CPU reference threads; BLAS remains single-threaded')
    parser.add_argument('--dependency-dir',type=Path,action='append',required=True)
    parser.add_argument('--cpu-build-receipt',type=Path,required=True)
    args = parser.parse_args()
    if args.cpu_threads < 1:
        parser.error('CPU threads must be positive')
    if args.cycles < 2:
        parser.error('At least two ABBA cycles are required')
    job = args.job.resolve()
    result = dict(accepted=False, release_accepted=False, formal_reference_accepted=False,
                  independent_review_pending=True, samples=[],
                  classification='Bounded same-host diagnostic measurement',
                  memory_limitations='Exact GPU owned buffers, available-memory endpoints and 50ms host-process RSS/private samples; sampled peaks can miss spikes, lifetime RSS includes earlier work; global/CUDA runtime allocation peaks unmeasured',
                  cold_limitations='Fresh owned Engine startup and first payload separately; CPU interpreter/import startup excluded')
    result['timing_boundary']=dict(
        user_api_total='CPU Python workflow or MATLAB solveHDivElements: FE, assembly, transfers where applicable, solve, physical result recovery and native destruction; Engine transport excluded',
        controller_elapsed='Outer Python clock including Engine call, result and diagnostic-stage RPC collection for GPU; direct workflow call for CPU; independent certification and idle sampling excluded',
        startup='Owned Engine startup and GPU initialization recorded separately; interpreter/import startup excluded',
        whole_seconds_compatibility='Alias of user_api_total_seconds, not controller elapsed')
    result['memory_sampling_timing_bias']='Sampler runs in the Python controller process: CPU workflow competes with its sampler thread, while MATLAB runs in a separate process. Ratios include this unquantified asymmetric overhead; memory was not measured in a separate run.'
    result['stage_comparison']='Backend-specific boundaries; CPU stages and GPU stages are not directly comparable'
    result['transfer_timing']='Included in GPU API total and combined stage times, not measured separately'
    engine = None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
        os.environ.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_THREADING_LAYER='SEQUENTIAL')
        dependency_dirs,dll_handles=dependency_directories(args.dependency_dir)
        result['dependency_directories']=list(map(str,dependency_dirs))
        result['threads']=thread_record(args.cpu_threads)
        manifest = json.loads((job/'manifest.json').read_text())
        cpu_manifest = json.loads(args.cpu_manifest.read_text())
        cpu_root = args.cpu_root.resolve()
        def verify():
            for name, digest in manifest.items():
                path=(job/name).resolve()
                if not path.is_relative_to(job) or sha(path) != digest:
                    raise RuntimeError('Payload hash mismatch: '+name)
            for name, digest in cpu_manifest.items():
                path=(cpu_root/name).resolve()
                if not path.is_relative_to(cpu_root) or sha(path) != digest:
                    raise RuntimeError('CPU reference hash mismatch: '+name)
        verify()
        sys.path.insert(0, str(args.cpu_root.resolve()))
        import numpy as np
        import ngsolve as ng
        import radia.vim._vim as cpu
        import radia._radia_pybind as native
        from hdiv_hdf5 import read_hdiv_h5
        import matlab.engine
        ng.SetNumThreads(args.cpu_threads)
        blas_limit,result['observed_blas_after_ngsolve_configuration']=single_thread_blas()
        result['reuse_policy']='Both routes rebuild FE, Gram, material operator and solver for every sample; no cross-sample operator/factor reuse'
        for module in (cpu,native):
            relative=Path(module.__file__).resolve().relative_to(cpu_root).as_posix()
            if relative not in cpu_manifest:
                raise RuntimeError('Loaded CPU module missing from reference manifest: '+relative)
        result['reference'] = dict(native_sha256=sha(native.__file__), frontend_sha256=sha(cpu.__file__),
                                   ngsolve=ng.__version__, python=sys.version,
                                   limitation='Loaded file identity; independent native build receipt must accompany release evidence')
        result['reference']['build_receipt']=cpu_receipt(args.cpu_build_receipt,native.__file__)
        result['reference']['build_receipt_sha256']=sha(args.cpu_build_receipt)
        result['payload_sha256'] = manifest
        result['cpu_reference_sha256'] = cpu_manifest
        matlab_path = job/'matlab'
        if args.wheel:
            import importlib.util
            import radia_simulink
            spec=importlib.util.spec_from_file_location('wheel_verifier',job/'verify_wheel.py')
            verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
            wheel=verifier.verify(args.wheel)
            verifier.verify_installed_identity(radia_simulink.verify(),wheel)
            result['wheel']=wheel
            matlab_path=radia_simulink.matlab_path()
        result['idle_before'] = idle()
        start = time.perf_counter()
        engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        result['engine_start_seconds'] = time.perf_counter()-start
        result['matlab'] = engine.version()
        result['owned_matlab_pid']=int(engine.feature('getpid'))
        result['matlab_max_num_comp_threads']=int(engine.maxNumCompThreads())
        result['matlab_process_affinity_mask']=process_affinity(result['owned_matlab_pid'])
        device_start=time.perf_counter()
        engine.eval("observedDevice=gpuDevice;observedGpu=struct('name',observedDevice.Name,'compute_capability',observedDevice.ComputeCapability,'selected_device',observedDevice.Index,'backend','MATLAB mexcuda');",nargout=0)
        result['gpu_initialization_seconds']=time.perf_counter()-device_start
        result['gpu']=json.loads(engine.eval('jsonencode(observedGpu)'))
        engine.restoredefaultpath(nargout=0)
        engine.addpath(str(job), str(matlab_path), nargout=0)
        if args.wheel:
            installed_manifest=json.loads((matlab_path.parent/'manifest.json').read_text())
            for name in ('radia_cuda_hmatrix','radia_cuda_projection'):
                loaded=Path(engine.which(name)).resolve()
                if loaded.parent!=matlab_path.resolve() or sha(loaded)!=installed_manifest['files']['matlab/'+loaded.name]:
                    raise RuntimeError('HDiv must load from the verified installed wheel')
        else:
            result['loaded_mex']={}
            for name in ('radia_cuda_hmatrix','radia_cuda_projection'):
                loaded=Path(engine.which(name)).resolve(strict=True)
                relative=loaded.relative_to(job).as_posix()
                if relative not in manifest or sha(loaded)!=manifest[relative]:
                    raise RuntimeError('HDiv MEX missing from verified job manifest')
                result['loaded_mex'][name]=dict(relative_path=relative,sha256=sha(loaded))
        for case in json.loads((job/'cases.json').read_text()):
            for name in ('mesh','contract','fixture'):
                if case[name] not in manifest:raise RuntimeError('Case input missing from payload manifest: '+name)
            checked_mesh(job/case['mesh'],job/case['contract'],manifest[case['mesh']],manifest[case['contract']],job/(case['name']+'-vol-check.json'))
            fixture,fixture_metadata = read_hdiv_h5(job/case['fixture'])
            expected,system,rhs=oracle_vectors(fixture)
            engine.workspace['mesh_path'] = str(job/case['mesh'])
            engine.workspace['degree'] = float(case['order'])
            engine.workspace['assembly_mode'] = case.get('assembly_mode','host-cache')
            # Independent HDF5 oracle is never passed to the GPU solver.
            def run(backend, phase):
                before = idle()
                blas_before=observed_blas()
                with ProcessMemory(os.getpid() if backend=='cpu' else result['owned_matlab_pid']) as memory:
                    start = time.perf_counter()
                    if backend == 'cpu':
                        stage_start=time.perf_counter()
                        mesh = ng.Mesh(str(job/case['mesh']))
                        with ng.TaskManager():
                            fes = ng.HDiv(mesh, order=case['order'])
                            setup_seconds=time.perf_counter()-stage_start;stage_start=time.perf_counter()
                            B,G,mass = cpu.build_charge_gram(fes,intorder=4,
                                curve_order=2 if case.get('geometry','curved')=='curved' else 0,
                                curve_gauss=8,ho_far_factor=1e30,gram_backend='exact-dense',exact_dense_memory_mb=64)
                            load=ng.LinearForm(fes)
                            load += ng.CoefficientFunction((1.,.2,.1))*fes.TestFunction()*ng.dx
                            load.Assemble()
                            assembly_seconds=time.perf_counter()-stage_start;stage_start=time.perf_counter()
                            answer=G.solve_configured_linear_material_mass_riesz(2.5,np.array(load.vec.FV().NumPy()),1e-11,1000)
                            solve_seconds=time.perf_counter()-stage_start;stage_start=time.perf_counter()
                            iterations=int(answer['iters'])
                            coefficients=np.asarray(answer['m']).copy()
                            mean=[]
                            for c in range(3):
                                form=ng.LinearForm(fes);form+=fes.TestFunction()[c]*ng.dx;form.Assemble()
                                mean.append(float(np.dot(form.vec.FV().NumPy(),coefficients)/ng.Integrate(1,mesh)))
                            recovery_seconds=time.perf_counter()-stage_start
                        stage_start=time.perf_counter()
                        del answer,B,G,mass,load,form,fes,mesh
                        stages=dict(fe_setup_seconds=setup_seconds,assembly_seconds=assembly_seconds,
                                    solve_seconds=solve_seconds,recovery_seconds=recovery_seconds,
                                    destroy_seconds=time.perf_counter()-stage_start)
                        elapsed=time.perf_counter()-start
                        owned=None
                    else:
                        engine.eval("timer=tic;solution=radia.gpu.solveHDivElements(mesh_path,degree,[1 .2 .1],2.5,1e-11,1000,assembly_mode);whole_seconds=toc(timer);assert(radia_cuda_hmatrix('live_count')==0);",nargout=0)
                        coefficients=np.asarray(engine.eval('solution.coefficients')).ravel()
                        mean=np.asarray(engine.eval('solution.mean_magnetization')).ravel()
                        owned=float(engine.eval('solution.combined_peak_owned_device_bytes'))
                        elapsed=float(engine.workspace['whole_seconds'])
                        stages={name:float(engine.eval('solution.'+name)) for name in
                                ('fe_setup_seconds','assembly_upload_seconds','solve_transfer_seconds','recovery_seconds','destroy_seconds')}
                        iterations=int(engine.eval('solution.iterations'))
                    controller_elapsed=time.perf_counter()-start
                expected_mean=fixture['mean_operator']@expected
                if coefficients.shape!=expected.shape or np.asarray(mean).shape!=expected_mean.shape:
                    raise RuntimeError('CPU/GPU coefficient or mean vector shape mismatch')
                error=float(np.linalg.norm(coefficients-expected)/np.linalg.norm(expected))
                residual=float(np.linalg.norm(rhs-system@coefficients)/np.linalg.norm(rhs))
                mean_error=float(np.linalg.norm(np.asarray(mean)-expected_mean)/np.linalg.norm(expected_mean))
                sample=dict(case=case['name'],backend=backend,phase=phase,whole_seconds=elapsed,
                            user_api_total_seconds=elapsed,
                            matlab_native_total_seconds=elapsed if backend=='gpu' else None,
                            outer_controller_elapsed_seconds=controller_elapsed,
                            iterations=iterations,
                            coefficient_error=error,mean_error=mean_error,true_relative_residual=residual,
                            stages=stages,reuse='none_across_samples',host_process_memory=memory.result(),
                            blas_before=blas_before,blas_after=observed_blas(),
                            exact_accounted_owned_peak_bytes=owned,idle_before=before,idle_after=idle())
                result['samples'].append(sample)
                if not all(math.isfinite(v) for v in (error,mean_error,residual,elapsed,controller_elapsed)) or min(elapsed,controller_elapsed)<=0 or max(error,mean_error)>1e-8 or residual>1e-10:
                    raise RuntimeError('Whole-solve parity/residual failed')
            run('cpu','first_payload');run('gpu','first_payload')
            for cycle in range(args.cycles):
                for backend in ('cpu','gpu','gpu','cpu'):
                    run(backend,'warm')
        result['idle_after']=idle();verify()
        result['warm_medians']={backend:statistics.median(s['whole_seconds'] for s in result['samples'] if s['phase']=='warm' and s['backend']==backend) for backend in ('cpu','gpu')} if len(json.loads((job/'cases.json').read_text()))==1 else None
        result['comparison_by_case']={}
        for name in {s['case'] for s in result['samples']}:
            comparison={}
            for boundary in ('user_api_total_seconds','outer_controller_elapsed_seconds'):
                medians={backend:statistics.median(s[boundary] for s in result['samples'] if s['case']==name and s['phase']=='warm' and s['backend']==backend) for backend in ('cpu','gpu')}
                first={backend:next(s[boundary] for s in result['samples'] if s['case']==name and s['phase']=='first_payload' and s['backend']==backend) for backend in ('cpu','gpu')}
                comparison[boundary]=dict(warm_medians=medians,cpu_over_gpu=medians['cpu']/medians['gpu'],first_payload_seconds=first,first_payload_cpu_over_gpu=first['cpu']/first['gpu'],ratio_scope='Bounded linear public case; backend-specific preconditioners; memory sampler asymmetry included; first payload excludes process startup')
            result['comparison_by_case'][name]=comparison
        result['accepted']=True
    except BaseException as exc:
        result['error']=f'{type(exc).__name__}: {exc}'
        raise
    finally:
        if engine is not None:
            try:
                engine.quit()
                result['owned_engine_quit']=True
            except BaseException as exc:
                result['release_error']=str(exc)
                result['accepted']=False
        (job/'formal-result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
