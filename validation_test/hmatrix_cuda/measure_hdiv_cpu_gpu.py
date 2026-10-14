"""Idle same-host independent CPU versus transfer-inclusive MATLAB HDiv timing.

Consumes a hash-verified job payload and an independently built CPU reference.
Python performs CPU validation and orchestration, never GPU solver work.
"""
from pathlib import Path
import argparse, ctypes, hashlib, json, os, statistics, sys, time, math
from validate_hdiv_artifact import idle


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
    args = parser.parse_args()
    if args.cpu_threads < 1:
        parser.error('CPU threads must be positive')
    if args.cycles < 2:
        parser.error('At least two ABBA cycles are required')
    job = args.job.resolve()
    result = dict(accepted=False, release_accepted=False, samples=[],
                  classification='Formal same-host CPU/GPU measurement',
                  memory_limitations='Exact GPU owned buffers and available-memory samples only; host/global/CUDA runtime allocation peaks unmeasured',
                  cold_limitations='Fresh owned Engine startup and first payload separately; CPU interpreter/import startup excluded')
    result['timing_boundary']='Public CPU whole solve and MATLAB solveHDivElements through result recovery and destruction; GPU Engine transport/report collection excluded'
    engine = None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
        os.environ.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', MKL_THREADING_LAYER='SEQUENTIAL')
        site=Path(sys.prefix)/'Lib/site-packages'
        dependency_dirs=[job/'runtime',Path(sys.prefix)/'bin',Path(sys.prefix)/'Library/bin',
                         site/'netgen',site/'netgen_mesher.libs',site/'ngsolve.libs',
                         site/'ngsolve_openblas',Path(sys.base_prefix)]
        os.environ['PATH']=os.pathsep.join(str(p) for p in dependency_dirs)+os.pathsep+os.environ['PATH']
        dll_handles=[os.add_dll_directory(str(p)) for p in dependency_dirs if p.is_dir()]
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
        from scipy.io import loadmat
        import matlab.engine
        ng.SetNumThreads(args.cpu_threads)
        result['threads']={'cpu_ngsolve':args.cpu_threads,'cpu_blas':1,'matlab_single_comp_thread':True}
        for module in (cpu,native):
            relative=Path(module.__file__).resolve().relative_to(cpu_root).as_posix()
            if relative not in cpu_manifest:
                raise RuntimeError('Loaded CPU module missing from reference manifest: '+relative)
        result['reference'] = dict(native_sha256=sha(native.__file__), frontend_sha256=sha(cpu.__file__),
                                   ngsolve=ng.__version__, python=sys.version,
                                   limitation='Loaded file identity; independent native build receipt must accompany release evidence')
        result['payload_sha256'] = manifest
        result['cpu_reference_sha256'] = cpu_manifest
        matlab_path = job/'matlab'
        if args.wheel:
            import importlib.util
            import radia_gpu
            spec=importlib.util.spec_from_file_location('wheel_verifier',job/'verify_wheel.py')
            verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
            wheel=verifier.verify(args.wheel)
            verifier.verify_installed_identity(radia_gpu.verify(),wheel)
            result['wheel']=wheel
            matlab_path=radia_gpu.matlab_path()
        result['idle_before'] = idle()
        start = time.perf_counter()
        engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        result['engine_start_seconds'] = time.perf_counter()-start
        result['matlab'] = engine.version()
        engine.restoredefaultpath(nargout=0)
        engine.addpath(str(job), str(matlab_path), nargout=0)
        if args.wheel:
            installed_manifest=json.loads((matlab_path.parent/'manifest.json').read_text())
            for name in ('radia_cuda_hmatrix','radia_cuda_projection'):
                loaded=Path(engine.which(name)).resolve()
                if loaded.parent!=matlab_path.resolve() or sha(loaded)!=installed_manifest['files']['matlab/'+loaded.name]:
                    raise RuntimeError('HDiv must load from the verified installed wheel')
        for case in json.loads((job/'cases.json').read_text()):
            fixture = loadmat(job/case['fixture'], simplify_cells=True)
            expected = np.asarray(fixture['expected']).ravel()
            system, rhs = fixture['system'], np.asarray(fixture['rhs']).ravel()
            engine.workspace['mesh_path'] = str(job/case['mesh'])
            engine.workspace['degree'] = float(case['order'])
            engine.workspace['assembly_mode'] = case.get('assembly_mode','host-cache')
            engine.workspace['fixture_path'] = str(job/case['fixture'])
            engine.eval('cpu=load(fixture_path);', nargout=0)
            def run(backend, phase):
                before = idle()
                start = time.perf_counter()
                if backend == 'cpu':
                    mesh = ng.Mesh(str(job/case['mesh']))
                    with ng.TaskManager():
                        fes = ng.HDiv(mesh, order=case['order'])
                        B,G,mass = cpu.build_charge_gram(fes,intorder=4,
                            curve_order=2 if case.get('geometry','curved')=='curved' else 0,
                            curve_gauss=8,ho_far_factor=1e30,gram_backend='exact-dense',exact_dense_memory_mb=64)
                        load=ng.LinearForm(fes)
                        load += ng.CoefficientFunction((1.,.2,.1))*fes.TestFunction()*ng.dx
                        load.Assemble()
                        answer=G.solve_configured_linear_material_mass_riesz(2.5,np.array(load.vec.FV().NumPy()),1e-11,1000)
                        coefficients=np.asarray(answer['m']).copy()
                        mean=[]
                        for c in range(3):
                            form=ng.LinearForm(fes);form+=fes.TestFunction()[c]*ng.dx;form.Assemble()
                            mean.append(float(np.dot(form.vec.FV().NumPy(),coefficients)/ng.Integrate(1,mesh)))
                    del answer,B,G,mass,load,fes,mesh
                    elapsed=time.perf_counter()-start
                    owned=None
                else:
                    engine.eval("timer=tic;solution=radia.gpu.solveHDivElements(mesh_path,degree,[1 .2 .1],2.5,1e-11,1000,assembly_mode);whole_seconds=toc(timer);assert(radia_cuda_hmatrix('live_count')==0);",nargout=0)
                    coefficients=np.asarray(engine.eval('solution.coefficients')).ravel()
                    mean=np.asarray(engine.eval('solution.mean_magnetization')).ravel()
                    owned=float(engine.eval('solution.combined_peak_owned_device_bytes'))
                    elapsed=float(engine.workspace['whole_seconds'])
                error=float(np.linalg.norm(coefficients-expected)/np.linalg.norm(expected))
                residual=float(np.linalg.norm(rhs-system@coefficients)/np.linalg.norm(rhs))
                expected_mean=fixture['mean_operator']@expected
                mean_error=float(np.linalg.norm(np.asarray(mean)-expected_mean)/np.linalg.norm(expected_mean))
                sample=dict(case=case['name'],backend=backend,phase=phase,whole_seconds=elapsed,
                            coefficient_error=error,mean_error=mean_error,true_relative_residual=residual,
                            exact_accounted_owned_peak_bytes=owned,idle_before=before,idle_after=idle())
                result['samples'].append(sample)
                if not all(math.isfinite(v) for v in (error,mean_error,residual,elapsed)) or elapsed<0 or max(error,mean_error)>1e-8 or residual>1e-10:
                    raise RuntimeError('Whole-solve parity/residual failed')
            run('cpu','first_payload');run('gpu','first_payload')
            for cycle in range(args.cycles):
                for backend in ('cpu','gpu','gpu','cpu'):
                    run(backend,'warm')
        result['idle_after']=idle();verify()
        result['warm_medians']={backend:statistics.median(s['whole_seconds'] for s in result['samples'] if s['phase']=='warm' and s['backend']==backend) for backend in ('cpu','gpu')} if len(json.loads((job/'cases.json').read_text()))==1 else None
        result['accepted']=True
    except BaseException as exc:
        result['error']=f'{type(exc).__name__}: {exc}'
        raise
    finally:
        if engine is not None:
            try:
                engine.quit()
            except BaseException as exc:
                result['release_error']=str(exc)
                result['accepted']=False
        (job/'formal-result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
