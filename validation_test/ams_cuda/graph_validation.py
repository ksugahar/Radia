"""Compare AMS Graph replay or CUDA matrix generation on identical HDF5 inputs.

Timing on RTX 5000/100 is diagnostic only. Python orchestrates validation;
it is not required by the MATLAB GPU runtime.
"""
from pathlib import Path
import argparse, hashlib, json, os, statistics, subprocess, sys, time
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[2]

def assert_idle_host():
    """Check between calls, outside timed intervals; never stop another job."""
    if os.name!='nt':
        raise RuntimeError('The idle-host probe currently requires Windows')
    for status in ('queued','in_progress'):
        runs=json.loads(subprocess.check_output(['gh','run','list','--repo','ksugahar/Radia',
            '--status',status,'--limit','1','--json','status'],text=True))
        if runs:
            raise RuntimeError('CI/preflight is active or queued; retry on an idle host')
    probe=r'''
$all=Get-CimInstance Win32_Process
$owned=@(OWNER_PID)
do {
    $previous=$owned.Count
    $owned+=@($all | Where-Object { $_.ParentProcessId -in $owned -and $_.ProcessId -notin $owned } | ForEach-Object ProcessId)
} while ($owned.Count -gt $previous)
$ids=@($all | Where-Object { $_.Name -match '^(python[\d.w]*|MATLAB|nvcc|ninja|cl|MSBuild|ngsolve)\.exe$' -and $_.ProcessId -notin $owned } | ForEach-Object ProcessId)
$before=@{}
foreach($pidValue in $ids) { $p=Get-Process -Id $pidValue -ErrorAction SilentlyContinue; if($p) { $before[$pidValue]=$p.CPU } }
Start-Sleep -Seconds 1
$busy=@(foreach($pidValue in $ids) { $p=Get-Process -Id $pidValue -ErrorAction SilentlyContinue; if($p -and $before.ContainsKey($pidValue) -and ($p.CPU-$before[$pidValue]) -gt 0.25) { $pidValue } })
ConvertTo-Json -InputObject $busy -Compress
'''.replace('OWNER_PID',str(os.getpid()))
    busy=json.loads(subprocess.check_output(['powershell','-NoProfile','-Command',probe],text=True))
    if busy:
        raise RuntimeError('Another numerical/build process is active; retry on an idle host')
    def gpu_util():
        values=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu',
            '--format=csv,noheader,nounits'],text=True)
        return max(float(v) for v in values.split())
    if gpu_util()>25:
        time.sleep(0.5)
        if gpu_util()>25:
            raise RuntimeError('GPU is busy outside this benchmark; retry on an idle host')

def cpu_oracle(path,oracle_dir=ROOT/'build-ams-cuda'):
    """Independent Python/NGSolve solve from the identical checked operators."""
    import ngsolve as ng
    from pyngcore import Array_I_S,Array_D_S
    sys.path.insert(0,str(ROOT/'packages/radia-simulink/src'))
    sys.path.insert(0,str(oracle_dir))
    from radia_simulink.ams_exchange import read_exchange
    import _radia_ams_snapshot as bridge
    ng.SetNumThreads(1)
    data=read_exchange(path);model=data['model']
    def matrix(csr):
        from radia.ngsolve_sparse import real_sparse_matrix_type
        rows=np.repeat(np.arange(int(csr['rows']),dtype=np.int32),np.diff(csr['rowptr']))
        return real_sparse_matrix_type().CreateFromCOO(Array_I_S(rows),Array_I_S(csr['columns']),
            Array_D_S(csr['values']),int(csr['rows']),int(csr['cols']))
    a=matrix(model['A']);g=matrix(model['G']);free=np.asarray(model['free'],dtype=bool)
    reference=bridge.snapshot(a,g,ng.BitArray(free.tolist()),*model['coordinates'].T,
        data['probes'],int(model['cycle']),1)
    probe_diff=float(np.linalg.norm(reference['cpu_probe_results']-data['expected'])/np.linalg.norm(data['expected']))
    if probe_diff>1e-12:raise RuntimeError('CPU oracle input/reference identity failed')
    pre=reference['cpu_preconditioner']
    rhs=a.CreateColVector();rhs.FV().NumPy()[:]=data['rhs']
    x=a.CreateRowVector();x.data=ng.CGSolver(a,pre=pre,tol=1e-9,maxiter=2000)*rhs
    residual=a.CreateColVector();residual.data=a*x-rhs
    true=float(np.linalg.norm(residual.FV().NumPy()[free])/np.linalg.norm(data['rhs'][free]))
    if true>1e-6:raise RuntimeError('CPU true residual failed')
    return dict(a=a,rhs=rhs,x=x,free=free,cpu=x.FV().NumPy().copy(),
        true_residual=true,probe_diff=probe_diff,dofs=int(model['A']['rows']),
        matrix_nonzeros=len(model['A']['values']))

def main(args):
    import netgen, matlab.engine
    net=Path(netgen.__file__).parent
    os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
    args.out.mkdir(parents=True,exist_ok=False)
    generation=args.comparison=='generation'
    record=dict(schema='radia.ams-cuda-matrix-generation.v1' if generation else 'radia.ams-cuda-graph-comparison.v1',accepted=False,
        timing_classification='diagnostic; not formal performance acceptance',a100_started=False,
        comparison='CPU versus CUDA Galerkin generation; identical remaining setup and GPU solve' if generation else 'direct CUDA kernels versus CUDA Graph; same host setup threads and GPU arithmetic',
        setup_threads=4,rounds=args.rounds,cases=[],
        idle_host_checks_enabled=args.require_idle_host,
        memory_measurement='setup: peak concurrently owned CUDA allocation bytes; solve: live allocation delta, not whole-process peak')
    eng=None
    try:
        if args.require_idle_host:assert_idle_host()
        startup=time.perf_counter()
        eng=matlab.engine.start_matlab('-nodesktop -nosplash')
        record['matlab_startup_seconds']=time.perf_counter()-startup
        eng.addpath(str(Path(__file__).parent),nargout=0)
        record['build']=json.loads(args.build_receipt.read_text())
        record['binary_sha256']=file_hash(args.binary_dir/'radia_cuda_ams.mexw64')
        record['oracle_sha256']=file_hash(next(args.oracle_dir.glob('_radia_ams_snapshot*.pyd')))
        record['validator_sha256']=file_hash(Path(__file__))
        record['benchmark_sha256']=file_hash(Path(__file__).with_name('benchmark.m'))
        declared=record['build'].get('binary_sha256',record['build'].get('sha256',{}).get('build-ams-cuda/radia_cuda_ams.mexw64'))
        if declared!=record['binary_sha256']:
            raise RuntimeError('CUDA binary does not match the selected build receipt')
        import ngsolve
        record['ngsolve_version']=ngsolve.__version__
        for path in args.inputs:
            case=dict(input=path.name,input_sha256=file_hash(path),input_bytes=path.stat().st_size,samples=[]);record['cases'].append(case)
            oracle=cpu_oracle(path,args.oracle_dir)
            if args.require_idle_host:assert_idle_host()
            case['cpu_true_relative_residual']=oracle['true_residual']
            case['cpu_probe_reference_relative_difference']=oracle['probe_diff']
            case.update(dofs=oracle['dofs'],matrix_nonzeros=oracle['matrix_nonzeros'])
            # Exclude first native/MATLAB/CUDA load and both backend warmups.
            case['warmups']=[]
            for mode in [False,True]:
                if args.require_idle_host:assert_idle_host()
                warm_path=args.out/f'{path.stem}-warm-{mode}.mat'
                eng.benchmark(str(path.resolve()),str(warm_path.resolve()),str(args.binary_dir),
                    str(args.matlab_dir),4.0,False if generation else mode,mode if generation else True,nargout=0)
                warm=loadmat(warm_path,simplify_cells=True)
                case['warmups'].append(dict(accelerated=mode,native_setup_seconds=float(warm['nativeSetupSeconds']),
                    first_solve_seconds=float(warm['coldSeconds']),pipeline_seconds=float(warm['pipelineSeconds']),
                    excluded_from_medians=True))
            ref=None
            for index,mode in enumerate([False,True,True,False]*args.rounds):
                if args.require_idle_host:assert_idle_host()
                graph=False if generation else mode
                target=args.out/f'{path.stem}-{index}.mat'
                eng.benchmark(str(path.resolve()),str(target.resolve()),str(args.binary_dir),
                    str(args.matlab_dir),4.0,graph,mode if generation else True,nargout=0)
                r=loadmat(target,simplify_cells=True)
                actual=np.asarray(r['pipelineSolution']).reshape(-1)
                if ref is None:ref=actual.copy()
                delta=float(np.linalg.norm(actual-ref)/max(np.linalg.norm(ref),np.finfo(float).tiny))
                cpu_delta=float(np.linalg.norm((actual-oracle['cpu'])[oracle['free']])/np.linalg.norm(oracle['cpu'][oracle['free']]))
                oracle['x'].FV().NumPy()[:]=actual
                residual=oracle['a'].CreateColVector();residual.data=oracle['a']*oracle['x']-oracle['rhs']
                cpu_checked_residual=float(np.linalg.norm(residual.FV().NumPy()[oracle['free']])/np.linalg.norm(oracle['rhs'].FV().NumPy()[oracle['free']]))
                native=r['nativeSetupInfo']
                sample=dict(graph=graph,accelerated=mode,setup_on_gpu=mode if generation else True,setup_seconds=float(r['nativeSetupSeconds']),
                    transfer_solve_download_seconds=float(r['inclusiveSeconds']),
                    total_seconds=float(r['pipelineSeconds']),
                    stage_sum_seconds=float(r['nativeSetupSeconds']+r['inclusiveSeconds']),
                    pipeline_gpu_live_allocation_delta_bytes=float(r['pipelineMemoryBytes']),
                    pipeline_setup=r['pipelineSetup'],
                    gpu_memory_admission_plan=r['memoryPlan'],
                    warmed_solve_median_seconds=float(np.median(r['times'])),
                    cpu_solution_relative_difference=cpu_delta,cpu_checked_true_relative_residual=cpu_checked_residual,
                    direct_solution_relative_difference=delta,probe_relative_difference=float(np.max(r['probeErrors'])),
                    hierarchy_relative_difference=float(r['hierarchyMaxRelativeDifference']),
                    true_relative_residual=float(r['info']['true_relative_residual']),
                    preconditioner_execution=r['info']['preconditioner_execution'],
                    device=r['device'],matlab_version=r['matlabVersion'],
                    gpu_live_allocation_delta_bytes=float(r['memoryBytes']),
                    gpu_setup_peak_owned_bytes=float(native['gpu_setup_peak_owned_bytes']),native_setup=native)
                case['samples'].append(sample)
                if generation and mode and native['gpu_products']<4:raise RuntimeError('Declared CUDA generation did not execute')
                if cpu_delta>1e-6 or cpu_checked_residual>1e-6 or delta>1e-10 or sample['true_relative_residual']>1e-6 or sample['probe_relative_difference']>1e-8 or sample['hierarchy_relative_difference']!=0:
                    raise RuntimeError('AMS accelerated numerical parity failed')
            case['before_after']={}
            for key in ['setup_seconds','transfer_solve_download_seconds','total_seconds','warmed_solve_median_seconds','pipeline_gpu_live_allocation_delta_bytes','gpu_live_allocation_delta_bytes','gpu_setup_peak_owned_bytes']:
                before=statistics.median(s[key] for s in case['samples'] if not s['accelerated'])
                after=statistics.median(s[key] for s in case['samples'] if s['accelerated'])
                case['before_after'][key]=(dict(before=before,after=after,delta=after-before) if key.endswith('_bytes') else dict(before=before,after=after,speedup=before/after))
        record['accepted']=True
    except BaseException as e:
        record['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        def convert(value):
            if isinstance(value,np.ndarray):return value.tolist()
            if isinstance(value,np.generic):return value.item()
            raise TypeError(type(value).__name__)
        (args.out/'result.json').write_text(json.dumps(record,indent=2,default=convert,allow_nan=False)+'\n')

def file_hash(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--inputs',type=Path,nargs='+',required=True)
    p.add_argument('--comparison',choices=['graph','generation'],default='graph')
    p.add_argument('--rounds',type=int,default=1)
    p.add_argument('--binary-dir',type=Path,default=ROOT/'build-ams-cuda')
    p.add_argument('--matlab-dir',type=Path,default=ROOT/'packages/radia-simulink/matlab')
    p.add_argument('--build-receipt',type=Path,default=ROOT/'build-ams-cuda/build.json')
    p.add_argument('--oracle-dir',type=Path,default=ROOT/'build-ams-cuda')
    p.add_argument('--require-idle-host',action='store_true',help='Reject CI/compute/GPU contention between measured calls')
    args=p.parse_args()
    if args.rounds<1:p.error('--rounds must be positive')
    main(args)
