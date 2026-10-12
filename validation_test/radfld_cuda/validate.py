"""ABBA CPU/CUDA retained-cloud RadFld validation; 100 timings are diagnostic.

The pre-CUDA CPU oracle is independently compiled from its recorded Git commit.
Python orchestrates validation only; the MATLAB GPU API has no Python adapter.
"""
from pathlib import Path
import argparse
import hashlib
import json
import statistics
import sys
import time
import importlib.util
import ctypes
import subprocess

import h5py
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def relative(value, reference):
    return float(np.linalg.norm(value-reference)/max(np.linalg.norm(reference), np.finfo(float).tiny))


def point_relative(value, reference):
    scale = np.linalg.norm(reference, axis=1)
    floor = max(float(scale.max())*1e-12, np.finfo(float).tiny)
    return float(np.max(np.linalg.norm(value-reference, axis=1)/np.maximum(scale, floor)))


def cpu_memory():
    """Windows process working set and private bytes, sampled outside timing."""
    class Counters(ctypes.Structure):
        _fields_ = [('cb', ctypes.c_ulong), ('PageFaultCount', ctypes.c_ulong)] + [
            (name, ctypes.c_size_t) for name in ('PeakWorkingSetSize', 'WorkingSetSize',
             'QuotaPeakPagedPoolUsage', 'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage',
             'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage', 'PrivateUsage')]
    counter = Counters()
    counter.cb = ctypes.sizeof(counter)
    ctypes.windll.kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    handle = ctypes.windll.kernel32.GetCurrentProcess()
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes = [ctypes.c_void_p, ctypes.POINTER(Counters), ctypes.c_ulong]
    if not query(handle, ctypes.byref(counter), counter.cb):
        raise ctypes.WinError()
    return dict(working_set_bytes=counter.WorkingSetSize, private_bytes=counter.PrivateUsage,
                process_lifetime_peak_working_set_bytes=counter.PeakWorkingSetSize)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--reference-dir', type=Path, required=True)
    p.add_argument('--sizes', type=int, nargs='+', default=[4096, 16384, 65536])
    p.add_argument('--targets', type=int, default=1024)
    p.add_argument('--rounds', type=int, default=2)
    p.add_argument('--threads', type=int, default=4)
    args = p.parse_args()
    if min(args.sizes+[args.targets, args.rounds, args.threads]) < 1:
        p.error('Sizes, rounds and threads must be positive')
    args.out.mkdir(parents=True, exist_ok=True)
    # Reuse the existing tested CI/process/GPU contention guard.
    spec = importlib.util.spec_from_file_location('ams_graph_validation', ROOT/'validation_test/ams_cuda/graph_validation.py')
    guard = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)
    import ngsolve as ng
    import netgen  # Initialize native DLL search paths before the oracle import.
    sys.path.insert(0, str(args.reference_dir))
    import _radia_field_reference as oracle
    import matlab.engine
    ng.SetNumThreads(args.threads)
    binary = ROOT/'build-field-cuda/radia_cuda_field.mexw64'
    receipt = json.loads((binary.parent/'build.json').read_text(encoding='utf-8'))
    reference_receipt = json.loads((args.reference_dir/'reference.json').read_text(encoding='utf-8'))
    if digest(binary) != receipt['binary_sha256'] or digest(oracle.__file__) != reference_receipt['binary_sha256']:
        raise RuntimeError('Build artifact hash mismatch')
    record = dict(schema='radia.radfld-cuda-validation.v1', accepted=False,
                  timing_classification='100 RTX 5000 diagnostic; not idle-mdx performance acceptance',
                  scope='existing retained point cloud only; no analytic/curved element substitution; no generic RadFld dispatch yet',
                  build=receipt, cpu_reference=reference_receipt, threads=args.threads,
                  script_sha256=digest(__file__), matlab_script_sha256=digest(Path(__file__).with_name('runFieldCase.m')),
                  matlab_wrapper_sha256=digest(ROOT/'packages/radia-gpu/matlab/+radia/+gpu/FieldCloud.m'),
                  idle_guard_sha256=digest(ROOT/'validation_test/ams_cuda/graph_validation.py'),
                  validation_commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
                  validation_dirty=bool(subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip()),
                  transfer_inclusive=True, a100_used=False, cases=[],
                  numerical_limits=dict(cuda_cpu_same_algorithm_relative=1e-12, tree_cpu_direct_relative=1e-5),
                  memory_limits='CUDA owned allocations are exact simultaneous buffer counts. cudaMemGetInfo and MATLAB available/host memory are samples. Python working set/private bytes are process samples; its peak is process-lifetime, not case-local. Host native per-case peak is not measured.',
                  timing_scope='CPU: construct original evaluator, evaluate, return host output. CUDA: public MATLAB constructor, source/tree upload, observation upload, evaluate, result download, native destroy and wait. HDF5 read and numerical validation excluded. MATLAB startup recorded separately.')
    eng = None
    try:
        guard.assert_idle_host()
        started = time.perf_counter()
        eng = matlab.engine.start_matlab('-nodesktop -nosplash')
        record['matlab_startup_seconds'] = time.perf_counter()-started
        eng.addpath(str(binary.parent), str(ROOT/'packages/radia-gpu/matlab'), str(Path(__file__).parent), nargout=0)
        record['matlab'] = eng.version()
        record['gpu'] = json.loads(eng.eval("jsonencode(struct('name',gpuDevice().Name,'total_memory',gpuDevice().TotalMemory,'available_memory',gpuDevice().AvailableMemory,'compute_capability',gpuDevice().ComputeCapability))"))
        guard.assert_idle_host()
        record['lifecycle'] = json.loads(eng.eval('jsonencode(testFieldCloud())'))
        for size in args.sizes:
            guard.assert_idle_host()
            rng = np.random.default_rng(20261012+size)
            xyz = rng.uniform(-0.5, 0.5, (size, 3))
            q = rng.normal(size=size)
            obs = rng.uniform(1.5, 3.0, (args.targets, 3))*rng.choice([-1., 1.], (args.targets, 3))
            # Include genuine near-field and coincident retained-cloud cases.
            obs[:min(16, args.targets)] = xyz[:min(16, args.targets)]+[0.001, -0.002, 0.003]
            obs[0] = xyz[0]
            masks = np.array([1., 5.]); signs = np.array([1., -1.]); angles = np.array([0.3, -0.2])
            leaf, theta = 32, 0.025
            input_file = args.out/f'input-{size}.h5'
            with h5py.File(input_file, 'w') as h:
                for name, value in dict(xyz=xyz.T, strength=q, observations=obs.T,
                                        masks=masks, signs=signs, angles=angles,
                                        leaf=float(leaf), theta=theta).items():
                    h.create_dataset(name, data=value)
            def cpu(algorithm):
                memory_before = cpu_memory()
                with ng.TaskManager():
                    start = time.perf_counter()
                    f = oracle.Field(xyz, q, masks, signs, angles, leaf, theta)
                    value = np.asarray(f.field(obs, algorithm))
                    del f
                    elapsed = time.perf_counter()-start
                return value, elapsed, dict(before=memory_before, after=cpu_memory())
            direct, _, _ = cpu('direct')
            references = {name: cpu(name)[0] for name in ('direct', 'tree')}
            tree_error = relative(references['tree'], direct)
            case = dict(sources=size, targets=args.targets, leaf_size=leaf, theta=theta,
                        image_count=2, input_sha256=digest(input_file), input_bytes=input_file.stat().st_size,
                        cpu_tree_vs_direct_relative=tree_error, samples=[])
            case['cpu_tree_vs_direct_max_point_relative'] = point_relative(references['tree'], direct)
            record['cases'].append(case)
            if tree_error > 1e-5 or case['cpu_tree_vs_direct_max_point_relative'] > 1e-5:
                raise RuntimeError('CPU tree accuracy exceeds declared tolerance')
            # Warm both native paths, including allocation/upload and destruction.
            for algorithm in ('direct', 'tree'):
                guard.assert_idle_host()
                eng.runFieldCase(str(input_file), algorithm, str(args.out/'warmup.mat'), nargout=0)
                for index, backend in enumerate(('cpu', 'cuda', 'cuda', 'cpu')*args.rounds):
                    guard.assert_idle_host()
                    if backend == 'cpu':
                        value, elapsed, host_memory = cpu(algorithm)
                        sample = dict(algorithm=algorithm, backend=backend, total_seconds=elapsed,
                                      cpu_host_memory=host_memory)
                    else:
                        raw = args.out/f'{size}-{algorithm}-{index}.mat'
                        eng.runFieldCase(str(input_file), algorithm, str(raw), nargout=0)
                        r = loadmat(raw, simplify_cells=True)
                        value = np.asarray(r['value']).reshape((-1, 3))
                        sample = dict(algorithm=algorithm, backend=backend, total_seconds=float(r['totalSeconds']),
                                      setup=r['setupInfo'], field=r['info'],
                                      matlab_memory_before_bytes=float(r['memoryBefore']),
                                      matlab_memory_after_bytes=float(r['memoryAfter']),
                                      matlab_host_used_before_bytes=float(r['hostBytesBefore']),
                                      matlab_host_used_after_bytes=float(r['hostBytesAfter']))
                    sample['relative_to_cpu_same_algorithm'] = relative(value, references[algorithm])
                    sample['relative_to_cpu_direct'] = relative(value, direct)
                    sample['max_point_relative_to_cpu_same_algorithm'] = point_relative(value, references[algorithm])
                    sample['max_point_relative_to_cpu_direct'] = point_relative(value, direct)
                    case['samples'].append(sample)
                    if (sample['relative_to_cpu_same_algorithm'] > 1e-12 or sample['relative_to_cpu_direct'] > 1e-5
                        or sample['max_point_relative_to_cpu_same_algorithm'] > 1e-12
                        or sample['max_point_relative_to_cpu_direct'] > 1e-5):
                        raise RuntimeError('CUDA/CPU field parity failed')
            case['comparison'] = {}
            for algorithm in ('direct', 'tree'):
                before = statistics.median(s['total_seconds'] for s in case['samples'] if s['algorithm']==algorithm and s['backend']=='cpu')
                after = statistics.median(s['total_seconds'] for s in case['samples'] if s['algorithm']==algorithm and s['backend']=='cuda')
                case['comparison'][algorithm] = dict(before_seconds=before, after_seconds=after, speedup=before/after)
        if digest(binary) != receipt['binary_sha256']:
            raise RuntimeError('CUDA binary changed during validation')
        record['accepted'] = True
    except BaseException as e:
        record['error'] = str(e)
        raise
    finally:
        if eng is not None:
            eng.quit()
        (args.out/'result.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
