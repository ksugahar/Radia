"""Idle same-host direct CPU versus installed-wheel FMM whole-API ABBA."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import statistics
import sys
import time


def main(args):
    job = args.job.resolve()
    sys.path.insert(0, str(job))
    from measure_native_cpu_gpu import Memory, idle, encode
    from fmm_numeric import direct
    import verify_wheel
    import radia_gpu
    import numpy as np
    from scipy.io import loadmat
    import matlab.engine
    expected_identity = verify_wheel.verify(args.wheel)
    verify_wheel.verify_installed_identity(radia_gpu.verify(), expected_identity)
    manifest = json.loads((job / 'manifest.json').read_text())

    def verify_inputs():
        for name, digest in manifest.items():
            path = (job / name).resolve()
            if not path.is_relative_to(job) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise RuntimeError('Measurement input identity changed')

    verify_inputs()
    data = loadmat(job / 'numeric-input.mat')
    if float(data['order'].item()) != 8 or float(data['separation'].item()) != 2:
        raise ValueError('This acceptance requires Order=8, Separation=2')
    record = dict(schema='radia.installed-fmm-performance.v1', accepted=False,
                  timing_accepted=False, release_accepted=False,
                  wheel=expected_identity, input_source_sha256=manifest, samples=[],
                  order=8, separation=2, cpu_threads=1,
                  cpu_algorithm='Independent NumPy direct sum, not CPU FMM',
                  cold_definition='First numerical call; interpreter/imports, MATLAB startup and GPU initialization are separate, not process-cold solver samples',
                  timing_boundary='Context construction, all uploads, apply, result download and context destruction; input-file read and output serialization excluded for both routes',
                  memory_limitations='50ms process samples may miss short spikes; lifetime RSS peak includes earlier work; native reported components do not certify global device, CUDA runtime or all MATLAB temporaries')
    engine = None
    oracle = None
    try:
        record['idle_before'] = idle()
        if 'A100' not in record['idle_before']['gpu']:
            raise RuntimeError('Formal FMM evaluation requires the selected A100 device')
        start = time.perf_counter()
        engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        record['matlab_start_seconds'] = time.perf_counter() - start
        payload = radia_gpu.matlab_path()
        engine.addpath(str(payload), str(job), nargout=0)
        loaded = Path(engine.which('radia_cuda_fmm'))
        installed_manifest = json.loads((payload.parent / 'manifest.json').read_text())
        if loaded.parent.resolve() != payload.resolve() or hashlib.sha256(loaded.read_bytes()).hexdigest() != installed_manifest['files']['matlab/radia_cuda_fmm.mexw64']:
            raise RuntimeError('FMM must load from the verified installed wheel')
        record['fmm_binary_sha256'] = hashlib.sha256(loaded.read_bytes()).hexdigest()
        record['matlab_version'] = engine.version()
        pid = int(engine.feature('getpid'))
        engine.eval('t=tic;d=gpuDevice;wait(d);init_seconds=toc(t);', nargout=0)
        record['gpu_initialization_seconds'] = float(engine.workspace['init_seconds'])
        for index, route in enumerate(['cpu', 'gpu'] + ['cpu', 'gpu', 'gpu', 'cpu'] * args.cycles):
            guard = idle()
            with Memory(os.getpid() if route == 'cpu' else pid) as memory:
                if route == 'cpu':
                    start = time.perf_counter()
                    values = direct(data['sources'], data['targets'], data['charge'], data['dipole'])
                    seconds = time.perf_counter() - start
                    native = {}
                    if oracle is None:
                        oracle = values
                else:
                    output = job / f'sample-{index}.mat'
                    engine.measure_native_sample(str(job), 'fmm-small', str(output), nargout=0)
                    sample = loadmat(output, simplify_cells=True)
                    values = sample['values']
                    native = sample['report']
                    seconds = float(native['whole_seconds'])
            if len(values) != 2 or any(not np.all(np.isfinite(x)) for x in values):
                raise RuntimeError('Invalid potential/field outputs')
            errors = [float(np.linalg.norm(x - y) / np.linalg.norm(y)) for x, y in zip(values, oracle)]
            if not np.all(np.isfinite(errors)) or max(errors) > 1e-8:
                raise RuntimeError(f'Independent direct-sum parity failed: {errors}')
            record['samples'].append(dict(backend=route, phase='first' if index < 2 else 'ABBA',
                whole_seconds=seconds, independent_relative_errors=errors,
                idle_before=guard, host_process_memory=memory.result(), native=native))
            (job / 'fmm-performance.json').write_text(json.dumps(record, indent=2, default=encode))
        warm = record['samples'][2:]
        cpu = statistics.median(x['whole_seconds'] for x in warm if x['backend'] == 'cpu')
        gpu = statistics.median(x['whole_seconds'] for x in warm if x['backend'] == 'gpu')
        record['summary'] = dict(cpu_median_seconds=cpu, gpu_transfer_inclusive_median_seconds=gpu,
                                 cpu_direct_over_gpu_speedup=cpu / gpu)
        record['idle_after'] = idle()
        verify_inputs()
        record.update(accepted=True, timing_accepted=True)
    except BaseException as error:
        record['error'] = repr(error)
        raise
    finally:
        if engine is not None:
            try:
                engine.quit()
            except BaseException as error:
                record.update(accepted=False, timing_accepted=False, release_error=repr(error))
        (job / 'fmm-performance.json').write_text(json.dumps(record, indent=2, default=encode))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', required=True, type=Path)
    parser.add_argument('--wheel', required=True, type=Path)
    parser.add_argument('--cycles', type=int, default=2)
    options = parser.parse_args()
    if options.cycles < 2:
        parser.error('At least two ABBA cycles required')
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    main(options)
