"""Compare direct CUDA launches and Graph replay on identical AMS HDF5 inputs.

Timing on RTX 5000/100 is diagnostic only. Python orchestrates validation;
it is not required by the MATLAB GPU runtime.
"""
from pathlib import Path
import argparse, json, os, statistics, sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[2]

def main(args):
    import netgen, matlab.engine
    net=Path(netgen.__file__).parent
    os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(net.parent/'ngsolve.libs'),str(net.parent/'netgen_mesher.libs'),str(Path(sys.base_prefix)),os.environ['PATH']])
    args.out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.ams-cuda-graph-comparison.v1',accepted=False,
        timing_classification='diagnostic; not formal performance acceptance',a100_started=False,
        comparison='direct CUDA kernels versus CUDA Graph; same host setup threads and GPU arithmetic',
        setup_threads=4,cases=[])
    eng=None
    try:
        eng=matlab.engine.start_matlab('-nodesktop -nosplash')
        eng.addpath(str(Path(__file__).parent),nargout=0)
        record['build']=json.loads((ROOT/'build-ams-cuda/build.json').read_text())
        for path in args.inputs:
            case=dict(input=path.name,samples=[]);record['cases'].append(case)
            ref=None
            for index,graph in enumerate([False,True,True,False]):
                target=args.out/f'{path.stem}-{index}.mat'
                eng.benchmark(str(path.resolve()),str(target.resolve()),str(ROOT/'build-ams-cuda'),
                    str(ROOT/'packages/radia-gpu/matlab'),4.0,graph,nargout=0)
                r=loadmat(target,simplify_cells=True)
                actual=np.asarray(r['inclusiveSolution']).reshape(-1)
                if ref is None:ref=actual.copy()
                delta=float(np.linalg.norm(actual-ref)/max(np.linalg.norm(ref),np.finfo(float).tiny))
                sample=dict(graph=graph,setup_seconds=float(r['nativeSetupSeconds']),
                    transfer_solve_download_seconds=float(r['inclusiveSeconds']),
                    total_seconds=float(r['nativeSetupSeconds']+r['inclusiveSeconds']),
                    warmed_solve_median_seconds=float(np.median(r['times'])),
                    direct_solution_relative_difference=delta,probe_relative_difference=float(np.max(r['probeErrors'])),
                    hierarchy_relative_difference=float(r['hierarchyMaxRelativeDifference']),
                    true_relative_residual=float(r['info']['true_relative_residual']),
                    preconditioner_execution=r['info']['preconditioner_execution'],
                    gpu_live_allocation_delta_bytes=float(r['memoryBytes']))
                case['samples'].append(sample)
                if delta>1e-10 or sample['true_relative_residual']>1e-6 or sample['probe_relative_difference']>1e-8:
                    raise RuntimeError('CUDA Graph numerical parity failed')
            case['before_after']={}
            for key in ['transfer_solve_download_seconds','total_seconds','warmed_solve_median_seconds']:
                before=statistics.median(s[key] for s in case['samples'] if not s['graph'])
                after=statistics.median(s[key] for s in case['samples'] if s['graph'])
                case['before_after'][key]=dict(before=before,after=after,speedup=before/after)
        record['accepted']=True
    except BaseException as e:
        record['error']=repr(e);raise
    finally:
        if eng is not None:eng.quit()
        (args.out/'result.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--inputs',type=Path,nargs='+',required=True)
    main(p.parse_args())
