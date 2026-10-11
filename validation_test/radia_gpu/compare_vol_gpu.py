"""Compare a checked .vol CPU reference with an installed CUDA MEX payload.

Timings are diagnostics for the declared mesh, not a GPU scaling claim.
Run with a MATLAB Engine capable Python and an independently installed wheel.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matlab.engine
import numpy as np
from scipy.io import loadmat

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bundle', required=True, type=Path)
parser.add_argument('--installed-python', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
parser.add_argument('--require-a100', action='store_true')
args = parser.parse_args()
base = args.bundle.resolve()
cpu = json.loads((base / 'cpu-result.json').read_text())
system = base / 'system.mat'
if not cpu['passed'] or hashlib.sha256(system.read_bytes()).hexdigest() != cpu['system_sha256']:
    raise RuntimeError('CPU reference provenance differs from exported system')
installed = json.loads(subprocess.check_output(
    [str(args.installed_python), '-I', '-c',
     'import json,radia_gpu; print(json.dumps(radia_gpu.verify()))'], text=True))
if not installed['passed']:
    raise RuntimeError('Installed payload verification failed')
args.out.mkdir(parents=True, exist_ok=False)
engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
try:
    engine.restoredefaultpath(nargout=0)
    engine.addpath(str(Path(__file__).resolve().parents[1] / 'ngsolve_cuda_mex'), nargout=0)
    engine.benchmark_gpu(str(system), str(args.out.resolve() / 'gpu.mat'),
                         installed['matlab_path'], nargout=0)
finally:
    engine.quit()
gpu = loadmat(args.out / 'gpu.mat', simplify_cells=True)
reference = loadmat(system, simplify_cells=True)['reference'].reshape(-1)
solution = gpu['solution'].reshape(-1)
delta = float(np.linalg.norm(solution-reference) / np.linalg.norm(reference))
gpu_report = gpu['result']
a100 = 'A100' in gpu_report['device']
if not np.isfinite(delta) or delta > 1e-8 or gpu_report['info']['true_relative_residual'] > 1e-6:
    raise RuntimeError('CUDA solve fails NGSolve parity or true residual acceptance')
if args.require_a100 and not a100:
    raise RuntimeError('Final A100 acceptance requires an actual A100 device')
report = dict(schema='radia.gpu.vol-cpu-gpu-comparison/v1', passed=True,
              source_commit=installed['source_commit'], cpu_reference=cpu,
              gpu_vs_ngsolve_direct_relative_l2=delta, gpu_report=gpu_report,
              a100_checked=a100, formal_performance_checked=False)
encoded = json.dumps(report, indent=2, default=lambda value: value.tolist())
(args.out / 'result.json').write_text(encoded+'\n')
print(encoded)
