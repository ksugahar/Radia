"""Validate an installed radia-gpu wheel through an owned MATLAB Engine session."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import matlab.engine

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--installed-python', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
parser.add_argument('--wheel', type=Path, help='Record verified exact wheel provenance for release')
parser.add_argument('--gpu-simulink-only', action='store_true',
                    help='Explicitly exclude external LTspice acceptance; record it as unchecked')
args = parser.parse_args()
report = json.loads(subprocess.check_output(
    [str(args.installed_python), '-I', '-c',
     'import json,radia_gpu; print(json.dumps(radia_gpu.verify()))'], text=True))
if not report['passed']:
    raise RuntimeError('Installed payload verification failed')
wheel_report = None
if args.wheel:
    verifier_path = Path(__file__).resolve().parents[2] / 'packages/radia-gpu/verify_wheel.py'
    spec = importlib.util.spec_from_file_location('radia_gpu_wheel_verifier', verifier_path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel_report = verifier.verify(args.wheel)
    verifier.verify_installed_identity(report, wheel_report)
args.out.resolve().mkdir(parents=True, exist_ok=False)
engine = matlab.engine.start_matlab('-nodesktop -nosplash')
try:
    # This session is owned by acceptance: exclude startup-added checkout paths.
    engine.restoredefaultpath(nargout=0)
    engine.addpath(str(Path(__file__).resolve().parent), nargout=0)
    engine.check_installed(report['matlab_path'], str(args.out.resolve()),
                           not args.gpu_simulink_only, nargout=0)
finally:
    engine.quit()
result_path = args.out / 'result.json'
result = json.loads(result_path.read_text())
result.update(source_commit=report['source_commit'], distribution_version=report['version'],
              manifest_sha256=report['manifest_sha256'])
if wheel_report:
    result['wheel_sha256'] = wheel_report['wheel_sha256']
result_path.write_text(json.dumps(result,indent=2)+'\n')
print(result_path.read_text())
