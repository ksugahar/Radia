"""Validate an installed radia-simulink wheel through an owned MATLAB Engine session."""
import argparse
import importlib.util
import json
import socket
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
     'import json,radia_simulink; print(json.dumps(radia_simulink.verify()))'], text=True))
if not report['passed']:
    raise RuntimeError('Installed payload verification failed')
wheel_report = None
if args.wheel:
    verifier_path = Path(__file__).resolve().parents[2] / 'packages/radia-simulink/verify_wheel.py'
    spec = importlib.util.spec_from_file_location('radia_simulink_wheel_verifier', verifier_path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel_report = verifier.verify(args.wheel)
    verifier.verify_installed_identity(report, wheel_report)
args.out.resolve().mkdir(parents=True, exist_ok=False)
engine = matlab.engine.start_matlab('-nodesktop -nosplash')
ic_acceptance = None
try:
    # This session is owned by acceptance: exclude startup-added checkout paths.
    engine.restoredefaultpath(nargout=0)
    engine.addpath(str(Path(__file__).resolve().parent), nargout=0)
    engine.check_installed(report['matlab_path'], str(args.out.resolve()),
                           not args.gpu_simulink_only, nargout=0)
    payload = Path(report['matlab_path'])
    if (payload/'radia_cuda_ic.mexw64').is_file():
        ic_script = Path(__file__).resolve().parents[1]/'ic_cuda/run.py'
        spec = importlib.util.spec_from_file_location('installed_ic_accuracy',ic_script)
        checker = importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)
        ic_acceptance = checker.check(engine,payload,payload,args.out/'ic')
        manifest = json.loads((payload.parent/'manifest.json').read_text())
        if ic_acceptance['binary_sha256'] != manifest['files']['matlab/radia_cuda_ic.mexw64']:
            raise RuntimeError('Installed IC differs from wheel manifest')
except BaseException as error:
    failed_path=args.out/'result.json'
    failed=json.loads(failed_path.read_text()) if failed_path.is_file() else {}
    failed.update(passed=False,ic_error=repr(error));failed_path.write_text(json.dumps(failed,indent=2)+'\n')
    raise
finally:
    engine.quit()
result_path = args.out / 'result.json'
result = json.loads(result_path.read_text())
if ic_acceptance is not None:
    result.update(ic_checked=True,ic_acceptance=ic_acceptance)
result.update(source_commit=report['source_commit'], distribution_version=report['version'],
              manifest_sha256=report['manifest_sha256'], os_hostname=socket.gethostname())
if wheel_report:
    result['wheel_sha256'] = wheel_report['wheel_sha256']
result_path.write_text(json.dumps(result,indent=2)+'\n')
print(result_path.read_text())
