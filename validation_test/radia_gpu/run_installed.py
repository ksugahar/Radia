"""Validate an installed radia-gpu wheel through an owned MATLAB Engine session."""
import argparse
import json
from pathlib import Path
import subprocess
import matlab.engine

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--installed-python', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
args = parser.parse_args()
report = json.loads(subprocess.check_output(
    [str(args.installed_python), '-c',
     'import json,radia_gpu; print(json.dumps(radia_gpu.verify()))'], text=True))
if not report['passed']:
    raise RuntimeError('Installed payload verification failed')
args.out.resolve().mkdir(parents=True, exist_ok=False)
engine = matlab.engine.start_matlab('-nodesktop -nosplash')
try:
    engine.addpath(str(Path(__file__).resolve().parent), nargout=0)
    engine.check_installed(report['matlab_path'], str(args.out.resolve()), nargout=0)
finally:
    engine.quit()
print((args.out / 'result.json').read_text())
