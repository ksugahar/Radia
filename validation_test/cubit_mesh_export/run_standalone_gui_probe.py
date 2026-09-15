"""Run the explicitly authorized GUI smoke with a wheel-only interpreter."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

assert importlib.util.find_spec('radia') is None, 'Use the isolated wheel venv'
import cubit_mesh_export
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=False)
env = os.environ.copy()
for key in ('PYTHONPATH', 'PYTHONHOME', 'CUBIT_PLUGIN_DIR'):
    env.pop(key, None)
env['CME_GUI_TEST_OUTPUT'] = str(out)
env['CME_GUI_TEST_SOURCE'] = str(Path(cubit_mesh_export.__file__).parent / 'cubit_gui')
env['CUBIT_MESH_EXPORT_PYTHON'] = sys.executable
probe = Path(__file__).with_name('standalone_gui_probe.py')
with (out / 'cubit.log').open('w', encoding='utf-8') as log:
    process = subprocess.Popen(['C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.exe',
                                '-noinitfile', '-nojournal', '-commandplugindir',
                                'C:/Program Files/Coreform Cubit 2025.12/bin/plugins', str(probe)],
                               cwd=out, env=env, stdin=subprocess.PIPE, stdout=log, stderr=log)
    print('Owned Cubit PID:', process.pid, flush=True)
    try:
        code = process.wait(timeout=60)
    except subprocess.TimeoutExpired:
        process.stdin.write(b'exit 0\n')
        process.stdin.flush()
        code = process.wait(timeout=15)
    print('Exit:', code)
result = out / 'gui-result.json'
print(result.read_text() if result.exists() else 'GUI result missing; inspect cubit.log')
payload = json.loads(result.read_text()) if result.exists() else {'passed': False}
payload['cubit_exit_code'] = code
payload['probe_passed'] = payload.get('passed', False)
payload['passed'] = payload['probe_passed'] and code == 0
(out / 'result.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
sys.exit(0 if payload['passed'] else 1)
