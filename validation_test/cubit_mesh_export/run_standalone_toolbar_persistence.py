"""Automatic registration and two real GUI starts, restoring exact preferences."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

assert importlib.util.find_spec('radia') is None
from cubit_mesh_export import toolbar_install
from cubit_mesh_export.install import _running_cubit_processes

if _running_cubit_processes():
    raise RuntimeError('Close existing Cubit sessions before this explicit GUI test')
root = Path(sys.argv[1]).resolve()
root.mkdir(parents=True, exist_ok=False)
ini = Path(os.environ['APPDATA']) / 'Coreform/Cubit.ini'
original = ini.read_bytes()
(root / 'Cubit.ini.before').write_bytes(original)
env = os.environ.copy()
env['CME_GUI_TEST_TOOLBAR_DESTINATION'] = str(root / 'imported')
runner = Path(__file__).with_name('run_standalone_gui_probe.py')
result = {'passed': False, 'preferences_restored': False}
try:
    # Exercise installer registration with a test-owned payload, preserving
    # other toolbar registrations and every user's installed toolbar payload.
    toolbar_install._profile_toolbar_dir = lambda _: root / 'imported'
    assert toolbar_install._register_official_toolbar(ini)
    assert toolbar_install._ini_has_official_toolbar(ini)
    assert not toolbar_install._register_official_toolbar(ini)
    for mode in ('restart-1', 'restart-2'):
        env['CME_GUI_TEST_TOOLBAR_MODE'] = 'restart'
        completed = subprocess.run([sys.executable, '-I', str(runner), str(root / mode)],
                                   env=env, timeout=220)
        (root / ('Cubit.ini.' + mode)).write_bytes(ini.read_bytes())
        if completed.returncode:
            raise RuntimeError(f'Automatically registered toolbar {mode} test failed')
    result['passed'] = True
finally:
    if _running_cubit_processes():
        result['restore_blocked'] = 'A Cubit process is still running; backup retained'
    else:
        ini.write_bytes(original)
        result['preferences_restored'] = ini.read_bytes() == original
    result['original_preferences_sha256'] = hashlib.sha256(original).hexdigest()
    result['passed'] = result['passed'] and result['preferences_restored']
    (root / 'result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result))
