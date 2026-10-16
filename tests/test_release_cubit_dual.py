"""Independent release receipts must not waive bytes, targets or GUI gates."""
import importlib.util
from pathlib import Path
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('dual_test', ROOT / 'tools/release_cubit_dual.py')
dual = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dual)


def test_worker_compiles_and_targets_only_dual():
    compile(dual.WORKER, '<worker>', 'exec')
    assert dual.TARGETS == ('lab', '100')
    # AGENTS.md: LAB runs verified wheels, 100 keeps the editable install.
    assert dual.TIERS == {'lab': 'wheel', '100': 'editable'}
    assert "'-e', str(package)" in dual.WORKER
    assert "'--force-reinstall'" in dual.WORKER
    assert 'pip\', \'uninstall' not in dual.WORKER
    # Cubit is closed only on explicit request, and then everything closed
    # is recorded in the receipt.
    head, forced = dual.WORKER.split("if cfg.get('force_close_cubit')", 1)
    assert 'taskkill' not in head
    assert "result['force_closed_cubit']" in forced


def test_force_and_run_from_are_explicit_options():
    args = dual.build_parser().parse_args([
        '--action', 'deploy', '--wheel', 'w.whl', '--source-sha', 'a' * 40,
        '--source-root-lab', r'C:\\release-dual\\x', '--source-root-100', r'W:\\x',
        '--evidence-lab', r'C:\\temp\\e', '--evidence-100', r'C:\\temp\\e'])
    assert args.run_from == 'lab' and args.force_close_cubit is False
    args = dual.build_parser().parse_args([
        '--action', 'deploy', '--wheel', 'w.whl', '--source-sha', 'a' * 40,
        '--source-root-lab', r'C:\\release-dual\\x', '--source-root-100', r'W:\\x',
        '--evidence-lab', r'C:\\temp\\e', '--evidence-100', r'C:\\temp\\e',
        '--run-from', '100', '--force-close-cubit'])
    assert args.run_from == '100' and args.force_close_cubit is True


def test_student_host_registers_every_profile():
    # 100 serves many student accounts; registering only the deploying
    # Administrator left every other profile loading a removed checkout.
    assert "profile_scope = ['--all-users'] if cfg['target'] == '100' else []" in dual.WORKER
    assert "'cubit_mesh_export.install', *profile_scope]" in dual.WORKER
    assert "'--verify-only', *profile_scope]" in dual.WORKER
    contract = dict(version='1.0.2', source_sha='a' * 40, wheel_sha256='b' * 64)
    receipt = dict(contract, schema=dual.SCHEMA, target='100', passed=True,
                   unrelated_packages_unchanged=True, smoke_test=True, toolbar_smoke=True,
                   mcp_selftest=True, mcp_cli_selftest=True, profile_scope=['--all-users'],
                   tier='editable')
    assert dual.check_receipt(receipt, contract, '100')
    assert not dual.check_receipt(dict(receipt, profile_scope=[]), contract, '100')


def test_dedicated_cli_owns_release_dual_entrypoint():
    args = dual.build_parser().parse_args([
        '--action', 'preflight',
        '--wheel', 'candidate.whl',
        '--source-sha', 'a' * 40,
        '--source-root-lab', r'S:\\synthetic-repo\\main-checkout',
        '--source-root-100', r'W:\\synthetic-repo\\main-checkout',
        '--evidence-lab', r'C:\\temp\\cubit-dual',
        '--evidence-100', r'C:\\temp\\cubit-dual',
    ])
    assert args.action == 'preflight'
    assert args.wheel == 'candidate.whl'


def test_receipt_requires_every_acceptance_field():
    contract = dict(version='1.0.2', source_sha='a' * 40, wheel_sha256='b' * 64)
    receipt = dict(contract, schema=dual.SCHEMA, target='lab', passed=True,
                   unrelated_packages_unchanged=True, smoke_test=True, toolbar_smoke=True,
                   mcp_selftest=True, mcp_cli_selftest=True, profile_scope=[],
                   tier='wheel', installed_payload_verified=True)
    assert dual.check_receipt(receipt, contract, 'lab')
    for key in receipt:
        damaged = {k: v for k, v in receipt.items() if k != key}
        assert not dual.check_receipt(damaged, contract, 'lab'), key
    assert not dual.check_receipt(receipt, contract, '100')
    assert not dual.check_receipt(receipt, dict(contract, wheel_sha256='c' * 64), 'lab')


def test_receipt_is_bound_to_the_target_tier():
    contract = dict(version='1.0.2', source_sha='a' * 40, wheel_sha256='b' * 64)
    common = dict(contract, schema=dual.SCHEMA, passed=True,
                  unrelated_packages_unchanged=True, smoke_test=True, toolbar_smoke=True,
                  mcp_selftest=True, mcp_cli_selftest=True)
    # An editable LAB receipt (the pre-v4 behaviour) no longer satisfies done.
    lab_editable = dict(common, target='lab', profile_scope=[], tier='editable',
                        installed_payload_verified=True)
    assert not dual.check_receipt(lab_editable, contract, 'lab')
    student_wheel = dict(common, target='100', profile_scope=['--all-users'], tier='wheel')
    assert not dual.check_receipt(student_wheel, contract, '100')


def test_wheel_transport_is_local_or_hash_rechecked_remote_copy(monkeypatch, tmp_path):
    wheel = tmp_path / 'cubit_mesh_export-1.0.2-cp312-cp312-win_amd64.whl'
    wheel.write_bytes(b'wheel')
    assert dual.transport_wheel(wheel, 'lab', 'lab', r'C:\temp\e') == str(wheel.resolve())
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs.get('input')))
        return None

    monkeypatch.setattr(dual.subprocess, 'run', fake_run)
    destination = dual.transport_wheel(wheel, 'lab', '100', r'C:\temp\e')
    assert destination == r'C:\temp\e\lab\wheel' + '\\' + wheel.name
    (mkdir, script), (copy, _) = calls
    assert mkdir[:5] == ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10']
    assert mkdir[5:8] == ['lab', 'python', '-'] and 'mkdir' in script
    assert copy[-1] == 'lab:' + destination.replace('\\', '/')


def _fake_cubit_wheel(path):
    contract_files = ('cubit_mesh_export.ccm', 'cubit_mesh_curver.pyd',
                      'toolbar_smoke.py', 'cubit_gui/toolbar_probe.py',
                      'cubit_gui/cubit_export_menu.py',
                      'cubit_gui/cubit_toolbar/toolbars/cubit_mesh_export_toolbar.ttb.tmpl',
                      'cubit_gui/solver_ready_sample.jou', 'native_payloads.json',
                      'mcp/server.py', 'mcp/_support/status.py')
    installer = ('from pathlib import Path\n'
                 'def _find_cubit_dir():\n    return Path(".")\n'
                 'def preflight(cubit_dir, verbose=False):\n    return True, []\n')
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('cubit_mesh_export-1.0.2.dist-info/METADATA',
                         'Name: cubit-mesh-export\nVersion: 1.0.2\n')
        for file in contract_files:
            archive.writestr('cubit_mesh_export/' + file, b'payload\n')
        archive.writestr('cubit_mesh_export/install.py', installer)


def _run_worker(cfg):
    import base64
    import json
    import subprocess
    encoded = base64.b64encode(json.dumps(cfg).encode()).decode()
    p = subprocess.run([sys.executable, '-', encoded], input=_mock_process_worker(), text=True,
                       capture_output=True, timeout=120)
    lines = [line for line in p.stdout.splitlines() if line.startswith('CUBIT_DUAL_RESULT=')]
    assert lines, p.stderr
    return p.returncode, json.loads(lines[-1].partition('=')[2])


@pytest.mark.parametrize('tamper', [False, True])
def test_wheel_tier_preflight_reads_only_hash_checked_wheel_bytes(tmp_path, tamper):
    from importlib import metadata
    try:
        pins = {name: metadata.version(name) for name in ('netgen-mesher', 'ngsolve')}
    except metadata.PackageNotFoundError:
        pytest.skip('worker dependency pins need netgen-mesher and ngsolve installed')
    wheel = tmp_path / 'cubit_mesh_export-1.0.2-cp312-cp312-win_amd64.whl'
    _fake_cubit_wheel(wheel)
    contract = dual.wheel_contract(wheel)
    if tamper:
        contract['wheel_sha256'] = '0' * 64
    out = tmp_path / 'evidence' / 'lab'
    cfg = dict(contract, schema=dual.SCHEMA, source_sha='a' * 40, target='lab',
               action='preflight', tier='wheel', source_root='', output=str(out),
               wheel_path=str(wheel), dependencies=pins, force_close_cubit=False)
    # A pre-existing tree under the evidence directory is never removed.
    foreign = out / 'wheel_payload'
    foreign.mkdir(parents=True)
    (foreign / 'keep.txt').write_text('not ours')
    returncode, result = _run_worker(cfg)
    assert result['tier'] == 'wheel'
    assert (foreign / 'keep.txt').read_text() == 'not ours'
    assert not list(out.glob('cubit-dual-wheel-*'))
    if tamper:
        assert returncode == 1 and result['passed'] is False
        assert 'Transported wheel differs' in result['error']
    else:
        assert returncode == 0 and result['passed'] is True, result.get('error')
        assert result['staging_removed'] is True


def _pins_or_skip():
    from importlib import metadata
    try:
        return {name: metadata.version(name) for name in ('netgen-mesher', 'ngsolve')}
    except metadata.PackageNotFoundError:
        pytest.skip('worker dependency pins need netgen-mesher and ngsolve installed')


def test_wheel_tier_rejects_a_member_that_leaves_the_staging_tree(tmp_path):
    pins = _pins_or_skip()
    wheel = tmp_path / 'cubit_mesh_export-1.0.2-cp312-cp312-win_amd64.whl'
    _fake_cubit_wheel(wheel)
    contract = dual.wheel_contract(wheel)
    with zipfile.ZipFile(wheel, 'a') as archive:
        archive.writestr('../escape.txt', b'x')
    contract['wheel_sha256'] = dual.digest(wheel.read_bytes())
    out = tmp_path / 'evidence' / 'lab'
    cfg = dict(contract, schema=dual.SCHEMA, source_sha='a' * 40, target='lab',
               action='preflight', tier='wheel', source_root='', output=str(out),
               wheel_path=str(wheel), dependencies=pins, force_close_cubit=False)
    returncode, result = _run_worker(cfg)
    assert returncode == 1 and 'Unsafe wheel member path' in result['error']
    # The owned directory was recorded before extraction failed, so it is
    # cleaned up and reported without replacing the original error.
    assert result['staging'] and result['staging_removed'] is True
    assert 'staging_cleanup_error' not in result
    assert not (tmp_path / 'evidence' / 'escape.txt').exists()
    assert not list(out.glob('cubit-dual-wheel-*'))


def _installable_fake_wheel(path):
    """A pip-installable wheel with the contract files (no Cubit, no MCP)."""
    import base64
    payload = {'cubit_mesh_export/__init__.py': b'__version__ = "1.0.2"\n',
               'cubit_mesh_export/install.py': (
                   b'from pathlib import Path\n'
                   b'def _find_cubit_dir():\n    return Path(".")\n'
                   b'def preflight(cubit_dir, verbose=False):\n    return True, []\n')}
    for file in ('cubit_mesh_export.ccm', 'cubit_mesh_curver.pyd', 'toolbar_smoke.py',
                 'cubit_gui/toolbar_probe.py', 'cubit_gui/cubit_export_menu.py',
                 'cubit_gui/cubit_toolbar/toolbars/cubit_mesh_export_toolbar.ttb.tmpl',
                 'cubit_gui/solver_ready_sample.jou', 'native_payloads.json',
                 'mcp/server.py', 'mcp/_support/status.py'):
        payload['cubit_mesh_export/' + file] = b'payload\n'
    info = 'cubit_mesh_export-1.0.2.dist-info/'
    payload[info + 'METADATA'] = b'Metadata-Version: 2.1\nName: cubit-mesh-export\nVersion: 1.0.2\n'
    payload[info + 'WHEEL'] = (b'Wheel-Version: 1.0\nGenerator: test\nRoot-Is-Purelib: false\n'
                               b'Tag: py3-none-any\n')
    record = []
    for name, data in payload.items():
        digest = base64.urlsafe_b64encode(dual.hashlib.sha256(data).digest()).rstrip(b'=')
        record.append(f'{name},sha256={digest.decode()},{len(data)}')
    record.append(info + 'RECORD,,')
    payload[info + 'RECORD'] = ('\n'.join(record) + '\n').encode()
    with zipfile.ZipFile(path, 'w') as archive:
        for name, data in payload.items():
            archive.writestr(name, data)


def test_wheel_tier_verify_accepts_only_the_exact_installed_wheel(tmp_path):
    """Install the wheel into an isolated venv and run the worker's verify phase."""
    import os
    import subprocess
    # The deployment worker runs in the target's own interpreter environment.
    # Test runners may export PYTHONPATH (the cubit-mesh-export workflow sets
    # it to the package source), which would put the real source package in
    # front of the installed wheel, so the isolated venv gets an explicit
    # environment without interpreter path overrides.
    isolated = {k: v for k, v in os.environ.items()
                if k not in ('PYTHONPATH', 'PYTHONHOME', 'PYTHONSTARTUP', 'PYTHONUSERBASE',
                             'PYTHONSAFEPATH', 'PYTHONNOUSERSITE')}
    venv = tmp_path / 'venv'
    subprocess.run([sys.executable, '-m', 'venv', '--system-site-packages', str(venv)],
                   check=True, timeout=300, env=isolated)
    python = venv / 'Scripts' / 'python.exe'
    wheel = tmp_path / 'cubit_mesh_export-1.0.2-py3-none-any.whl'
    _installable_fake_wheel(wheel)
    subprocess.run([str(python), '-m', 'pip', 'install', '--no-deps', '--no-index', str(wheel)],
                   check=True, capture_output=True, timeout=300, env=isolated)
    import base64
    import json
    contract = dual.wheel_contract(wheel)
    pins = {}
    for name in ('netgen-mesher', 'ngsolve'):
        probe = subprocess.run([str(python), '-c',
                                'import importlib.metadata as m,sys;print(m.version(sys.argv[1]))', name],
                               capture_output=True, text=True, env=isolated)
        if probe.returncode:
            pytest.skip('worker dependency pins need netgen-mesher and ngsolve installed')
        pins[name] = probe.stdout.strip()

    def verify(installed_wheel_name, env):
        out = tmp_path / ('evidence-' + installed_wheel_name)
        cfg = dict(contract, schema=dual.SCHEMA, source_sha='a' * 40, target='lab',
                   action='verify', tier='wheel', source_root='', output=str(out),
                   wheel_path=str(wheel), dependencies=pins, force_close_cubit=False)
        encoded = base64.b64encode(json.dumps(cfg).encode()).decode()
        p = subprocess.run([str(python), '-', encoded], input=_mock_process_worker(), text=True,
                           capture_output=True, timeout=300, env=env)
        line = [l for l in p.stdout.splitlines() if l.startswith('CUBIT_DUAL_RESULT=')][-1]
        return json.loads(line.partition('=')[2])

    # The exact wheel passes identity and payload verification; the run then
    # stops at the MCP self-test that this synthetic wheel cannot provide.
    result = verify('exact', isolated)
    assert result.get('installed_payload_verified') is True, result.get('error')
    assert result['installed']['direct_url']['archive_info'] is not None
    assert 'mcp.server' in result['error']

    # A source tree shadowing the installed wheel (the CI PYTHONPATH case) is
    # refused by the fresh-interpreter identity check, not reported as the wheel.
    shadow = dict(isolated, PYTHONPATH=str(ROOT / 'packages' / 'cubit-mesh-export' / 'src'))
    result = verify('shadowed', shadow)
    assert 'installed_payload_verified' not in result
    assert 'Wrong installed version' in result.get('error', ''), result.get('error')

    # A modified installed file is caught by the payload check.
    site = venv / 'Lib' / 'site-packages' / 'cubit_mesh_export'
    (site / 'toolbar_smoke.py').write_text('changed\n')
    result = verify('modified', isolated)
    assert 'installed_payload_verified' not in result
    assert 'Published wheel/source mismatch: cubit_mesh_export/toolbar_smoke.py' in result['error']


def test_wheel_contract_includes_embedded_probe_and_bytes(tmp_path):
    wheel = tmp_path / 'candidate.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('cubit_mesh_export-1.0.2.dist-info/METADATA',
                         'Name: cubit-mesh-export\nVersion: 1.0.2\n')
        for file in ('cubit_mesh_export.ccm', 'cubit_mesh_curver.pyd',
                     'toolbar_smoke.py', 'cubit_gui/toolbar_probe.py',
                     'cubit_gui/cubit_export_menu.py',
                     'cubit_gui/cubit_toolbar/toolbars/cubit_mesh_export_toolbar.ttb.tmpl',
                     'cubit_gui/solver_ready_sample.jou', 'native_payloads.json',
                     'mcp/server.py',
                     'mcp/_support/status.py', 'mcp/_support/LICENSE-BSD-3-Clause.txt'):
            archive.writestr('cubit_mesh_export/' + file, b'test\r\n')
    result = dual.wheel_contract(wheel)
    assert result['version'] == '1.0.2'
    assert result['wheel_sha256'] == dual.digest(wheel.read_bytes())
    assert result['files']['cubit_mesh_export/toolbar_smoke.py']['sha256'] == dual.digest(b'test\n')
    license_file = result['files']['cubit_mesh_export/mcp/_support/LICENSE-BSD-3-Clause.txt']
    assert license_file == {'sha256': dual.digest(b'test\n'), 'text': True}
    assert result['files']['cubit_mesh_export/cubit_mesh_export.ccm']['sha256'] == dual.digest(b'test\r\n')


def test_wheel_contract_rejects_other_distribution(tmp_path):
    wheel = tmp_path / 'core.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('unrelated-1.0.dist-info/METADATA',
                         'Name: unrelated\nVersion: 1.0\n')
    with pytest.raises(ValueError, match='Not a cubit'):
        dual.wheel_contract(wheel)


def test_standalone_preflight_never_checks_or_imports_radia(monkeypatch, tmp_path):
    sys.path.insert(0, str(ROOT / 'packages/cubit-mesh-export/src'))
    from cubit_mesh_export import install
    def forbidden():
        raise AssertionError('Standalone deployment consulted Radia')
    assert not hasattr(install, '_check_radia_compat')
    monkeypatch.setattr(install, '_running_cubit_processes', lambda: [])
    monkeypatch.setattr(install, '_critical_plugin_files', lambda root: [])
    assert install.preflight(tmp_path, verbose=False) == (True, [])


def test_standalone_installer_rejects_reverse_radia_integration_flag(monkeypatch):
    from cubit_mesh_export import install
    monkeypatch.setattr(sys, 'argv', ['cubit-plugin-install', '--check-radia-compat'])
    with pytest.raises(SystemExit) as error:
        install.main()
    assert error.value.code == 2


def _mock_process_worker():
    # Mock only process commands; wheel installs/identity probes keep their real isolation.
    return """import subprocess
_actual_run = subprocess.run
def _process_mock(args, **kwargs):
    if args[0] == 'tasklist':
        return subprocess.CompletedProcess(args, 0, '', '')
    if args[0] == 'taskkill':
        raise AssertionError('Unexpected process termination in fixture test')
    return _actual_run(args, **kwargs)
subprocess.run = _process_mock
""" + dual.WORKER


def _process_guard(monkeypatch, images, force, fail_pid=None):
    """Execute the actual worker guard with synthetic tasklist/taskkill only."""
    import ast
    import csv
    import io
    import subprocess
    from types import SimpleNamespace
    active = [(name, 100 + i) for i, name in enumerate(images)]
    calls = []
    def fake_run(args, **kwargs):
        calls.append(args)
        if args[0] == 'tasklist':
            out = io.StringIO()
            writer = csv.writer(out)
            for name, pid in active:
                writer.writerow([name, pid, 'Console', 1, '10 K', 'Running', 'owner', '0:00', 'Cubit'])
            return subprocess.CompletedProcess(args, 0, out.getvalue(), '')
        assert args[:3] == ['taskkill', '/F', '/PID']
        pid = int(args[3])
        if pid == fail_pid:
            return subprocess.CompletedProcess(args, 1, '', 'termination denied')
        active[:] = [(name, value) for name, value in active if value != pid]
        return subprocess.CompletedProcess(args, 0, '', '')
    monkeypatch.setattr(subprocess, 'run', fake_run)
    tree = ast.parse(dual.WORKER)
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                   and node.name in ('command', 'running_cubit_processes')]
    ns = dict(CUBIT_PROCESS_IMAGES=dual.CUBIT_PROCESS_IMAGES, events=[], subprocess=subprocess,
              cfg=dict(force_close_cubit=force, action='preflight'), result={}, cubit_dir=Path('.'),
              installer=SimpleNamespace(preflight=lambda *a, **k: (True, [])))
    exec(compile(ast.Module(body=definitions, type_ignores=[]), '<worker-functions>', 'exec'), ns)
    block = next(node for node in tree.body if isinstance(node, ast.Try))
    start = next(i for i, node in enumerate(block.body) if isinstance(node, ast.If)
                 and "cfg.get('force_close_cubit')" in ast.unparse(node.test))
    # Guard ends immediately before the deployment/verification branches.
    end = next(i for i in range(start + 1, len(block.body))
               if isinstance(block.body[i], ast.If)
               and ast.unparse(block.body[i].test) == "cfg['action'] == 'deploy'")
    guard = ast.Module(body=block.body[start:end], type_ignores=[])
    error = None
    try:
        exec(compile(guard, '<worker-guard>', 'exec'), ns)
    except RuntimeError as exc:
        error = str(exc)
    return ns, calls, error


@pytest.mark.parametrize('image', ['coreform_cubit.exe', 'clarox.exe', 'cubitx.exe'])
def test_worker_preflight_refuses_each_cubit_image(monkeypatch, image):
    _, calls, error = _process_guard(monkeypatch, [image], False)
    assert error and image in error
    assert not any(call[0] == 'taskkill' for call in calls)


def test_force_close_targets_exact_cubit_images_and_records_receipt(monkeypatch):
    images = ['coreform_cubit.exe', 'CLAROX.EXE', 'cubitx.exe',
              'python.exe', 'cubit-plugin-install.exe', 'not-clarox.exe']
    ns, calls, error = _process_guard(monkeypatch, images, True)
    assert error is None
    assert [call for call in calls if call[0] == 'taskkill'] == [
        ['taskkill', '/F', '/PID', str(pid)] for pid in (100, 101, 102)]
    assert [p['image'] for p in ns['result']['force_closed_cubit']] == images[:3]
    assert all(p['user'] == 'owner' and p['session'] == 'Console'
               for p in ns['result']['force_closed_cubit'])
    sys.path.insert(0, str(ROOT / 'packages/cubit-mesh-export/src'))
    from cubit_mesh_export import install
    assert set(dual.CUBIT_PROCESS_IMAGES) == {'coreform_cubit.exe', 'clarox.exe', 'cubitx.exe'}
    assert all(install._is_cubit_process_name(name) for name in dual.CUBIT_PROCESS_IMAGES)
    assert all(not install._is_cubit_process_name(name) for name in images[3:])


def test_force_close_failure_stops_and_retains_completed_closures(monkeypatch):
    ns, calls, error = _process_guard(monkeypatch, list(dual.CUBIT_PROCESS_IMAGES),
                                      True, fail_pid=101)
    assert error and 'termination denied' in error
    assert [p['pid'] for p in ns['result']['force_closed_cubit']] == [100]
    assert ns['result']['force_closed_cubit'][0]['image'] == 'coreform_cubit.exe'
    assert not any(call == ['taskkill', '/F', '/PID', '102'] for call in calls)


@pytest.mark.parametrize('tier', ['wheel', 'editable'])
def test_preflight_installer_imports_selected_helpers_before_stale_package(tmp_path, monkeypatch, tier):
    """Exercise the worker in a fresh process with a hostile old PYTHONPATH."""
    import subprocess
    stale = tmp_path / 'stale' / 'cubit_mesh_export'
    stale.mkdir(parents=True)
    (stale / '__init__.py').write_text('')
    (stale / '_layout.py').write_text('raise RuntimeError("STALE layout imported")\n')
    monkeypatch.setenv('PYTHONPATH', str(stale.parent))
    wheel = tmp_path / 'candidate.whl'
    _fake_cubit_wheel(wheel)
    with zipfile.ZipFile(wheel, 'a') as archive:
        archive.writestr('cubit_mesh_export/__init__.py', '')
        archive.writestr('cubit_mesh_export/_layout.py', 'selected = True\n')
    # Replace the synthetic installer without duplicate zip members.
    with zipfile.ZipFile(wheel) as archive:
        payload = {name: archive.read(name) for name in archive.namelist()}
    payload['cubit_mesh_export/install.py'] = (
        b'from pathlib import Path\nfrom cubit_mesh_export._layout import selected\n'
        b'assert selected\ndef _find_cubit_dir():\n    return Path(".")\n'
        b'def preflight(cubit_dir, verbose=False):\n    return True, []\n')
    with zipfile.ZipFile(wheel, 'w') as archive:
        for name, data in payload.items():
            archive.writestr(name, data)
    source = tmp_path / 'release'
    sha = 'a' * 40
    if tier == 'editable':
        src = source / 'packages/cubit-mesh-export/src'
        src.mkdir(parents=True)
        with zipfile.ZipFile(wheel) as archive:
            archive.extractall(src)
        for args in (['init'], ['add', '.'], ['-c', 'user.name=Test', '-c',
                     'user.email=test@example.invalid', 'commit', '-m', 'fixture']):
            subprocess.run(['git', '-C', str(source), *args], check=True, capture_output=True)
        sha = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    cfg = dict(dual.wheel_contract(wheel), schema=dual.SCHEMA, source_sha=sha,
               target='100' if tier == 'editable' else 'lab', tier=tier,
               action='preflight', source_root=str(source) if tier == 'editable' else '',
               output=str(tmp_path / 'evidence'), wheel_path=str(wheel),
               dependencies=_pins_or_skip(), force_close_cubit=False)
    code, result = _run_worker(cfg)
    assert code == 0 and result['passed'], result.get('error')


def test_configured_ssh_alias_is_used_for_both_transports(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(dual.subprocess, 'run', lambda cmd, **kw: calls.append(cmd))
    wheel = tmp_path / 'candidate.whl'
    wheel.write_bytes(b'wheel')
    dual.transport_wheel(wheel, 'lab', '100', r'C:\temp\e', '102')
    assert calls[0][5] == '102'
    assert calls[1][-1].startswith('102:')
    args = dual.build_parser().parse_args([
        '--action', 'preflight', '--wheel', str(wheel), '--source-sha', 'a' * 40,
        '--source-root-100', 'release', '--evidence-lab', 'e', '--evidence-100', 'e',
        '--ssh-alias-lab', '102'])
    assert args.ssh_alias_lab == '102' and args.ssh_alias_100 == '100'


@pytest.mark.parametrize('failure', ['dns', 'config', 'missing'])
def test_ssh_alias_resolution_fails_clearly(monkeypatch, failure):
    from types import SimpleNamespace
    def fake_run(cmd, **kwargs):
        assert cmd == ['ssh', '-G', 'configured-alias']
        if failure == 'config':
            raise dual.subprocess.CalledProcessError(255, cmd, stderr='Bad configuration')
        return SimpleNamespace(stdout='' if failure == 'missing' else 'hostname remote.invalid\nport 22\n')
    def dns(*args, **kwargs):
        raise dual.socket.gaierror('host not found')
    monkeypatch.setattr(dual.subprocess, 'run', fake_run)
    monkeypatch.setattr(dual.socket, 'getaddrinfo', dns)
    with pytest.raises(RuntimeError, match="lab: SSH alias 'configured-alias'.*--ssh-alias-lab"):
        dual.check_ssh_alias('lab', 'configured-alias')


def test_ssh_alias_checks_effective_hostname_not_alias(monkeypatch):
    from types import SimpleNamespace
    monkeypatch.setattr(dual.subprocess, 'run', lambda *a, **k: SimpleNamespace(
        stdout='hostname remote.example\nport 2222\n'))
    calls = []
    monkeypatch.setattr(dual.socket, 'getaddrinfo', lambda *a, **k: calls.append(a))
    assert dual.check_ssh_alias('lab', '102') == 'remote.example'
    assert calls == [('remote.example', '2222')]

@pytest.mark.parametrize('run_from', ['100', 'lab'])
def test_orchestrator_checks_remote_alias_before_transport_and_uses_it_for_worker(monkeypatch, tmp_path, run_from):
    import json
    from types import SimpleNamespace
    events = []
    remote = next(t for t in dual.TARGETS if t != run_from)
    alias = 'custom-' + remote
    monkeypatch.setattr(dual, 'wheel_contract', lambda _: dict(version='1.0.2', files={}, requires=[]))
    monkeypatch.setattr(dual, 'check_ssh_alias', lambda target, name: events.append(('check', target, name)))
    def transport(wheel, target, local, output, name):
        events.append(('transport', target, name))
        return 'wheel.whl'
    monkeypatch.setattr(dual, 'transport_wheel', transport)
    def worker(cmd, **kwargs):
        events.append(('worker', cmd))
        return SimpleNamespace(returncode=0, stdout='CUBIT_DUAL_RESULT=' + json.dumps({'passed': True}), stderr='')
    monkeypatch.setattr(dual.subprocess, 'run', worker)
    args = dual.build_parser().parse_args([
        '--action', 'preflight', '--wheel', 'wheel.whl', '--source-sha', 'a' * 40,
        '--source-root-100', str(ROOT), '--source-root-lab', str(ROOT),
        '--evidence-lab', str(tmp_path), '--evidence-100', str(tmp_path),
        '--run-from', run_from, '--ssh-alias-' + remote, alias])
    assert dual.run(args) == 0
    assert events[0] == ('check', remote, alias)
    commands = [event[1] for event in events if event[0] == 'worker']
    ssh = next(cmd for cmd in commands if cmd[0] == 'ssh')
    assert ssh[5] == alias
    assert len(commands) == 2


def test_unresolved_alias_aborts_orchestrator_before_transfer_or_workers(monkeypatch):
    monkeypatch.setattr(dual, 'wheel_contract', lambda _: dict(version='1.0.2'))
    def fail(*args):
        raise RuntimeError('unresolved configured SSH alias')
    monkeypatch.setattr(dual, 'check_ssh_alias', fail)
    monkeypatch.setattr(dual.subprocess, 'run', lambda *a, **k: pytest.fail('started worker'))
    monkeypatch.setattr(dual, 'transport_wheel', lambda *a: pytest.fail('transferred wheel'))
    args = dual.build_parser().parse_args([
        '--action', 'deploy', '--wheel', 'wheel.whl', '--source-sha', 'a' * 40,
        '--source-root-100', 'missing', '--evidence-lab', 'e', '--evidence-100', 'e',
        '--run-from', '100'])
    with pytest.raises(RuntimeError, match='unresolved configured SSH alias'):
        dual.run(args)
