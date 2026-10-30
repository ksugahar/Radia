"""Guards: release receipts/target/gates. Assumes: isolated synthetic artifacts."""
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
    assert dual.TARGETS == ('100',)
    # AGENTS.md: validation runs verified wheels, 100 keeps the editable install.
    assert dual.TIERS == {'100': 'editable'}
    assert "'-e', str(package)" in dual.WORKER
    assert "'--force-reinstall'" not in dual.WORKER
    assert 'pip\', \'uninstall' not in dual.WORKER
    # Cubit is closed only on explicit request, and then everything closed
    # is recorded in the receipt.
    head, forced = dual.WORKER.split("if cfg.get('force_close_cubit')", 1)
    assert 'taskkill' not in head
    assert "result['force_closed_cubit']" in forced


def test_done_reports_only_the_build_target(tmp_path, monkeypatch, capsys):
    import json
    from types import SimpleNamespace
    root = tmp_path / 'source'
    metadata = root / 'packages/cubit-mesh-export/pyproject.toml'
    metadata.parent.mkdir(parents=True)
    metadata.write_text('[project]\ndependencies=["netgen-mesher==6.2.2608","ngsolve==6.2.2608"]\n')
    evidence = tmp_path / 'evidence'
    (evidence / '100').mkdir(parents=True)
    contract = dict(version='2.2.8', source_sha='a'*40, wheel_sha256='b'*64)
    receipt = dict(contract, schema=dual.SCHEMA, target='100', tier='editable',
        profile_scope=['--all-users'], passed=True, unrelated_packages_unchanged=True,
        smoke_test=True, toolbar_smoke=True, mcp_selftest=True, mcp_cli_selftest=True)
    (evidence/'100/deploy.json').write_text(json.dumps(receipt))
    monkeypatch.setattr(dual, 'wheel_contract', lambda path: dict(contract))
    monkeypatch.setattr(dual, 'verify_published', lambda name, record: None)
    def process(argv, **kwargs):
        output = contract['source_sha']+'\n' if argv[0]=='git' else 'CUBIT_DUAL_RESULT={"passed":true}\n'
        return SimpleNamespace(returncode=0,stdout=output,stderr='')
    monkeypatch.setattr(dual.subprocess, 'run', process)
    args = SimpleNamespace(action='done',run_from='100',source_root_100=str(root),
        evidence_100=str(evidence),wheel=tmp_path/'synthetic.whl',source_sha=contract['source_sha'])
    assert dual.run(args)==0
    report=json.loads((evidence/'done.json').read_text())
    assert report['targets']==['100'] and report['passed'] is True, report
    output=capsys.readouterr().out
    assert 'PASS cubit-mesh-export:' in output, output
    assert 'PASS release-dual:' not in output, output


def test_force_and_run_from_are_explicit_options():
    args = dual.build_parser().parse_args([
        '--action', 'deploy', '--wheel', 'w.whl', '--source-sha', 'a' * 40, '--source-root-100', r'W:\\x', '--evidence-100', r'C:\\temp\\e'])
    assert args.run_from == '100' and args.force_close_cubit is False
    args = dual.build_parser().parse_args([
        '--action', 'deploy', '--wheel', 'w.whl', '--source-sha', 'a' * 40, '--source-root-100', r'W:\\x', '--evidence-100', r'C:\\temp\\e',
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
        '--source-root-100', r'W:\\synthetic-repo\\main-checkout',
        '--evidence-100', r'C:\\temp\\cubit-dual',
    ])
    assert args.action == 'preflight'
    assert args.wheel == 'candidate.whl'


def test_receipt_requires_every_acceptance_field():
    contract = dict(version='1.0.2', source_sha='a' * 40, wheel_sha256='b' * 64)
    receipt = dict(contract, schema=dual.SCHEMA, target='100', passed=True,
                   unrelated_packages_unchanged=True, smoke_test=True, toolbar_smoke=True,
                   mcp_selftest=True, mcp_cli_selftest=True, profile_scope=['--all-users'],
                   tier='editable')
    assert dual.check_receipt(receipt, contract, '100')
    for key in receipt:
        damaged = {k: v for k, v in receipt.items() if k != key}
        assert not dual.check_receipt(damaged, contract, '100'), key
    assert not dual.check_receipt(receipt, contract, 'validation')
    assert not dual.check_receipt(receipt, dict(contract, wheel_sha256='c' * 64), '100')


def test_receipt_is_bound_to_the_target_tier():
    contract = dict(version='1.0.2', source_sha='a' * 40, wheel_sha256='b' * 64)
    common = dict(contract, schema=dual.SCHEMA, passed=True,
                  unrelated_packages_unchanged=True, smoke_test=True, toolbar_smoke=True,
                  mcp_selftest=True, mcp_cli_selftest=True)
    # An editable validation receipt (the pre-v4 behaviour) no longer satisfies done.
    lab_editable = dict(common, target='validation', profile_scope=[], tier='editable',
                        installed_payload_verified=True)
    assert not dual.check_receipt(lab_editable, contract, 'validation')
    student_wheel = dict(common, target='100', profile_scope=['--all-users'], tier='wheel')
    assert not dual.check_receipt(student_wheel, contract, '100')




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
    import os
    import subprocess
    encoded = base64.b64encode(json.dumps(cfg).encode()).decode()
    p = subprocess.run([sys.executable, '-', encoded], input=_mock_process_worker(), text=True,
                       capture_output=True, timeout=120, env=os.environ.copy())
    lines = [line for line in p.stdout.splitlines() if line.startswith('CUBIT_DUAL_RESULT=')]
    assert lines, p.stderr
    return p.returncode, json.loads(lines[-1].partition('=')[2])




def _pins_or_skip():
    from importlib import metadata
    try:
        return {name: metadata.version(name) for name in ('netgen-mesher', 'ngsolve')}
    except metadata.PackageNotFoundError:
        pytest.skip('worker dependency pins need netgen-mesher and ngsolve installed')




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


@pytest.mark.parametrize('tier', ['editable'])
def test_preflight_installer_imports_selected_helpers_before_stale_package(tmp_path, monkeypatch, tier):
    """Exercise the worker in a fresh process with a hostile old PYTHONPATH."""
    import shutil
    import subprocess
    git = shutil.which('git')
    assert git is not None, 'Git is required to create the editable-source fixture'
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
            subprocess.run([git, '-C', str(source), *args], check=True, capture_output=True)
        sha = subprocess.check_output([git, '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    cfg = dict(dual.wheel_contract(wheel), schema=dual.SCHEMA, source_sha=sha,
               target='100' if tier == 'editable' else 'validation', tier=tier,
               action='preflight', source_root=str(source) if tier == 'editable' else '',
               output=str(tmp_path / 'evidence'), wheel_path=str(wheel),
               dependencies=_pins_or_skip(), force_close_cubit=False)
    code, result = _run_worker(cfg)
    assert code == 0 and result['passed'], result.get('error')











@pytest.mark.parametrize("target", ["validation", "compute"])
def test_compute_deployment_rejected_before_wheel_access(target):
    from types import SimpleNamespace
    with pytest.raises(ValueError, match="forbidden"):
        dual.run(SimpleNamespace(run_from=target))
