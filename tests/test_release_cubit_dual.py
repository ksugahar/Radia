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
    assert "'-e', str(package)" in dual.WORKER
    assert 'pip\', \'uninstall' not in dual.WORKER
    assert 'taskkill' not in dual.WORKER


def test_dedicated_cli_owns_release_dual_entrypoint():
    args = dual.build_parser().parse_args([
        '--action', 'preflight',
        '--wheel', 'candidate.whl',
        '--source-sha', 'a' * 40,
        '--source-root-lab', r'S:\\Radia\\01_GitHub',
        '--source-root-100', r'W:\\00_CAE\\Radia\\01_GitHub',
        '--evidence-lab', r'C:\\temp\\cubit-dual',
        '--evidence-100', r'C:\\temp\\cubit-dual',
    ])
    assert args.action == 'preflight'
    assert args.wheel == 'candidate.whl'


def test_receipt_requires_every_acceptance_field():
    contract = dict(version='1.0.2', source_sha='a' * 40, wheel_sha256='b' * 64)
    receipt = dict(contract, schema=dual.SCHEMA, target='lab', passed=True,
                   unrelated_packages_unchanged=True, smoke_test=True, toolbar_smoke=True,
                   mcp_selftest=True, mcp_cli_selftest=True)
    assert dual.check_receipt(receipt, contract, 'lab')
    for key in receipt:
        damaged = {k: v for k, v in receipt.items() if k != key}
        assert not dual.check_receipt(damaged, contract, 'lab'), key
    assert not dual.check_receipt(receipt, contract, '100')
    assert not dual.check_receipt(receipt, dict(contract, wheel_sha256='c' * 64), 'lab')


def test_wheel_contract_includes_embedded_probe_and_bytes(tmp_path):
    wheel = tmp_path / 'candidate.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('cubit_mesh_export-1.0.2.dist-info/METADATA',
                         'Name: cubit-mesh-export\nVersion: 1.0.2\n')
        for file in ('cubit_mesh_export.ccm', 'cubit_mesh_curver.pyd',
                     'toolbar_smoke.py', 'cubit_gui/toolbar_probe.py', 'mcp/server.py',
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
    monkeypatch.setattr(install, '_check_radia_compat', forbidden)
    monkeypatch.setattr(install, '_running_cubit_processes', lambda: [])
    monkeypatch.setattr(install, '_critical_plugin_files', lambda root: [])
    assert install.preflight(tmp_path, verbose=False) == (True, [])


def test_optional_integration_check_is_explicit(monkeypatch):
    from cubit_mesh_export import install
    monkeypatch.setattr(sys, 'argv', ['cubit-plugin-install', '--check-radia-compat'])
    monkeypatch.setattr(install, '_check_radia_compat', lambda: (False, 'integration not accepted'))
    with pytest.raises(SystemExit) as error:
        install.main()
    assert error.value.code == 4
