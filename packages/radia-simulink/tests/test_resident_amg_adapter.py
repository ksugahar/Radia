"""Language lifecycle boundaries; fake native owner is not GPU acceptance."""
import importlib.util
import hashlib
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import pytest


@pytest.fixture
def adapter(monkeypatch):
    path = Path(__file__).parents[1] / 'src/radia_simulink/resident_amg.py'
    before = set(sys.modules)
    spec = importlib.util.spec_from_file_location('resident_adapter_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert not {'numpy', 'ngsolve', 'matlab.engine', 'radia_simulink._cuda_runtime'} & (set(sys.modules) - before)
    return module


def test_wrong_backend_fails_before_runtime_loading(adapter, monkeypatch):
    def forbidden():
        raise AssertionError('Numerical runtime loaded before backend validation')
    monkeypatch.setattr(adapter, '_qualified_native', forbidden)
    for backend in ('cpu', None, True, 'GPU'):
        with pytest.raises(ValueError):
            adapter.ResidentAMG([], None, backend=backend)


def test_owner_reuse_errors_and_idempotent_close(adapter, monkeypatch):
    class Native:
        def __init__(self, *_):
            self.solves = self.closes = 0
        def solve(self, rhs, tolerance, limit):
            self.solves += 1
            return rhs, {'backend': 'gpu', 'solve_count': self.solves}
        def close(self):
            self.closes += 1
    monkeypatch.setattr(adapter, '_qualified_native', lambda: SimpleNamespace(Resident=Native))
    owner = adapter.ResidentAMG([], None, backend='gpu')
    native = owner._owner
    for tolerance in (True, 0., -1., float('nan'), float('inf'), 1e-5):
        with pytest.raises(ValueError):
            owner.solve([], tolerance=tolerance)
    for limit in (True, 0, 1.5, 100001):
        with pytest.raises(ValueError):
            owner.solve([], max_iterations=limit)
    assert native.solves == 0
    with owner:
        assert owner.solve([1.])[1]['solve_count'] == 1
        assert owner.solve([2.])[1]['solve_count'] == 2
    owner.close()
    assert native.closes == 1
    with pytest.raises(RuntimeError, match='closed'):
        owner.solve([1.])


@pytest.mark.parametrize('mutation', [
    'none', 'source', 'binary', 'dirty', 'fault-missing', 'fault-true',
    'device-fault', 'diagnostic', 'engine', 'runtime', 'redistribution', 'proof',
    'receipt-missing', 'binary-missing',
])
def test_receipt_admission_precedes_cuda_loading(adapter, monkeypatch, tmp_path, mutation):
    package = ModuleType('resident_receipt_contract')
    package.__path__ = []
    package.verify = lambda: {'source_commit': '1' * 40}
    monkeypatch.setitem(sys.modules, package.__name__, package)
    adapter.__package__ = package.__name__
    adapter.__file__ = str(tmp_path / 'resident_amg.py')
    name = '_resident_amg' + adapter.sysconfig.get_config_var('EXT_SUFFIX')
    (tmp_path / name).write_bytes(b'contract-only-not-native')
    receipt = dict(schema='radia.resident-amg-python-build.v1', source_commit='1' * 40,
                   source_dirty=False, binary=name,
                   extension_suffix=adapter.sysconfig.get_config_var('EXT_SUFFIX'),
                   binary_sha256=hashlib.sha256((tmp_path / name).read_bytes()).hexdigest(),
                   fault_injection=False, device_fault_injection=False,
                   lifecycle_diagnostic=False, matlab_engine_user_dependency=False,
                   runtime_import_checked=True, redistribution_accepted=True,
                   runtime_validation_sha256='2' * 64, loaded_runtime_sha256='3' * 64,
                   redistribution_reference_sha256='4' * 64)
    changes = {
        'source': ('source_commit', '5' * 40), 'binary': ('binary_sha256', '0' * 64),
        'dirty': ('source_dirty', True), 'fault-true': ('fault_injection', True),
        'device-fault': ('device_fault_injection', True),
        'diagnostic': ('lifecycle_diagnostic', True),
        'engine': ('matlab_engine_user_dependency', True),
        'runtime': ('runtime_import_checked', False),
        'redistribution': ('redistribution_accepted', False),
        'proof': ('redistribution_reference_sha256', ''),
    }
    if mutation == 'fault-missing':
        receipt.pop('fault_injection')
    elif mutation in changes:
        key, value = changes[mutation]
        receipt[key] = value
    (tmp_path / 'resident-amg-build.json').write_text(json.dumps(receipt))
    if mutation == 'receipt-missing':
        (tmp_path / 'resident-amg-build.json').unlink()
    elif mutation == 'binary-missing':
        (tmp_path / name).unlink()
    calls = []
    runtime = ModuleType(package.__name__ + '._cuda_runtime')
    runtime.prepare_cuda_runtime = lambda: calls.append('runtime')
    monkeypatch.setitem(sys.modules, runtime.__name__, runtime)
    marker = object()
    monkeypatch.setattr(adapter, 'import_module', lambda *_: marker)
    if mutation == 'none':
        assert adapter._qualified_native() is marker
        assert calls == ['runtime']
    else:
        with pytest.raises(ImportError, match='receipt'):
            adapter._qualified_native()
        assert calls == []
