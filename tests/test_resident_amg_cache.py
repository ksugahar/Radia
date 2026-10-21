"""CPU-only contracts for explicit resident-AMG caches and selection."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.sparse import csr_matrix, diags

PROBES = Path(__file__).parents[1] / 'validation_test' / 'magnetostatics'


@pytest.fixture
def cache(monkeypatch):
    monkeypatch.syspath_prepend(str(PROBES))
    import mixed_omega_amg_cache
    return mixed_omega_amg_cache


def test_linear_cache_contract_probe(monkeypatch, tmp_path, cache):
    """One executable contract owns linear cache acceptance/rejection signals."""
    spec = importlib.util.spec_from_file_location('cache_contract', PROBES / 'resident_amg_cache_probe.py')
    probe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probe)
    out = tmp_path / 'result.json'
    monkeypatch.setattr(sys, 'argv', ['probe', '--out', str(out)])
    probe.main()
    result = json.loads(out.read_text())
    assert result['passed']
    assert set(result['rejections']) == {'changed_operator', 'changed_structure', 'changed_candidate', 'changed_settings', 'overwrite', 'corrupted_payload'}


def test_nonlinear_cache_transfers_only_and_invalidation(cache, tmp_path):
    import h5py
    p = csr_matrix((np.ones(9), (np.arange(9), np.arange(9)//3)), shape=(9, 3))
    a = diags([-np.ones(8), 4*np.ones(9), -np.ones(8)], [-1, 0, 1]).tocsr()
    b = np.ones((9, 1))
    parameters = {'damping': .8}
    layout = cache.nonlinear_layout(a.indptr, a.indices)
    path = tmp_path / 'transfers.h5'
    cache.save_nonlinear_transfers(path, [p], b, parameters=parameters, layout=layout)
    loaded, receipt = cache.load_nonlinear_transfers(path, b, parameters=parameters, layout=layout)
    np.testing.assert_array_equal(loaded[0].toarray(), p.toarray())
    assert receipt['transfers_only'] and receipt['galerkin_operators_rebuilt']
    with h5py.File(path, 'r') as f:
        assert set(f) == {'candidate', 'p0'}  # No stale operator values stored.
    for candidate, settings, structure in [(b*2, parameters, layout), (b, {}, layout), (b, parameters, 'wrong')]:
        with pytest.raises(ValueError):
            cache.load_nonlinear_transfers(path, candidate, parameters=settings, layout=structure)
    with pytest.raises(FileExistsError):
        cache.save_nonlinear_transfers(path, [p], b, parameters=parameters, layout=layout)
    with h5py.File(path, 'r+') as f:
        f['p0/val'][0] += .1
    with pytest.raises(ValueError, match='payload'):
        cache.load_nonlinear_transfers(path, b, parameters=parameters, layout=layout)


@pytest.mark.parametrize('settings', [
    {'nonlinear_cache_mode': 'reuse'},
    {'nonlinear_cache_mode': 'unknown'},
    {'cache_mode': 'reuse'},
    {'cache_mode': 'unknown'},
    {'cache_path': 'unused', 'updated_transfers': []},
])
def test_invalid_selection_fails_before_optional_runtime_import(cache, monkeypatch, settings):
    from mixed_omega_amg_resident import ResidentAMG
    import builtins
    original = builtins.__import__
    def guarded(name, *args, **kwargs):
        if name in {'pyamg', 'ngsolve', 'radia_simulink', 'matlab.engine', '_resident_amg'}:
            raise AssertionError('Optional runtime imported before selection validation')
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', guarded)
    with pytest.raises(ValueError):
        ResidentAMG(None, None, None, directory='unused', candidate=None, **settings)

@pytest.mark.parametrize('bad_solution', [np.array([0., 0.]), np.array([np.nan, 0.])])
def test_host_original_residual_rejects_native_false_pass(cache, bad_solution):
    from mixed_omega_amg_resident import ResidentAMG
    class Native:
        def solve(self, rhs, tolerance, limit):
            return bad_solution, {'true_relative_residual': 0., 'certificate_gate': 1e-11}
    solver = ResidentAMG.__new__(ResidentAMG)
    solver.operator = diags([2., 4.]).tocsr()
    solver.native = Native()
    with pytest.raises(RuntimeError, match='Independent host original residual'):
        solver.solve(np.ones(2))


@pytest.mark.parametrize('rhs', [np.array([2., 8.]), np.zeros(2)])
def test_host_original_residual_records_independent_certificate(cache, rhs):
    from mixed_omega_amg_resident import ResidentAMG
    class Native:
        def solve(self, values, tolerance, limit):
            return values / np.array([2., 4.]), {'certificate_gate': 1e-11}
    solver = ResidentAMG.__new__(ResidentAMG)
    solver.operator = diags([2., 4.]).tocsr()
    solver.native = Native()
    solution, info = solver.solve(rhs)
    assert info['host_true_relative_residual'] == 0.
    assert info['host_certificate_gate'] == 1e-6
    assert info['certificate_gate'] == 1e-11
    assert info['host_residual_seconds'] >= 0.
