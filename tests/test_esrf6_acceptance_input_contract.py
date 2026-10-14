"""Fail-closed acceptance and restaging checks, without native imports."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'case6_contract', ROOT / 'validation_test/esrf_three_engine/case6_acceptance_contract.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.fixture
def accepted():
    path = ROOT / 'validation_test/esrf_three_engine/results/candidate_4d85e72cc/case6_bdm1_mass_riesz.json'
    return json.loads(path.read_bytes())


def test_original_accepted_record_passes_unchanged(accepted):
    original = copy.deepcopy(accepted)
    MODULE.require_accepted_hdiv(accepted)
    assert accepted == original


@pytest.mark.parametrize('value', [-1e-9, -float('inf'), float('inf'), float('nan'),
                                 True, False, '0', None, 1.0])
def test_invalid_residual_is_rejected(accepted, value):
    accepted['nonlinear_stats']['nonlinear_final_relative_residual'] = value
    with pytest.raises(ValueError):
        MODULE.require_accepted_hdiv(accepted)


@pytest.mark.parametrize('value', [0, -1, float('inf'), -float('inf'), float('nan'),
                                 True, False, '2e-5', None])
def test_invalid_tolerance_is_rejected(accepted, value):
    accepted['options']['nl_tol'] = value
    with pytest.raises(ValueError):
        MODULE.require_accepted_hdiv(accepted)


@pytest.mark.parametrize('name', ['completed', 'accepted', 'nonlinear_converged_final_stage'])
@pytest.mark.parametrize('value', [1, 'true', [], False, None])
def test_flags_require_literal_true(accepted, name, value):
    target = accepted['nonlinear_stats'] if name.startswith('nonlinear_') else accepted
    target[name] = value
    with pytest.raises(ValueError):
        MODULE.require_accepted_hdiv(accepted)


def test_zero_and_boundary_residual_are_valid(accepted):
    for value in (0, accepted['options']['nl_tol']):
        accepted['nonlinear_stats']['nonlinear_final_relative_residual'] = value
        MODULE.require_accepted_hdiv(accepted)


def test_missing_residual_is_not_accepted(accepted):
    del accepted['nonlinear_stats']['nonlinear_final_relative_residual']
    with pytest.raises(ValueError):
        MODULE.require_accepted_hdiv(accepted)


@pytest.fixture
def restaged(tmp_path):
    package = tmp_path / 'new-venv/radia'
    (package / 'vim').mkdir(parents=True)
    (package / '__init__.py').write_bytes(b'# fixture package\n')
    (package / 'vim/_solve.py').write_bytes(b'# fixture solver\n')
    paths = {'native': package / '_radia_pybind.pyd',
             'iron': tmp_path / 'new-iron.vol', 'wheel': tmp_path / 'new-candidate.whl'}
    for name, path in paths.items():
        path.write_bytes(('fixture-' + name).encode())
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    implementation = {
        'identities': {name: {'path': 'Z:/deleted/old-run/' + path.name, 'sha256': digest(path)}
                       for name, path in paths.items()},
        'python_sources': {str(path.relative_to(package)).replace('/', '\\'): digest(path)
                           for path in package.rglob('*.py')},
    }
    return implementation, package, paths


def verify(restaged):
    implementation, package, paths = restaged
    return MODULE.verify_relocated_identity(implementation, package, paths['iron'], paths['wheel'])


def test_restaging_ignores_deleted_historical_paths(restaged):
    implementation, package, paths = restaged
    original = copy.deepcopy(implementation)
    record = verify(restaged)
    assert implementation == original
    assert record['package_path'] == str(package.resolve())
    for name, path in paths.items():
        assert record['identities'][name]['path'] == str(path.resolve())
        assert record['identities'][name]['sha256'] == implementation['identities'][name]['sha256']


@pytest.mark.parametrize('name', ['native', 'iron', 'wheel'])
def test_restaged_changed_bytes_are_rejected(restaged, name):
    restaged[2][name].write_bytes(b'incorrect bytes')
    with pytest.raises(RuntimeError, match=name):
        verify(restaged)


@pytest.mark.parametrize('name', ['native', 'iron', 'wheel'])
def test_missing_explicit_artifact_is_rejected(restaged, name):
    restaged[2][name].unlink()
    with pytest.raises(FileNotFoundError):
        verify(restaged)


def test_installed_python_drift_is_rejected(restaged):
    (restaged[1] / 'vim/_solve.py').write_bytes(b'changed solver')
    with pytest.raises(RuntimeError, match='Python sources'):
        verify(restaged)


@pytest.mark.parametrize('name', ['native', 'iron', 'wheel'])
def test_missing_identity_is_rejected(restaged, name):
    del restaged[0]['identities'][name]
    with pytest.raises(ValueError, match='exactly'):
        verify(restaged)


def test_unknown_identity_is_not_silently_ignored(restaged):
    restaged[0]['identities']['unknown'] = {'path': 'old'}
    with pytest.raises(ValueError, match='exactly'):
        verify(restaged)
