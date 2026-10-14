"""Fast gates: a completed diagnostic is not numerical release acceptance."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest


@pytest.fixture
def lane():
    directory = Path(__file__).resolve().parents[1] / 'validation_test/omega_quadrature'
    sys.path.insert(0, str(directory))
    try:
        spec = importlib.util.spec_from_file_location('_omega_validation', directory / 'run.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.path.remove(str(directory))


def sample():
    return {'linear_residual': {'free_dofs': {'relative': 1e-12},
            'blocks': {'interface_constraint': {'relative': 1e-15}}},
            'assembled_energy': {'energy': 2.0, 'default_reconstruction_difference': 0.0,
            'fixed_rule_identity': {'W': 3.0, 'J': 2.0, 'source_offset': 1.0,
                                    'W_minus_J_minus_offset': 0.0}}}


def test_consistent_diagnostics(lane):
    assert all(lane.gates(sample()).values())


@pytest.mark.parametrize('value', [None, float('nan'), 1e-4, -1e-12])
def test_missing_or_bad_constraint_residual_fails(lane, value):
    data = sample()
    data['linear_residual']['blocks']['interface_constraint']['relative'] = value
    assert not lane.gates(data)['residual_interface_constraint']


def test_energy_mismatch_fails(lane):
    data = sample()
    data['assembled_energy']['fixed_rule_identity']['W_minus_J_minus_offset'] = 0.01
    assert not lane.gates(data)['same_rule_energy']


@pytest.mark.parametrize('value', [float('inf'), float('-inf'), float('nan')])
@pytest.mark.parametrize('key', ['W', 'J', 'source_offset', 'W_minus_J_minus_offset'])
def test_nonfinite_identity_fails(lane, key, value):
    data = sample()
    data['assembled_energy']['fixed_rule_identity'][key] = value
    assert not lane.gates(data)['same_rule_energy']


@pytest.mark.parametrize('key', ['energy', 'default_reconstruction_difference'])
def test_infinite_reconstruction_fails(lane, key):
    data = sample()
    data['assembled_energy'][key] = float('inf')
    assert not lane.gates(data)['assembly_reconstruction']


def test_negative_free_residual_fails(lane):
    data = sample()
    data['linear_residual']['free_dofs']['relative'] = -1e-12
    assert not lane.gates(data)['free_residual']


def test_norm_nonfinite_rejected(lane):
    from diagnostics import norm_report
    with pytest.raises(ValueError):
        norm_report(float('nan'), 1.0)


def test_nearly_zero_block_rhs_does_not_silently_pass(lane):
    data = sample()
    data['linear_residual']['blocks']['phi_total'] = {
        'rhs_l2': 1e-15, 'residual_l2': 2e-16, 'relative': 0.2}
    assert not lane.gates(data)['residual_phi_total']


def test_wheel_rejects_editable(lane, monkeypatch):
    from types import SimpleNamespace
    monkeypatch.setitem(sys.modules, 'radia', SimpleNamespace(__file__=__file__))
    dist = SimpleNamespace(read_text=lambda name: '{"dir_info":{"editable":true}}')
    monkeypatch.setattr(lane.importlib.metadata, 'distribution', lambda name: dist)
    with pytest.raises(RuntimeError, match='editable'):
        lane.check_runtime('wheel')


def test_preserved_candidate_is_hold_and_retains_failed_gate(lane):
    directory = Path(lane.__file__).parent / 'candidate_de7feea0'
    result = json.loads((directory / 'threads8/result.json').read_text())
    preflight = json.loads((directory / 'preflight.json').read_text())
    telemetry = json.loads((directory / 'threads8/telemetry.json').read_text())
    assert result['completed'] and result['source_unchanged'] and result['mesh_unchanged']
    assert result['acceptance'].startswith('HOLD:')
    assert preflight['candidate_unpublished']
    assert result['mesh_sha256'] == preflight['mesh_sha256']
    assert preflight['native_sha256'] in result['implementation']['native'].values()
    assert result['runtime']['direct_url']['archive_info']['hashes']['sha256'] == preflight['wheel_sha256']
    assert telemetry['exit_code'] == 2
    assert result['controls']['orders'] == [1]
    assert result['controls']['bonuses'] == [4, 8]
    assert result['controls']['threads'] == 8
    assert len(result['rows']) == 2
    for row in result['rows']:
        assert lane.gates({'linear_residual': row['linear_residual'],
                           'assembled_energy': row['energy']}) == row['gates']
    assert not result['rows'][0]['gates']['residual_phi_total']
    assert all(result['rows'][1]['gates'].values())
