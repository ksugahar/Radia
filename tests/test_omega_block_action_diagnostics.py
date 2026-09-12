"""Check diagnostic arithmetic with explicit matrices, without native imports."""
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


def function(name, filename='diagnostics.py'):
    path = Path(__file__).resolve().parents[1] / 'validation_test/omega_quadrature' / filename
    tree = ast.parse(path.read_text(encoding='utf-8'))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    namespace = {'np': np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace[name]


class Vector:
    def __init__(self, values):
        self.values = np.array(values, dtype=float)

    def CreateVector(self):
        return Vector(np.zeros_like(self.values))

    def FV(self):
        return self

    def NumPy(self):
        return self.values

    def __setitem__(self, key, value):
        self.values[key] = value

    @property
    def data(self):
        return self.values

    @data.setter
    def data(self, value):
        self.values[:] = value

    def __sub__(self, other):
        return self.values - other


def result(matrix, x, rhs, free=(True, True, True)):
    class Matrix:
        def __mul__(self, vector):
            return np.array(matrix) @ vector.values

        def Inverse(self, mask, inverse):
            assert inverse == 'pardiso'

            class Inverse:
                def __mul__(self, vector):
                    value = np.zeros(len(mask))
                    selected = np.array(mask, dtype=bool)
                    value[selected] = np.linalg.solve(np.array(matrix)[np.ix_(selected, selected)],
                                                      vector.values[selected])
                    return value
            return Inverse()

    fes = SimpleNamespace(ndof=3, FreeDofs=lambda: free,
                          Range=lambda i: SimpleNamespace(start=i, stop=i+1))
    return {'fes': fes, 'solution': SimpleNamespace(vec=Vector(x), components=(0, 1, 2)),
            'system': {'bilinear_form': SimpleNamespace(mat=Matrix()),
                       'linear_form': SimpleNamespace(vec=Vector(rhs))}}


def test_cancellation_is_reported_without_granting_acceptance():
    report = function('block_action_residual')(result(
        [[1, 0, 0], [0, 1, -1], [0, 0, 1]], [1, 1, 1], [1, 1e-15, 1]))
    row = report['blocks']['phi_total']
    assert row['residual_l2'] == 1e-15
    assert row['rhs_l2'] == 1e-15
    assert row['column_action_l2']['phi_total'] == 1
    assert row['column_action_l2']['interface_constraint'] == 1
    assert row['action_relative'] == pytest.approx(5e-16)
    assert report['acceptance_evidence'] is False


def test_bad_solution_remains_visible():
    report = function('block_action_residual')(result(np.eye(3), [1, 2, 1], [1, 1, 1]))
    assert report['blocks']['phi_total']['action_relative'] == pytest.approx(1/3)


def test_nonfree_rows_are_excluded_but_columns_are_retained():
    report = function('block_action_residual')(result(
        [[1, 0, 0], [3, 1, 0], [0, 0, 1]], [1, 1, 1], [10, 4, 1], (False, True, True)))
    assert report['blocks']['phi_reduced']['action_relative'] is None
    assert report['blocks']['phi_total']['column_action_l2']['phi_reduced'] == 3
    assert report['blocks']['phi_total']['residual_l2'] == 0


@pytest.mark.parametrize('bad', [np.nan, np.inf])
def test_nonfinite_is_rejected(bad):
    with pytest.raises(ValueError, match='nonfinite'):
        function('block_action_residual')(result(np.eye(3), [bad, 1, 1], [1, 1, 1]))


def test_zero_scale_is_undecidable():
    report = function('block_action_residual')(result(np.zeros((3, 3)), [0, 0, 0], [0, 0, 0]))
    assert all(row['action_relative'] is None for row in report['blocks'].values())


def test_shared_volume_average_and_nonfinite_field():
    observe = function('field_observations')
    case = {'observation_samples': np.zeros((8, 3)), 'observation_centres': [[0, 0, 0]]}
    value = observe({'B_cf': lambda point: (1, 2, 3)}, lambda *point: point, case)
    assert value['B_average_T'] == [[1, 2, 3]]
    assert value['acceptance_evidence'] is False
    with pytest.raises(ValueError, match='invalid magnetic'):
        observe({'B_cf': lambda point: (1, np.nan, 3)}, lambda *point: point, case)


def test_algebraic_cli_preserves_hold_and_never_runs_energy_audit(tmp_path, monkeypatch):
    import argparse
    from contextlib import nullcontext
    import json
    import platform
    import sys
    import time

    path = Path(__file__).resolve().parents[1] / 'validation_test/omega_quadrature/run.py'
    tree = ast.parse(path.read_text(encoding='utf-8'))
    main_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    calls = []
    solver = SimpleNamespace(
        project_source_total_hodge=lambda *a, **kw: {
            'potential': 0, 'harmonic_field': 0, 'relative_harmonic_norm': 0},
        solve_magnetostatic_mixed_total_reduced_omega_kelvin=lambda *a, **kw: (
            calls.append(kw) or {'fes': SimpleNamespace(ndof=3), 'linear_residual': {'raw': 'retained'}}))
    case = dict(mesh=SimpleNamespace(ne=1), H_s=0, H_ext=0, controls={}, radius=1, center=(0, 0, 0), mu_r=1)

    def forbidden(*args, **kwargs):
        raise AssertionError('algebraic mode must not run an energy or embedding audit')

    namespace = dict(argparse=argparse, Path=Path, json=json, platform=platform, time=time,
                     check_runtime=lambda mode: {'mode': mode},
                     identity=lambda *a: {'source': 'fixed'}, digest=lambda p: 'hash',
                     load_module=lambda p, name: solver if 'solver' in name else SimpleNamespace(create_case=lambda p: case),
                     ng=SimpleNamespace(SetNumThreads=lambda n: None, TaskManager=nullcontext),
                     block_action_residual=lambda r: {'acceptance_evidence': False},
                     field_observations=lambda *a: None, audit_energy=forbidden,
                     constraint_violation=forbidden)
    exec(compile(ast.Module(body=[main_node], type_ignores=[]), str(path), 'exec'), namespace)
    output = tmp_path / 'result.json'
    monkeypatch.setattr(sys, 'argv', ['run.py', '--mode', 'research', '--research-solver', 'solver.py',
                                    '--factory', 'factory.py', '--mesh', 'mesh.vol', '--output', str(output),
                                    '--orders', '1', '--bonuses', '4', '--algebraic-only'])
    assert namespace['main']() == 2
    report = json.loads(output.read_text())
    assert len(calls) == 1 and calls[0]['return_system'] is True
    assert report['completed'] and report['acceptance'].startswith('HOLD:')
    assert report['nesting'] == []
    row = report['rows'][0]
    assert row['linear_residual'] == {'raw': 'retained'}
    assert row['gates'] == {'energy_audit_completed': False, 'three_engine_acceptance': False}
    assert 'energy' not in row


def field_payload(bonus):
    return dict(completed=True, source_unchanged=True, mesh_unchanged=True,
                implementation={'native': {'test.pyd': 'fixture'}}, mesh_sha256='fixture',
                case={'mu_r': 1000}, runtime={'mode': 'wheel'},
                controls=dict(source_order=3, threads=8, evaluation_order=16, algebraic_only=True),
                rows=[dict(order=1, ndof=3, bonus=bonus, field_observations={
                    'centres_m': [[0, 0, 0]], 'samples_m': [[0, 0, 0]]*8,
                    'B_samples_T': [[1, 2, 3]]*8, 'B_average_T': [[1, 2, 3]],
                    'observable': 'shared eight-point volume average'})])


def test_field_delta_is_not_a_convergence_certificate():
    a, b = field_payload(4), field_payload(8)
    observe = b['rows'][0]['field_observations']
    observe['B_samples_T'] = [[2, 2, 3]]*8
    observe['B_average_T'] = [[2, 2, 3]]
    report = function('compare_fields', 'compare_fields.py')(a, b)
    assert report['difference_rms_T'] == 1
    assert report['difference_relative_rms'] == pytest.approx(1/np.sqrt(17))
    assert report['acceptance'].startswith('HOLD:')


@pytest.mark.parametrize('damage', ['identity', 'unfinished', 'average', 'points', 'order', 'nan'])
def test_field_delta_rejects_incomparable_inputs(damage):
    a, b = field_payload(4), field_payload(8)
    row = b['rows'][0]
    if damage == 'identity':
        b['implementation']['native'] = {'test.pyd': 'different'}
    elif damage == 'unfinished':
        b['completed'] = False
    elif damage == 'average':
        row['field_observations']['B_average_T'] = [[3, 2, 1]]
    elif damage == 'points':
        row['field_observations']['samples_m'] = [[1, 0, 0]]*8
    elif damage == 'order':
        row['order'] = 2
    else:
        row['field_observations']['B_samples_T'][0] = [np.nan, 0, 0]
    with pytest.raises(ValueError):
        function('compare_fields', 'compare_fields.py')(a, b)


@pytest.mark.parametrize('fail_after_update', [False, True])
def test_residual_correction_measures_and_restores_original(fail_after_update):
    data = result(np.eye(3), [1, 2, 1], [1, 1, 1])
    original = data['solution'].vec.values.copy()
    corrected = function('residual_correction_observation')
    count = []

    def fields(*args):
        count.append(1)
        if len(count) == 2 and fail_after_update:
            raise ValueError('field probe failure')
        return data['solution'].vec.values.tolist()

    corrected.__globals__.update(field_observations=fields,
                                 block_action_residual=function('block_action_residual'))
    if fail_after_update:
        with pytest.raises(ValueError, match='field probe failure'):
            corrected(data, None, {'observation_samples': [1]})
    else:
        report = corrected(data, None, {'observation_samples': [1]})
        assert report['field_before'] == [1, 2, 1]
        assert report['field_after'] == [1, 1, 1]
        assert report['correction_coefficient_l2'] == 1
        assert report['corrected_block_residual']['blocks']['phi_total']['residual_l2'] == 0
        assert report['acceptance_evidence'] is False
    np.testing.assert_array_equal(data['solution'].vec.values, original)
