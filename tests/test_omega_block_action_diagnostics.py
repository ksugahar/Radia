"""Check diagnostic arithmetic with explicit matrices, without native imports."""
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


def function(name):
    path = Path(__file__).resolve().parents[1] / 'validation_test/omega_quadrature/diagnostics.py'
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
