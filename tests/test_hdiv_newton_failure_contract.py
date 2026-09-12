"""Execute the production inner-solve closure without constructing a FE system.

Only the native calls and surrounding closure variables are test doubles.
This tests exception metadata, not native convergence or numerical correctness.
"""
import ast
from pathlib import Path

import numpy as np
import pytest


def _inner_solve(preconditioner, result):
    source = Path(__file__).resolve().parents[1] / 'src/radia/vim/_solve.py'
    tree = ast.parse(source.read_text(encoding='utf-8'))
    outer = next(node for node in tree.body
                 if isinstance(node, ast.FunctionDef)
                 and node.name == '_solve_nonlinear_energy_cpp')
    inner = next(node for node in outer.body
                 if isinstance(node, ast.FunctionDef) and node.name == '_solve_W')
    calls = []

    def native(*args, **kwargs):
        calls.append((args, kwargs))
        return result

    captured = []
    namespace = dict(np=np, cg_tol=1e-8, cg_maxit=7, n_face=2, H=object(),
                     inner_preconditioner=preconditioner,
                     _h_solve_auto_prec=native, _h_solve_mass_riesz=native,
                     _capture_cpp_solve_timings=captured.append)
    exec(compile(ast.Module(body=[inner], type_ignores=[]), str(source), 'exec'), namespace)
    return namespace['_solve_W'], calls, captured


@pytest.mark.parametrize('preconditioner', ['jacobi', 'mass-riesz'])
def test_failure_records_the_requested_inner_solve(preconditioner):
    result = {'iters': 7, 'm': [1., 2.]}
    solve, calls, captured = _inner_solve(preconditioner, result)
    with pytest.raises(RuntimeError, match='context=newton_stage_2_iteration_4') as exc:
        solve(object(), np.ones(2), tol_override=1e-5,
              context='newton_stage_2_iteration_4')
    error = exc.value
    assert error.linear_context == 'newton_stage_2_iteration_4'
    assert error.linear_iterations == 7
    assert error.linear_requested_tolerance == 1e-5
    assert error.linear_preconditioner == preconditioner
    assert 'operator is SPD' not in str(error)
    assert captured == [result]
    assert len(calls) == 1


@pytest.mark.parametrize('preconditioner', ['jacobi', 'mass-riesz'])
def test_success_and_tolerance_floor_are_unchanged(preconditioner):
    result = {'iters': 3, 'm': [1., 2.]}
    solve, calls, captured = _inner_solve(preconditioner, result)
    warm = np.array([0.5, 1.])
    values, iterations = solve(object(), np.ones(2), tol_override=1e-10, x0=warm)
    np.testing.assert_array_equal(values, result['m'])
    assert iterations == 3
    assert calls[0][0][5] == 1e-8
    assert calls[0][1]['x0'] is warm
    assert captured == [result]


def test_all_inner_call_sites_supply_a_context():
    source = Path(__file__).resolve().parents[1] / 'src/radia/vim/_solve.py'
    tree = ast.parse(source.read_text(encoding='utf-8'))
    outer = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name == '_solve_nonlinear_energy_cpp')
    calls = [node for node in ast.walk(outer) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name) and node.func.id == '_solve_W']
    assert len(calls) == 3
    assert all(any(keyword.arg == 'context' for keyword in call.keywords) for call in calls)
