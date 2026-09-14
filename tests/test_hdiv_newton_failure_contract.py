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


@pytest.mark.parametrize('residual,mode,exhausted,expected', [
    (1e-6, 'tolerance', False, True), (3e-5, 'tolerance', False, False),
    (None, 'tolerance', False, False), (float('nan'), 'tolerance', False, False),
    (float('inf'), 'tolerance', False, False), (1e-6, 'settled', False, False),
    (1e-6, 'line-search-exhausted', True, False), (1e-6, 'tolerance', True, False),
])
def test_shared_hdiv_runner_checks_actual_residual(residual, mode, exhausted, expected):
    import contextlib
    import time
    from types import SimpleNamespace
    source = Path(__file__).resolve().parents[1] / 'validation_test/c_type_three_engine/run_three_engine.py'
    tree = ast.parse(source.read_text(encoding='utf-8'))
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'solve_hdiv')
    result = {'ndof': 3, 'iters': 1, 'nonlinear_solve_stats': {
        'nonlinear_converged_final_stage': True, 'nonlinear_convergence_mode': mode,
        'nonlinear_line_search_exhausted': exhausted, 'nonlinear_final_relative_residual': residual}}
    ns = dict(np=np, time=time, MU0=4e-7*np.pi,
              ng=SimpleNamespace(Mesh=object, TaskManager=contextlib.nullcontext),
              rad=SimpleNamespace(RadiaField=lambda *a: None, Fld=lambda *a: np.zeros((1, 3))),
              vim=SimpleNamespace(Solve=lambda *a, **kw: result,
                                  FieldFromSolution=lambda *a, **kw: np.zeros((1, 3))))
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), ns)
    _, diagnostics = ns['solve_hdiv'](SimpleNamespace(ne=1, nv=8), 1, [], nonlinear=True,
        order=1, gram_eps=1e-12, nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80,
        points=np.zeros((1, 3)))
    assert diagnostics['nonlinear_stats']['converged'] is expected
