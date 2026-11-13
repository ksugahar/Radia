"""Energy-Newton step acceptance: a residual certificate first, then Armijo on a directly formed energy change.

Guards: the `_line_search` closure of `_solve_nonlinear_energy_cpp` and `_energy_armijo_line_search`.
Assumes: synthetic one-unknown energies; the line-search code is executed from the source, not the solver."""
import ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np

SOURCE=Path(__file__).resolve().parents[1]/"src/radia/vim/_solve.py"


def _line_search(env):
    tree=ast.parse(SOURCE.read_text(encoding="utf-8"))
    outer=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="_solve_nonlinear_energy_cpp")
    search=next(n for n in outer.body if isinstance(n,ast.FunctionDef) and n.name=="_line_search")
    helper=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="_energy_armijo_line_search")
    env.update(np=np)
    exec(compile(ast.Module(body=[helper,search],type_ignores=[]),str(SOURCE),"exec"),env)
    return env["_line_search"]


def test_root_certificate_precedes_rounded_energy_armijo():
    stats={}
    def update(x):
        material.field=x.copy()
        return 1e20  # Every Armijo decrement is lost: only the certificate can accept.
    material=SimpleNamespace(field=np.zeros(1),update=update)
    search=_line_search(dict(material=material,stats=stats,_bH=lambda x:x,
                             _N_apply=lambda x:np.zeros_like(x),nl_tol=5e-14,residual_scale=1.))
    alpha,_,backtracks=search(np.zeros(1),np.ones(1),np.ones(1),1.,0.,np.zeros(1))
    assert (alpha,backtracks)==(1.,0), (alpha,backtracks)
    assert stats["nonlinear_residual_certified_trials"]==1, stats
    # Armijo may accept a non-root, but only a root produces a residual
    # certificate, and a nonfinite direction produces neither.
    stats.clear()
    search(np.zeros(1),np.array([.5]),np.ones(1),.5,0.,np.zeros(1))
    assert not stats, stats
    assert search(np.zeros(1),np.array([np.nan]),np.ones(1),1.,0.,np.zeros(1))==(None,None,0)
    assert not stats, stats


def _jittered_problem():
    """E(x) = 1/2 (x-c)^2 + 1/2 x.Nx - rhs.x with N = 1 and root 1.5e4.

    The N-apply carries jitter of -1e-7 at the current point and +1e-7 at every
    trial, the sign pattern of a failing H-matrix realization: the totals then
    differ by about +1.5e-3 while a Newton step lowers E by about 1e-6."""
    c,rhs=1e4,np.array([2e4]);m=np.array([1.5e4+1e-3])
    def update(x):
        material.field=x-c
        return 0.5*float((x-c)@(x-c))
    material=SimpleNamespace(field=None,update=update)
    def apply(x):
        return x+(-1e-7 if np.array_equal(x,m) else 1e-7)
    stats={}
    search=_line_search(dict(material=material,stats=stats,_bH=lambda f:f,_N_apply=apply,
                             nl_tol=1e-12,residual_scale=1.))
    material_energy=update(m);N_m=apply(m)
    R=material.field+N_m-rhs
    return search,m,rhs,R,material_energy,N_m,stats


def test_armijo_accepts_a_newton_step_below_the_energy_total_rounding():
    search,m,rhs,R,material_energy,N_m,stats=_jittered_problem()
    dm=-R/2.0                                    # exact Newton step, J = 2
    alpha,trial_energy,backtracks=search(m,dm,rhs,float(-dm@R),material_energy,N_m)
    assert (alpha,backtracks)==(1.,0), (alpha,backtracks)
    assert np.isfinite(trial_energy) and not stats, (trial_energy,stats)


def test_armijo_uses_the_curvature_of_the_step():
    search,m,rhs,R,material_energy,N_m,stats=_jittered_problem()
    dm=-1.5*R                                    # three Newton steps: E rises at alpha 1, falls at 1/2
    alpha,_,backtracks=search(m,dm,rhs,float(-dm@R),material_energy,N_m)
    assert (alpha,backtracks)==(.5,1), (alpha,backtracks)


def test_compensated_budget_rejects_implicit_conversion_and_overflow():
    import pytest
    source=Path(__file__).resolve().parents[1]/"src/radia/vim/_solve.py"
    tree=ast.parse(source.read_text(encoding="utf-8"));function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="_checked_compensated_memory_mb")
    env=dict(np=np);exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),"exec"),env)
    check=env["_checked_compensated_memory_mb"]
    for value in (0,1,16,512,np.int64(16)):
        assert check(value)==int(value)
    for value in (True,False,-1,513,2**64,16.,float("nan"),"16",None):
        with pytest.raises(ValueError,match="budget"):check(value)
