"""A rounded merit must not reject an independently certified Newton root."""
import ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np


def test_root_certificate_precedes_rounded_energy_armijo():
    source=Path(__file__).resolve().parents[1]/"src/radia/vim/_solve.py"
    tree=ast.parse(source.read_text(encoding="utf-8"))
    outer=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="_solve_nonlinear_energy_cpp")
    search=next(n for n in outer.body if isinstance(n,ast.FunctionDef) and n.name=="_line_search")
    material=SimpleNamespace(field=np.zeros(1));stats={}
    def energy(x,rhs):
        material.field=x.copy()
        return 1e20  # Every Armijo decrement is lost to rounding.
    env=dict(np=np,material=material,stats=stats,_energy=energy,_bH=lambda x:x,
             _N_apply=lambda x:np.zeros_like(x),nl_tol=5e-14,residual_scale=1.)
    exec(compile(ast.Module(body=[search],type_ignores=[]),str(source),"exec"),env)
    alpha,_,backtracks=env["_line_search"](np.zeros(1),np.ones(1),np.ones(1),1e20,1.)
    assert alpha==1. and backtracks==0
    assert stats["nonlinear_residual_certified_trials"]==1
    # A rounded energy may satisfy non-strict Armijo, but cannot produce a
    # residual certificate for a non-root or a nonfinite trial.
    stats.clear()
    env["_line_search"](np.zeros(1),np.array([.5]),np.ones(1),1e20,1.)
    assert not stats
    env["_line_search"](np.zeros(1),np.array([np.nan]),np.ones(1),1e20,1.)
    assert not stats


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
