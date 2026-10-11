"""Finite large-loop route boundaries, without allocating large matrices."""
from types import SimpleNamespace
import pytest
from radia.bem_loop_work import _check_loop_work_size
from radia.panels.calc_inductance import _resolve_workpiece_backend


class Handle:
    def __init__(self, route='p1-entry-on-demand', valid=True):
        self.route,self.valid=route,valid
    def IsValid(self):return self.valid
    def GetStats(self):return dict(construction_route=self.route)


def body(ndof=10000, route='p1-entry-on-demand'):
    return SimpleNamespace(ndof=ndof,_loop_body_backend='hacapk',
        _SL_hacapk=Handle(route),_DL_hacapk=Handle(route))


def test_large_work_requires_actual_both_routes_and_bounded_full_dofs():
    _check_loop_work_size(body(),20000,'fmm')
    for count,backend,solver,match in (
        (20001,'fmm',body(),'20000'),
        (20000,'fmm',body(10001),'10000'),
        (14001,'dense',body(),'FMM'),
        (14001,'auto',body(),'FMM'),
        (14001,'fmm',body(route='dense-entry'),'actual on-demand')):
        with pytest.raises(ValueError,match=match):
            _check_loop_work_size(solver,count,backend)
    solver=body();solver._DL_hacapk=Handle(valid=False)
    with pytest.raises(ValueError,match='actual on-demand'):
        _check_loop_work_size(solver,14001,'fmm')
    # The existing <=14000-face routes keep their prior contract.
    _check_loop_work_size(body(route='dense-entry'),14000,'dense')


def args(**overrides):
    values=dict(wp_bem_backend='hacapk',impedance_model='sibc',h1_order=1,
                coupling_mode='weak',wp_loop_work_backend='fmm')
    values.update(overrides)
    return SimpleNamespace(**values)


def test_large_weak_cli_is_explicit_bounded_p1_and_fmm():
    assert _resolve_workpiece_backend(args(),1,10000)=='hacapk'
    for count,settings,match in (
        (10001,{},'10000'),
        (7001,dict(wp_bem_backend='intree-dense'),'7000'),
        (7001,dict(wp_bem_backend='auto'),'7000'),
        (7001,dict(coupling_mode='strong'),'7000'),
        (7001,dict(h1_order=2),'P1'),
        (7001,dict(wp_loop_work_backend='dense'),'FMM')):
        with pytest.raises(ValueError,match=match):
            _resolve_workpiece_backend(args(**settings),1,count)
    with pytest.raises(ValueError,match='multiple holes'):
        _resolve_workpiece_backend(args(),2,100)


@pytest.mark.parametrize('module,classname',[
    ('radia.bem_coupled_solver','CoupledBEMSolver'),
    ('radia.peec_coupled_bem_solver','CoupledPEECBEMSolver')])
def test_public_strong_api_refuses_large_loop_before_iterate_allocation(module,classname):
    from importlib import import_module
    cls=getattr(import_module(module),classname)
    solver=object.__new__(cls)
    solver.wp_solver=SimpleNamespace(ndof=7001)
    with pytest.raises(ValueError,match='7000-nodal'):
        solver.solve(1+1j,2.,loop_dof=True)
