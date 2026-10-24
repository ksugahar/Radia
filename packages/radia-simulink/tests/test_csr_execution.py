"""Guards: NGSolve CPU independence and strict optional CUDA dispatch/residuals.
Assumes: NGSolve is installed; these CSR contracts need no Radia distribution.
"""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('csr_execution',ROOT/'src/radia_simulink/csr.py')
csr=importlib.util.module_from_spec(spec);spec.loader.exec_module(csr)


@pytest.fixture
def spd_system():
    ng=pytest.importorskip('ngsolve')
    from netgen.geom2d import unit_square
    mesh=ng.Mesh(unit_square.GenerateMesh(maxh=1))
    space=ng.H1(mesh,order=1)
    u,v=space.TnT()
    form=ng.BilinearForm(space,symmetric=False)
    form+=(ng.grad(u)*ng.grad(v)+u*v)*ng.dx
    with ng.TaskManager():
        form.Assemble()
    a=form.mat
    b=a.CreateColVector();b.FV().NumPy()[:]=np.arange(1,a.height+1,dtype=np.float64)
    return ng,a,b


def test_ngsolve_cpu_solve_is_independent_of_cuda_binding(monkeypatch,spd_system):
    ng,a,b=spd_system
    # An accidental GPU import is an error even if CUDA happens to be installed.
    monkeypatch.setitem(sys.modules,'radia_simulink._csr_cuda',None)
    monkeypatch.setitem(sys.modules,'radia_simulink._cuda_runtime',None)
    with ng.TaskManager():
        x,info=csr.solve_spd(a,b,backend='cpu',tolerance=1e-12)
    rows,columns,values=a.COO()
    dense=np.zeros((a.height,a.width))
    np.add.at(dense,(np.asarray(rows),np.asarray(columns)),np.asarray(values))
    expected=np.linalg.solve(dense,b.FV().NumPy())
    np.testing.assert_allclose(x,expected,rtol=1e-13,atol=1e-14)
    assert info['backend']=='cpu' and info['true_relative_residual']<1e-13


@pytest.mark.parametrize('backend',[None,'auto','cuda',''])
def test_backend_must_be_explicit(backend):
    with pytest.raises(ValueError,match='backend'):
        csr.solve_spd(None,None,backend=backend)


def test_backend_is_required():
    with pytest.raises(TypeError):csr.solve_spd(None,None)


def test_gpu_missing_binding_propagates_without_cpu_fallback(monkeypatch):
    import types
    package=types.ModuleType('radia_simulink');package.__path__=[]
    monkeypatch.setitem(sys.modules,'radia_simulink',package)
    monkeypatch.setitem(sys.modules,'radia_simulink._cuda_runtime',types.SimpleNamespace(prepare_cuda_runtime=lambda:None))
    monkeypatch.setitem(sys.modules,'radia_simulink._csr_cuda',None)
    original=csr.__package__;csr.__package__='radia_simulink'
    class Vector:
        def FV(self):return self
        def NumPy(self):return np.ones(1)
    class Matrix:
        height=width=1
        def Inverse(self,**kwargs):raise AssertionError('CPU fallback forbidden')
    try:
        with pytest.raises(ImportError):csr.solve_spd(Matrix(),Vector(),backend='gpu')
    finally:csr.__package__=original



def test_gpu_result_cannot_hide_a_bad_original_equation_residual(monkeypatch,spd_system):
    import types
    ng,a,b=spd_system
    fake=types.SimpleNamespace(solve=lambda *args:(np.zeros(a.height),dict(backend='gpu',true_relative_residual=0.)))
    package=types.ModuleType('radia_simulink');package.__path__=[];package._csr_cuda=fake
    monkeypatch.setitem(sys.modules,'radia_simulink',package)
    monkeypatch.setitem(sys.modules,'radia_simulink._cuda_runtime',types.SimpleNamespace(prepare_cuda_runtime=lambda:None))
    original=csr.__package__;csr.__package__='radia_simulink'
    try:
        with pytest.raises(RuntimeError,match='original-system residual'):
            csr.solve_spd(a,b,backend='gpu')
    finally:csr.__package__=original
