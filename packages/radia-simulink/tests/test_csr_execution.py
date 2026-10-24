"""Independent NGSolve CPU solve and strict optional CUDA dispatch contracts."""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('csr_execution',ROOT/'src/radia_simulink/csr.py')
csr=importlib.util.module_from_spec(spec);spec.loader.exec_module(csr)


def test_ngsolve_cpu_solve_is_independent_of_cuda_binding(monkeypatch):
    from radia.ngsolve_sparse import real_sparse_matrix_type
    ng=pytest.importorskip('ngsolve')
    # An accidental GPU import is an error even if CUDA happens to be installed.
    monkeypatch.setitem(sys.modules,'radia_simulink._csr_cuda',None)
    monkeypatch.setitem(sys.modules,'radia_simulink._cuda_runtime',None)
    a=real_sparse_matrix_type().CreateFromCOO([0,0,1,1,1,2,2],[0,1,0,1,2,1,2],
                                          [4.,-1.,-1.,4.,-1.,-1.,3.],3,3)
    b=a.CreateColVector();b.FV().NumPy()[:]=[1.,2.,3.]
    with ng.TaskManager():
        x,info=csr.solve_spd(a,b,backend='cpu',tolerance=1e-12)
    expected=np.linalg.solve(np.array([[4.,-1.,0.],[-1.,4.,-1.],[0.,-1.,3.]]),[1.,2.,3.])
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



def test_gpu_result_cannot_hide_a_bad_original_equation_residual(monkeypatch):
    from radia.ngsolve_sparse import real_sparse_matrix_type
    import types
    ng=pytest.importorskip('ngsolve')
    a=real_sparse_matrix_type().CreateFromCOO([0],[0],[2.],1,1)
    b=a.CreateColVector();b.FV().NumPy()[:]=[1.]
    fake=types.SimpleNamespace(solve=lambda *args:(np.zeros(1),dict(backend='gpu',true_relative_residual=0.)))
    package=types.ModuleType('radia_simulink');package.__path__=[];package._csr_cuda=fake
    monkeypatch.setitem(sys.modules,'radia_simulink',package)
    monkeypatch.setitem(sys.modules,'radia_simulink._cuda_runtime',types.SimpleNamespace(prepare_cuda_runtime=lambda:None))
    original=csr.__package__;csr.__package__='radia_simulink'
    try:
        with pytest.raises(RuntimeError,match='original-system residual'):
            csr.solve_spd(a,b,backend='gpu')
    finally:csr.__package__=original
