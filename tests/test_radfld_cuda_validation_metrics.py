"""The numerical acceptance gate must not hide error through norm overflow."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

@pytest.fixture
def metrics(monkeypatch):
    folder=Path(__file__).resolve().parents[1]/'validation_test/radfld_cuda'
    monkeypatch.syspath_prepend(str(folder))
    spec=importlib.util.spec_from_file_location('radfld_metrics',folder/'validate.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

@pytest.mark.parametrize('magnitude',[1e-250,1.,1e250])
def test_vector_error_is_scale_invariant(metrics,magnitude):
    reference=np.array([[1.,-2.,3.],[4.,1.,-2.]])*magnitude
    value=reference*1.25
    assert metrics.relative(value,reference)==pytest.approx(.25)
    assert metrics.point_relative(value,reference)==pytest.approx(.25)

@pytest.mark.parametrize('bad',[np.nan,np.inf,-np.inf])
def test_nonfinite_results_cannot_pass(metrics,bad):
    value=np.array([[bad,0.,0.]])
    for check in (metrics.relative,metrics.point_relative):
        with pytest.raises(RuntimeError,match='Nonfinite'):
            check(value,np.ones((1,3)))

def test_extreme_opposite_vectors_and_zero_fields(metrics):
    reference=np.full((2,3),1e308)
    for check in (metrics.relative,metrics.point_relative):
        assert check(-reference,reference)==pytest.approx(2.)
        assert check(np.zeros((2,3)),np.zeros((2,3)))==0.
        with pytest.raises(RuntimeError,match='shape mismatch'):
            check(np.ones((1,3)),np.ones((2,3)))
