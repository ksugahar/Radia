"""Independent AMS numerical acceptance must reject nonfinite oracle data."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

path=Path(__file__).resolve().parents[3]/'validation_test/ams_cuda/ams_tuning.py'
spec=importlib.util.spec_from_file_location('ams_acceptance',path)
acceptance=importlib.util.module_from_spec(spec)
spec.loader.exec_module(acceptance)

@pytest.mark.parametrize('index',range(4))
@pytest.mark.parametrize('invalid',[np.nan,np.inf,-np.inf])
def test_nonfinite_solution_reference_operator_or_rhs_is_rejected(index,invalid):
    vectors=[np.array([1.,2.]) for _ in range(4)]
    vectors[index][0]=invalid
    with pytest.raises(RuntimeError,match='nonfinite'):
        acceptance.checked_errors(*vectors,np.array([True,True]),1e-10)

def test_zero_rhs_and_original_equation_residual():
    zero=np.zeros(2);free=np.array([True,True])
    assert acceptance.checked_errors(zero,zero,zero,zero,free,1e-10)==dict(solution_relative=0.,true_relative_residual=0.)
    with pytest.raises(RuntimeError,match='true residual failed'):
        acceptance.checked_errors(np.ones(2),np.ones(2),np.ones(2),np.full(2,2.),free,1e-10)

@pytest.mark.parametrize('index',[1,3])
def test_overflowing_reference_norm_is_rejected(index):
    vectors=[np.ones(2) for _ in range(4)];vectors[index][:]=1e308
    with np.errstate(over='ignore'),pytest.raises(RuntimeError,match='nonfinite norm'):
        acceptance.checked_errors(*vectors,np.array([True,True]),1e-10)
