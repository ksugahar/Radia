"""A candidate IC binary must be tied to its clean GPU source receipt."""
from pathlib import Path
import importlib.util
import pytest

path=Path(__file__).resolve().parents[1]/'verify_wheel.py'
spec=importlib.util.spec_from_file_location('ic_wheel_verifier',path)
verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)

def records():
    manifest=dict(source_commit='a'*40,files={'matlab/radia_cuda_ic.mexw64':'b'*64})
    receipt=dict(schema='radia.gpu.ic-payload.v1',source_commit='a'*40,source_dirty=False,
                 binary_sha256='b'*64,gpu_numerical_factorization=True,no_cpu_solver_fallback=True)
    return manifest,receipt

def test_exact_ic_receipt():
    verifier.verify_ic_receipt(*records())

@pytest.mark.parametrize('field,value', [('source_commit','c'*40),('source_dirty',True),
    ('binary_sha256','d'*64),('gpu_numerical_factorization',False),('no_cpu_solver_fallback',False)])
def test_reject_drift_or_cpu_substitution(field,value):
    manifest,receipt=records();receipt[field]=value
    with pytest.raises(ValueError,match='IC CUDA binary receipt mismatch'):
        verifier.verify_ic_receipt(manifest,receipt)
