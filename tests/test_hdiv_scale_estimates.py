"""Small contract checks for scale preflight arithmetic, without FE/GPU execution."""
import importlib.util
from pathlib import Path
import pytest
path=Path(__file__).parents[1]/'validation_test/hmatrix_cuda/prepare_hdiv_scale.py'
spec=importlib.util.spec_from_file_location('hdiv_scale',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def test_dense_scale_estimates_and_guards():
    result=mod.estimate(10000,500,200,2,45,16)
    assert result['charges']==3200 and result['mass_dense_bytes']==800000000
    assert result['potential_cache_bytes']==657920000
    assert result['native_dimension_gate'] and result['dense_array_gate'] and result['cache_workspace_lower_bound_gate']
    assert not mod.estimate(17000,500,200,2,45,16)['native_dimension_gate']
    large=mod.estimate(10000,2000,1000,2,45,16)
    assert large['estimated_resident_potential_cache_bytes']<=256*2**20
    assert large['potential_cache_bytes']>2*2**30
    assert not large['native_dimension_gate']
    with pytest.raises(ValueError):mod.estimate(-1,1,1,2,45,16)


def test_sparse_staging_cannot_certify_large_dense_solver_capacity():
    # A dimension-only hypothetical case, not a generated or accepted FE mesh.
    result=mod.estimate(100000,5000,1000,2,45,16)
    assert result['mass_dense_bytes']==80_000_000_000
    assert not result['mass_cholesky_array_gate']
    assert not result['dense_gram_array_gate']
    assert not result['native_dimension_gate']
    # Source tiling can fit while both dense numerical allocations still fail.
    assert result['cache_workspace_lower_bound_gate']
