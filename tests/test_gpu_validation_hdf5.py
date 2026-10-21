"""Current GPU validation helpers use HDF5; historical receipts stay immutable."""
import ast
import importlib.util
from pathlib import Path
import re
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('relative',[
    'tests/matlab/test_prepared_hdiv_shared_workflow.m',
    'tests/matlab/test_prepared_hdiv_host_setup.m',
    'validation_test/ams_cuda/checkSharedAMS.m',
    'validation_test/ams_cuda/check_shared_ams.py',
    'validation_test/ngsolve_cuda_mex/checkCsrBackends.m',
    'validation_test/ngsolve_cuda_mex/check_csr_backends.py',
])
def test_active_owned_helpers_do_not_use_mat_storage(relative):
    source=(ROOT/relative).read_text(encoding='utf-8')
    if relative.endswith('.py'):
        tree=ast.parse(source)
        assert not any(isinstance(node,ast.ImportFrom) and node.module=='scipy.io'
                       and any(alias.name in ('loadmat','savemat') for alias in node.names)
                       for node in ast.walk(tree))
    else:
        assert not re.search(r'\b(?:save|load|matfile)\s*\(',source)
    assert not re.search(r"\.mat['\"]",source)


def test_csr_cpu_oracle_fixtures_roundtrip_hdf5(tmp_path,monkeypatch):
    pytest.importorskip('ngsolve')
    pytest.importorskip('h5py')
    spec=importlib.util.spec_from_file_location('csr_validation_fixture',
        ROOT/'validation_test/ngsolve_cuda_mex/check_csr_backends.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    prepared=module.fixtures()
    monkeypatch.setattr(module,'fixtures',lambda:prepared)
    output=tmp_path/'fixtures.h5'
    monkeypatch.setattr(sys,'argv',['check_csr_backends.py','--phase','fixtures',
        '--build-dir',str(tmp_path),'--output',str(output)])
    original_path=list(sys.path)
    try:
        module.main()
        from radia_simulink.ams_exchange import read_exchange
        decoded=read_exchange(output)
    finally:
        sys.path[:]=original_path
    assert len(decoded['cases'])==2
    for (_,_,expected),actual in zip(module.fixtures(),decoded['cases']):
        assert set(actual)==set(expected)
        for key in expected:
            np.testing.assert_array_equal(actual[key],expected[key])
            assert actual[key].dtype==expected[key].dtype


def test_csr_fixture_storage_survives_ngsolve_matrix_owner_drop(tmp_path,monkeypatch):
    import gc
    pytest.importorskip('ngsolve')
    pytest.importorskip('h5py')
    spec=importlib.util.spec_from_file_location('csr_owned_fixture',
        ROOT/'validation_test/ngsolve_cuda_mex/check_csr_backends.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    owners=module.fixtures()
    expected=[{name:value.copy() for name,value in fields.items()} for _,_,fields in owners]
    snapshots=[fields for _,_,fields in owners]
    del owners
    gc.collect()
    # A standalone exported fixture must not retain borrowed CSR memory.
    assert all(value.flags.owndata for fields in snapshots for value in fields.values())
    monkeypatch.syspath_prepend(str(ROOT/'packages/radia-simulink/src'))
    from radia_simulink.ams_exchange import write_exchange,read_exchange
    path=tmp_path/'owned-fixtures.h5'
    write_exchange(path,dict(cases=snapshots))
    for actual,wanted in zip(read_exchange(path)['cases'],expected):
        for name in wanted:np.testing.assert_array_equal(actual[name],wanted[name])
