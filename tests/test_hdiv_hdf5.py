"""Typed independent CPU interchange contracts; no GPU or MATLAB launch."""
import importlib.util
from pathlib import Path
import h5py
import numpy as np
import pytest
from scipy import sparse

spec = importlib.util.spec_from_file_location('hdiv_hdf5', Path(__file__).parents[1] / 'validation_test/hmatrix_cuda/hdiv_hdf5.py')
io = importlib.util.module_from_spec(spec)
spec.loader.exec_module(io)


def test_nonsquare_complex_csr_and_utf8_roundtrip(tmp_path):
    fields = {'matrix': np.arange(6, dtype=np.float64).reshape(2, 3),
              'operator': sparse.csr_matrix([[1+2j, 0, 3], [0, 4-1j, 0]]),
              'mask': np.array([True, False])}
    metadata = {'units': 'SI', 'labels': ['boundary', '温度']}
    path = tmp_path / 'fixture.h5'
    io.write_hdiv_h5(path, fields, metadata)
    result, actual = io.read_hdiv_h5(path)
    assert actual == metadata
    for name in ('matrix', 'mask'):
        np.testing.assert_array_equal(result[name], fields[name])
    np.testing.assert_array_equal(result['operator'].toarray(), fields['operator'].toarray())
    with pytest.raises(FileExistsError):
        io.write_hdiv_h5(path, fields, metadata)


def test_corruption_budget_and_format_fail_loudly(tmp_path):
    path = tmp_path / 'fixture.hdf5'
    io.write_hdiv_h5(path, {'values': np.arange(8, dtype=np.float64)}, {'units': 'SI', 'labels': []})
    with pytest.raises(ValueError, match='budget'):
        io.read_hdiv_h5(path, max_bytes=8)
    with h5py.File(path, 'r+') as f:
        f['fields/values/real'][0] = 99
    with pytest.raises(ValueError, match='SHA256'):
        io.read_hdiv_h5(path)
    with pytest.raises(ValueError, match='HDF5 only'):
        io.read_hdiv_h5(tmp_path / 'fixture.mat')


def test_rank_three_is_explicitly_unsupported(tmp_path):
    with pytest.raises(ValueError, match='rank-one/two'):
        io.write_hdiv_h5(tmp_path / 'rank3.h5',
                         {'values': np.zeros((2, 3, 4), dtype=np.float64)},
                         {'units': 'SI', 'labels': []})
