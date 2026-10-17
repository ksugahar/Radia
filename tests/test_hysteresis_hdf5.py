"""HDF5-only branch contract; no native solver or private data needed."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import types

import h5py
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'validation_test/hysteresis/fixtures/potter_schmulian/B_input.hdf5'


def load_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_payload_and_direct_readers(monkeypatch):
    reader = load_source('hysteresis_hdf5_reader', ROOT / 'src/radia/hysteresis_io.py')
    proof = json.loads((ROOT / 'validation_test/hysteresis/hdf5_fixture_migration.json').read_text())
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == proof['hdf5_sha256']
    with h5py.File(FIXTURE, 'r') as data:
        for record in proof['numeric_arrays']:
            array = data[record['dataset']][...]
            assert list(array.shape) == record['shape']
            assert array.dtype.str == record['dtype']
            assert hashlib.sha256(array.tobytes(order='C')).hexdigest() == record['raw_sha256']
    loops, db, levels = reader.load_hdf5(FIXTURE)
    assert len(loops) == 20 and db == 0.095 and levels.shape == (20,)
    # Inject only the explicitly loaded Python reader; no native solver is tested.
    monkeypatch.setitem(sys.modules, 'radia', types.ModuleType('radia'))
    monkeypatch.setitem(sys.modules, 'radia.hysteresis_io', reader)
    bqm = load_source('hysteresis_bqm_reader', ROOT / 'validation_test/hysteresis/bqm_model.py')
    direct = bqm.load_binput_data(FIXTURE)
    with h5py.File(FIXTURE, 'r') as data:
        for i, (h, b) in enumerate(direct['BH']):
            np.testing.assert_array_equal(h, data[f'BH/{i:04d}/H'][...].flatten())
            np.testing.assert_array_equal(b, data[f'BH/{i:04d}/B'][...].flatten())
    assert not hasattr(reader, 'load_mat')


def test_rejects_unsupported_schema(tmp_path):
    reader = load_source('hysteresis_hdf5_reader_invalid', ROOT / 'src/radia/hysteresis_io.py')
    path = tmp_path / 'unsupported.hdf5'
    with h5py.File(path, 'w') as data:
        data.attrs['schema'] = 'unknown'
    with pytest.raises(ValueError, match='schema'):
        reader.load_hdf5(path)
