"""Protect AMS numeric layout/type and the versioned HDF5 boundary."""
import importlib.util
from pathlib import Path
import json
import h5py
import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("ams_exchange", Path(__file__).resolve().parents[1] / "src/radia_simulink/ams_exchange.py")
exchange = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exchange)


def test_nested_hdf5_preserves_indices_layout_and_empty_arrays(tmp_path):
    payload = {"snapshot": {"schema": "radia.ams-accelerator.v1", "beta_zero": False,
        "spaces": [{"levels": [{"rows": 3, "columns": np.array([0, 2], dtype=np.int32)}]}]},
        "probes": np.arange(15, dtype=np.float64).reshape(3, 5),
        "empty": np.empty(0, dtype=np.int32)}
    path = tmp_path / "input.hdf5"
    receipt = exchange.write_exchange(path, payload)
    actual = exchange.read_exchange(path)
    assert receipt["format"] == "HDF5"
    assert actual["snapshot"]["beta_zero"] is False
    indices = actual["snapshot"]["spaces"][0]["levels"][0]["columns"]
    assert indices.dtype == np.int32
    np.testing.assert_array_equal(indices, [0, 2])
    np.testing.assert_array_equal(actual["probes"], payload["probes"])
    assert actual["empty"].shape == (0,)
    assert len(exchange.file_sha256(path)) == 64
    with pytest.raises(FileExistsError):
        exchange.write_exchange(path, payload)


def test_hdf5_contract_rejects_unsupported_arrays_and_corrupt_shapes(tmp_path):
    with pytest.raises(ValueError, match="FP64/int32/uint8"):
        exchange.write_exchange(tmp_path / "bad.hdf5", {"values": np.ones(3, dtype=np.float32)})
    path = tmp_path / "checked.hdf5"
    exchange.write_exchange(path, {"values": np.ones((2, 3))})
    with h5py.File(path, "r+") as handle:
        manifest = json.loads(handle["manifest"][...].tobytes())
        manifest["tree"]["fields"]["values"]["shape"] = [3, 2]
        del handle["manifest"]
        handle.create_dataset("manifest", data=np.frombuffer(json.dumps(manifest).encode(), dtype=np.uint8))
    with pytest.raises(ValueError, match="shape/type mismatch"):
        exchange.read_exchange(path)
