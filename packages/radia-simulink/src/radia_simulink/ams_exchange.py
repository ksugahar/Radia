"""Versioned plain-HDF5 numeric exchange for CPU-owned AMS snapshots.

This is HDF5, not MATLAB's proprietary MAT v7.3 struct encoding. The matching
MATLAB reader is radia.gpu.readAMSExchange. Arrays keep their NumPy shapes,
FP64 values and zero-based int32 indices; the reader restores dimension order.
All sizes use this format, including artifacts exceeding 2 GB.
"""
from __future__ import annotations

import json
from pathlib import Path

SCHEMA = "radia.ams-hdf5.v1"


def write_exchange(path, payload):
    """Write without collecting array data into a serialized in-memory blob."""
    import h5py
    import numpy as np

    path = Path(path)
    if path.suffix.lower() not in (".h5", ".hdf5"):
        raise ValueError("AMS exchange requires an .h5 or .hdf5 file")
    count = 0
    numeric_bytes = 0
    with h5py.File(path, "x") as handle:
        arrays = handle.create_group("arrays")

        def node(value):
            nonlocal count, numeric_bytes
            if isinstance(value, dict):
                if any(not isinstance(k, str) or not k.isidentifier() for k in value):
                    raise ValueError("Exchange fields must be valid identifiers")
                return {"kind": "struct", "fields": {k: node(v) for k, v in value.items()}}
            if isinstance(value, (list, tuple)):
                return {"kind": "cell", "items": [node(v) for v in value]}
            if isinstance(value, np.ndarray):
                if value.dtype not in (np.dtype("float64"), np.dtype("int32"), np.dtype("uint8")) or not 1 <= value.ndim <= 2:
                    raise ValueError("AMS exchange supports FP64/int32/uint8 vectors and matrices")
                name = f"a{count:06d}"; count += 1
                # h5py streams directly from contiguous native snapshot buffers.
                # No compression, whole-payload copy, or MAT 2 GB length field.
                arrays.create_dataset(name, data=value)
                numeric_bytes += value.nbytes
                return {"kind": "array", "path": "/arrays/" + name,
                        "shape": list(value.shape), "dtype": value.dtype.name}
            if isinstance(value, (bool, np.bool_)):
                return {"kind": "logical", "value": bool(value)}
            if isinstance(value, (int, float, np.integer, np.floating)):
                return {"kind": "number", "value": float(value)}
            if isinstance(value, str):
                return {"kind": "string", "value": value}
            raise TypeError(f"Unsupported exchange type: {type(value).__name__}")

        manifest = {"schema": SCHEMA, "index_base": 0, "precision": "FP64",
                    "array_order": "numpy", "units": "SI", "tree": node(payload)}
        encoded = json.dumps(manifest, allow_nan=False, separators=(",", ":")).encode("utf-8")
        handle.create_dataset("manifest", data=np.frombuffer(encoded, dtype=np.uint8))
    return {"schema": SCHEMA, "format": "HDF5", "numeric_bytes": numeric_bytes,
            "file_bytes": path.stat().st_size, "arrays": count}


def read_exchange(path):
    """Read the same schema for diagnostics and paired timing input identity."""
    import h5py
    import numpy as np

    with h5py.File(path, "r") as handle:
        manifest = json.loads(handle["manifest"][...].tobytes().decode("utf-8"))
        if (manifest.get("schema") != SCHEMA or manifest.get("index_base") != 0
                or manifest.get("precision") != "FP64" or manifest.get("array_order") != "numpy"):
            raise ValueError("Unsupported AMS HDF5 exchange contract")

        def node(spec):
            kind = spec["kind"]
            if kind == "struct":
                return {k: node(v) for k, v in spec["fields"].items()}
            if kind == "cell":
                return [node(v) for v in spec["items"]]
            if kind == "array":
                name = spec["path"]
                if not name.startswith("/arrays/a") or not name[9:].isdigit():
                    raise ValueError("Invalid AMS dataset path")
                dataset = handle[name]
                if list(dataset.shape) != spec["shape"] or dataset.dtype != np.dtype(spec["dtype"]):
                    raise ValueError("AMS dataset shape/type mismatch")
                return dataset[...]
            if kind in ("number", "logical", "string"):
                return spec["value"]
            raise ValueError("Unsupported AMS HDF5 node")
        return node(manifest["tree"])


def file_sha256(path):
    """Hash even multi-GB inputs using bounded memory."""
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
