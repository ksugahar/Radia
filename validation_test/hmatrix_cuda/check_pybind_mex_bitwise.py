"""Compare separately produced same-kernel pybind11/MEX HDF5 records; never run GPU."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

CASES = {"device-flat-tet": {"apply", "linear_solution", "nonlinear_solution", "nonlinear_residual_history"}}
# Build and run identity both records must share; the kernel object hash binds both binaries to one nvcc object.
SHARED = ("source_sha", "kernel_object_sha256", "cuda_compiler", "python_binary_sha256", "mex_binary_sha256",
          "input_sha256", "device_uuid", "architecture", "case", "output_names")


def check_receipts(a: dict, b: dict) -> set:
    """Return the declared outputs after checking frontends, shared identity and loaded binaries."""
    if a.get("frontend") != "pybind11" or b.get("frontend") != "mex":
        raise ValueError("Expected independently produced pybind11 and mex records")
    for key in SHARED:
        if not a.get(key) or a[key] != b.get(key):
            raise ValueError(f"Missing or unequal receipt field: {key}")
    if a.get("artifact_sha256") != a["python_binary_sha256"] or b.get("artifact_sha256") != b["mex_binary_sha256"]:
        raise ValueError("A loaded frontend binary differs from the build receipt of the shared kernel object")
    if a["case"] not in CASES:
        raise ValueError("Unknown parity case")
    expected = CASES[a["case"]]
    if set(a["output_names"]) != expected or len(a["output_names"]) != len(expected):
        raise ValueError(f"Expected outputs {sorted(expected)}, got {a['output_names']}")
    return expected


def compare(python_path: Path, mex_path: Path) -> dict:
    import h5py

    with h5py.File(python_path, "r") as python, h5py.File(mex_path, "r") as mex:
        a = json.loads(python.attrs["receipt_json"])
        b = json.loads(mex.attrs["receipt_json"])
        expected = check_receipts(a, b)
        for name in sorted(expected):
            x, y = python["outputs"][name][...], mex["outputs"][name][...]
            if x.size == 0 or x.dtype.kind != "f" or x.dtype.itemsize != 8:
                raise ValueError(f"Expected nonempty FP64 output: {name}")
            if x.shape != y.shape or x.dtype != y.dtype or x.tobytes() != y.tobytes():
                raise ValueError(f"Bitwise mismatch: {name}; dtype {x.dtype}/{y.dtype}, shape {x.shape}/{y.shape}")
        return {"passed": True, "case": a["case"], "outputs": sorted(expected),
                "identity": {key: a[key] for key in SHARED if key not in ("case", "output_names")},
                "scope": "Recorded same-kernel frontend parity only; not CPU/GPU accuracy or performance"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-record", type=Path, required=True)
    parser.add_argument("--mex-record", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.python_record, args.mex_record)
    with args.report.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
