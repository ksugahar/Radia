"""Compare separately admitted same-kernel Python/MEX HDF5 records; never run GPU."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def validate_deterministic_flags(flags) -> None:
    required = {"--fmad=false", "--ftz=false", "--prec-div=true", "--prec-sqrt=true"}
    # This field contains only the four canonical controls, not the full command.
    # Reject duplicates and extra aliases/overrides rather than guessing precedence.
    if not isinstance(flags, list) or any(not isinstance(flag, str) for flag in flags):
        raise ValueError("Deterministic flags must be a list of canonical controls")
    if len(flags) != len(required) or set(flags) != required:
        raise ValueError("Exact deterministic compiler controls required; overrides forbidden")


def compare(python_path: Path, mex_path: Path) -> dict:
    import h5py

    keys = ("source_sha", "input_sha256", "cuda_compiler", "architecture",
            "device_uuid", "deterministic_flags", "case", "output_names")
    with h5py.File(python_path, "r") as python, h5py.File(mex_path, "r") as mex:
        a = json.loads(python.attrs["receipt_json"])
        b = json.loads(mex.attrs["receipt_json"])
        if a.get("frontend") != "pybind11" or b.get("frontend") != "mexcuda":
            raise ValueError("Expected independently produced pybind11 and mexcuda records")
        for key in keys:
            if not a.get(key) or a[key] != b.get(key):
                raise ValueError(f"Missing or unequal receipt field: {key}")
        validate_deterministic_flags(a["deterministic_flags"])
        if a["case"] == "normalized-leaf-resident-boundary":
            expected = {"apply"}
        elif a["case"] == "explicit-ng-host-hdiv-bh":
            expected = {"linear_solution", "cubic_solution", "cubic_residual_history"}
        else:
            raise ValueError("Unknown parity case")
        if set(a["output_names"]) != expected:
            raise ValueError(f"Expected outputs {sorted(expected)}, got {a['output_names']}")
        for name in expected:
            x, y = python["outputs"][name][...], mex["outputs"][name][...]
            if x.size == 0 or x.dtype.kind != "f" or x.dtype.itemsize != 8:
                raise ValueError(f"Expected nonempty FP64 output: {name}")
            if x.shape != y.shape or x.dtype != y.dtype or x.tobytes() != y.tobytes():
                raise ValueError(f"Bitwise mismatch: {name}; dtype {x.dtype}/{y.dtype}, shape {x.shape}/{y.shape}")
        return {"passed": True, "case": a["case"], "outputs": sorted(expected),
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
