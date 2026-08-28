"""Generate the native MATLAB Sobol direction-number payload.

The input is SciPy 1.17.1's checked Joe--Kuo criterion-6 table.  The output is
a deterministic little-endian binary file so QMCSampler can use all 21,201
dimensions without importing Python or SciPy at runtime.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np
import scipy
import scipy.stats._sobol as scipy_sobol

EXPECTED_SCIPY_VERSION = "1.17.1"
MAGIC = b"RADIA_SOBOL_JK6"
SCHEMA_VERSION = 1
MAX_DIMENSION = 21_201
MAX_DEGREE = 18


def build_payload() -> bytes:
    if scipy.__version__ != EXPECTED_SCIPY_VERSION:
        raise RuntimeError(
            f"Expected scipy=={EXPECTED_SCIPY_VERSION}, found {scipy.__version__}."
        )
    source = Path(scipy_sobol.__file__).with_name("_sobol_direction_numbers.npz")
    with np.load(source) as archive:
        poly = np.asarray(archive["poly"], dtype="<u4")
        vinit = np.asarray(archive["vinit"], dtype="<u4")
    if poly.shape != (MAX_DIMENSION,) or vinit.shape != (
        MAX_DIMENSION,
        MAX_DEGREE,
    ):
        raise RuntimeError(
            f"Unexpected SciPy Sobol table shapes: poly={poly.shape}, "
            f"vinit={vinit.shape}."
        )
    header = struct.pack(
        "<16sIII", MAGIC, SCHEMA_VERSION, MAX_DIMENSION, MAX_DEGREE
    )
    return header + poly.tobytes(order="C") + vinit.tobytes(order="C")


def main() -> None:
    default = (
        Path(__file__).resolve().parents[3]
        / "matlab"
        / "+radia"
        / "+optuna"
        / "+internal"
        / "sobol_direction_numbers.bin"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=default)
    args = parser.parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(payload)
    print(f"{args.output} ({len(payload)} bytes)")


if __name__ == "__main__":
    main()
