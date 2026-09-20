"""Compare two standalone builds of the real C ACA routine, without pybind11."""

import argparse
import ctypes
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np


def evaluate(library, matrix, eps=1e-10):
    pointer = ctypes.POINTER(ctypes.c_double)
    function = library.aca_probe
    function.argtypes = [pointer, ctypes.c_int, ctypes.c_int, ctypes.c_double,
                         pointer, ctypes.POINTER(ctypes.c_int)]
    function.restype = ctypes.c_int
    matrix = np.ascontiguousarray(matrix, dtype=float)
    output = np.zeros_like(matrix)
    calls = ctypes.c_int()
    started = time.perf_counter()
    rank = function(matrix.ctypes.data_as(pointer), *matrix.shape, eps,
                    output.ctypes.data_as(pointer), ctypes.byref(calls))
    elapsed = time.perf_counter() - started
    if rank < 0:
        raise RuntimeError("native ACA probe failed")
    scale = max(float(np.max(np.abs(matrix))), np.finfo(float).tiny)
    # Normalize before norms to keep scale tests clear of square underflow.
    denominator = np.linalg.norm(matrix / scale)
    error = np.linalg.norm((output - matrix) / scale)
    return {"rank": rank, "relative_error": float(error / denominator if denominator else error),
            "entry_evaluations": calls.value, "elapsed_seconds": elapsed}


def cases():
    hidden = np.array([[1., 1., 1.], [1., 1., 2.], [1., 1., 1.]])
    for scale in (1e-80, 1., 1e80):
        for permutation in itertools.permutations(range(3)):
            yield f"hidden-{scale:g}-{permutation}", hidden[list(permutation)] * scale
    yield "zero", np.zeros((7, 11))
    initial_zero = np.zeros((5, 9))
    initial_zero[2] = np.arange(1., 10.)
    yield "initial-zero-row", initial_zero
    rng = np.random.default_rng(782)
    for m, n in ((8, 13), (13, 8), (64, 48)):
        yield f"signed-rank2-{m}-{n}", rng.normal(size=(m, 2)) @ rng.normal(size=(2, n))
        yield f"full-rank-{m}-{n}", rng.normal(size=(m, n))
    # Smooth separated 3-D Laplace leaf, with and without regularization.
    t = np.linspace(-.5, .5, 6)
    points = np.array(list(itertools.product(t, repeat=3)))
    delta = points[:, None, :] - (points[None, :, :] + [3., 0., 0.])
    for epsilon in (0., .03, .3):
        yield f"laplace-leaf-{epsilon}", 1 / np.sqrt(np.sum(delta**2, axis=2) + epsilon**2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before, after = ctypes.CDLL(str(args.before)), ctypes.CDLL(str(args.after))
    records = []
    for name, matrix in cases():
        records.append({"case": name, "shape": list(matrix.shape),
                        "before": evaluate(before, matrix), "after": evaluate(after, matrix)})
    accepted = all(r["after"]["relative_error"] < 1e-9 for r in records)
    sources = [Path(__file__), Path(__file__).with_name("aca_probe_bridge.c"),
               Path(__file__).parents[2] / "src/ext/HACApK/cHACApK_base.c"]
    result = {"schema": "radia.aca-residual-guard.v1", "hostname": platform.node(),
              "python": platform.python_version(), "numpy": np.__version__,
              "aca_tolerance": 1e-10, "acceptance_budget": 1e-9,
              "production_hmatrix_rebuilt": False,
              "provenance": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in [args.before, args.after, *sources]},
              "records": records, "all_candidate_cases_pass": accepted}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for r in records:
        print(f"{r['case']}: {r['before']['relative_error']:.3e} -> "
              f"{r['after']['relative_error']:.3e}; entries "
              f"{r['before']['entry_evaluations']} -> {r['after']['entry_evaluations']}")
    if not accepted:
        raise SystemExit("candidate accuracy gate failed")


if __name__ == "__main__":
    main()
