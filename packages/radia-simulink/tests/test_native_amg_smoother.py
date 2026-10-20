"""Small explicit-GPU ABI/numerical check for the shared AMG smoothing phase.

Select through RADIA_CSR_NATIVE_TEST_DIR; absence is an explicit hardware-lane
skip. A configured but missing/broken binary is a failure, never a CPU fallback.
"""
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pytest


def test_shared_amg_smoothing_native_gpu():
    configured = os.environ.get("RADIA_CSR_NATIVE_TEST_DIR")
    if not configured:
        pytest.skip("Optional native GPU lane requires RADIA_CSR_NATIVE_TEST_DIR")
    build = Path(configured)
    receipt = json.loads((build / "python-build.json").read_text())
    binary = build / receipt["binary"]
    import hashlib
    assert hashlib.sha256(binary.read_bytes()).hexdigest() == receipt["binary_sha256"]
    runtime = Path(os.environ["RADIA_VALIDATION_MATLAB_ROOT"])
    handle = os.add_dll_directory(str(runtime / "bin/win64")) if sys.platform == "win32" else None
    try:
        spec = importlib.util.spec_from_file_location("_csr_cuda", binary)
        module = importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        for n in (3, 64):
            dense = np.diag(np.full(n, 4.)) + np.diag(np.full(n-1, -1.), 1) + np.diag(np.full(n-1, -1.), -1)
            row, col = np.nonzero(dense)
            values = np.ascontiguousarray(dense[row, col])
            pointers = np.concatenate(([0], np.cumsum(np.bincount(row, minlength=n)))).astype(np.int32)
            columns = col.astype(np.int32)
            rhs = np.linspace(1., 2., n)
            weights = np.abs(dense).sum(axis=1)
            for sweeps in (1, 2, 5):
                expected = rhs / weights
                for _ in range(1, sweeps):
                    expected = expected + (rhs - dense @ expected) / weights
                actual = module._amg_l1_smooth(pointers, columns, values, rhs, weights, sweeps)
                np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-14)
            zero = module._amg_l1_smooth(pointers, columns, values, np.zeros(n), weights, 5)
            assert np.count_nonzero(zero) == 0
            with pytest.raises(TypeError):
                module._amg_l1_smooth(pointers.astype(np.int64), columns, values, rhs, weights, 1)
            with pytest.raises(RuntimeError):
                module._amg_l1_smooth(pointers, columns, values, rhs, np.zeros(n), 1)
            with pytest.raises(RuntimeError):
                module._amg_l1_smooth(pointers, columns, values, rhs, weights, 0)
            bad = columns.copy();bad[0] = n
            with pytest.raises(RuntimeError):
                module._amg_l1_smooth(pointers, bad, values, rhs, weights, 1)
        assert "matlab.engine" not in sys.modules
    finally:
        if handle is not None:
            handle.close()
