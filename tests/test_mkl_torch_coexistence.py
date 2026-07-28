"""Windows oneMKL/PyTorch runtime coexistence regression tests."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO_ROOT / "src"), env.get("PYTHONPATH", "")]
    )
    env.pop("KMP_DUPLICATE_LIB_OK", None)
    return env


@pytest.mark.skipif(sys.platform != "win32", reason="Windows DLL-runtime test")
def test_torch_first_peec_solve_uses_tbb_without_openmp_abort():
    if importlib.util.find_spec("torch") is None:
        pytest.skip("PyTorch is not installed")

    script = r"""
import os
import torch
import numpy as np

from radia.coil_from_cad import build_peec_from_path
from radia.peec_topology import PEECCircuitSolver

assert os.environ["MKL_THREADING_LAYER"] == "TBB"
path = np.array([
    [0.0, 0.0, 0.0],
    [0.01, 0.0, 0.0],
    [0.01, 0.01, 0.0],
])
topology = build_peec_from_path(
    path,
    np.full(2, 1.0e-3),
    np.full(2, 1.0e-3),
)
impedance = PEECCircuitSolver(topology).compute_port_impedance(1.0e5)
assert np.isfinite([impedance.real, impedance.imag]).all()
"""
    env = _subprocess_env()
    env.pop("MKL_THREADING_LAYER", None)
    proc = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert proc.returncode == 0, (
        f"PyTorch-first PEEC subprocess exited {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )


@pytest.mark.skipif(sys.platform != "win32", reason="Windows DLL-runtime test")
def test_explicit_mkl_threading_layer_is_preserved():
    env = _subprocess_env()
    env["MKL_THREADING_LAYER"] = "SEQUENTIAL"
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "import os, radia; "
            "assert os.environ['MKL_THREADING_LAYER'] == 'SEQUENTIAL'",
        ],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stderr
