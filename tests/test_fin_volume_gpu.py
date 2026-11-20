"""CPU/GPU parity of the device-resident fin-volume PEEC solve.

Guards: the GPU route (L assembled, energy-gated, factorised and solved on
the device) reproduces the CPU route on identical bars: the same shared
oriented-box kernel, the same nodal formulation and terminals.
Assumes: the verified CPU route is the reference; device math differs from
the host CRT by a few ulp per entry, so agreement is to 1e-9 relative.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).resolve().parent / "coil_from_cad" / "fixtures"
KW = dict(frequency=222e3, sigma=5.8e7, current_A=2.0, n_layers=2, n_lanes=16,
          n_stations=7, lane_grading="uniform")


@pytest.fixture(scope="module")
def pair():
    if os.environ.get("RADIA_PEEC_FIN_GPU_TEST") != "1":
        pytest.skip("Optional GPU lane requires RADIA_PEEC_FIN_GPU_TEST=1")
    from radia.fin_volume import solve_fin_volume

    step = FIXTURES / "beak_fin_short.step"
    return solve_fin_volume(step, backend="cpu", **KW), solve_fin_volume(step, backend="gpu", **KW)


def test_gpu_matches_cpu_on_identical_bars(pair):
    cpu, gpu = pair
    meta = gpu["gpu"]
    assert gpu["backend"] == "gpu" and meta["backend"] == "gpu", meta
    assert meta["true_relative_residual"] <= 1e-6, meta
    for key in ("R_ohm", "L_H"):
        assert abs(gpu[key] - cpu[key]) <= 1e-9 * abs(cpu[key]), (key, gpu[key], cpu[key])
    di = np.linalg.norm(gpu["current"] - cpu["current"]) / np.linalg.norm(cpu["current"])
    assert di <= 1e-9, di
    assert meta["cholesky"] == cpu["inductance_energy"]["cholesky"], (meta, cpu["inductance_energy"])


def test_gpu_route_moves_no_matrix(pair):
    _, gpu = pair
    m, n = gpu["n_branches"], gpu["n_nodes"]
    meta = gpu["gpu"]
    # inputs: 13 doubles + 1 byte + 2 ints per bar, one int per node, rules, injection
    assert meta["bytes_to_device"] < 200 * (m + n) + 64 * 1024, meta
    assert meta["bytes_to_host"] < 64 * (m + n) + 64 * 1024, meta


def test_backend_is_explicit():
    from radia.fin_volume import solve_fin_volume

    with pytest.raises(ValueError, match="backend must be"):
        solve_fin_volume(FIXTURES / "beak_fin_short.step", backend="auto", **KW)
