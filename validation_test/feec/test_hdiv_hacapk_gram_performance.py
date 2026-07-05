"""Validation gate for the HDiv-VIM HACApK charge-Gram build.

Correctness tests already lock entries, spectra, and symmetric matvecs.  This
file adds the missing production guard: the build path must expose HACApK stats,
enter the compressed/low-rank regime on a modest hex grid, and stay within a
generous wall-time envelope so accidental dense rebuilds fail loudly.
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import Solve  # noqa: E402


def _hex_cube(n):
    return MakeStructured3DMesh(
        hexes=True, nx=n, ny=n, nz=n,
        mapping=lambda x, y, z: (0.02 * x, 0.02 * y, 0.02 * z),
    )


@pytest.mark.flaky(reruns=2, reruns_delay=1)
def test_hex_chargegram_hacapk_build_enters_lowrank_regime():
    """A 4^3 hex RT1 cube is large enough for HACApK ACA leaves and still small enough for validation."""
    with ng.TaskManager():
        res = Solve(_hex_cube(4), mu_r=1000.0,
                    H_ext=ng.CoefficientFunction((0.0, 0.0, 1.0e4)),
                    gram_eps=1.0e-4)
    st = res["hmat_stats"]
    assert res["n_charge"] == st["n_dof"]
    assert st["n_lowrank"] > 0, f"HACApK never entered low-rank mode: {st}"
    assert st["compression"] < 0.98, f"HACApK compression unexpectedly weak: {st}"
    assert st["memory_mb"] <= st["dense_memory_mb"], f"HACApK memory stats inverted: {st}"
    assert st["build_time"] < 30.0, f"charge-Gram build time looks like a dense/uncached regression: {st}"
    assert np.isfinite(res["demag"]) and 0.30 < res["demag"] < 0.36


@pytest.mark.flaky(reruns=2, reruns_delay=1)
def test_hex_chargegram_build_stats_scale_monotonically():
    """The validation lane records build stats on two sizes without using a brittle exact timing ratio."""
    rows = []
    with ng.TaskManager():
        for n in (2, 4):
            res = Solve(_hex_cube(n), mu_r=1000.0,
                        H_ext=ng.CoefficientFunction((0.0, 0.0, 1.0e4)),
                        gram_eps=1.0e-4)
            rows.append((res["n_charge"], res["hmat_stats"]))
    assert rows[1][0] > rows[0][0]
    assert rows[1][1]["dense_memory_mb"] > rows[0][1]["dense_memory_mb"]
    assert rows[1][1]["build_time"] < 40.0
