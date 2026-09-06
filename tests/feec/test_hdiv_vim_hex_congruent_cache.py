"""Translation-congruent shared cache of the HEX charge-Gram near blocks.

On a swept (extruded) HEX mesh every host is repeated once per layer, so each expensive near block
(pair-domain Duffy / near-band product, ~265 ms on hibino) recurs once per layer.  GetHexBlock keys
the instance-shared cache by (host templates, quantized centre offset) instead of host ids, so the
copies are served from one evaluation.  This test locks (a) the cache engages on a NON-uniform
swept mesh (the older uniform-lattice cache cannot: the x-spacing is graded), (b) the Gram it
produces equals the uncached Gram to round-off, and (c) the environment latch turns it off.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

import radia  # noqa: E402
from radia.vim import _vim as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def _swept_mesh(nz=4):
    """Graded x-spacing (non-uniform, so no global lattice), uniform layers along z."""
    return MakeStructured3DMesh(
        hexes=True, nx=3, ny=2, nz=nz,
        mapping=lambda x, y, z: (x + 0.35 * x * x, 0.8 * y, 0.5 * z))


def _dense_gram(mesh):
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M, _ = V._build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False)
        import scipy.sparse as sp
        n = sp.csr_matrix(B).shape[0]
        Gd = np.array([[G.entry(i, j) for j in range(n)] for i in range(n)])
        stats = dict(G.stats())
    return Gd, stats


def test_congruent_cache_engages_and_reproduces_the_uncached_gram(tmp_path):
    Gd, stats = _dense_gram(_swept_mesh())
    assert stats["hex_congruent_cache_enabled"] == 1.0
    assert stats["hex_congruent_ready"] == 1.0
    n_hosts = stats["hex_congruent_hosts"]
    # 4 layers of 6 cells: the 6 cross-section cells give 6 cell templates (plus the end/side faces),
    # far fewer than the 24 cells + boundary faces.
    assert 0 < stats["hex_congruent_templates"] < 0.6 * n_hosts
    assert stats["hex_general_shared_lookups"] > 0
    assert stats["hex_general_shared_hits"] > 0.3 * stats["hex_general_shared_lookups"]
    assert stats["nondefault_performance_override_active"] == 0.0

    # Same Gram with the cache latched off in a fresh interpreter (the latch is read once per process).
    out = tmp_path / "uncached.npy"
    script = (
        "import sys, numpy as np, ngsolve as ng\n"
        "from ngsolve.meshes import MakeStructured3DMesh\n"
        "from radia.vim import _vim as V\n"
        "import scipy.sparse as sp\n"
        "mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=2, nz=4, "
        "mapping=lambda x, y, z: (x + 0.35 * x * x, 0.8 * y, 0.5 * z))\n"
        "with ng.TaskManager():\n"
        "    fes = ng.HDiv(mesh, order=1)\n"
        "    B, G, M, _ = V._build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False)\n"
        "    n = sp.csr_matrix(B).shape[0]\n"
        "    Gd = np.array([[G.entry(i, j) for j in range(n)] for i in range(n)])\n"
        "    st = dict(G.stats())\n"
        "assert st['hex_congruent_cache_enabled'] == 0.0 and st['nondefault_performance_override_active'] == 1.0, st\n"
        f"np.save(r'{out}', Gd)\n"
    )
    env = dict(os.environ)
    env["RADIA_HDIV_DISABLE_CONGRUENT_CACHE"] = "1"
    env["PYTHONPATH"] = str(Path(radia.__file__).resolve().parents[1]) + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run([sys.executable, "-c", script], check=True, env=env, cwd=str(REPO))
    Gu = np.load(out)
    scale = np.abs(np.diag(Gd)).max()
    assert np.abs(Gd - Gu).max() < 1e-12 * scale, np.abs(Gd - Gu).max() / scale
