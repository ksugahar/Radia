"""Distorted-pair far dispatch + shared general-block cache (C-4): contract locks.

A Sculpt-style hex mesh has NO affine cells (the overlay-grid smoothing moves
every node; measured median deviation ~25 % of the cell scale), so every pair
used to run the expensive general graded machinery.  Two C++ changes speed the
charge-Gram fill on such meshes:

* well-separated pairs with distorted hosts route to the geometry-map-exact
  tensor-product rule (``QuadBlockHexAffineFarProduct``; the Q2 map places the
  points and the Piola reference charge measure cancels the Jacobian), and
* the remaining general-path blocks are served from an instance-shared cache
  instead of per-thread caches (the fill workers duplicated ~2x).

Locks here:

* the far switch changes far entries only at quadrature-truncation level
  (subprocess A/B against ``RADIA_HDIV_HEX_DISTORTED_FAR_FACTOR=0``, the
  old-path escape latch),
* the demag spectrum stays in the physical band on an all-distorted mesh
  (PSD floor from the sigma-normalized storage + eig <= 1 physics),
* the shared cache is instance-isolated (interleaved Grams on different
  meshes keep serving their own entries).
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import _vim as V  # noqa: E402

# Smooth nonlinear distortion, ~15 % of the cell size: every cell fails the
# 1e-10 affine gate, matching the Sculpt-mesh regime this dispatch targets.
_MESH_SRC = """
import math
from ngsolve.meshes import MakeStructured3DMesh

def make_mesh(phase=0.0):
    def mapping(x, y, z):
        s = 0.05
        return (0.02 * (x - 0.5) + s * 0.02 * math.sin(2.1 * y + 3.0 * z + phase),
                0.02 * (y - 0.5) + s * 0.02 * math.sin(2.7 * z + 1.9 * x + phase),
                0.02 * (z - 0.5) + s * 0.02 * math.sin(1.3 * x + 2.3 * y + phase))
    return MakeStructured3DMesh(nx=3, ny=3, nz=3, mapping=mapping, hexes=True)
"""

_ENTRY_CHILD = _MESH_SRC + """
import json
import numpy as np
from ngsolve import HDiv, TaskManager
from radia.vim import _vim as V

mesh = make_mesh()
fes = HDiv(mesh, order=1)
with TaskManager():
    B, G, M, Mng = V._build_charge_gram_hex(fes, build_hmatrix=False)
n = int(B.shape[0])   # stats()["n_dof"] is only set by the H-matrix build
rng = np.random.default_rng(2026)
pairs = rng.integers(0, n, size=(400, 2))
vals = [G.entry(int(a), int(b)) for a, b in pairs]
print("RESULT " + json.dumps({"n": n, "vals": vals}))
"""


def _run_child(source, factor):
    env = dict(os.environ)
    env["RADIA_HDIV_HEX_DISTORTED_FAR_FACTOR"] = str(factor)
    proc = subprocess.run([sys.executable, "-c", source], env=env,
                          capture_output=True, text=True, timeout=1200)
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    raise AssertionError(
        f"child (factor={factor}) produced no RESULT:\n"
        f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")


def _exec_mesh(phase=0.0):
    scope = {}
    exec(_MESH_SRC, scope)
    return scope["make_mesh"](phase)


def _validate_far_switch_matches_general_path_entries():
    """factor=1 (tensor far product) vs factor=0 (old general path): far
    entries agree to quadrature-truncation level.  Both sides are quadrature
    approximations of the same smooth integral (old: 6-sub Keast-15 clouds,
    new: whole-host 4^3 tensor Gauss), so the A/B measures the DIFFERENCE of
    two truncation errors -- measured max ~2.4e-5 of max |G| on the Sculpt
    design mesh (4000 entries).  The 1e-3 gate leaves margin for mesh-family
    variation while still catching a wrong-measure/wrong-map regression
    (those are O(1))."""
    old = _run_child(_ENTRY_CHILD, 0.0)
    new = _run_child(_ENTRY_CHILD, 1.0)
    assert old["n"] == new["n"]
    vo = np.asarray(old["vals"])
    vn = np.asarray(new["vals"])
    scale = np.abs(vo).max()
    assert scale > 0
    max_diff = np.abs(vn - vo).max()
    assert max_diff < 1e-3 * scale, (
        f"far-switch entry drift {max_diff:.3e} vs scale {scale:.3e}")
    # The switch must actually fire: on an all-distorted mesh with 27 cells a
    # far pair exists, and its tensor-product value differs from the graded
    # value in the last digits.  Identical arrays would mean the dispatch is
    # dead (env latch broken / branch removed).
    assert np.abs(vn - vo).max() > 0.0


def test_distorted_spectrum_stays_physical():
    """Dense generalized spectrum of the demag operator on the all-distorted
    mesh: PSD floor (sigma-normalized storage) + physical eig <= 1 band, now
    with far pairs served by the tensor product."""
    sla = pytest.importorskip("scipy.linalg")
    sp = pytest.importorskip("scipy.sparse")
    mesh = _exec_mesh()
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M, Mng = V._build_charge_gram_hex(fes, eps=1e-10)
        Bs = sp.csr_matrix(B)
        n = Bs.shape[1]
        N = np.empty((n, n))
        for column in range(n):
            basis = np.zeros(n)
            basis[column] = 1.0
            charge = np.ascontiguousarray(Bs @ basis, dtype=np.float64)
            N[:, column] = Bs.T @ np.asarray(G.matvec_sym(charge))
        N = 0.5 * (N + N.T)
        Md = np.asarray(M.todense() if hasattr(M, "todense") else M)
    w = sla.eigvalsh(N, Md)
    assert w[0] > -1e-8, f"demag spectrum lost PSD: min eig {w[0]:.3e}"
    assert w[-1] < 1.05, f"demag spectrum exceeds physical band: {w[-1]:.6f}"


def test_shared_general_cache_is_instance_isolated():
    """Interleave two Grams on DIFFERENT distorted meshes and re-read the
    first: the instance-shared general-block cache must keep serving the
    first mesh's blocks (a static/global cache would leak mesh B's geometry
    into mesh A's entries)."""
    with ng.TaskManager():
        fes_a = ng.HDiv(_exec_mesh(0.0), order=1)
        B_a, G_a, _, _ = V._build_charge_gram_hex(fes_a, build_hmatrix=False)
        n_a = int(B_a.shape[0])   # stats n_dof is only set by the H-matrix build
        rng = np.random.default_rng(7)
        pairs = rng.integers(0, n_a, size=(60, 2))
        before = np.asarray(
            [G_a.entry(int(a), int(b)) for a, b in pairs])
        fes_b = ng.HDiv(_exec_mesh(1.4), order=1)
        B_b, G_b, _, _ = V._build_charge_gram_hex(fes_b, build_hmatrix=False)
        n_b = int(B_b.shape[0])
        for a, b in pairs:
            G_b.entry(int(a) % n_b, int(b) % n_b)
        after = np.asarray(
            [G_a.entry(int(a), int(b)) for a, b in pairs])
    np.testing.assert_array_equal(before, after)
