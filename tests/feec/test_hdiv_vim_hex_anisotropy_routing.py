"""Flat HEX pairs must leave the exponentially convergent BDM1 pair rules.

The pair-domain Duffy rule and the plain pair product rule converge
exponentially on well-shaped hosts, but each sees a feature of relative width
1/AR on a host of aspect ratio AR -- the gap in the product rule, the angular
1/X of the cone in the Duffy rule -- so their rate is divided by AR.  Measured
2026-09-20 against analytic-reduced references (constant-monomial entries of a
single a x a x a/AR hex, a = 0.1 m), relative error:

    AR                         1      3      8     20     100
    side-face self (Duffy)   4e-9   8e-6   3e-4   1e-2    9e-2
    opposite faces (product) 2e-9   4e-5   9e-3   1e-1    1e+0
    graded near tensor       2e-6   6e-5   8e-5   2e-5    5e-5

Those block errors cost the assembled Gram its definiteness: with the Duffy
family on every near pair the smallest eigenvalue of one flattened unit hex ran
+2.5e-7 at AR 12 and -9.2e-5 at AR 14, and a 100 x 100 x 1 mm shield cell
reached -3.46 against a largest eigenvalue of 26.8, which broke the shield solve
down in CG at iteration 11.  Above the measured boundary the graded near tensor
rule is used instead; it restores definiteness and is flat in AR.

The threshold sits at the definiteness boundary rather than at the accuracy
crossover (~AR 3) on purpose: between 3 and 12 the Duffy family is the less
accurate one but still definite, and it is the only family with a consistent
shape derivative, so meshes in that band stay differentiable.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh

DEFINITE_LIMIT = 12.0


def _gram(aspect_ratio, nx=1, ny=1):
    from radia.vim._vim import _charge_basis_hex, build_charge_gram

    mesh = MakeStructured3DMesh(
        hexes=True, nx=nx, ny=ny, nz=1,
        mapping=lambda x, y, z: (x, y, z / aspect_ratio))
    fes = ng.HDiv(mesh, order=1)
    with ng.TaskManager():
        _charge_basis_hex(fes, cob_quad=3)
        _, gram, _ = build_charge_gram(fes, eps=1e-10, leafsize=256, eta=2.0)
    n = int(gram.stats()["n_dof"])
    dense = np.array([[gram.entry(i, j) for j in range(n)] for i in range(n)])
    return gram, 0.5 * (dense + dense.T)


@pytest.mark.parametrize("aspect_ratio", [1.0, 8.0, 12.0])
def test_up_to_the_definiteness_limit_nothing_is_routed(aspect_ratio):
    gram, dense = _gram(aspect_ratio)
    stats = gram.stats()
    assert stats["hex_exponential_rule_max_anisotropy"] == DEFINITE_LIMIT
    assert stats["hex_pair_anisotropy_graded"] == 0.0, (
        "a pair at or below the measured definiteness limit must stay on the "
        "Duffy family, which is the only one with a consistent shape "
        "derivative")
    assert np.linalg.eigvalsh(dense)[0] > 0.0


@pytest.mark.parametrize("aspect_ratio", [14.0, 20.0, 100.0])
def test_beyond_the_limit_routing_restores_definiteness(aspect_ratio):
    gram, dense = _gram(aspect_ratio)
    assert gram.stats()["hex_pair_anisotropy_graded"] > 0.0
    smallest = np.linalg.eigvalsh(dense)[0]
    assert smallest > 0.0, (
        f"flat hex AR {aspect_ratio} still indefinite: lambda_min {smallest:.3e}")


def test_the_duffy_family_is_what_loses_definiteness(monkeypatch):
    """The routing is load-bearing, not incidental.

    Raising the threshold puts the flat pair back on the Duffy family, which is
    the configuration that produced the negative eigenvalues; the knob exists
    so that comparison runs on one binary.
    """
    monkeypatch.setenv("RADIA_HDIV_HEX_MAX_ANISOTROPY", "1e9")
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    script = (
        "import sys, numpy as np, ngsolve as ng\n"
        f"sys.path.insert(0, r'{root / 'src'}')\n"
        "from ngsolve.meshes import MakeStructured3DMesh\n"
        "from radia.vim._vim import _charge_basis_hex, build_charge_gram\n"
        "m = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1,"
        " mapping=lambda x, y, z: (x, y, z/20.0))\n"
        "f = ng.HDiv(m, order=1)\n"
        "with ng.TaskManager():\n"
        "    _charge_basis_hex(f, cob_quad=3)\n"
        "    _, g, _ = build_charge_gram(f, eps=1e-10, leafsize=256, eta=2.0)\n"
        "n = int(g.stats()['n_dof'])\n"
        "G = np.array([[g.entry(i, j) for j in range(n)] for i in range(n)])\n"
        "print(np.linalg.eigvalsh(0.5*(G+G.T))[0])\n")
    out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                         text=True, encoding="utf-8", errors="replace",
                         cwd=str(root), check=False)
    assert out.returncode == 0, out.stderr
    assert float(out.stdout.strip().splitlines()[-1]) < 0.0, (
        "the Duffy family no longer loses definiteness on a flat hex, so this "
        "routing may no longer be needed")


def test_a_routed_pair_refuses_its_shape_derivative():
    """A Gram block integrated on the graded rule has no consistent
    directional derivative, so it must raise rather than return a 2 %
    inconsistent gradient."""
    from radia.vim._vim import _charge_basis_hex, build_charge_gram

    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1,
                                mapping=lambda x, y, z: (x, y, z / 20.0))
    fes = ng.HDiv(mesh, order=1)
    with ng.TaskManager():
        cb = _charge_basis_hex(fes, cob_quad=3)
        _, gram, _ = build_charge_gram(fes, eps=1e-10, leafsize=256, eta=2.0)
    cells = np.asarray(cb["cell_nodes"]).reshape(-1, 27, 3)
    faces = np.asarray(cb["face_nodes"]).reshape(-1, 9, 3)
    velocity = np.array([0.021, -0.017, 0.013])
    with pytest.raises(RuntimeError, match="length-scale ratio"):
        gram.hex_charge_gram_directional_derivative(
            np.zeros_like(cells) + velocity, np.zeros_like(faces) + velocity)


def test_the_dilation_identity_still_holds_below_the_limit():
    """Euler homogeneity of degree -1 on a pair that stays on the Duffy
    family, which is what the routing must not disturb."""
    from radia.vim._vim import _charge_basis_hex, build_charge_gram

    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1,
                                mapping=lambda x, y, z: (x, y, z / 8.0))
    fes = ng.HDiv(mesh, order=1)
    with ng.TaskManager():
        cb = _charge_basis_hex(fes, cob_quad=3)
        _, gram, _ = build_charge_gram(fes, eps=1e-10, leafsize=256, eta=2.0)
    cells = np.asarray(cb["cell_nodes"]).reshape(-1, 27, 3)
    faces = np.asarray(cb["face_nodes"]).reshape(-1, 9, 3)
    scaling = np.asarray(
        gram.hex_charge_gram_directional_derivative(cells, faces))
    n = scaling.shape[0]
    value = np.array([[gram.entry(i, j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(scaling, -value, rtol=3e-10, atol=3e-13)
