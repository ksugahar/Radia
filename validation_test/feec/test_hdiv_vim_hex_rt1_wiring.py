"""Golden: the PRODUCTION vim wiring of the pure-hex RT1 charge Gram.

`radia.vim._vim.build_charge_gram(HDiv(hexmesh, order=1))` must route to the hex-mode C++
`_ChargeGramHMatrix` (Q1 volume charge on the L2-order-1 space + 27-node Q2 geometry) and produce a demag
operator N = B^T G B whose generalized spectrum vs the HDiv mass is inside the physical bound [0, 1] (max
~0.998 on a structured cube), with a cube demag factor ~1/3.

This locks THREE things at once so they cannot drift silently:
  * the block-memo build (bit-identical to the per-entry callback -- a regression would move the spectrum),
  * the hex-mode default near_grade=0.6 (a looser/tighter grade shifts the near-quadrature -> the spectrum),
  * the build_charge_gram hex ROUTING (a pure-hex mesh must NOT fall through to the TET-only raise).

Self-contained (MakeStructured3DMesh, no fixture); ~10 s.  See memory hdiv-tet-hex-coupling-pyramid-gated
(the block-memo perf fix + the near_grade=0.6 eig de-risk).
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
import scipy.linalg as sla  # noqa: E402

from radia.vim._vim import build_charge_gram  # noqa: E402


def _materialize_N(B, G):
    n = B.shape[1]
    N = np.zeros((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        N[:, j] = B.T @ np.asarray(G.matvec_sym((B @ e).tolist()), float)
        e[j] = 0.0
    return 0.5 * (N + N.T)


@pytest.mark.parametrize(
    "name,mapping",
    [
        ("affine", lambda x, y, z: (0.02 * x, 0.02 * y, 0.02 * z)),
        ("distorted", lambda x, y, z: (0.02 * (x + 0.15 * y * z),
                                       0.02 * (y + 0.1 * x * z),
                                       0.02 * (z + 0.12 * x * y))),
    ],
)
def test_hex_rt1_wiring_spectrum_in_bound(name, mapping):
    mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=mapping)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = build_charge_gram(fes)          # <-- production wiring, auto hex branch
        N = _materialize_N(B, G)
        w = sla.eigh(N, M_mass.toarray(), eigvals_only=True)
        # cube demag factor via the uniform-Mz Rayleigh quotient
        gf = ng.GridFunction(fes)
        gf.Set(ng.CoefficientFunction((0, 0, 1.0)))
        mu = gf.vec.FV().NumPy().copy()
        demag = float((mu @ (N @ mu)) / (mu @ (M_mass @ mu)))

    # physical demag bound: 0 <= eig(M^-1 N) <= 1 (PSD floor + the self-energy cap)
    assert w.min() > -1e-8, f"{name}: N not PSD (min eig {w.min():.2e})"
    assert w.max() < 1.005, f"{name}: demag spectrum escaped [0,1] (max eig {w.max():.6f})"
    assert w.max() > 0.99, f"{name}: max eig {w.max():.6f} unexpectedly low (near-quadrature regression?)"
    # a solid cube's demag factor is ~1/3 (uniform magnetisation)
    assert 0.31 < demag < 0.35, f"{name}: cube demag factor {demag:.4f} off ~1/3"


def test_pure_hex_routes_to_hex_branch():
    """A pure-hex HDiv(order=1) must reach the hex Gram, NOT the TET-only ValueError."""
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1,
                                mapping=lambda x, y, z: (0.02 * x, 0.02 * y, 0.02 * z))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = build_charge_gram(fes)
        # 1 hex -> 8 Q1 volume charges + 6 faces * 4 = 24 surface charges = 32 charges; ndof = HDiv order-1 dofs
        assert B.shape[0] == 8 + 6 * 4
        assert G.ndof() == B.shape[0]


# The worst hex of a REAL Cubit cylinder mesh (|J|min/max = 0.376 -- genuine trilinear warp).  The
# near-affine boxes above CANNOT catch warped-hex regressions: the shipped pre-Piola model passed them at
# 0.998 while CONVERGING to eig 1.0105 > 1 on this hex (the non-Coulomb-orthogonal L2 charge projection
# leaked the physical bound by O(warp^2)).  This fixture locks the Piola-exact charge model on real warp.
_WORST_CYL_HEX_CORNERS = [   # x-fastest corner order (ix + 2*iy + 4*iz)
    (-0.5, 1e-16, 0.5),
    (-0.4619397662556434, 0.1913417161825449, 0.5),
    (-0.3799374424886703, -5.15188134101e-05, 0.5),
    (-0.3135968748399616, 0.0996211631527812, 0.5),
    (-0.5, -1e-16, 0.3),
    (-0.461939766255643, 0.1913417161825455, 0.3),
    (-0.3799374424886702, -5.151881341e-05, 0.2999999999999999),
    (-0.3135968748399616, 0.0996211631527814, 0.2999999999999999),
]


def test_hex_rt1_strongly_warped_real_hex_spectrum():
    """Piola-exact charge model on GENUINE warp: the worst cylinder hex's demag spectrum must respect the
    physical bound [0, 1].  (Pre-Piola this converged to 1.0105 -- a real-mesh-only failure mode.)"""
    corners = np.array(_WORST_CYL_HEX_CORNERS)

    def trilerp(x, y, z):
        out = np.zeros(3)
        for i in range(8):
            rx, ry, rz = i & 1, (i >> 1) & 1, (i >> 2) & 1
            out += (x if rx else 1 - x) * (y if ry else 1 - y) * (z if rz else 1 - z) * corners[i]
        return tuple(out)

    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1, mapping=trilerp)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = build_charge_gram(fes)
        N = _materialize_N(B, G)
        w = sla.eigh(N, M_mass.toarray(), eigvals_only=True)
    assert w.min() > -1e-8, f"warped hex: N not PSD (min eig {w.min():.2e})"
    assert w.max() < 1.005, f"warped hex: demag spectrum escaped [0,1] (max eig {w.max():.6f})"
