"""Golden test: analytic VOLUME Gram (#B) -- the tet Newtonian potential phi_tet + the full analytic
charge Gram that closes the nonlinear sharp-body non-uniform residual.

The surface Wilton Gram (test_hdiv_vim_wilton_gram.py) makes the demag FACTOR exact (uniform M = pure
surface charge).  For NON-uniform M (nonlinear, div M != 0) the VOLUME-involving Gram blocks (cell-cell,
cell-face) matter, and the centroid-monopole + sub-point near-correction left a residual (cube
nonlinear ~8.7% vs Radia).  The analytic tet potential

    phi_tet(V, P) = INT_tet 1/|P-r'| dV' = (1/2) sum_{4 faces} d_face * tri_potential(face, P)

(divergence theorem nabla'^2 R = 2/R; reuses the exact Wilton triangle potential -> NO new singular
quadrature) gives the full analytic charge Gram (build_demag(analytic_gram=True)).  It substantially
closes the residual (cube nonlinear vs Radia: 13% -> 6.2% at H0=2e5, 5.1% -> 1.5% at H0=5e5) while
keeping the linear demag factor exact.

Locks: phi_tet vs fine numerical; the analytic-Gram demag factor -> 1/3 on sphere AND cube.
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _core as tet  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402


def test_phi_tet_matches_numerical():
    """The tet Newtonian potential (Wilton-face divergence-theorem sum) matches fine volume quadrature."""
    V = np.array([[0, 0, 0], [1.0, 0, 0], [0.2, 0.9, 0], [0.3, 0.2, 1.1]])

    def numref(P, nsub=60):
        lam = []
        for i in range(nsub):
            for j in range(nsub - i):
                for k in range(nsub - i - j):
                    l = nsub - 1 - i - j - k
                    lam.append([(i + .25) / nsub, (j + .25) / nsub, (k + .25) / nsub, (l + .25) / nsub])
        lam = np.array(lam)
        C = lam @ V
        vol = abs(np.linalg.det(V[1:] - V[0])) / 6.0
        return float(np.sum((vol / len(C)) / np.linalg.norm(C - P, axis=1)))

    for P in ([0.3, 0.3, 1.6], [1.6, 0.3, 0.4], [0.3, 0.3, -0.5]):
        a = tet.phi_tet(V, np.array(P, float))
        n = numref(np.array(P, float))
        assert abs(a - n) < 1e-2 * abs(n), f"phi_tet {a:.5f} vs numref {n:.5f} at {P}"


def _demag_analytic(geo_kind):
    g = CSGeometry()
    if geo_kind == "sphere":
        g.Add(Sphere(Pnt(0, 0, 0), 1.0)); h = 0.4
    else:
        g.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))); h = 0.3
    with ng.TaskManager():
        mesh = ng.Mesh(g.GenerateMesh(maxh=h))
    return tet.demag_factor(tet.build_demag(mesh, 4, analytic_gram=True))


def test_analytic_gram_demag_factor_sphere():
    D = _demag_analytic("sphere")
    assert abs(D - 1.0 / 3) < 5e-3, f"analytic-Gram sphere demag {D:.4f} not within 0.5% of 1/3"


def test_analytic_gram_demag_factor_cube():
    D = _demag_analytic("cube")
    assert abs(D - 1.0 / 3) < 5e-3, f"analytic-Gram cube demag {D:.4f} not within 0.5% of 1/3"
