"""Golden test: the exact analytic flat-triangle potential (radia.vim.tri_potential).

tri_potential is the exact Newtonian potential of a uniformly-charged flat triangle (the Wilton/Graglia
analytic integral) -- the building block of the polytope / face self-energy diagonal in the HDiv-VIM
charge Gram.  This locks it against a fine numerical reference.  (The demag-factor tests that used the
now-removed dense Wilton/monopole Gram path live elsewhere -- demag ~1/3 is locked via the production
C++ charge-Gram path in test_hdiv_vim_volume_gram / test_hdiv_vim_demag_solve.)

Reference: Wilton et al., IEEE TAP 32(3):276 (1984); Graglia, IEEE TAP 41(10):1448 (1993).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

from radia.vim import _core as tet  # noqa: E402


def test_tri_potential_matches_numerical():
    """The analytic triangle potential matches a fine numerical reference (off-plane points)."""
    V = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.3, 0.9, 0]])

    def numref(r, nsub=300):
        lam = []
        for i in range(nsub):
            for j in range(nsub - i):
                k = nsub - 1 - i - j
                lam.append([(i + 1 / 3) / nsub, (j + 1 / 3) / nsub, (k + 1 / 3) / nsub])
        lam = np.array(lam)
        C = lam @ V
        A = 0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0]))
        return float(np.sum((A / len(C)) / np.linalg.norm(C - r, axis=1)))

    for r in ([0.4, 0.3, 0.5], [0.4, 0.3, 1.0], [1.2, 0.2, 0.3], [0.4, 0.3, -0.4]):
        a = tet.tri_potential(V, np.array(r, float))
        n = numref(np.array(r, float))
        assert abs(a - n) < 5e-3 * abs(n), f"tri_potential {a:.5f} vs numref {n:.5f} at {r}"
