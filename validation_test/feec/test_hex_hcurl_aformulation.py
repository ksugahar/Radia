"""
A-formulation (curl-curl) weak form on high-order HEX with Hcurl edge elements.

The magnetostatic A-formulation and eddy-current solvers rest on the curl-curl
operator int(nu curl(A) . curl(v)) on Hcurl (edge) elements.
validation_test/cubit/test_hex_highorder_fem.py validates the H1 (scalar)
high-order hex; this validates the VECTOR (Hcurl) curl-curl operator at high
order on hex -- the differentiator extended to edge elements.

Manufactured solution (definite curl-curl + mass, well-posed on Hcurl):
    A_exact = (0, 0, sin(pi x) sin(pi y))
    curl curl A_exact = 2 pi^2 A_exact   (verified by hand)
    => f = curl curl A + A = (2 pi^2 + 1) A_exact
    tangential trace of A_exact = 0 on the unit cube  => homogeneous Hcurl Dirichlet

Solve int(curl(u).curl(v) + u.v) = int f.v and check ||A - A_exact||_L2 drops at the
high-order (spectral) rate on genuine hex edge elements.
"""
import math
import pytest

pytestmark = pytest.mark.slow

ng = pytest.importorskip("ngsolve")
from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction, Integrate,
                     curl, dx, x, y, z, sin, pi, CF)
from ngsolve.meshes import MakeStructured3DMesh

_AEX = CF((0, 0, sin(pi * x) * sin(pi * y)))
_F = (2 * pi * pi + 1) * _AEX


def _hex_mesh(n):
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n)


def _solve_l2(mesh, order):
    fes = HCurl(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(curl(u) * curl(v) * dx + u * v * dx).Assemble()
    f = LinearForm(_F * v * dx).Assemble()
    gf = GridFunction(fes)
    gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return math.sqrt(Integrate((gf - _AEX) * (gf - _AEX), mesh))


def test_hcurl_p_refinement_monotone():
    """Hcurl curl-curl L2 error drops monotonically and steeply with order on hex."""
    mesh = _hex_mesh(8)
    errs = [_solve_l2(mesh, p) for p in (0, 1, 2, 3)]
    for lo, hi in zip(errs[1:], errs[:-1]):
        assert lo < hi, "Hcurl p-refinement not monotone: %s" % errs
    assert errs[0] / errs[1] > 10.0
    assert errs[2] / errs[3] > 10.0


def test_hcurl_high_order_accuracy():
    """order=3 Hcurl on an 8^3 hex mesh reaches spectral accuracy (< 1e-6)."""
    assert _solve_l2(_hex_mesh(8), 3) < 1e-6


def test_hcurl_h_convergence_order1():
    """h-refinement at Hcurl order=1 converges at >= the optimal O(h^2) rate
    (observed ~3 for this smooth manufactured field on hex)."""
    e4 = _solve_l2(_hex_mesh(4), 1)
    e8 = _solve_l2(_hex_mesh(8), 1)
    rate = math.log(e4 / e8) / math.log(2.0)
    assert 1.8 < rate < 3.4, "Hcurl order=1 h-rate %.2f not in [2,3] band" % rate


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
