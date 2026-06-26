"""
High-order HEX finite-element convergence in NGSolve (the radia differentiator).

The flagship cubit-mesh-export path feeds Cubit hex meshes into NGSolve at high polynomial
order. validation_test/cubit/test_ngsolve_volume_hex_sphere.py checks the GEOMETRIC accuracy of the
curved hex export (needs live Cubit). This test checks the complementary, more fundamental
property -- that NGSolve actually SOLVES a PDE on high-order hex elements at the OPTIMAL rate --
using a structured hex mesh, so it runs with no Cubit dependency (pure NGSolve, CI-friendly).

Manufactured solution (method of manufactured solutions):
    u_exact = sin(pi x) sin(pi y) sin(pi z)   on the unit cube, u = 0 on the boundary
    -lap u  = 3 pi^2 u                          (the assembled right-hand side)

Verifies on genuine 8-vertex hexes:
  * p-refinement: L2 error drops monotonically and steeply with polynomial order,
  * h-refinement: convergence at the optimal O(h^{p+1}) rate,
  * high-order accuracy: p=4 on an 8^3 hex mesh reaches < 1e-6.
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, Integrate,
                     grad, dx, x, y, z, sin, pi, VOL)
from ngsolve.meshes import MakeStructured3DMesh

_UEX = sin(pi * x) * sin(pi * y) * sin(pi * z)
_RHS = 3 * pi * pi * _UEX


def _hex_mesh(n):
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n)


def _solve_l2(mesh, p):
    fes = H1(mesh, order=p, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    f = LinearForm(_RHS * v * dx).Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return math.sqrt(Integrate((gfu - _UEX) ** 2, mesh))


def test_mesh_is_hex():
    """MakeStructured3DMesh(hexes=True) produces 8-vertex hex volume elements."""
    mesh = _hex_mesh(4)
    el0 = list(mesh.Elements(VOL))[0]
    assert len(el0.vertices) == 8, "expected 8-vertex hexes, got %d" % len(el0.vertices)
    assert mesh.ne == 64


def test_p_refinement_monotone():
    """L2 error drops monotonically and steeply with polynomial order on hex."""
    mesh = _hex_mesh(8)
    errs = [_solve_l2(mesh, p) for p in (1, 2, 3, 4)]
    for lo, hi in zip(errs[1:], errs[:-1]):
        assert lo < hi, "p-refinement not monotone: %s" % errs
    # each order should buy at least ~10x (high-order, not merely algebraic)
    assert errs[0] / errs[1] > 10.0
    assert errs[2] / errs[3] > 10.0


def test_high_order_accuracy():
    """p=4 on an 8^3 hex mesh reaches spectral accuracy (< 1e-6)."""
    assert _solve_l2(_hex_mesh(8), 4) < 1e-6


def test_optimal_h_rate_order2():
    """h-refinement at p=2 converges at the optimal O(h^{p+1}) = O(h^3) rate."""
    e4 = _solve_l2(_hex_mesh(4), 2)
    e8 = _solve_l2(_hex_mesh(8), 2)
    rate = math.log(e4 / e8) / math.log(2.0)
    assert 2.5 < rate < 3.6, "p=2 h-rate %.2f not ~3 (optimal)" % rate


# --- variable-coefficient magnetostatic-form weak form on hex -----------------
# The core magnetostatic operator is int mu(x) grad(u) . grad(w) with spatially
# varying mu. This checks that the VARIABLE-coefficient weak form (not just the
# constant-coefficient Poisson above) keeps optimal high-order accuracy on hexes,
# via a manufactured solution with a smooth mu(x) = 1 + 0.5 r^2.
from ngsolve import CF, cos  # noqa: E402

_MU = 1 + 0.5 * (x * x + y * y + z * z)
_GMU = CF((x, y, z))
_GUE = CF((pi * cos(pi * x) * sin(pi * y) * sin(pi * z),
           pi * sin(pi * x) * cos(pi * y) * sin(pi * z),
           pi * sin(pi * x) * sin(pi * y) * cos(pi * z)))
_F_VAR = 3 * pi * pi * _MU * _UEX - (_GMU * _GUE)   # -div(mu grad u)


def _solve_l2_varcoeff(mesh, p):
    fes = H1(mesh, order=p, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(_MU * grad(u) * grad(v) * dx).Assemble()
    f = LinearForm(_F_VAR * v * dx).Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return math.sqrt(Integrate((gfu - _UEX) ** 2, mesh))


def test_variable_coeff_high_order():
    """int mu(x) grad(u).grad(w) (variable mu) keeps spectral accuracy on hex:
    monotone p-refinement and p=4 on 8^3 reaching < 1e-6."""
    mesh = _hex_mesh(8)
    errs = [_solve_l2_varcoeff(mesh, p) for p in (1, 2, 3, 4)]
    for lo, hi in zip(errs[1:], errs[:-1]):
        assert lo < hi, "variable-coeff p-refinement not monotone: %s" % errs
    assert errs[-1] < 1e-6, "variable-coeff p=4 error %.2e too high" % errs[-1]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
