"""
Hexahedral-element magnetostatic CONVERGENCE validation.

This complements test_tetrahedron_hexahedron.py (API smoke tests) with a *quantitative*
correctness guarantee for the hex element's solved magnetization, and it pins down a
subtlety in examples/cube_uniform_field/hexahedron/benchmark_hex.py:

  The benchmark reports "error %" of a hex cube vs the analytic average magnetization
      M = chi*H / (1 + chi*N),   N = 1/3   (cube magnetometric demag factor)
  but that analytic value is the chi -> 0 (uniform-M) limit. At HIGH chi (mu_r=1000)
  the magnetization concentrates at edges/corners and is strongly non-uniform, so the
  true cube average legitimately exceeds the N=1/3 value (~19% at mu_r=1000). That gap
  is PHYSICS, not solver error.

  In the LOW-chi regime the N=1/3 average is essentially exact, so the solved hex cube
  must converge to it -- that is what these tests assert, isolating the element's
  discretization accuracy from the high-chi shape effect.
"""
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../examples/cube_uniform_field'))
import radia as rad
from benchmark_common import generate_hex_mesh, CUBE_SIZE, H_EXT, MU_0


def _solve_hex_cube(mu_r, n_div):
    """Solve an n_div^3 hex cube of linear material mu_r in a uniform axial H, return M_avg_z."""
    rad.UtiDelAll()
    cnt = rad.ObjCnt([rad.ObjHexahedron(v, [0, 0, 0])
                      for v in generate_hex_mesh(n_div, CUBE_SIZE)])
    rad.MatApl(cnt, rad.MatLin(mu_r))
    grp = rad.ObjCnt([cnt, rad.ObjBckg(lambda p: [0, 0, MU_0 * H_EXT])])
    rad.Solve(grp, 1e-6, 2000, 0)          # LU, linear
    Ms = rad.ObjM(cnt)
    return sum(m[1][2] for m in Ms) / len(Ms)


def _analytic_cube(mu_r):
    chi = mu_r - 1.0
    return chi * H_EXT / (1.0 + chi / 3.0)   # N = 1/3 cube magnetometric (exact as chi->0)


class TestHexDemagConvergence:
    """Hex cube converges to the exact N=1/3 cube magnetization in the low-chi regime."""

    def test_lowchi_accurate(self):
        """At mu_r=1.05 a coarse N=6 hex cube already matches N=1/3 to < 0.05%."""
        M_an = _analytic_cube(1.05)
        err = abs(_solve_hex_cube(1.05, 6) - M_an) / M_an * 100
        assert err < 0.05, "low-chi N=6 hex error %.4f%% too high" % err

    def test_convergence_with_refinement(self):
        """Refining the hex mesh reduces the low-chi error (N=8 better than N=4)."""
        M_an = _analytic_cube(1.05)
        err4 = abs(_solve_hex_cube(1.05, 4) - M_an) / M_an * 100
        err8 = abs(_solve_hex_cube(1.05, 8) - M_an) / M_an * 100
        assert err8 < err4, "no convergence: err(N=8)=%.4f%% !< err(N=4)=%.4f%%" % (err8, err4)
        assert err8 < 0.02, "N=8 low-chi error %.4f%% too high" % err8

    def test_highchi_exceeds_uniform_demag(self):
        """At mu_r=1000 the cube average legitimately EXCEEDS the N=1/3 (uniform-M) value
        because M is non-uniform -- documents that benchmark_hex's high-chi 'error' is physics."""
        M_an = _analytic_cube(1000.0)
        M_hi = _solve_hex_cube(1000.0, 8)
        assert M_hi > M_an, "high-chi cube average should exceed the uniform-M N=1/3 value"
        dev = (M_hi - M_an) / M_an * 100
        assert 10.0 < dev < 30.0, "high-chi deviation %.1f%% outside expected band" % dev


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
