"""
Periodic / ANTI-periodic boundary conditions (the motor sector-model 数式).

PM machines are modelled on ONE pole pitch and tiled to the whole machine with
boundary conditions on the two cut edges: PERIODIC when an even number of poles is
spanned, ANTI-PERIODIC (field repeats with opposite sign) across a single pole pitch.
Getting these right is what lets radia-motor model a sector instead of the full machine.

Forged on the unit square [0,1]^2, periodic in x, Dirichlet u=0 on y=0,1, with the
method of manufactured solutions:
  * PERIODIC:      u = sin(2 pi x) sin(pi y)   (u(0,y)=u(1,y)),     -lap u = 5 pi^2 u
  * ANTI-PERIODIC: u = cos(pi x)  sin(pi y)    (u(1,y) = -u(0,y)),  -lap u = 2 pi^2 u

We solve int grad(u).grad(v) = int f v on a Periodic FE space (phase=+1 / phase=-1) and
check the FEM solution recovers the exact manufactured field, and -- for the anti-periodic
case -- that the solution genuinely satisfies u(1,y) = -u(0,y). Pure NGSolve (Bash-robust).
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, Periodic, BilinearForm, LinearForm, GridFunction,
                     Integrate, Norm, grad, dx, x, y, sin, cos, pi)
from ngsolve.meshes import MakeStructured2DMesh


def _mesh(n):
    return MakeStructured2DMesh(quads=False, nx=n, ny=n, periodic_x=True)


def _solve(fes, f_cf):
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    f = LinearForm(f_cf * v * dx).Assemble()
    gf = GridFunction(fes)
    gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec     # default inverse (real or complex)
    return gf


def _l2(gf, ref, mesh):
    return math.sqrt(Integrate(Norm(gf - ref) ** 2, mesh).real)


def test_periodic_recovers_manufactured():
    """Periodic (phase +1) FE space recovers u = sin(2 pi x) sin(pi y)."""
    mesh = _mesh(32)
    fes = Periodic(H1(mesh, order=3, dirichlet="bottom|top"))
    uex = sin(2 * pi * x) * sin(pi * y)
    gf = _solve(fes, 5 * pi * pi * uex)
    assert _l2(gf, uex, mesh) < 1e-4, "periodic L2 error %.2e" % _l2(gf, uex, mesh)


def test_antiperiodic_recovers_manufactured():
    """Anti-periodic (phase -1) FE space recovers u = cos(pi x) sin(pi y)."""
    mesh = _mesh(32)
    fes = Periodic(H1(mesh, order=3, dirichlet="bottom|top", complex=True), phase=[-1])
    uex = cos(pi * x) * sin(pi * y)
    gf = _solve(fes, 2 * pi * pi * uex)
    assert _l2(gf, uex, mesh) < 1e-4, "anti-periodic L2 error %.2e" % _l2(gf, uex, mesh)


def test_antiperiodic_sign_flip_across_pitch():
    """The anti-periodic solution genuinely satisfies u(1,y) = -u(0,y)."""
    mesh = _mesh(32)
    fes = Periodic(H1(mesh, order=3, dirichlet="bottom|top", complex=True), phase=[-1])
    uex = cos(pi * x) * sin(pi * y)
    gf = _solve(fes, 2 * pi * pi * uex)
    for yy in (0.3, 0.5, 0.7):
        left = gf(mesh(0.0, yy))
        right = gf(mesh(1.0, yy))
        assert abs(right + left) < 1e-3, "not anti-periodic at y=%.1f: u0=%.4f u1=%.4f" % (yy, left, right)


def test_periodic_h_convergence():
    """Periodic BC keeps optimal high-order convergence (order 3 -> rate ~4)."""
    uex = sin(2 * pi * x) * sin(pi * y)
    e1 = _l2(_solve(Periodic(H1(_mesh(8), order=3, dirichlet="bottom|top")), 5 * pi * pi * uex), uex, _mesh(8))
    e2 = _l2(_solve(Periodic(H1(_mesh(16), order=3, dirichlet="bottom|top")), 5 * pi * pi * uex), uex, _mesh(16))
    rate = math.log(e1 / e2) / math.log(2.0)
    assert rate > 3.0, "periodic h-rate %.2f below optimal" % rate


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
