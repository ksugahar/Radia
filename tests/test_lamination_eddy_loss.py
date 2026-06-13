"""
Classical lamination EDDY LOSS (the iron-loss eddy 数式 for motor coreloss).

A lamination of thickness d (conductivity sigma, permeability mu) carrying a sinusoidal
in-plane flux density B_y(t)=B_p sin(wt) develops eddy currents J_z across its thickness.
The 1-D magnetic diffusion across [-d/2, d/2] is

    d^2 B_y/dx^2 = k^2 B_y ,   k=(1+j)/delta,  delta=sqrt(2/(w mu sigma)),  B_y(+-d/2)=B_p,

exact field B_y(x)=B_p cosh(kx)/cosh(kd/2), eddy current J_z=-(1/mu) dB_y/dx, and the
time-average loss per unit volume

    P = (1/d) int_{-d/2}^{d/2} |J_z|^2/(2 sigma) dx = (1/(2 sigma mu^2 d)) int |dB_y/dx|^2 .

In the THIN limit d << delta this collapses to the textbook lamination law (mu cancels):

    P_thin = sigma d^2 w^2 B_p^2 / 24 = pi^2 sigma d^2 f^2 B_p^2 / 6 .

This forges the eddy-loss weak form, gating the FEM loss against (i) quadrature of the
EXACT field (no hand-derived loss constant) and (ii) the thin-limit d^2-law, and verifies
the classic d^2 and f^2 scalings. Pure NGSolve + math (Bash-robust).
"""
import math
import cmath
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, GridFunction, Integrate, Norm,
                     grad, dx, x, BND)
from ngsolve.meshes import Make1DMesh

MU0 = 4e-7 * math.pi
MUR = 1000.0
MU = MUR * MU0
SIGMA = 2.0e6          # electrical steel
BP = 1.0               # surface flux density [T]


def _params(d, f):
    w = 2 * math.pi * f
    delta = math.sqrt(2.0 / (w * MU * SIGMA))
    k = (1 + 1j) / delta
    return w, delta, k


def _fem_loss(d, f, order=4, n=400):
    w, delta, k = _params(d, f)
    mesh = Make1DMesh(n, mapping=lambda t: -d / 2 + t * d)
    fes = H1(mesh, order=order, complex=True, dirichlet="left|right")
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx + (k * k) * u * v * dx).Assemble()
    g = GridFunction(fes)
    g.Set(BP, BND)
    r = g.vec.CreateVector()
    r.data = -a.mat * g.vec
    g.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    dB = grad(g)[0]
    return (1.0 / (2 * SIGMA * MU * MU * d)) * Integrate(Norm(dB) ** 2, mesh)


def _exact_loss(d, f, nq=6000):
    """Loss from quadrature of the EXACT field derivative (independent closed form)."""
    w, delta, k = _params(d, f)
    cd = cmath.cosh(k * d / 2)
    acc = 0.0
    for i in range(nq):
        xx = -d / 2 + (i + 0.5) / nq * d
        dB = BP * k * cmath.sinh(k * xx) / cd
        acc += abs(dB) ** 2
    integral = acc * (d / nq)
    return (1.0 / (2 * SIGMA * MU * MU * d)) * integral


def _thin_law(d, f):
    return math.pi ** 2 * SIGMA * d ** 2 * f ** 2 * BP ** 2 / 6.0


@pytest.mark.parametrize("d,f", [(0.5e-3, 200.0), (0.35e-3, 400.0)])
def test_fem_loss_matches_exact_field(d, f):
    """FEM eddy loss equals the loss from quadrature of the exact field."""
    fem = _fem_loss(d, f)
    ref = _exact_loss(d, f)
    assert abs(fem - ref) / ref < 1e-3, "d=%.1emm f=%g: FEM=%.4e exact=%.4e" % (d * 1e3, f, fem, ref)


def test_thin_limit_d2_law():
    """In the thin regime d<<delta the loss matches the textbook pi^2 sigma d^2 f^2 B^2/6."""
    d, f = 0.2e-3, 50.0
    w, delta, k = _params(d, f)
    assert d / delta < 0.2, "not in the thin regime (d/delta=%.2f)" % (d / delta)
    fem = _fem_loss(d, f)
    law = _thin_law(d, f)
    assert abs(fem - law) / law < 0.03, "thin loss %.4e vs d^2-law %.4e" % (fem, law)


def test_d_squared_scaling():
    """At fixed frequency (thin regime), halving the thickness quarters the loss."""
    f = 50.0
    p1 = _fem_loss(0.2e-3, f)
    p2 = _fem_loss(0.1e-3, f)
    assert abs((p1 / p2) - 4.0) / 4.0 < 0.05, "d^2 scaling broken: ratio %.3f" % (p1 / p2)


def test_f_squared_scaling():
    """In the thin regime the loss scales as f^2 (doubling f quadruples the loss)."""
    d = 0.15e-3
    p1 = _fem_loss(d, 50.0)
    p2 = _fem_loss(d, 100.0)
    assert abs((p2 / p1) - 4.0) / 4.0 < 0.05, "f^2 scaling broken: ratio %.3f" % (p2 / p1)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
