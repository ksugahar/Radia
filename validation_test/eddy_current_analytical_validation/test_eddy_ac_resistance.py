"""
Round-conductor AC resistance (skin effect) vs the Bessel/Kelvin closed form.

The radia-motor winding-loss 数式: a solid round conductor carrying AC has
    R_ac/R_dc = Re[ (k a / 2) I0(k a) / I1(k a) ],  k = (1+j)/delta,
    delta = sqrt(2/(omega mu sigma)).
This solves the radial complex magnetic-diffusion equation for the axial current
density J_z(r) with the axisymmetric (r-weighted) weak form on a 1-D radial mesh
[0, a] (regularity at r=0 natural, J_z(a) prescribed), computes R_ac/R_dc from the
solved profile, and checks it against the exact Bessel ratio -- cross-validating the
jmag_converter cap round_conductor_ac_resistance_bessel. Pure NGSolve + scipy.
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
sp = pytest.importorskip("scipy.special")
from scipy.special import iv
from ngsolve import (H1, BilinearForm, GridFunction, Integrate, Norm,
                     grad, dx, x, BND)
from ngsolve.meshes import Make1DMesh

MU0 = 4e-7 * math.pi
SIGMA = 5.8e7   # copper


@pytest.fixture(autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _fem_rac_over_rdc(a_over_delta, a=0.005, order=4, n=600):
    delta = a / a_over_delta
    omega = 2.0 / (MU0 * SIGMA * delta ** 2)
    mesh = Make1DMesh(n, mapping=lambda t: t * a)
    fes = H1(mesh, order=order, complex=True, dirichlet="right")
    u, v = fes.TnT()
    k2 = 1j * omega * MU0 * SIGMA
    A = BilinearForm(grad(u) * grad(v) * x * dx + k2 * u * v * x * dx).Assemble()
    g = GridFunction(fes)
    g.Set(1.0, BND)
    r = g.vec.CreateVector()
    r.data = -A.mat * g.vec
    g.vec.data += A.mat.Inverse(fes.FreeDofs()) * r
    I = 2 * math.pi * complex(Integrate(g * x, mesh))           # total current
    Jint = 2 * math.pi * Integrate(Norm(g) ** 2 * x, mesh)      # int |J|^2 dA
    return (math.pi * a * a) * Jint / abs(I) ** 2              # R_ac/R_dc


def _bessel_rac_over_rdc(a_over_delta):
    ka = (1 + 1j) * a_over_delta
    return (ka / 2 * iv(0, ka) / iv(1, ka)).real


@pytest.mark.parametrize("a_over_delta", [1.0, 3.0, 5.0])
def test_round_wire_ac_resistance(a_over_delta):
    """FEM R_ac/R_dc matches the Bessel/Kelvin closed form to high precision."""
    fem = _fem_rac_over_rdc(a_over_delta)
    an = _bessel_rac_over_rdc(a_over_delta)
    assert abs(fem - an) / an < 1e-6, \
        "a/delta=%.1f: R_ac/R_dc FEM=%.6f != Bessel=%.6f" % (a_over_delta, fem, an)


def test_dc_limit_unity():
    """As a/delta -> 0 the AC resistance ratio -> 1 (uniform current)."""
    assert abs(_fem_rac_over_rdc(0.05) - 1.0) < 1e-3


def test_high_freq_asymptote():
    """For a >> delta, R_ac/R_dc -> a/(2 delta) + 1/4 (classic skin-effect asymptote)."""
    ad = 20.0
    fem = _fem_rac_over_rdc(ad, n=2000)
    asymp = ad / 2.0 + 0.25
    assert abs(fem - asymp) / asymp < 0.02, "R_ac/R_dc=%.4f vs asymptote %.4f" % (fem, asymp)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
