"""
Eddy-current field penetration / shielding in a conducting cylinder (axial AC field).

The eddy-current diffusion 数式 (induction heating, transformer lamination, motor
rotor shielding): an infinite conducting cylinder (radius a, conductivity sigma,
permeability mu) in a uniform AXIAL time-harmonic field H0 e^{jwt}. By symmetry the
interior field is H_z(r) only, with azimuthal eddy currents J_phi = -dH_z/dr, and

    (1/r) d/dr( r dH_z/dr ) = k^2 H_z,   k = (1+j)/delta,  delta = sqrt(2/(w mu sigma)),

so the exact field is the Kelvin/Bessel profile

    H_z(r) = H0 * I0(k r) / I0(k a),      shielding factor  H_z(0)/H0 = 1/I0(k a).

This forges the eddy weak form TWO ways, both gated by that closed form (no fabrication):
  (A) a 1-D radial (r-weighted) complex FEM on [0, a]  -> matches the Bessel profile to
      machine precision (this IS the operator, discretized);
  (B) a genuine 2-D scalar Helmholtz solve on an UNSTRUCTURED disk mesh
      (int grad(H).grad(v) + k^2 int H v = 0, Dirichlet H_z=H0 on r=a) -> reproduces the
      same Bessel interior to FEM accuracy. Outside a cylinder the induced azimuthal
      currents are solenoid-like (≈0 external field), so H_z=H0 on r=a is the exact BC.
  (C) low-frequency limit ka->0 : H_z -> H0 (no shielding);
  (D) eddy power loss per unit length from the solved profile vs the analytic Bessel value.

Pure NGSolve + scipy.special (Bash-robust, no MCP / Cubit / COMSOL dependency).
"""
import math
import cmath
import pytest

ng = pytest.importorskip("ngsolve")
sp = pytest.importorskip("scipy.special")
from scipy.special import iv
from ngsolve import (H1, BilinearForm, GridFunction, Integrate, Norm,
                     grad, dx, x, y, BND, CF, sqrt as ngsqrt)
from ngsolve.meshes import Make1DMesh

MU0 = 4e-7 * math.pi
SIGMA = 3.5e7      # aluminium-ish


def _k_omega(a, a_over_delta, mu_r=1.0):
    """Return (k, omega, mu) for a given a/delta."""
    delta = a / a_over_delta
    mu = mu_r * MU0
    omega = 2.0 / (mu * SIGMA * delta ** 2)
    k = (1 + 1j) / delta
    return k, omega, mu


def _bessel_Hz(r, a, k, H0=1.0):
    return H0 * iv(0, k * r) / iv(0, k * a)


# ---------------------------------------------------------------------------
# (A) 1-D radial complex FEM:  int (H' v') r dr + k^2 int H v r dr = 0
# ---------------------------------------------------------------------------
def _fem1d_profile(a, a_over_delta, order=4, n=500, H0=1.0):
    k, _, _ = _k_omega(a, a_over_delta)
    mesh = Make1DMesh(n, mapping=lambda t: t * a)
    fes = H1(mesh, order=order, complex=True, dirichlet="right")  # r=a is "right"
    u, v = fes.TnT()
    A = BilinearForm(grad(u) * grad(v) * x * dx + (k * k) * u * v * x * dx).Assemble()
    g = GridFunction(fes)
    g.Set(H0, BND)
    r = g.vec.CreateVector()
    r.data = -A.mat * g.vec
    g.vec.data += A.mat.Inverse(fes.FreeDofs()) * r
    return mesh, g, k


@pytest.mark.parametrize("a_over_delta", [0.5, 2.0, 5.0])
def test_radial_matches_bessel(a_over_delta):
    """1-D radial FEM field equals the Kelvin/Bessel profile at sample radii."""
    a = 0.01
    mesh, g, k = _fem1d_profile(a, a_over_delta)
    for frac in (0.0, 0.25, 0.5, 0.75):
        r = frac * a
        fem = complex(g(mesh(r)))
        an = _bessel_Hz(r, a, k)
        assert abs(fem - an) / abs(an) < 1e-5, \
            "a/d=%.1f r/a=%.2f: FEM=%s Bessel=%s" % (a_over_delta, frac, fem, an)


@pytest.mark.parametrize("a_over_delta", [1.0, 3.0, 6.0])
def test_shielding_factor(a_over_delta):
    """Center shielding factor |H_z(0)/H0| matches |1/I0(ka)| and decreases with frequency."""
    a = 0.01
    mesh, g, k = _fem1d_profile(a, a_over_delta)
    fem_center = abs(complex(g(mesh(0.0))))
    an_center = abs(1.0 / iv(0, k * a))
    assert abs(fem_center - an_center) / an_center < 1e-4
    assert fem_center < 1.0   # field is shielded from the core


def test_low_frequency_no_shielding():
    """ka -> 0 : the core field approaches the applied field (no shielding)."""
    a = 0.01
    mesh, g, k = _fem1d_profile(a, 0.05)
    assert abs(abs(complex(g(mesh(0.0)))) - 1.0) < 1e-3


def test_high_frequency_skin_confined():
    """a >> delta : the core is strongly shielded, |H_z(0)/H0| ~ |1/I0(ka)| << 1."""
    a = 0.01
    mesh, g, k = _fem1d_profile(a, 8.0, n=1500)
    fem_center = abs(complex(g(mesh(0.0))))
    an_center = abs(1.0 / iv(0, k * a))
    assert abs(fem_center - an_center) / an_center < 1e-3
    assert fem_center < 0.05   # deeply shielded core


# ---------------------------------------------------------------------------
# (D) eddy power loss per unit length from the 1-D profile vs analytic Bessel.
#     P' = (1/(2 sigma)) int |J_phi|^2 dA,  J_phi = -dH_z/dr.
#     Analytic: P' = (pi a H0^2 /(sigma delta)) * Re[ k a I1(ka) conj... ] -- instead of a
#     hand-derived constant (fabrication risk), gate the FEM loss against the loss computed
#     by RADIAL QUADRATURE of the EXACT Bessel derivative (independent closed-form integrand).
# ---------------------------------------------------------------------------
def test_eddy_loss_matches_bessel_quadrature():
    """FEM eddy loss equals the loss from quadrature of the exact Bessel current."""
    a, a_over_delta, H0 = 0.01, 3.0, 1.0
    k, omega, mu = _k_omega(a, a_over_delta)
    mesh, g, _ = _fem1d_profile(a, a_over_delta, H0=H0)

    # FEM loss: J_phi = -dH/dr ; P' = 1/(2 sigma) int |J|^2 * 2 pi r dr
    Jfem = -grad(g)[0]
    P_fem = (1.0 / (2 * SIGMA)) * 2 * math.pi * Integrate(Norm(Jfem) ** 2 * x, mesh)

    # Reference: same integral with the EXACT Bessel current J = -H0 k I1(kr)/I0(ka)
    nq = 4000
    I0ka = iv(0, k * a)
    P_ref = 0.0
    for i in range(nq):
        r = (i + 0.5) / nq * a
        J = -H0 * k * iv(1, k * r) / I0ka
        P_ref += (abs(J) ** 2) * r
    P_ref *= (1.0 / (2 * SIGMA)) * 2 * math.pi * (a / nq)

    assert abs(P_fem - P_ref) / P_ref < 1e-3, "P_fem=%.6e P_ref=%.6e" % (P_fem, P_ref)


# ---------------------------------------------------------------------------
# (B) genuine 2-D scalar Helmholtz eddy solve on an UNSTRUCTURED disk mesh.
# ---------------------------------------------------------------------------
def _disk_mesh(a, maxh):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    geo = SplineGeometry()
    geo.AddCircle((0, 0), a, bc="outer")
    return Mesh(geo.GenerateMesh(maxh=maxh))


def test_2d_disk_reproduces_bessel():
    """A real 2-D unstructured scalar-Helmholtz eddy solve on a disk reproduces the
    Bessel interior field and center shielding factor to FEM accuracy."""
    a, a_over_delta, H0 = 0.01, 3.0, 1.0
    k, _, _ = _k_omega(a, a_over_delta)
    mesh = _disk_mesh(a, a / 30.0)
    fes = H1(mesh, order=4, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    A = BilinearForm(grad(u) * grad(v) * dx + (k * k) * u * v * dx).Assemble()
    g = GridFunction(fes)
    g.Set(CF(H0), BND)
    r = g.vec.CreateVector()
    r.data = -A.mat * g.vec
    g.vec.data += A.mat.Inverse(fes.FreeDofs()) * r

    # center shielding factor
    cen = abs(complex(g(mesh(0.0, 0.0))))
    an_cen = abs(1.0 / iv(0, k * a))
    assert abs(cen - an_cen) / an_cen < 1e-2, "2D center %.5f vs Bessel %.5f" % (cen, an_cen)

    # a mid-radius sample point (radially symmetric -> any angle)
    rs = 0.5 * a
    fem = complex(g(mesh(rs, 0.0)))
    an = _bessel_Hz(rs, a, k)
    assert abs(fem - an) / abs(an) < 1e-2, "2D r/a=0.5 %s vs Bessel %s" % (fem, an)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
