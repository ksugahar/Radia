"""
Back-EMF / flux-linkage of a rotating PM rotor (the motor voltage 数式).

The generated-voltage half of the machine: a rotating permanent-magnet rotor links a
flux Phi(theta) through a stator search coil, and the back-EMF is e = -dPhi/dt =
-omega dPhi/dtheta. This composes the PM source (test_pm_cylinder_aform) with a flux-
linkage probe, and exploits linearity: with the magnet at angle theta the magnetization
is M0(cos theta, sin theta), so the field is

    A(theta) = cos(theta) A_x + sin(theta) A_y          (A_x, A_y solved ONCE),

where A_x, A_y are the vector potentials for M along x / y. The flux linked by a diametric
coil with sides at (+-d, 0) is the 2-D flux-linkage Phi = A_z(+d,0) - A_z(-d,0).

For the transversely magnetized cylinder (moment m = M0 pi a^2) the exterior is a 2-D
dipole, A_z(r,phi) = (mu0 m / (2 pi r)) sin(phi - theta), so the coil on the x-axis links

    Phi(theta) = -(mu0 M0 a^2 / d) sin(theta)
               (peak when M is ALONG the coil flux-axis y, i.e. theta=90 deg; ~0 at theta=0),
    e(theta)   = -dPhi/dtheta = (mu0 M0 a^2 / d) cos(theta)   (back-EMF, per unit omega,
               90 deg out of phase with the flux).

Tests: the peak flux linkage (M||y) matches mu0 M0 a^2/d, a quadrature magnet (M||x) links
~zero, Phi(theta) is a pure sinusoid, and -dPhi/dtheta is its (cosine) quadrature = the
back-EMF. Pure NGSolve (Bash-robust).
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx, CF)

pytestmark = pytest.mark.slow

MU0 = 4e-7 * math.pi
A_RAD = 0.01
D_COIL = 0.05         # coil side radius [m]
R = 0.40
M0 = 1.0e6
PHI_PEAK = MU0 * M0 * A_RAD ** 2 / D_COIL     # |Phi|_max = mu0 M0 a^2 / d


def _setup(maxh=0.0025):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    geo = SplineGeometry()
    geo.AddRectangle((-R, -R), (R, R), bc="outer", leftdomain=1)
    geo.AddCircle((0, 0), A_RAD, leftdomain=2, rightdomain=1, bc="magsurf")
    geo.SetMaterial(1, "air"); geo.SetMaterial(2, "mag")
    geo.SetDomainMaxH(2, A_RAD / 10)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    fes = H1(mesh, order=4, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm((1.0 / MU0) * grad(u) * grad(v) * dx).Assemble()
    ainv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    one = CF([1.0 if m == "mag" else 0.0 for m in mesh.GetMaterials()])
    # RHS for M=(M0,0) and M=(0,M0):  int (M_x dv/dy - M_y dv/dx)
    fx = LinearForm(M0 * one * grad(v)[1] * dx).Assemble()      # M along x
    fy = LinearForm(-M0 * one * grad(v)[0] * dx).Assemble()     # M along y
    Ax = GridFunction(fes); Ax.vec.data = ainv * fx.vec
    Ay = GridFunction(fes); Ay.vec.data = ainv * fy.vec
    phi1 = Ax(mesh(D_COIL, 0.0)) - Ax(mesh(-D_COIL, 0.0))       # linkage, M along x
    phi2 = Ay(mesh(D_COIL, 0.0)) - Ay(mesh(-D_COIL, 0.0))       # linkage, M along y
    return phi1, phi2


@pytest.fixture(scope="module")
def _phi():
    return _setup()


def _Phi(theta, phi1, phi2):
    return math.cos(theta) * phi1 + math.sin(theta) * phi2


def test_peak_flux_linkage(_phi):
    """Peak flux linkage (magnet ALONG the coil flux-axis, M||y) = mu0 M0 a^2 / d."""
    phi1, phi2 = _phi
    rel = abs(abs(phi2) - PHI_PEAK) / PHI_PEAK
    assert rel < 0.08, "peak linkage |phi2|=%.4e vs analytic %.4e (rel %.3f)" % (abs(phi2), PHI_PEAK, rel)


def test_quadrature_magnet_links_zero(_phi):
    """A magnet along x (perpendicular to the coil flux-axis) links ~zero flux."""
    phi1, phi2 = _phi
    assert abs(phi1) / PHI_PEAK < 0.02, "M||x linkage not ~0: phi1/peak=%.3e" % (phi1 / PHI_PEAK)


def test_flux_linkage_sinusoid(_phi):
    """Phi(theta) is a pure sinusoid: ~0 at theta=0, +-peak at +-90 deg (= phi2 sin theta)."""
    phi1, phi2 = _phi
    for deg in (0, 45, 90, 135, 180, 270):
        th = math.radians(deg)
        phi = _Phi(th, phi1, phi2)
        assert abs(phi - phi2 * math.sin(th)) < 0.02 * PHI_PEAK, \
            "theta=%d: Phi=%.4e not ~phi2 sin=%.4e" % (deg, phi, phi2 * math.sin(th))


def test_back_emf_is_quadrature_sinusoid(_phi):
    """Back-EMF e(theta) = -dPhi/dtheta traces the cosine quadrature of the flux,
    amplitude mu0 M0 a^2/d, 90 deg out of phase with Phi(theta)."""
    phi1, phi2 = _phi
    h = 1e-4
    worst = 0.0
    for deg in range(0, 360, 15):
        th = math.radians(deg)
        emf = -(_Phi(th + h, phi1, phi2) - _Phi(th - h, phi1, phi2)) / (2 * h)   # -dPhi/dtheta
        worst = max(worst, abs(emf - (-phi2 * math.cos(th))))                    # ~ -phi2 cos theta
    assert worst < 0.02 * PHI_PEAK, "back-EMF deviates from the cosine quadrature: worst=%.3e" % worst


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
