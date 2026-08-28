"""
CAPSTONE: PM rotor torque-vs-angle curve (PM source + applied field + MST torque).

This composes three separately-forged blocks into the motor torque signature:
  * the PM magnetization source in the A-formulation   (test_pm_cylinder_aform),
  * a uniform applied field via the Dirichlet lift       (test_maxwell_stress_torque),
  * torque extraction via the Maxwell stress tensor      (test_maxwell_stress_torque).

A rigid magnetized cylinder (moment per length m = M pi a^2 along x) sits in a uniform
applied field B0 rotated to angle phi. Rotating the FIELD is equivalent to rotating the
rotor, so no moving mesh is needed -- and the operator is unchanged across phi, so a single
factorization is reused for the whole sweep. The torque on the magnet is the elementary
dipole-in-field couple

    tau_z(phi) = m x B0 = M pi a^2 B0 sin(phi),

so the MST contour torque must trace a clean sine: zero when aligned (phi=0, pi), peak at
phi=pi/2, restoring sign, and magnitude M pi a^2 B0 sin(phi) at every angle. Pure NGSolve.
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, Integrate,
                     VectorH1, grad, dx, x, y, CF)

pytestmark = pytest.mark.slow

MU0 = 4e-7 * math.pi
A_RAD = 0.01
RC = 0.020
R = 0.30
M0 = 1.0e6           # magnetization [A/m]
B0 = 0.5             # applied flux density [T]
M_MOMENT = M0 * math.pi * A_RAD ** 2     # dipole moment per length
TAU_PEAK = M_MOMENT * B0                  # = M pi a^2 B0


def _setup(maxh=0.0025):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    geo = SplineGeometry()
    geo.AddRectangle((-R, -R), (R, R), bc="outer", leftdomain=1)
    geo.AddCircle((0, 0), RC, leftdomain=2, rightdomain=1, bc="contour")
    geo.AddCircle((0, 0), A_RAD, leftdomain=3, rightdomain=2, bc="magsurf")
    geo.SetMaterial(1, "air"); geo.SetMaterial(2, "core"); geo.SetMaterial(3, "mag")
    geo.SetDomainMaxH(3, A_RAD / 10); geo.SetDomainMaxH(2, 0.0015)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    fes = H1(mesh, order=4, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm((1.0 / MU0) * grad(u) * grad(v) * dx).Assemble()
    Mx = CF([M0 if m == "mag" else 0.0 for m in mesh.GetMaterials()])     # M = M0 x_hat
    f = LinearForm(Mx * grad(v)[1] * dx).Assemble()                        # PM source
    ainv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")        # reused across phi
    return mesh, fes, a, f, ainv


def _torque_at(mesh, fes, a, f, ainv, phi):
    B0x, B0y = B0 * math.cos(phi), B0 * math.sin(phi)
    A = GridFunction(fes)
    A.Set(CF(B0x * y - B0y * x), definedon=mesh.Boundaries("outer"))      # uniform-field lift
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * A.vec
    A.vec.data += ainv * res
    # MST torque (project B to VectorH1 for a reliable internal-interface trace)
    Bgf = GridFunction(VectorH1(mesh, order=3))
    Bgf.Set(CF((grad(A)[1], -grad(A)[0])))
    Bx, By = Bgf[0], Bgf[1]
    B2 = Bx * Bx + By * By
    nx, ny = x / RC, y / RC
    Bn = Bx * nx + By * ny
    tnx = (1.0 / MU0) * (Bn * Bx - 0.5 * B2 * nx)
    tny = (1.0 / MU0) * (Bn * By - 0.5 * B2 * ny)
    return Integrate(x * tny - y * tnx, mesh, definedon=mesh.Boundaries("contour"))


@pytest.fixture(scope="module")
def _sweep():
    mesh, fes, a, f, ainv = _setup()
    degs = [0, 30, 90, 150, 180]
    return {d: _torque_at(mesh, fes, a, f, ainv, math.radians(d)) for d in degs}


def test_sine_law(_sweep):
    """tau(phi) matches the dipole couple M pi a^2 B0 sin(phi) at every angle."""
    for d, tau in _sweep.items():
        an = TAU_PEAK * math.sin(math.radians(d))
        assert abs(tau - an) <= 0.08 * TAU_PEAK, \
            "phi=%d deg: tau=%.4f vs M pi a^2 B0 sin=%.4f" % (d, tau, an)


def test_zero_torque_when_aligned(_sweep):
    """Aligned (phi=0) and anti-aligned (phi=180) give ~zero torque."""
    assert abs(_sweep[0]) / TAU_PEAK < 0.03, "phi=0 torque not ~0: %.4f" % _sweep[0]
    assert abs(_sweep[180]) / TAU_PEAK < 0.03, "phi=180 torque not ~0: %.4f" % _sweep[180]


def test_peak_and_restoring_sign(_sweep):
    """Peak (aligning) torque at phi=90, and the torque is restoring for 0<phi<180."""
    assert abs(_sweep[90] - TAU_PEAK) / TAU_PEAK < 0.08, "phi=90 peak %.4f vs %.4f" % (_sweep[90], TAU_PEAK)
    assert _sweep[30] > 0 and _sweep[90] > 0 and _sweep[150] > 0, \
        "torque not restoring (should align M toward B0): %s" % _sweep
    # symmetry of the sine about 90 deg: tau(30) ~ tau(150) ~ 0.5 * peak
    assert abs(_sweep[30] - _sweep[150]) / TAU_PEAK < 0.03, "sine not symmetric about 90 deg"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
