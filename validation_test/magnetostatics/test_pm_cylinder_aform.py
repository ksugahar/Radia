"""
Permanent-magnet source in the A-formulation (the PM-machine field 数式).

A rigid permanent magnet (relative permeability ~1, magnetization M) enters the 2-D
vector-potential magnetostatic problem through the magnetization source term: from
curl H = 0 with H = (1/mu0) curl A - M, the weak form is

    int (1/mu0) grad(A_z) . grad(v) dx = int ( M_x dv/dy - M_y dv/dx ) dx .

Forged on a transversely magnetized circular cylinder (radius a, M = M x_hat, infinite in
z), which has the exact closed-form solution
  * INTERIOR: uniform B = mu0 M / 2 x_hat  (the 2-D transverse demag factor is 1/2),
  * EXTERIOR: a 2-D line dipole; on the magnetization axis  B_x(d,0) = mu0 M a^2 / (2 d^2)
    (which equals mu0 M/2 at the pole d=a, i.e. continuous).

Tests: the interior field equals mu0 M/2 and is uniform, the exterior on-axis field matches
the 2-D dipole law, and B scales linearly with M. Pure NGSolve (Bash-robust).
"""
import math
import pytest

pytestmark = [pytest.mark.slow, pytest.mark.usefixtures("ngsolve_taskmanager")]

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, Integrate,
                     grad, dx, CF)

MU0 = 4e-7 * math.pi
A_RAD = 0.01          # cylinder radius [m]
R = 0.30             # far boundary [m]


def _solve(M, maxh=0.0025):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    geo = SplineGeometry()
    geo.AddRectangle((-R, -R), (R, R), bc="outer", leftdomain=1)
    geo.AddCircle((0, 0), A_RAD, leftdomain=2, rightdomain=1, bc="magsurf")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "mag")
    geo.SetDomainMaxH(2, A_RAD / 12)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    fes = H1(mesh, order=4, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm((1.0 / MU0) * grad(u) * grad(v) * dx).Assemble()
    Mx = CF([M if m == "mag" else 0.0 for m in mesh.GetMaterials()])   # M = M x_hat
    f = LinearForm(Mx * grad(v)[1] * dx).Assemble()                    # int M_x dv/dy
    A = GridFunction(fes)
    A.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return mesh, A


def _B(A):
    return grad(A)[1], -grad(A)[0]      # (Bx, By) = (dA/dy, -dA/dx)


def test_internal_field_is_half_mu0_M():
    """Interior B is uniform and equals mu0 M / 2 x_hat."""
    M = 1.0e6
    mesh, A = _solve(M)
    Bx, By = _B(A)
    dom = mesh.Materials("mag")
    area = Integrate(CF(1), mesh, definedon=dom)
    bx = Integrate(Bx, mesh, definedon=dom) / area
    by = Integrate(By, mesh, definedon=dom) / area
    target = MU0 * M / 2.0
    assert abs(bx - target) / target < 0.03, "interior Bx=%.4f vs mu0 M/2=%.4f" % (bx, target)
    assert abs(by) / target < 0.02, "interior By not ~0: %.4e" % by
    # uniformity: spread of Bx over the magnet is small
    std = math.sqrt(Integrate((Bx - bx) ** 2, mesh, definedon=dom) / area)
    assert std / target < 0.05, "interior field not uniform: std/mean=%.3f" % (std / target)


def test_external_dipole_on_axis():
    """Exterior on-axis field matches the 2-D dipole law B_x = mu0 M a^2 / (2 d^2)."""
    M = 1.0e6
    mesh, A = _solve(M)
    Bx, _ = _B(A)
    for d in (0.02, 0.03):
        bx = Bx(mesh(d, 0.0))
        an = MU0 * M * A_RAD ** 2 / (2.0 * d ** 2)
        assert abs(bx - an) / an < 0.08, "d=%.3f: Bx=%.4e vs dipole %.4e" % (d, bx, an)


def test_linearity_in_M():
    """B scales linearly with the magnetization."""
    m1, A1 = _solve(5.0e5)
    m2, A2 = _solve(1.0e6)
    dom = m1.Materials("mag")
    b1 = Integrate(_B(A1)[0], m1, definedon=dom) / Integrate(CF(1), m1, definedon=dom)
    b2 = Integrate(_B(A2)[0], m2, definedon=dom) / Integrate(CF(1), m2, definedon=dom)
    assert abs((b2 / b1) - 2.0) / 2.0 < 0.02, "not linear in M: ratio %.3f" % (b2 / b1)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
