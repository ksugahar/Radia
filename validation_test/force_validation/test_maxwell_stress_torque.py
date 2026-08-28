"""
TORQUE extraction via the Maxwell stress tensor (the motor torque 数式).

Torque is the rotational complement of test_maxwell_stress_force: the motor torque
(cogging, Arkkio, GalFer) is the contour integral of r x (T . n),

    tau_z = closed_int ( x (T.n)_y - y (T.n)_x ) dl ,
    T_ij = (1/mu0)( B_i B_j - 1/2 delta_ij |B|^2 ).

Gated against the elementary LORENTZ COUPLE (no hand-derived constant -> fabrication-proof):
a 2-D magnetic line-dipole = two anti-parallel line currents +-I at +-(s/2) x_hat, placed in
a uniform applied field B0=(B0x,B0y) imposed by the Dirichlet lift A_z = B0x y - B0y x on
the outer boundary. The couple on the pair about its centre is
    tau_z = sum_i r_i x (I_i z_hat x B0) = s I B0x
(the mutual wire force is central -> contributes ZERO torque about the centre). So the MST
contour torque around the dipole must equal s I B0x: it scales with B0x, flips sign with
B0x, and VANISHES for a purely transverse field (B0x=0, dipole aligned).

Uses the same internal-interface fix as the force test (project B to VectorH1 before the
contour integral). Pure NGSolve (Bash-robust).
"""
import math
import pytest

pytestmark = pytest.mark.slow

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, Integrate,
                     VectorH1, grad, dx, x, y, CF)

MU0 = 4e-7 * math.pi
I = 100.0                # A (wire +I / -I)
S = 0.02                 # dipole arm (wire separation) [m]
RW = 0.002               # wire radius [m]
RC = 0.020               # MST contour radius (encloses both wires) [m]
R = 0.40                 # outer boundary [m]


def _solve(B0x, B0y, maxh=0.003):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    geo = SplineGeometry()
    geo.AddRectangle((-R, -R), (R, R), bc="outer", leftdomain=1)
    geo.AddCircle((0, 0), RC, leftdomain=2, rightdomain=1, bc="contour")
    geo.AddCircle(( S / 2, 0), RW, leftdomain=3, rightdomain=2, bc="wps")
    geo.AddCircle((-S / 2, 0), RW, leftdomain=4, rightdomain=2, bc="wns")
    geo.SetMaterial(1, "air"); geo.SetMaterial(2, "core")
    geo.SetMaterial(3, "wp");  geo.SetMaterial(4, "wn")
    geo.SetDomainMaxH(3, RW / 3); geo.SetDomainMaxH(4, RW / 3)
    geo.SetDomainMaxH(2, 0.0015)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    fes = H1(mesh, order=4, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm((1.0 / MU0) * grad(u) * grad(v) * dx).Assemble()
    Jdens = I / (math.pi * RW ** 2)
    Jz = CF([Jdens if m == "wp" else (-Jdens if m == "wn" else 0.0)
             for m in mesh.GetMaterials()])
    f = LinearForm(Jz * v * dx).Assemble()

    A = GridFunction(fes)
    A.Set(CF(B0x * y - B0y * x), definedon=mesh.Boundaries("outer"))   # uniform-field lift
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * A.vec
    A.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
    return mesh, A


def _torque_z(mesh, A):
    """tau_z on everything inside the contour, via the MST (B projected to VectorH1)."""
    Bgf = GridFunction(VectorH1(mesh, order=3))
    Bgf.Set(CF((grad(A)[1], -grad(A)[0])))
    Bx, By = Bgf[0], Bgf[1]
    B2 = Bx * Bx + By * By
    nx, ny = x / RC, y / RC
    Bn = Bx * nx + By * ny
    tnx = (1.0 / MU0) * (Bn * Bx - 0.5 * B2 * nx)
    tny = (1.0 / MU0) * (Bn * By - 0.5 * B2 * ny)
    integrand = x * tny - y * tnx
    return Integrate(integrand, mesh, definedon=mesh.Boundaries("contour"))


def test_torque_matches_lorentz_couple():
    """MST torque equals the Lorentz couple s I B0x."""
    B0x, B0y = 1.0e-3, 5.0e-4
    mesh, A = _solve(B0x, B0y)
    tau = _torque_z(mesh, A)
    tau_an = S * I * B0x
    rel = abs(tau - tau_an) / abs(tau_an)
    assert rel < 0.08, "MST tau=%.5e vs Lorentz %.5e (rel %.3f)" % (tau, tau_an, rel)


def test_torque_flips_sign_with_field():
    """Reversing B0x reverses the torque (restoring couple)."""
    mesh1, A1 = _solve(+1.0e-3, 0.0)
    mesh2, A2 = _solve(-1.0e-3, 0.0)
    t1, t2 = _torque_z(mesh1, A1), _torque_z(mesh2, A2)
    assert t1 * t2 < 0, "torque did not flip sign: %.4e, %.4e" % (t1, t2)
    assert abs(t1 + t2) / abs(t1) < 0.05, "magnitudes not symmetric: %.4e, %.4e" % (t1, t2)


def test_transverse_field_zero_torque():
    """A purely transverse field (B0x=0, dipole aligned) produces ~zero torque."""
    mesh, A = _solve(0.0, 1.0e-3)
    tau = _torque_z(mesh, A)
    ref = S * I * 1.0e-3            # the torque scale if the field were along x
    assert abs(tau) / ref < 0.05, "aligned-dipole torque not ~zero: %.4e (scale %.4e)" % (tau, ref)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
