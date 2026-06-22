"""Regression guard + root-cause record for the H1Henrotte order=1 defect.

solve_axi_magnetostatic rejects order=1 because the order-1 (P1 triangle / Q1
quad) basis {1, r^2, z} CANNOT represent a uniform axial field:

    uniform B_z = B0   <=>   A_phi = B0 * r / 2   (odd in r)

The reconstructed B_z = dA/dr + A/r maps {1, r^2, z} -> {1/r, r, z/r}, which has
no constant term, so a uniform B_z is structurally unrepresentable.  This is a
REPRESENTATION defect, not a solver issue: projecting the EXACT A_phi into the
order-1 space and reconstructing B_z still gives ~74 % RMS error that does NOT
decrease under refinement, while order=2 converges O(h^2).

If a future radia.axifem build fixes the order-1 element (e.g. the documented
r_i*psi/r V-DOF transform), `test_p1_cannot_represent_uniform_field` will start
FAILING -- that is the signal to drop the order>=2 guard in
solve_axi_magnetostatic.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, GridFunction, grad, Integrate, CoefficientFunction,
                     TaskManager)
from ngsolve import x as r_cf
from netgen.occ import OCCGeometry, MoveTo
from radia.axifem import H1Henrotte
from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic

MU0 = 4e-7 * math.pi


def _box(h):
    f = MoveTo(0.2, -0.4).Rectangle(0.8, 0.8).Face()
    f.edges.name = "bnd"
    return Mesh(OCCGeometry(f, dim=2).GenerateMesh(maxh=h))


def _bz_repr_error(order, h, B0=0.5):
    """RMS error of reconstructed B_z when the EXACT A_phi=B0*r/2 is projected
    into the order-`order` H1Henrotte space (no solve)."""
    mesh = _box(h)
    fes = H1Henrotte(mesh, order=order)
    gfu = GridFunction(fes)
    gfu.Set(B0 * r_cf / 2.0)
    Bz = grad(gfu)[0] + gfu / r_cf
    num = Integrate((Bz - B0) ** 2 * r_cf, mesh)
    den = B0 ** 2 * Integrate(r_cf, mesh)
    return math.sqrt(num / den)


def test_order1_rejected():
    """solve_axi_magnetostatic must reject order=1 with an informative error."""
    mesh = _box(0.2)
    nu = CoefficientFunction(1.0 / MU0)
    raised = False
    try:
        solve_axi_magnetostatic(mesh, nu, order=1)
    except ValueError as e:
        raised = True
        assert "order >= 2" in str(e)
    assert raised, "order=1 should raise ValueError"
    print("[OK] order=1 rejected with informative ValueError")


def test_p1_cannot_represent_uniform_field():
    """ROOT CAUSE: order=1 B_z reconstruction is ~74% off and non-convergent;
    order=2 converges.  (If this fails, the C++ P1 element was fixed.)"""
    e1_coarse = _bz_repr_error(1, 0.2)
    e1_fine = _bz_repr_error(1, 0.05)
    e2_coarse = _bz_repr_error(2, 0.2)
    e2_fine = _bz_repr_error(2, 0.05)
    print(f"  order1 B_z repr err: h=0.2 {e1_coarse:.3f}  h=0.05 {e1_fine:.3f}")
    print(f"  order2 B_z repr err: h=0.2 {e2_coarse:.3e}  h=0.05 {e2_fine:.3e}")
    # P1: huge and NON-convergent (refining does not help)
    assert e1_coarse > 0.5, f"expected P1 coarse >50%, got {e1_coarse}"
    assert e1_fine > 0.5, f"expected P1 fine still >50%, got {e1_fine}"
    assert e1_fine > 0.5 * e1_coarse, "P1 should NOT converge (no >2x improvement)"
    # P2: small and convergent
    assert e2_fine < 1e-3, f"expected P2 fine <1e-3, got {e2_fine}"
    assert e2_fine < e2_coarse, "P2 should converge"
    print("[OK] root cause confirmed: P1 non-convergent, P2 converges")


def test_order2_sphere_baseline():
    """order=2 still validates on the magnetized sphere (unchanged baseline)."""
    A, MU_R, HC, RFAR = 1.0, 2.0, 3.0e5, 50.0
    B_exact = 2.0 * MU0 * MU_R * HC / (MU_R + 2.0)
    from netgen.occ import WorkPlane, Glue, X, Y
    disk = WorkPlane().Circle(0, 0, A).Face()
    half = MoveTo(0, -RFAR).Rectangle(RFAR, 2 * RFAR).Face()
    magnet = disk * half
    magnet.faces.name = "magnet"; magnet.maxh = 0.05
    magnet.edges.Min(X).name = "axis"
    air = MoveTo(0, -RFAR).Rectangle(RFAR, 2 * RFAR).Face() - magnet
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "outer"; air.edges.Max(Y).name = "outer"
    air.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, magnet]), dim=2).GenerateMesh(maxh=0.5))
    nu = CoefficientFunction(1.0 / (MU0 * mesh.MaterialCF({"magnet": MU_R}, default=1.0)))
    with TaskManager():
        gfu = solve_axi_magnetostatic(mesh, nu, magnets={"magnet": (HC, 90.0)}, order=2)
    Bz = grad(gfu)[0] + gfu / r_cf
    vol = Integrate(r_cf, mesh, definedon=mesh.Materials("magnet"))
    bz = Integrate(r_cf * Bz, mesh, definedon=mesh.Materials("magnet")) / vol
    err = abs((bz - B_exact) / B_exact)
    print(f"  order2 sphere <B_z>={bz:.6f} exact {B_exact:.6f} err {100*err:.3f}%")
    assert err < 5e-3, f"order2 sphere off by {100*err:.3f}%"
    print("[OK] order=2 sphere baseline preserved")


if __name__ == "__main__":
    test_order1_rejected()
    test_p1_cannot_represent_uniform_field()
    test_order2_sphere_baseline()
    print("\nAll order=1 guard / root-cause tests passed.")
