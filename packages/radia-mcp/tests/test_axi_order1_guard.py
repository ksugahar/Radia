"""order=1 axisymmetric magnetostatics: the A=psi defect, and its V-DOF fix.

The order-1 (P1 triangle / Q1 quad) {1, r^2, z} basis used with the *symbolic
A=psi* reconstruction B_z = dA/dr + A/r CANNOT represent a uniform axial field:

    uniform B_z = B0   <=>   A_phi = B0 * r / 2   (odd in r)

maps {1, r^2, z} -> {1/r, r, z/r}, which has no constant term (root-cause test
below: ~74 % RMS error, non-convergent).  solve_axi_magnetostatic therefore does
NOT use the A=psi form at order 1; it dispatches to the V-DOF custom-BFI path
(K_ij = 2pi/mu r_i r_j INT grad(psi_i).grad(psi_j)/r, == the flux-function / FEMM
linear-element form), where  A = sum V_i r_i psi_i / r  DOES represent a uniform
B_z exactly (sum_i r_i^2 psi_i == r^2).  So order=1 now converges on the sphere
(FEMM-P1-like), validated below.
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
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia.axifem import H1Henrotte
from radia_mcp.radia_ngsolve.solve import (solve_axi_magnetostatic,
                                           axi_vdof_magnet_bz_average)

MU0 = 4e-7 * math.pi
_SPHERE = dict(A=1.0, MU_R=2.0, HC=3.0e5, RFAR=50.0)
_B_EXACT = 2.0 * MU0 * _SPHERE["MU_R"] * _SPHERE["HC"] / (_SPHERE["MU_R"] + 2.0)


def _box(h):
    f = MoveTo(0.2, -0.4).Rectangle(0.8, 0.8).Face()
    f.edges.name = "bnd"
    return Mesh(OCCGeometry(f, dim=2).GenerateMesh(maxh=h))


def _sphere_mesh(hmag, hfar=0.5):
    A, RFAR = _SPHERE["A"], _SPHERE["RFAR"]
    disk = WorkPlane().Circle(0, 0, A).Face()
    half = MoveTo(0, -RFAR).Rectangle(RFAR, 2 * RFAR).Face()
    magnet = disk * half
    magnet.faces.name = "magnet"; magnet.maxh = hmag; magnet.edges.Min(X).name = "axis"
    air = MoveTo(0, -RFAR).Rectangle(RFAR, 2 * RFAR).Face() - magnet
    air.faces.name = "air"; air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "outer"; air.edges.Max(Y).name = "outer"
    air.edges.Min(Y).name = "outer"
    return Mesh(OCCGeometry(Glue([air, magnet]), dim=2).GenerateMesh(maxh=hfar))


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


def test_order1_vdof_converges():
    """order=1 (P1) now SOLVES via the V-DOF custom-BFI path: the magnetized
    sphere <B_z> converges toward exact 2 mu0 mu_r Hc/(mu_r+2) under refinement
    (FEMM-P1-like), where the symbolic A=psi order-1 cannot represent it at all."""
    HC = _SPHERE["HC"]
    errs = []
    for hmag in (0.10, 0.05):
        mesh = _sphere_mesh(hmag)
        nu = CoefficientFunction(1.0 / (MU0 * mesh.MaterialCF({"magnet": _SPHERE["MU_R"]},
                                                              default=1.0)))
        with TaskManager():
            gfu = solve_axi_magnetostatic(mesh, nu, magnets={"magnet": (HC, 90.0)}, order=1)
            bz = axi_vdof_magnet_bz_average(mesh, nu, {"magnet": (HC, 90.0)}, gfu)["magnet"]
        err = abs((bz - _B_EXACT) / _B_EXACT)
        errs.append(err)
        print(f"  order1 sphere hmag={hmag}: <B_z>={bz:.6f} exact {_B_EXACT:.6f} err {100*err:.3f}%")
    assert errs[0] < 0.03, f"order1 coarse err {100*errs[0]:.2f}% too large"
    assert errs[1] < errs[0] + 1e-9, "order1 should converge (finer <= coarser)"
    assert errs[1] < 0.015, f"order1 fine err {100*errs[1]:.2f}% too large"
    print("[OK] order=1 V-DOF magnetostatic converges on the sphere")


def test_p1_cannot_represent_uniform_field():
    """ROOT CAUSE the V-DOF path exists to avoid: the *symbolic A=psi* order-1
    B_z reconstruction is ~74% off and non-convergent; order=2 converges."""
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
    test_order1_vdof_converges()
    test_p1_cannot_represent_uniform_field()
    test_order2_sphere_baseline()
    print("\nAll order=1 V-DOF / root-cause tests passed.")
