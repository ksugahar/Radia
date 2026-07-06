"""order=1 axisymmetric magnetostatics: V-DOF path and uniform-field guards.

The historical order-1 (P1 triangle / Q1 quad) failure was the symbolic A=psi
reconstruction for B_z = dA/dr + A/r.  The production solver avoids that path:
solve_axi_magnetostatic dispatches order=1 to the V-DOF custom-BFI form
(K_ij = 2pi/mu r_i r_j INT grad(psi_i).grad(psi_j)/r, i.e. the flux-function /
FEMM linear-element form).

The current contract is therefore positive: order=1 must converge on the
magnetized sphere and its uniform-field representation must remain small and
refinement-improving, while order=2 remains the more accurate baseline.
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


def test_p1_uniform_field_representation_is_convergent():
    """The current order=1 path must not regress to the old non-convergent
    symbolic A=psi behaviour; order=2 should still be clearly sharper."""
    e1_coarse = _bz_repr_error(1, 0.2)
    e1_fine = _bz_repr_error(1, 0.05)
    e2_coarse = _bz_repr_error(2, 0.2)
    e2_fine = _bz_repr_error(2, 0.05)
    print(f"  order1 B_z repr err: h=0.2 {e1_coarse:.3f}  h=0.05 {e1_fine:.3f}")
    print(f"  order2 B_z repr err: h=0.2 {e2_coarse:.3e}  h=0.05 {e2_fine:.3e}")
    # P1: small enough for the order=1 production path and refinement-improving.
    assert e1_coarse < 0.10, f"expected P1 coarse <10%, got {e1_coarse}"
    assert e1_fine < 0.03, f"expected P1 fine <3%, got {e1_fine}"
    assert e1_fine < 0.5 * e1_coarse, "P1 should improve by at least 2x"
    # P2: smaller and convergent.
    assert e2_fine < 1e-3, f"expected P2 fine <1e-3, got {e2_fine}"
    assert e2_fine < e2_coarse, "P2 should converge"
    assert e2_coarse < e1_coarse, "P2 should beat P1 on the coarse mesh"
    assert e2_fine < 0.2 * e1_fine, "P2 fine error should be much smaller than P1"
    print("[OK] P1 uniform-field representation is convergent; P2 remains sharper")


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
    test_p1_uniform_field_representation_is_convergent()
    test_order2_sphere_baseline()
    print("\nAll order=1 V-DOF / uniform-field guard tests passed.")
