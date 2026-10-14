"""2D PLANAR time-harmonic eddy currents + current-driven circuit coupling
(FEMM prob2big analog), validated against the exact round-wire skin-effect
AC resistance (Kelvin ber/bei functions).

Exercises radia_mcp.radia_ngsolve.solve.solve_planar_eddy (driven conductor via
a NumberSpace axial-E unknown) and force.ohmic_loss_2d. For a round wire of
radius a, conductivity sigma, at q = a*sqrt(w*mu0*sigma):
    Rac/Rdc = (q/2)*(ber*bei' - bei*ber')/(ber'^2 + bei'^2)
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from scipy.special import kelvin
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_planar_eddy
from radia_mcp.radia_ngsolve.force import ohmic_loss_2d

MU0 = 4e-7 * math.pi
A, SIGMA, I_TOT, Q = 1e-3, 5.8e7, 1.0, 4.0
OMEGA = Q * Q / (A * A * MU0 * SIGMA)
R_FAR = 12 * A


def rac_over_rdc_exact(q):
    be, ke, bep, kep = kelvin(q)
    ber, bei, berp, beip = be.real, be.imag, bep.real, bep.imag
    return (q / 2.0) * (ber * beip - bei * berp) / (berp ** 2 + beip ** 2)


def build_mesh():
    delta = math.sqrt(2.0 / (OMEGA * MU0 * SIGMA))
    wire = WorkPlane().Circle(0, 0, A).Face(); wire.faces.name = "wire"
    wire.maxh = delta / 4.0
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - wire; air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air, wire]), dim=2).GenerateMesh(maxh=R_FAR / 8))


def main():
    mesh = build_mesh()
    sigma = mesh.MaterialCF({"wire": SIGMA}, default=0.0)
    nu = CoefficientFunction(1.0 / MU0)
    with TaskManager():
        gfu = solve_planar_eddy(mesh, nu, sigma, OMEGA,
                                driven_region="wire", total_current=I_TOT, order=3)
        Az, Vc = gfu.components
        Ez = -1j * OMEGA * Az + Vc
        P = ohmic_loss_2d(Ez, mesh, sigma, region="wire")

    Rac = 2.0 * P / (I_TOT ** 2)
    Rdc = 1.0 / (SIGMA * math.pi * A * A)
    ratio_fem = Rac / Rdc
    ratio_exact = rac_over_rdc_exact(Q)
    err = (ratio_fem - ratio_exact) / ratio_exact
    print(f"  q={Q}  f={OMEGA/2/math.pi:.0f} Hz  delta/a={math.sqrt(2/(OMEGA*MU0*SIGMA))/A:.3f}")
    print(f"  Rac/Rdc  FEM={ratio_fem:.5f}  exact(Kelvin)={ratio_exact:.5f}  err={100*err:+.2f}%")

    assert abs(err) < 0.01, f"Rac/Rdc off by {100*err:.2f}% (>1%)"
    print(f"\n[OK] 2D planar eddy + current-driven circuit validated vs round-wire "
          f"Kelvin Rac/Rdc ({100*abs(err):.2f}%).")


def test_skin_slab_driven_surface():
    # #130: solve_planar_eddy DRIVEN-SURFACE mode (dirichlet_value=A0) reproduces the slab skin-effect
    # A_z(0)/A0 = 1/cosh(gamma a), gamma=(1+j)/delta -- the complex analytic, to machine precision.
    import cmath
    a, b, sigma, f, A0 = 0.03, 0.02, 1.0e7, 60.0, 1.0
    omega = 2 * math.pi * f
    delta = math.sqrt(2.0 / (omega * MU0 * sigma))
    gamma = (1 + 1j) / delta
    slab = MoveTo(-a, 0).Rectangle(2 * a, b).Face()
    slab.faces.name = "c"
    slab.edges.Nearest((-a, b / 2)).name = "driven"
    slab.edges.Nearest((a, b / 2)).name = "driven"
    mesh = Mesh(OCCGeometry(slab, dim=2).GenerateMesh(maxh=(2 * a) / 30))
    A = solve_planar_eddy(mesh, CoefficientFunction(1.0 / MU0), sigma, omega,
                          dirichlet="driven", dirichlet_value=A0, order=3)
    atten = complex(A(mesh(0.0, b / 2))) / A0
    atten_ref = 1.0 / cmath.cosh(gamma * a)
    assert abs(atten / atten_ref - 1.0) < 1e-3, (atten, atten_ref)


if __name__ == "__main__":
    main()
    test_skin_slab_driven_surface()
