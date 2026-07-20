"""FEMM scalar-solver extensions: convection (Robin) BC for heat, and AC
current flow (complex conductivity). Validated against closed forms.

  Robin/convection: 1D slab, Dirichlet T0 at x=0, convection (h, T_inf) at x=L,
    insulated sides. Series thermal resistance ->
        T_L = T_inf + (T0 - T_inf)/(1 + h L / k),  q = (T0-T_inf)/(L/k + 1/h).

  AC current flow: coaxial lossy dielectric, complex admittance per length
        Y = 2 pi (sigma + j w eps) / ln(b/a) = G + j w C.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from scipy.optimize import brentq
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, X, Y
from radia_mcp.radia_ngsolve.scalar_fem2d import (
    EPS0, SIGMA_SB, solve_thermal, solve_current_flow_ac, admittance_ac)


def check_robin_convection():
    k, L, Wy, T0, h, Tinf = 50.0, 0.10, 0.02, 100.0, 200.0, 20.0
    slab = MoveTo(0, 0).Rectangle(L, Wy).Face()
    slab.faces.name = "wall"
    slab.edges.Min(X).name = "hot"
    slab.edges.Max(X).name = "conv"
    # top/bottom (Min/Max Y) left unnamed -> natural Neumann (insulated)
    mesh = Mesh(OCCGeometry(slab, dim=2).GenerateMesh(maxh=Wy / 3))

    T = solve_thermal(mesh, CoefficientFunction(k), {"hot": T0},
                      convection={"conv": (h, Tinf)})

    TL_exact = Tinf + (T0 - Tinf) / (1.0 + h * L / k)
    q_exact = (T0 - Tinf) / (L / k + 1.0 / h)
    TL_fem = T(mesh(0.9999 * L, Wy / 2))
    Tmid_fem = T(mesh(0.5 * L, Wy / 2))
    Tmid_exact = T0 - q_exact * (0.5 * L) / k
    eL = (TL_fem - TL_exact) / TL_exact
    em = (Tmid_fem - Tmid_exact) / Tmid_exact
    print(f"  [Robin] T_L FEM={TL_fem:.4f} exact={TL_exact:.4f} ({100*eL:+.3f}%); "
          f"T_mid FEM={Tmid_fem:.4f} exact={Tmid_exact:.4f} ({100*em:+.3f}%)")
    assert abs(eL) < 5e-3 and abs(em) < 5e-3, "convection BC off >0.5%"


def check_ac_current_flow():
    a, b = 1e-3, 5e-3
    sigma, eps = 1e-3, 4.0 * EPS0
    omega = 2 * math.pi * 1e6

    outer = WorkPlane().Circle(0, 0, b).Face()
    hole = WorkPlane().Circle(0, 0, a).Face()
    ann = outer - hole
    ann.faces.name = "diel"
    ann.edges.Nearest((a, 0)).name = "inner"
    ann.edges.Nearest((b, 0)).name = "outer"
    mesh = Mesh(OCCGeometry(ann, dim=2).GenerateMesh(maxh=b / 18))

    V = solve_current_flow_ac(mesh, sigma, eps, omega,
                              {"inner": 1.0, "outer": 0.0})
    Y = admittance_ac(V, mesh, sigma, eps, omega, 1.0)

    G_ex = 2 * math.pi * sigma / math.log(b / a)
    C_ex = 2 * math.pi * eps / math.log(b / a)
    Y_ex = complex(G_ex, omega * C_ex)
    eG = (Y.real - Y_ex.real) / Y_ex.real
    eB = (Y.imag - Y_ex.imag) / Y_ex.imag
    print(f"  [AC I-flow] G FEM={Y.real:.6e} exact={Y_ex.real:.6e} ({100*eG:+.3f}%); "
          f"wC FEM={Y.imag:.6e} exact={Y_ex.imag:.6e} ({100*eB:+.3f}%)")
    assert abs(eG) < 1e-2 and abs(eB) < 1e-2, "AC admittance off >1%"


def check_radiation():
    k, L, Wy, T0, eps, Tinf = 20.0, 0.05, 0.02, 800.0, 0.8, 300.0  # T in Kelvin
    slab = MoveTo(0, 0).Rectangle(L, Wy).Face()
    slab.faces.name = "wall"
    slab.edges.Min(X).name = "hot"
    slab.edges.Max(X).name = "rad"
    mesh = Mesh(OCCGeometry(slab, dim=2).GenerateMesh(maxh=Wy / 3))

    T = solve_thermal(mesh, CoefficientFunction(k), {"hot": T0},
                      radiation={"rad": (eps, Tinf)})

    # conduction flux k(T0-T_L)/L  ==  radiation eps*sigma*(T_L^4 - Tinf^4)
    TL_exact = brentq(lambda t: k * (T0 - t) / L - eps * SIGMA_SB * (t**4 - Tinf**4),
                      Tinf, T0)
    TL_fem = T(mesh(0.9999 * L, Wy / 2))
    e = (TL_fem - TL_exact) / TL_exact
    print(f"  [Radiation] T_L FEM={TL_fem:.3f} K  exact(brentq)={TL_exact:.3f} K "
          f"({100*e:+.3f}%)")
    assert abs(e) < 5e-3, f"radiation BC off {100*e:.3f}% (>0.5%)"


def main():
    with TaskManager():
        check_robin_convection()
        check_radiation()
        check_ac_current_flow()
    print("\n[OK] scalar extensions validated: convection (Robin) + radiation (T^4) "
          "BC + AC current flow (complex admittance) vs closed forms.")


if __name__ == "__main__":
    main()
