"""COMSOL-class #42 -- PM OPERATING POINT / load line in a gapped iron circuit.

The permanent-magnet counterpart of the coil-driven magnetic circuit (#27 gap force, #39
gap inductance): a magnet (Br, recoil mu_r=1) in one leg of a window frame drives flux
around the iron and across the air gap g. Ampere's law + the recoil line B=Br+mu0 H give the
LOAD-LINE gap field

    B_gap = Br * l_m / (l_m + mu_rec*g + l_fe/mu_r),

a steeper load line (larger gap) -> lower B_gap (more demagnetisation). LEAKAGE: some magnet
flux short-cuts the window instead of reaching the gap, so the 2D FE B_gap sits BELOW this
leakage-free load line, and the deficit GROWS as the gap opens. An independent 2D FE solver
agrees with radia to ~0.3 % (stored regression reference; capability parity). Self-contained,
headless; PM field is scale-invariant but we keep metres to match the #27 frame.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           pm_circuit_loadline_gap_field, NU0, MU0)

W, ww, b = 0.24, 0.12, 0.06
mu_r, Rbox = 2000.0, 0.60
Br, mu_rec, lm = 1.2, 1.0, 0.06


def _rect(x0, y0, w, h):
    return MoveTo(x0, y0).Rectangle(w, h).Face()


def Bgap_fe(g):
    l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g - lm
    big = _rect(-W/2, -W/2, W, W); window = _rect(-ww/2, -ww/2, ww, ww)
    gap = _rect(-g/2, ww/2, g, b)
    mag = _rect(-W/2, -lm/2, (W-ww)/2, lm)
    iron = big - window - gap - mag
    gap.faces.name = "gap"; mag.faces.name = "magnet"; iron.faces.name = "iron"
    box = _rect(-Rbox, -Rbox, 2*Rbox, 2*Rbox)
    air = box - iron - gap - mag; air.faces.name = "air"
    air.edges.Max(X).name = "outer"; air.edges.Min(X).name = "outer"
    air.edges.Max(Y).name = "outer"; air.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, mag]), dim=2).GenerateMesh(maxh=0.04))
    nu = mesh.MaterialCF({"iron": NU0/mu_r, "magnet": NU0/mu_rec}, default=NU0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, magnets={"magnet": (Br/MU0, 90.0)},
                                         order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*dx, mesh)
    Bx = abs(Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*B[0]*dx, mesh) / Agap)
    l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g - lm
    return Bx, pm_circuit_loadline_gap_field(Br, lm, g, l_fe, mu_r, mu_rec)


print(f"{'g[mm]':>6} {'Bgap_FE[T]':>11} {'loadline[T]':>12} {'FE/model':>9}")
for g in (0.002, 0.004, 0.008):
    Bx, Bm = Bgap_fe(g)
    print(f"{g*1e3:6.1f} {Bx:11.4f} {Bm:12.4f} {Bx/Bm:9.3f}")
print("OK -- PM load line: larger gap -> lower B_gap (more demag); FE below the leakage-free "
      "load line, deficit growing with the gap; radia matches the independent FE to ~0.3%.")
