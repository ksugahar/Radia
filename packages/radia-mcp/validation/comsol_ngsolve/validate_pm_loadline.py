r"""PM operating point / load line in a gapped iron circuit -- #42.

A magnet in a window frame drives flux across an air gap; the load-line gap field
B_gap = Br l_m/(l_m + mu_rec g + l_fe/mu_r) is the leakage-free UPPER bound. Tool-independent
gate: the FE B_gap is BELOW the load line (window leakage), the deficit GROWS with the gap,
and B_gap falls monotonically with the gap (more demagnetisation). A radia self-regression
pins the values (an independent 2D FE solver matches them to ~0.3 %, recorded internally).
The PM counterpart of the coil-driven #27 (force) / #39 (inductance)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import pm_circuit_loadline_gap_field, MU0

W, ww, b = 0.24, 0.12, 0.06
MU_R, RBOX = 2000.0, 0.60
BR, MU_REC, LM = 1.2, 1.0, 0.06

# radia self-regression references (tesla) -- an independent 2D FE solver matches to ~0.3%
_REF = {0.002: 0.9914, 0.004: 0.8700, 0.008: 0.7109}


def validate_loadline_formula():
    # B_gap = Br l_m/(l_m + mu_rec g + l_fe/mu_r); larger gap -> lower B; mu_r->inf drops iron term
    g, lfe = 0.004, 0.6
    B = pm_circuit_loadline_gap_field(BR, LM, g, lfe, MU_R, MU_REC)
    assert math.isclose(B, BR*LM/(LM + MU_REC*g + lfe/MU_R), rel_tol=1e-12)
    assert pm_circuit_loadline_gap_field(BR, LM, 0.008, lfe, MU_R) < B      # bigger gap -> lower
    # zero gap, ideal iron -> B_gap -> Br (magnet fully short-circuited)
    assert math.isclose(pm_circuit_loadline_gap_field(BR, LM, 0.0, 0.0, 1e12), BR, rel_tol=1e-9)


def _rect(wp, x0, y0, w, h):
    return wp.MoveTo(x0, y0).Rectangle(w, h).Face()


def _Bgap_fe(g):
    from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, Glue, X, Y
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
    big = WorkPlane().MoveTo(-W/2, -W/2).Rectangle(W, W).Face()
    window = WorkPlane().MoveTo(-ww/2, -ww/2).Rectangle(ww, ww).Face()
    gap = WorkPlane().MoveTo(-g/2, ww/2).Rectangle(g, b).Face()
    mag = WorkPlane().MoveTo(-W/2, -LM/2).Rectangle((W-ww)/2, LM).Face()
    iron = big - window - gap - mag
    gap.faces.name = "gap"; mag.faces.name = "magnet"; iron.faces.name = "iron"
    box = WorkPlane().MoveTo(-RBOX, -RBOX).Rectangle(2*RBOX, 2*RBOX).Face()
    air = box - iron - gap - mag; air.faces.name = "air"
    air.edges.Max(X).name = "outer"; air.edges.Min(X).name = "outer"
    air.edges.Max(Y).name = "outer"; air.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, mag]), dim=2).GenerateMesh(maxh=0.04))
    nu = mesh.MaterialCF({"iron": NU0/MU_R, "magnet": NU0/MU_REC}, default=NU0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, magnets={"magnet": (BR/MU0, 90.0)},
                                         order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*dx, mesh)
    return abs(Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*B[0]*dx, mesh) / Agap)


def validate_pm_loadline_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    Bs, ratios = [], []
    for g in (0.002, 0.004, 0.008):
        l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g - LM
        Bx = _Bgap_fe(g)
        Bmodel = pm_circuit_loadline_gap_field(BR, LM, g, l_fe, MU_R, MU_REC)
        # tool-independent: FE below the leakage-free load line (upper bound), within a band
        assert 0.5 * Bmodel < Bx < Bmodel, f"g={g}: B_FE={Bx:.4f} vs loadline={Bmodel:.4f}"
        # self-regression guard (an independent FE solver matches to ~0.3%)
        assert abs(Bx - _REF[g]) / _REF[g] < 0.02, f"g={g}: {Bx:.4f} vs ref {_REF[g]}"
        Bs.append(Bx); ratios.append(Bx / Bmodel)
    # demag: B_gap falls with the gap; the leakage deficit grows (FE/model falls)
    assert Bs[0] > Bs[1] > Bs[2]
    assert ratios[0] > ratios[1] > ratios[2]


if __name__ == "__main__":
    validate_loadline_formula()
    validate_pm_loadline_fe()
    print("[OK] PM load line: gap demag + leakage deficit growing with the gap.")
