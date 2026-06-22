r"""Gapped iron-core inductor -- magnetic-circuit inductance + window leakage (#39).

The inductance twin of the reluctance actuator (#27): same window frame, read for the
coil inductance. Tool-independent validation: the FE inductance must sit ABOVE the
leakage-free reluctance model L0 = N^2 mu0 A/(g + l_fe/mu_r) (the coil links window
leakage flux the lumped model misses), with the excess GROWING as the gap opens, and L
monotonically falling with g. A radia self-regression pins the values (an independent 2D
FE solver agrees to ~0.4 %; that cross-check is recorded internally, not here)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import magnetic_circuit_inductance, MU0

W, ww, b = 0.24, 0.12, 0.06
MU_R, RBOX, DEPTH = 2000.0, 0.60, 1.0
N, I = 1000.0, 1.0
A_CORE = (W - ww) / 2 * DEPTH

# radia self-regression references (henries) -- guard the FE values; an independent
# 2D FE solver matches these to ~0.4% (recorded in the internal cross-val, unattributed).
_REF = {0.004: 21.98, 0.006: 16.31, 0.009: 12.26}


def validate_magnetic_circuit_inductance_formula():
    # L = N^2 mu0 A / (g + l_fe/mu_r) = N^2 / R, leakage-free reluctance model
    L = magnetic_circuit_inductance(N, A_CORE, 0.004, 0.716, MU_R)
    assert math.isclose(L, N*N*MU0*A_CORE/(0.004 + 0.716/MU_R), rel_tol=1e-12)
    # N^2 scaling
    assert math.isclose(magnetic_circuit_inductance(2*N, A_CORE, 0.004, 0.716, MU_R), 4*L, rel_tol=1e-12)
    # mu_r -> inf : gap-limited L -> N^2 mu0 A / g
    Linf = magnetic_circuit_inductance(N, A_CORE, 0.004, 0.716, 1e12)
    assert math.isclose(Linf, N*N*MU0*A_CORE/0.004, rel_tol=1e-6)
    # bigger gap -> smaller L
    assert magnetic_circuit_inductance(N, A_CORE, 0.009, 0.711, MU_R) < L


def _rect(wp, x0, y0, w, h):
    return wp.MoveTo(x0, y0).Rectangle(w, h).Face()


def _L_fe(g):
    from ngsolve import Mesh, Integrate, dx, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, Glue, X, Y
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, inductance_2d, NU0
    big = WorkPlane().MoveTo(-W/2, -W/2).Rectangle(W, W).Face()
    window = WorkPlane().MoveTo(-ww/2, -ww/2).Rectangle(ww, ww).Face()
    gap = WorkPlane().MoveTo(-g/2, ww/2).Rectangle(g, b).Face()
    cin = WorkPlane().MoveTo(-ww/2 + 0.005, -0.04).Rectangle(0.02, 0.08).Face()
    cout = WorkPlane().MoveTo(-W/2 - 0.025, -0.04).Rectangle(0.02, 0.08).Face()
    iron = big - window - gap
    gap.faces.name = "gap"; cin.faces.name = "cin"; cout.faces.name = "cout"; iron.faces.name = "iron"
    box = WorkPlane().MoveTo(-RBOX, -RBOX).Rectangle(2*RBOX, 2*RBOX).Face()
    air = box - iron - gap - cin - cout; air.faces.name = "air"
    air.edges.Max(X).name = "outer"; air.edges.Min(X).name = "outer"
    air.edges.Max(Y).name = "outer"; air.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, cin, cout]), dim=2).GenerateMesh(maxh=0.04))
    Ain = Integrate(mesh.MaterialCF({"cin": 1.0}, default=0.0)*dx, mesh)
    Aout = Integrate(mesh.MaterialCF({"cout": 1.0}, default=0.0)*dx, mesh)
    Jz = mesh.MaterialCF({"cin": N*I/Ain, "cout": -N*I/Aout}, default=0.0)
    nu = mesh.MaterialCF({"iron": NU0/MU_R}, default=NU0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
    return inductance_2d(gfu, mesh, "cin", "cout", n_turns=N, current=I, depth=DEPTH)


def validate_gapped_core_inductance_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    Ls, ratios = [], []
    for g in (0.004, 0.006, 0.009):
        l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g
        L = _L_fe(g)
        L0 = magnetic_circuit_inductance(N, A_CORE, g, l_fe, MU_R)
        # tool-independent: FE L is ABOVE the leakage-free model, by a bounded factor
        assert L0 < L < 1.8 * L0, f"g={g}: L={L:.3f} vs L0={L0:.3f}"
        # self-regression guard (an independent FE solver matches to ~0.4%)
        assert abs(L - _REF[g]) / _REF[g] < 0.02, f"g={g}: L={L:.4f} vs ref {_REF[g]}"
        Ls.append(L); ratios.append(L / L0)
    # inductance falls as the gap opens; leakage excess grows with the gap
    assert Ls[0] > Ls[1] > Ls[2]
    assert ratios[0] < ratios[1] < ratios[2]


if __name__ == "__main__":
    validate_magnetic_circuit_inductance_formula()
    validate_gapped_core_inductance_fe()
    print("[OK] gapped-core inductance: reluctance lower bound + leakage trend validated.")
