r"""AC internal-inductance roll-off of a round wire (skin effect) -- regression test (#45).

L_int(omega)/L_int_dc = (4/q)[ber ber'+bei bei']/[ber'^2+bei'^2] (Kelvin, q=sqrt(2)a/delta),
the inductive twin of the skin-effect Rac/Rdc. radia's `solve_planar_eddy` on the WIRE ALONE
(A_z=0 on the surface -> Z=Vc/I is the pure internal impedance) reproduces both Rac=Re(Z) and
L_int=Im(Z)/omega. L_int -> mu0/8pi (#36 DC limit) at low freq, rolls off ~4/q high. Pure
Kelvin closed form -> tool-independent."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (internal_inductance_round_wire,
                                           skin_effect_resistance_ratio,
                                           skin_effect_internal_inductance_ratio)

MU0 = 4e-7 * math.pi
A, SIGMA, I = 1e-3, 5.8e7, 1.0


def validate_kelvin_ratio_limits():
    pytest.importorskip("scipy")
    # low frequency: both ratios -> 1 (DC: uniform current, L_int = mu0/8pi)
    assert abs(skin_effect_resistance_ratio(0.01) - 1.0) < 1e-3
    assert abs(skin_effect_internal_inductance_ratio(0.01) - 1.0) < 1e-3
    # high frequency asymptotes (q = sqrt(2) a/delta): Rac/Rdc -> q/(2 sqrt2) + 1/4,
    # L_int/L_dc -> 2 sqrt2 / q  (current confined to the ~delta skin)
    q = 30.0
    assert abs(skin_effect_resistance_ratio(q) - (q / (2 * math.sqrt(2)) + 0.25)) / (q / (2 * math.sqrt(2))) < 0.02
    assert abs(skin_effect_internal_inductance_ratio(q) - 2 * math.sqrt(2) / q) / (2 * math.sqrt(2) / q) < 0.05
    # monotone: Rac rises, L_int falls with q
    assert skin_effect_resistance_ratio(1) < skin_effect_resistance_ratio(2) < skin_effect_resistance_ratio(4)
    assert skin_effect_internal_inductance_ratio(1) > skin_effect_internal_inductance_ratio(2) > skin_effect_internal_inductance_ratio(4)


def validate_ac_internal_inductance_fe():
    pytest.importorskip("scipy")
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.solve import solve_planar_eddy

    LINT_DC = internal_inductance_round_wire()
    RDC = 1.0 / (SIGMA * math.pi * A * A)

    def radia_RL(q):
        omega = q * q / (A * A * MU0 * SIGMA)
        delta = math.sqrt(2.0 / (omega * MU0 * SIGMA))
        wire = WorkPlane().Circle(0, 0, A).Face(); wire.faces.name = "wire"; wire.edges.name = "surf"
        wire.maxh = delta / 4.0
        mesh = Mesh(OCCGeometry(wire, dim=2).GenerateMesh(maxh=A / 8))
        sigma = mesh.MaterialCF({"wire": SIGMA}, default=0.0)
        with TaskManager():
            gfu = solve_planar_eddy(mesh, CoefficientFunction(1.0 / MU0), sigma, omega,
                                    driven_region="wire", total_current=I, order=4, dirichlet="surf")
        Az, Vc = gfu.components
        Z = complex(Vc(mesh(0.0, 0.0))) / I
        return Z.real / RDC, (Z.imag / omega) / LINT_DC

    ratios = []
    for q in (0.5, 2.0, 4.0, 8.0):
        rac_fe, lint_fe = radia_RL(q)
        lint_k = skin_effect_internal_inductance_ratio(q)
        rac_k = skin_effect_resistance_ratio(q)
        assert abs(lint_fe - lint_k) / lint_k < 5e-3, f"q={q}: L_int {lint_fe} vs {lint_k}"
        assert abs(rac_fe - rac_k) / rac_k < 1e-2, f"q={q}: Rac {rac_fe} vs {rac_k}"  # cross-check
        ratios.append(lint_fe)
    # low-q -> DC limit mu0/8pi (ratio ~1); roll-off monotone
    assert abs(ratios[0] - 1.0) < 0.01
    assert ratios[0] > ratios[1] > ratios[2] > ratios[3]


if __name__ == "__main__":
    validate_kelvin_ratio_limits()
    validate_ac_internal_inductance_fe()
    print("[OK] AC internal-inductance roll-off validated vs Kelvin (DC limit mu0/8pi + roll-off).")
