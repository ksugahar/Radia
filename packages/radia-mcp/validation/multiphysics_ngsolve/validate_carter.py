r"""Carter coefficient (slotting -> effective air-gap increase) -- regression test (#57).

A slot opening adds reluctance to the air gap, so it behaves magnetically larger: g_eff = k_C g,
k_C = tau_s/(tau_s - gamma g), gamma = (b_o/g)^2/(5 + b_o/g). Closed-form helper (tool-independent)
+ a scalar-magnetic-potential gap-permeance FE (k_C = P_smooth/P_slot)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import carter_coefficient, effective_air_gap


def validate_carter_closed_form():
    tau_s, g = 0.012, 0.001
    # no slot opening -> k_C = 1 (no correction)
    assert math.isclose(carter_coefficient(tau_s, g, 0.0), 1.0, rel_tol=1e-12)
    # k_C >= 1 and grows with the opening
    kcs = [carter_coefficient(tau_s, g, bo) for bo in (0.001, 0.002, 0.004, 0.006)]
    assert all(k >= 1.0 for k in kcs)
    assert all(b > a for a, b in zip(kcs, kcs[1:]))
    # effective gap = k_C * g
    bo = 0.004
    assert math.isclose(effective_air_gap(tau_s, g, bo), carter_coefficient(tau_s, g, bo) * g, rel_tol=1e-12)
    # textbook value (b_o/g = 4): gamma = 16/9, k_C = 12/(12 - 16/9) ~ 1.174
    assert math.isclose(carter_coefficient(0.012, 0.001, 0.004), 1.1739, abs_tol=2e-3)


def validate_carter_fe_permeance():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, grad, InnerProduct, Integrate, dx, CoefficientFunction, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.scalar_fem2d import solve_poisson_2d

    def fe_kc(tau_s, g, bo, ds=None, order=3):
        ds = ds or 12 * g
        h = tau_s / 2.0
        xw = h - bo / 2.0
        wp = WorkPlane()
        wp.MoveTo(0, 0).LineTo(h, 0).LineTo(h, g + ds).LineTo(xw, g + ds).LineTo(xw, g).LineTo(0, g).Close()
        face = wp.Face()
        face.edges.Nearest((h / 2, 0)).name = "rotor"
        face.edges.Nearest(((h + xw) / 2, g + ds)).name = "stator"
        face.edges.Nearest((xw, g + ds / 2)).name = "stator"
        face.edges.Nearest((xw / 2, g)).name = "stator"
        mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=g / 2))
        phi = solve_poisson_2d(mesh, CoefficientFunction(1.0), {"stator": 1.0, "rotor": 0.0}, order=order)
        return (h / g) / Integrate(InnerProduct(grad(phi), grad(phi)) * dx, mesh)

    with TaskManager():
        for (tau_s, g, bo) in [(0.012, 0.001, 0.002), (0.012, 0.001, 0.004), (0.020, 0.002, 0.008)]:
            kc_fe = fe_kc(tau_s, g, bo)
            kc_cf = carter_coefficient(tau_s, g, bo)
            assert math.isclose(kc_fe, kc_cf, rel_tol=5e-3), f"{bo}: FE {kc_fe} vs {kc_cf}"
            assert kc_fe > 1.0


if __name__ == "__main__":
    validate_carter_closed_form()
    validate_carter_fe_permeance()
    print("[OK] Carter coefficient: scalar-potential gap permeance == k_C formula to <0.5%.")
