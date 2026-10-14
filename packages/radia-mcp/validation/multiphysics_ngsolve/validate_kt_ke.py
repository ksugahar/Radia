r"""Kt = Ke (torque constant == back-EMF constant) -- the PM-machine identity. Regression test.

Energy conservation: a current I in a winding with PM flux linkage lambda_m(theta) feels the
co-energy torque T = I dlambda_m/dtheta, so Kt=(T(I)-T(0))/I = dlambda_m/dtheta = Ke. radia:
2-pole PM rotor + search coil + iron ring (the back-EMF geometry), at theta=90 deg. Ke from
pm_machine_constant (PM-only flux-linkage FD); Kt from the eggshell torque (current minus
cogging)/I. MESHED IN METRES (current-driven B is not scale-invariant); depth cancels in Kt=Ke.
Tool-independent (self-consistency), no commercial data.
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (reluctivity, solve_planar_magnetostatic,
                                           coil_flux_linkage_2d, pm_machine_constant)
from radia_mcp.radia_ngsolve.force import eggshell_torque_2d

MU0 = 4e-7 * math.pi
HC = 1.0 / MU0
Rr, Rbore, Rsout, Rout = 0.010, 0.014, 0.024, 0.070
rc, rw, N = 0.012, 0.001, 100
RB_IN, RB_OUT = 0.0102, 0.0108
I_COIL, DTH = 10.0, 5.0


def _mesh():
    rotor = WorkPlane().Circle(0, 0, Rr).Face(); rotor.faces.name = "rotor"
    cp = WorkPlane().Circle(0, rc, rw).Face(); cp.faces.name = "coilpos"
    cn = WorkPlane().Circle(0, -rc, rw).Face(); cn.faces.name = "coilneg"
    bore = WorkPlane().Circle(0, 0, Rbore).Face()
    gap = bore - rotor - cp - cn; gap.faces.name = "air"
    iron = WorkPlane().Circle(0, 0, Rsout).Face() - bore; iron.faces.name = "iron"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    outer = box - WorkPlane().Circle(0, 0, Rsout).Face(); outer.faces.name = "air"
    rotor.faces.maxh = 0.8e-3; gap.faces.maxh = 0.4e-3
    cp.faces.maxh = 0.15e-3; cn.faces.maxh = 0.15e-3; iron.faces.maxh = 2e-3
    return Mesh(OCCGeometry(Glue([rotor, cp, cn, gap, iron, outer]), dim=2).GenerateMesh(maxh=8e-3))


def validate_kt_equals_ke():
    mesh = _mesh()
    nu = reluctivity(mesh, {"iron": 1000.0})
    area = Integrate(mesh.MaterialCF({"coilpos": 1.0}, default=0.0) * dx, mesh)
    Jz = mesh.MaterialCF({"coilpos": N*I_COIL/area, "coilneg": -N*I_COIL/area}, default=0.0)

    def lam_tau(th, cur):
        A = solve_planar_magnetostatic(mesh, nu, magnets={"rotor": (HC, th)},
                                       Jz=(Jz if cur else None), order=3, dirichlet="outer")
        B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
        return (coil_flux_linkage_2d(A, mesh, "coilpos", "coilneg", N, depth=1.0),
                eggshell_torque_2d(B, mesh, (0.0, 0.0), RB_IN, RB_OUT))

    with TaskManager():
        _, T0 = lam_tau(90.0, False)
        lam_p, _ = lam_tau(90.0 + DTH, False)
        lam_n, _ = lam_tau(90.0 - DTH, False)
        _, T_I = lam_tau(90.0, True)
    Ke = pm_machine_constant(lam_p, lam_n, math.radians(DTH))
    Kt = (T_I - T0) / I_COIL
    rel = abs(abs(Kt) - abs(Ke)) / abs(Ke)
    print(f"[Kt=Ke] Ke={Ke:.4e}  Kt={Kt:.4e}  ({100*rel:.2f}%)")
    assert rel < 0.04, f"Kt {Kt:.3e} vs Ke {Ke:.3e} ({100*rel:.1f}%)"


if __name__ == "__main__":
    validate_kt_equals_ke()
    print("[OK] Kt = Ke (torque constant = back-EMF constant) validated.")
