r"""Synchronous-machine Ld / Lq (inductance saliency) -- regression test.

A SALIENT iron rotor (rectangle a x b, mu_r=1000, NO magnet) in an iron stator ring with a
stator coil (N=100, I=1 A).  Coil self-inductance L = lambda/I depends on rotor orientation:
long axis along the coil-flux axis -> Ld (high); rotated 90 deg -> Lq (low).  Ld != Lq is
the saliency (basis of reluctance torque + MTPA).

Meshed in METRES -- a current-driven B is NOT scale-invariant (unlike the PM cases), so
L = solve.inductance_2d(...) is in henries directly.  radia-ngsolve vs a stored regression
reference (independently cross-validated): Ld 0.20 %, Lq 0.05 %, saliency Ld/Lq 0.26 %.
Reference values stored below -> CI-portable.

LESSON: inductance is sensitive to the IN-COIL field (coil self-field / internal inductance)
-- a coarse coil mesh undershoots L by ~2 %; refine the coil region (~7 elements across the
wire) to converge.
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

from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue
from radia_mcp.radia_ngsolve.solve import (reluctivity, solve_planar_magnetostatic,
                                           inductance_2d)

Rbore, Rsout, Rout = 14e-3, 24e-3, 70e-3
rc, rw = 12e-3, 1.0e-3
A_LONG, A_SHORT = 18e-3, 6e-3
DEPTH, NTURNS, CURRENT, MUR = 1e-3, 100, 1.0, 1000.0
LD_REF, LQ_REF = 3.60516e-5, 2.94840e-5        # stored regression references [H], 1 mm stack


def _mesh(a, b):
    rotor = MoveTo(-a / 2, -b / 2).Rectangle(a, b).Face(); rotor.faces.name = "rotor"
    coilpos = WorkPlane().Circle(0, rc, rw).Face(); coilpos.faces.name = "coilpos"
    coilneg = WorkPlane().Circle(0, -rc, rw).Face(); coilneg.faces.name = "coilneg"
    disk_bore = WorkPlane().Circle(0, 0, Rbore).Face()
    gap = disk_bore - rotor - coilpos - coilneg; gap.faces.name = "air"
    disk_so = WorkPlane().Circle(0, 0, Rsout).Face()
    stator = disk_so - disk_bore; stator.faces.name = "stator"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    outer = box - disk_so; outer.faces.name = "air"
    rotor.faces.maxh = 1.0e-3; gap.faces.maxh = 0.8e-3
    coilpos.faces.maxh = 0.15e-3; coilneg.faces.maxh = 0.15e-3; stator.faces.maxh = 2.5e-3
    return Mesh(OCCGeometry(Glue([rotor, coilpos, coilneg, gap, stator, outer]),
                            dim=2).GenerateMesh(maxh=8e-3))


def _inductance(a, b):
    mesh = _mesh(a, b)
    nu = reluctivity(mesh, {"rotor": MUR, "stator": MUR})
    j0 = NTURNS * CURRENT / (math.pi * rw * rw)
    Jz = CoefficientFunction([{"coilpos": j0, "coilneg": -j0}.get(m, 0.0)
                              for m in mesh.GetMaterials()])
    with TaskManager():
        A = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
        return inductance_2d(A, mesh, "coilpos", "coilneg", NTURNS, CURRENT, depth=DEPTH)


def validate_ldlq_saliency():
    Ld = _inductance(A_LONG, A_SHORT)
    Lq = _inductance(A_SHORT, A_LONG)
    print(f"[Ld/Lq] radia Ld={Ld:.4e} Lq={Lq:.4e} (Ld/Lq={Ld/Lq:.3f}) | "
          f"ref Ld={LD_REF:.4e} Lq={LQ_REF:.4e} (Ld/Lq={LD_REF/LQ_REF:.3f})")
    assert abs(Ld - LD_REF) / LD_REF < 0.03, f"Ld {Ld:.3e} vs ref {LD_REF:.3e}"
    assert abs(Lq - LQ_REF) / LQ_REF < 0.03, f"Lq {Lq:.3e} vs ref {LQ_REF:.3e}"
    assert Ld > Lq, "salient rotor must give Ld > Lq"
    sal, salref = Ld / Lq, LD_REF / LQ_REF
    assert abs(sal - salref) / salref < 0.015, f"saliency {sal:.3f} vs ref {salref:.3f}"


if __name__ == "__main__":
    validate_ldlq_saliency()
    print("[OK] Ld/Lq saliency validated against the reference.")
