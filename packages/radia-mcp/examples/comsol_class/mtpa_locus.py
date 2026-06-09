"""COMSOL-class #26 -- MTPA (Maximum-Torque-Per-Ampere) locus of a salient PM machine.

The dq control capstone: from FE-extracted machine parameters -- the PM flux linkage
lambda_m, and the d-/q-axis inductances Ld/Lq (saliency) -- find the stator-current
angle that maximises torque per ampere. dq torque (motor convention):

    T = (3/2) p [ lambda_m iq + (Ld - Lq) id iq ],   id = -I sin g, iq = I cos g.

This wires three radia-ngsolve FE solves (lambda_m from a PM-rotor flux linkage; Ld and
Lq from a salient-iron-rotor coil inductance along the d-/q-axis) into
``solve.mtpa_operating_point`` and sweeps the current magnitude.  Self-contained,
headless, tool-independent (the MTPA optimality is checked analytically in
tests/test_motor_mtpa.py against a numeric argmax; here we just exercise the pipeline).
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue
from radia_mcp.radia_ngsolve.solve import (reluctivity, solve_planar_magnetostatic,
                                           inductance_2d, coil_flux_linkage_2d,
                                           dq_torque, mtpa_operating_point)

Rbore, Rsout, Rout = 14e-3, 24e-3, 70e-3
rc, rw = 12e-3, 1.0e-3
A_LONG, A_SHORT = 18e-3, 6e-3
DEPTH, NTURNS, CURRENT, MUR = 1e-3, 100, 1.0, 1000.0
HC_MAG, MUR_REC, P = 4.0e5, 1.05, 2          # rotor PM coercivity / recoil mu_r / pole pairs


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
    Jz = CoefficientFunction([{"coilpos": j0, "coilneg": -j0}.get(m, 0.0) for m in mesh.GetMaterials()])
    with TaskManager():
        A = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
    return inductance_2d(A, mesh, "coilpos", "coilneg", NTURNS, CURRENT, depth=DEPTH)


def _lambda_m():
    # PM rotor magnetized along the d-axis = +x (the coil-flux axis: the (0,+/-rc) coil
    # links X-directed flux, so Ld uses the rotor long-axis along x). NO coil current
    # -> coil flux linkage = lambda_m. (A +y magnet gives A_z ~ x/r^2 = 0 on the x=0
    #  coil axis -> would read lambda_m ~ 0; the easy axis must be the d-axis.)
    mesh = _mesh(A_LONG, A_SHORT)
    nu = reluctivity(mesh, {"rotor": MUR_REC, "stator": MUR})
    with TaskManager():
        A = solve_planar_magnetostatic(mesh, nu, magnets={"rotor": (HC_MAG, 0.0)},
                                       order=3, dirichlet="outer")
    return abs(coil_flux_linkage_2d(A, mesh, "coilpos", "coilneg", NTURNS, depth=DEPTH))


Ld = _inductance(A_LONG, A_SHORT)
Lq = _inductance(A_SHORT, A_LONG)
lam = _lambda_m()
print(f"FE-extracted: lambda_m={lam:.4e} Wb   Ld={Ld:.4e} H   Lq={Lq:.4e} H   (Ld/Lq={Ld/Lq:.3f}, p={P})")
print(f"{'I[A]':>6} {'gamma[deg]':>10} {'id[A]':>8} {'iq[A]':>8} {'T[N.m]':>11} {'T/Tq':>7}")
for I in (2.0, 5.0, 10.0, 20.0, 40.0):
    g, idc, iqc, T = mtpa_operating_point(lam, Ld, Lq, I, P)
    Tq = dq_torque(lam, Ld, Lq, 0.0, I, P)            # pure q-axis (no saliency exploitation)
    print(f"{I:6.1f} {math.degrees(g):10.2f} {idc:8.3f} {iqc:8.3f} {T:11.4e} {T/Tq:7.4f}")
print("OK -- FE-extracted (lambda_m, Ld, Lq) -> MTPA locus via solve.mtpa_operating_point")
