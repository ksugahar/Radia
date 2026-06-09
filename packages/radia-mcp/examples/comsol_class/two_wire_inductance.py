"""COMSOL-class #30 -- TWO-WIRE LINE external inductance (transmission-line parameter).

The magnetic twin of the two-wire CAPACITANCE (#12) and the parallel of the coax line (#25).
Two parallel wires (radius a, centre separation D) carrying +-I form a 2D line-dipole whose
field energy CONVERGES (B ~ 1/r^2 far field), so the external inductance per unit length is

    L_ext = 2 W_ext/I^2 = (mu0/pi) acosh(D/2a)   [H/m]   (-> (mu0/pi) ln(D/a) for D >> a).

With the two-wire C = pi eps0/acosh(D/2a) this gives v = 1/sqrt(LC) = c and
Z0 = sqrt(L/C) = 120 acosh(D/2a) [ohm]. radia gets L from the FE field energy
(force.magnetic_energy_2d). Self-contained, headless, tool-independent (exact analytic ref).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, reluctivity,
                                           two_wire_external_inductance)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d, MU0

EPS0 = 8.8541878128e-12
a, D, Rout, I = 0.001, 0.010, 0.12, 1.0


def build():
    w1 = WorkPlane().Circle(-D/2, 0, a).Face(); w1.faces.name = "w1"
    w2 = WorkPlane().Circle(D/2, 0, a).Face(); w2.faces.name = "w2"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    air = box - w1 - w2; air.faces.name = "air"
    w1.faces.maxh = a/4; w2.faces.maxh = a/4; air.faces.maxh = Rout/8
    for w in (w1, w2):
        for e in w.edges:
            e.maxh = a/4                           # refine near the wire surface (1/r energy peak)
    return Mesh(OCCGeometry(Glue([w1, w2, air]), dim=2).GenerateMesh(maxh=Rout/8))


mesh = build()
Jz = CoefficientFunction([{"w1": I/(math.pi*a*a), "w2": -I/(math.pi*a*a)}.get(m, 0.0)
                         for m in mesh.GetMaterials()])
A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
L_fe = 2.0 * magnetic_energy_2d(B, mesh, "air") / (I*I)
L_cf = two_wire_external_inductance(D, a)
C_cf = math.pi*EPS0/math.acosh(D/(2*a))
errL = abs(L_fe-L_cf)/L_cf*100
print(f"L_ext: FE={L_fe:.6e} H/m   exact (mu0/pi)acosh(D/2a)={L_cf:.6e} H/m   err={errL:.2f}%")
print(f"Z0 (FE L, analytic C): {math.sqrt(L_fe/C_cf):.2f} ohm  vs 120 acosh={120*math.acosh(D/(2*a)):.2f} ohm")
print(f"wave speed v=1/sqrt(LC)={1/math.sqrt(L_fe*C_cf):.4e} m/s  vs c={1/math.sqrt(MU0*EPS0):.4e} m/s")
print("OK" if errL < 3.0 else f"OFF ({errL:.2f}%)")
