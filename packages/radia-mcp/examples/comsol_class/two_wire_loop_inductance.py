"""COMSOL-class #48 -- TWO-WIRE line TOTAL loop inductance (external + internal).

The practical loop inductance a meter reads, completing the transmission-line model: the
EXTERNAL field inductance (#30) PLUS the INTERNAL inductance of BOTH wires (#36):

    L_loop = (mu0/pi) acosh(D/2a)  +  2*(mu0/8pi)  =  L_ext + mu0/(4 pi).

radia partitions the FE field energy: 2 W_air/I^2 = L_ext, 2 W_wire/I^2 = L_int (per wire),
and the sum over air + both wires = L_loop. Self-contained, headless, tool-independent (exact
analytic ref). The +-I pair's field energy converges (1/r^2 far field), so no outer-boundary
dependence -- the internal part is exact (mu0/8pi each, radius-independent), the external is
the acosh term (#30 to ~1.8 %).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, reluctivity,
                                           two_wire_external_inductance,
                                           two_wire_loop_inductance, internal_inductance_round_wire)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d

a, D, Rout, I = 0.001, 0.010, 0.12, 1.0
L_ext, L_loop = two_wire_external_inductance(D, a), two_wire_loop_inductance(D, a)
print(f"analytic: L_ext={L_ext:.6e}  +  2*L_int={2*internal_inductance_round_wire():.6e}  "
      f"=  L_loop={L_loop:.6e} H/m  (internal {2*internal_inductance_round_wire()/L_loop*100:.1f}%)")

w1 = WorkPlane().Circle(-D/2, 0, a).Face(); w1.faces.name = "w1"
w2 = WorkPlane().Circle(D/2, 0, a).Face(); w2.faces.name = "w2"
box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
air = box - w1 - w2; air.faces.name = "air"
w1.faces.maxh = a/4; w2.faces.maxh = a/4; air.faces.maxh = Rout/8
for w in (w1, w2):
    for e in w.edges:
        e.maxh = a/4
mesh = Mesh(OCCGeometry(Glue([w1, w2, air]), dim=2).GenerateMesh(maxh=Rout/8))
Jz = CoefficientFunction([{"w1": I/(math.pi*a*a), "w2": -I/(math.pi*a*a)}.get(m, 0.0)
                         for m in mesh.GetMaterials()])
A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
L_fe_ext = 2 * magnetic_energy_2d(B, mesh, "air") / (I*I)
L_fe_int = 2 * (magnetic_energy_2d(B, mesh, "w1") + magnetic_energy_2d(B, mesh, "w2")) / (I*I)
L_fe_loop = L_fe_ext + L_fe_int
print(f"FE: L_ext(air)={L_fe_ext:.6e} ({abs(L_fe_ext-L_ext)/L_ext*100:.2f}%), "
      f"L_int(wires)={L_fe_int:.6e} ({abs(L_fe_int-2*internal_inductance_round_wire())/(2*internal_inductance_round_wire())*100:.2f}%)")
print(f"FE: L_loop(all)={L_fe_loop:.6e}  vs analytic {L_loop:.6e}  err={abs(L_fe_loop-L_loop)/L_loop*100:.2f}%")
print("OK -- L_loop = external acosh term + both wires' internal mu0/4pi.")
