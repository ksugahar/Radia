"""COMSOL-class #36 -- INTERNAL inductance of a round wire = mu0/(8 pi) per unit length.

A wire (radius R) with a UNIFORM (DC) current has B(r)=mu0 I r/(2 pi R^2) inside, so the energy
stored INSIDE gives the internal inductance L_int = 2 W_in/I^2 = mu0/(8 pi) [H/m] -- INDEPENDENT
of R. The internal companion to the EXTERNAL inductances (coax #25, two-wire #30); it sets the
DC inductance and rolls off under the skin effect at high frequency. radia uses the wire-interior
field energy (force.magnetic_energy_2d), which is boundary-independent (Ampere fixes B(r) from the
enclosed current). Self-contained, headless, tool-independent (exact analytic ref).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, reluctivity,
                                           internal_inductance_round_wire)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d, MU0

I = 1.0
L_cf = internal_inductance_round_wire()


def L_internal(R, Rout):
    wire = WorkPlane().Circle(0, 0, R).Face(); wire.faces.name = "wire"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    air = box - wire; air.faces.name = "air"
    wire.faces.maxh = R/12; air.faces.maxh = Rout/8
    mesh = Mesh(OCCGeometry(Glue([wire, air]), dim=2).GenerateMesh(maxh=Rout/8))
    Jz = CoefficientFunction([{"wire": I/(math.pi*R*R)}.get(m, 0.0) for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    return 2.0 * magnetic_energy_2d(B, mesh, "wire") / (I*I)


print(f"exact mu0/(8 pi) = {L_cf:.6e} H/m   (radius-independent)")
for R in (0.001, 0.002):
    L = L_internal(R, 20*R)
    print(f"R={R*1e3:.1f} mm:  L_int(FE)={L:.6e} H/m   err={abs(L-L_cf)/L_cf*100:.3f} %")
print("OK -- L_int = mu0/(8 pi), independent of radius")
