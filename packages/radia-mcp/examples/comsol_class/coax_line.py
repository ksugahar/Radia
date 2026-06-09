"""COMSOL-class #25 -- COAXIAL transmission line: inductance + characteristic impedance.

The magnetic companion to the coax capacitance (#6/#7).  Inner conductor (r<a) carries +I,
shield (b<r<c) carries -I; in the dielectric a<r<b the field is the exact Ampere field
B = mu0 I/(2 pi r), so the external inductance per unit length is

    L_ext = 2W/I^2 = (mu0/2pi) ln(b/a),   W = int_dielectric |B|^2/(2 mu0) dA.

With the coax capacitance C = 2 pi eps0/ln(b/a) this gives the textbook results

    Z0 = sqrt(L/C) = (1/2pi) sqrt(mu0/eps0) ln(b/a) = 60 ln(b/a) [ohm],   v = 1/sqrt(LC) = c.

radia-ngsolve gets L from the FE field energy (force.magnetic_energy_2d).  Self-contained,
headless, tool-independent (analytic reference).
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, reluctivity
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d, MU0

EPS0 = 8.8541878128e-12
a, b, c, Rout, I = 0.002, 0.008, 0.010, 0.030, 1.0      # metres, amperes


def build():
    da = WorkPlane().Circle(0, 0, a).Face(); da.faces.name = "inner"
    db = WorkPlane().Circle(0, 0, b).Face()
    diel = db - da; diel.faces.name = "diel"
    dc = WorkPlane().Circle(0, 0, c).Face()
    shield = dc - db; shield.faces.name = "shield"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    air = box - dc; air.faces.name = "air"
    da.faces.maxh = a/8; diel.faces.maxh = (b-a)/24       # fine near r=a (1/r energy peaks there)
    shield.faces.maxh = (c-b)/4; air.faces.maxh = Rout/8
    return Mesh(OCCGeometry(Glue([da, diel, shield, air]), dim=2).GenerateMesh(maxh=Rout/6))


mesh = build()
Jz = CoefficientFunction([{"inner": I/(math.pi*a*a),
                           "shield": -I/(math.pi*(c*c-b*b))}.get(m, 0.0)
                          for m in mesh.GetMaterials()])
A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
L_fe = 2.0 * magnetic_energy_2d(B, mesh, "diel") / (I * I)       # L = 2W/I^2 (dielectric energy)
L_cf = MU0/(2*math.pi)*math.log(b/a)
C_cf = 2*math.pi*EPS0/math.log(b/a)
Z0_fe, Z0_cf = math.sqrt(L_fe/C_cf), 60.0*math.log(b/a)
v_fe = 1.0/math.sqrt(L_fe*C_cf)

errL = abs(L_fe-L_cf)/L_cf*100
print(f"L_ext: FE={L_fe:.6e} H/m   closed={L_cf:.6e} H/m   err={errL:.3f} %")
print(f"Z0 (FE L, analytic C): {Z0_fe:.3f} ohm   vs 60 ln(b/a)={Z0_cf:.3f} ohm   ({100*abs(Z0_fe-Z0_cf)/Z0_cf:.2f} %)")
print(f"wave speed v=1/sqrt(LC)={v_fe:.4e} m/s   vs c={1/math.sqrt(MU0*EPS0):.4e} m/s")
print("OK" if errL < 1.0 else f"OFF ({errL:.2f}%)")
