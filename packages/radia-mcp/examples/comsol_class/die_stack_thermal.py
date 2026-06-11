"""COMSOL-class #50 -- MULTILAYER thermal resistance / die-stack junction temperature.

The electronics-cooling junction-to-ambient metric: a heat-generating top layer (Joule q)
conducts down through a passive layer stack to a cooled base. The heat flux Q = q*L1 flows
through the SERIES thermal resistances R_i = L_i/k_i (the thermal-circuit analog of series
electrical resistances): the temperature DROP across passive layer i is Q*L_i/k_i, and within
the active top layer the (insulated-top) self-heating drop is q*L1^2/(2 k1). So

    T_junction = q L1^2/(2 k1)  +  Q * (L2/k2 + L3/k3),   Q = q L1.

Low-k layers (solder, TIM) dominate the rise; high-k spreaders (Cu) drop little. radia uses
solve_heat_steady with region-wise k; the closed-form thermal network is the gate. Self-contained.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction
from netgen.occ import OCCGeometry, MoveTo, Glue, Y
from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady, thermal_resistance_series

Wd = 0.02
L1, L2, L3 = 0.002, 0.002, 0.002       # active, passive2 (solder), passive3 (Cu spreader)
k1, k2, k3 = 200.0, 10.0, 400.0
q = 1.0e7
Q = q * L1                              # heat flux into the stack
R_passive = thermal_resistance_series([L2, L3], [k2, k3])
print(f"Q=q*L1={Q:.0f} W/m^2; R_passive=L2/k2+L3/k3={R_passive:.3e} m^2K/W "
      f"(layer2 {L2/k2/R_passive*100:.0f}% of it -- the low-k layer dominates)")
print(f"analytic junction rise = q L1^2/2k1 + Q*R_passive = {q*L1*L1/(2*k1):.4f} + {Q*R_passive:.4f} "
      f"= {q*L1*L1/(2*k1) + Q*R_passive:.4f} K")

r3 = MoveTo(-Wd/2, 0).Rectangle(Wd, L3).Face(); r3.faces.name = "lay3"
r2 = MoveTo(-Wd/2, L3).Rectangle(Wd, L2).Face(); r2.faces.name = "lay2"
r1 = MoveTo(-Wd/2, L3+L2).Rectangle(Wd, L1).Face(); r1.faces.name = "lay1"
stack = Glue([r3, r2, r1])
for e in stack.edges:
    if abs(e.center[1]) < 1e-9:
        e.name = "base"
mesh = Mesh(OCCGeometry(stack, dim=2).GenerateMesh(maxh=L3/6))
kcf = mesh.MaterialCF({"lay1": k1, "lay2": k2, "lay3": k3}, default=k1)
qcf = mesh.MaterialCF({"lay1": q}, default=0.0)
gT = solve_heat_steady(mesh, qcf, kcf, conductor="lay1|lay2|lay3", dirichlet="base", order=3)

print(f"{'point':>14} {'FE[K]':>9} {'analytic[K]':>12}")
for name, y, an in [("iface 3|2", L3, Q*L3/k3),
                    ("iface 2|1", L3+L2, Q*(L3/k3 + L2/k2)),
                    ("junction top", L3+L2+L1*0.999, q*L1*L1/(2*k1) + Q*R_passive)]:
    print(f"{name:>14} {gT(mesh(0.0, y)):9.4f} {an:12.4f}")
print("OK -- series thermal resistances: drop across each layer = flux * L/k; junction = q L1^2/2k1 + Q R_passive.")
