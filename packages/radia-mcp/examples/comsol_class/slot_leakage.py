"""COMSOL-class #54 -- slot-leakage inductance of a machine slot (slotted-iron FE).

The leakage flux that crosses a stator/rotor slot from tooth to tooth and links the conductor
WITHOUT reaching the air gap -- the leakage reactance every winding carries. For a rectangular
slot the specific permeance is

    lambda_s = h_c/(3 w_s) + h_o/w_o ,   L = mu0 N^2 l_stk lambda_s.

The **1/3 factor** is exact: a uniform-current conductor builds the cross-slot MMF LINEARLY (0
at the slot bottom to NI at the conductor top), and integrating B^2 over that triangular field
gives 1/3. Modelled as the well-posed 1-D-equivalent BVP -- the conductor fills the slot with
uniform J, the tooth walls (mu_r->inf) become NEUMANN edges (B enters perpendicular) and the
yoke flux line at the slot top is DIRICHLET A=0. Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, slot_leakage_permeance,
                                           slot_leakage_inductance, NU0, MU0)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d

I = 1.0


def slot_L_fe(w_s, h_c, h_o=0.0, w_o=None, order=3):
    """FE slot-leakage inductance per axial metre (single turn): L' = 2W/I^2."""
    cond = MoveTo(0, 0).Rectangle(w_s, h_c).Face(); cond.faces.name = "cond"
    if h_o > 0:
        w_o = w_o or w_s
        opening = MoveTo((w_s - w_o) / 2, h_c).Rectangle(w_o, h_o).Face(); opening.faces.name = "air"
        shape = Glue([cond, opening])
    else:
        shape = cond
    shape.edges.Max(Y).name = "top"          # yoke flux line: Dirichlet A=0
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=min(w_s, w_o or w_s) / 6))
    Jz = CoefficientFunction([I / (w_s * h_c) if m == "cond" else 0.0 for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=order, dirichlet="top")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    return 2 * sum(magnetic_energy_2d(B, mesh, m) for m in set(mesh.GetMaterials())) / I**2


print("slot-leakage permeance:  lambda_s = h_c/(3 w_s) + h_o/w_o   (L' = mu0 lambda_s per metre)")
print("case                         FE L' [uH/m]  closed form    err")
with TaskManager():
    # 1) closed slot (conductor fills it): exact 1/3 factor
    w_s, h_c = 0.004, 0.020
    L = slot_L_fe(w_s, h_c); cf = MU0 * slot_leakage_permeance(h_c, w_s)
    print(f"closed slot w={w_s*1e3:.0f}/h={h_c*1e3:.0f}mm   {L*1e6:8.4f}     {cf*1e6:8.4f}    {(L/cf-1)*100:+.3f}%")
    # 2) partial fill, same-width opening (no step): still exact
    h_o = 0.010
    L = slot_L_fe(w_s, h_c, h_o, w_s); cf = MU0 * slot_leakage_permeance(h_c, w_s, h_o, w_s)
    print(f"+ same-width opening h_o={h_o*1e3:.0f}mm  {L*1e6:8.4f}     {cf*1e6:8.4f}    {(L/cf-1)*100:+.3f}%")
    # 3) practical semi-closed slot: narrowed tooth-tip opening -> shoulder fringing
    w_s, h_c, h_o, w_o = 0.006, 0.024, 0.004, 0.002
    L = slot_L_fe(w_s, h_c, h_o, w_o); cf = MU0 * slot_leakage_permeance(h_c, w_s, h_o, w_o)
    print(f"semi-closed w_o={w_o*1e3:.0f}<w_s={w_s*1e3:.0f}mm  {L*1e6:8.4f}     {cf*1e6:8.4f}    {(L/cf-1)*100:+.3f}%  (FE>cf: tooth-shoulder fringing)")

# inductance of a real slot (turns, stack length)
L_slot = slot_leakage_inductance(0.020, 0.004, axial_length=0.10, turns=12)
print(f"\nexample slot: 12 turns, 100 mm stack, 4x20 mm -> L_leak = {L_slot*1e6:.2f} uH")
print("OK -- 1/3 conductor term + same-width opening EXACT; narrowed opening = lower bound (fringing).")
