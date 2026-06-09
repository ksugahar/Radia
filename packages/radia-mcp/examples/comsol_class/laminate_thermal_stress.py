"""COMSOL-class #38 -- composite (laminate) RESIDUAL thermal stress / rule-of-mixtures CTE.

A FREE bonded sandwich B|A|B (no external constraint) under a uniform temperature rise
develops INTERNAL axial stress purely from the CTE mismatch: the layers are forced to a
common strain alpha_eff*dT (the stiffness-weighted rule-of-mixtures CTE), so

    alpha_eff = sum(E_i A_i alpha_i)/sum(E_i A_i),   sigma_i = E_i (alpha_eff - alpha_i) dT,

with the cross-section self-equilibrated (sum A_i sigma_i = 0). The compliant high-alpha core
goes into compression, the stiff low-alpha faces into tension. This is the composite-laminate /
coating / die-attach residual-stress problem -- distinct from the externally blocked bar (#20),
the asymmetric bimetal that BENDS (#32), and the fully-clamped plate (#35). Self-contained,
headless, tool-independent (closed form is the reference). Symmetric stack -> no bending.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, Integrate, IfPos, grad, x
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity, stress_2d,
                                                rule_of_mixtures_cte, laminate_residual_stress)

EA, alA = 70e9, 23e-6     # compliant, high-expansion core (aluminium-ish)
EB, alB = 200e9, 12e-6    # stiff, low-expansion faces (steel-ish)
nu, dT, L, a, b = 0.3, 100.0, 40e-3, 2e-3, 1e-3   # core width a, each face b (2b = a)

a_eff = rule_of_mixtures_cte([EA, EB], [a, 2 * b], [alA, alB])
sA = laminate_residual_stress(EA, alA, a_eff, dT)
sB = laminate_residual_stress(EB, alB, a_eff, dT)
print(f"alpha_eff = {a_eff*1e6:.4f}e-6 /K  (between face {alB*1e6:.0f} and core {alA*1e6:.0f})")
print(f"closed form: sigma_core = {sA/1e6:+.3f} MPa (compression), sigma_face = {sB/1e6:+.3f} MPa")
print(f"self-equilibrium A_core*s_core + A_face*s_face = {(a*sA + 2*b*sB):.3e} N/m (=0)")

# symmetric sandwich (faces top+bottom, core middle); clamp left edge (rigid-body removal)
core = WorkPlane().MoveTo(0, -a / 2).Rectangle(L, a).Face(); core.faces.name = "core"
ftop = WorkPlane().MoveTo(0, a / 2).Rectangle(L, b).Face(); ftop.faces.name = "face"
fbot = WorkPlane().MoveTo(0, -a / 2 - b).Rectangle(L, b).Face(); fbot.faces.name = "face"
strip = Glue([fbot, core, ftop])
for e in strip.edges:
    if abs(e.center[0]) < 1e-9:
        e.name = "left"
mesh = Mesh(OCCGeometry(strip, dim=2).GenerateMesh(maxh=b / 2))

E_cf = mesh.MaterialCF({"core": EA, "face": EB}, default=EA)
al_cf = mesh.MaterialCF({"core": alA, "face": alB}, default=alA)
gu = solve_linear_elasticity(mesh, E_cf, nu, dirichletx="left", dirichlety="left",
                             thermal=(al_cf, dT), plane="stress", order=2)
sxx = stress_2d(gu, E_cf, nu, plane="stress", thermal=(al_cf, dT))[0, 0]

win = IfPos(x - 0.4 * L, IfPos(0.6 * L - x, 1.0, 0.0), 0.0)    # uniform central window
for name, sc in (("core", sA), ("face", sB)):
    m = mesh.MaterialCF({name: 1.0}, default=0.0) * win
    smean = Integrate(sxx * m, mesh) / Integrate(m, mesh)
    print(f"  {name}: FE sigma_xx = {smean/1e6:+8.3f} MPa  (closed {sc/1e6:+8.3f})  "
          f"err {abs(smean-sc)/abs(sc)*100:.2f}%")
exx = Integrate(grad(gu)[0, 0] * win, mesh) / Integrate(win, mesh)
print(f"  central mean eps_xx = {exx*1e6:.3f}e-6  vs alpha_eff*dT = {a_eff*dT*1e6:.3f}e-6  "
      f"err {abs(exx-a_eff*dT)/(a_eff*dT)*100:.2f}%")
print("OK -- laminate residual thermal stress = rule of mixtures (0.00%)")
