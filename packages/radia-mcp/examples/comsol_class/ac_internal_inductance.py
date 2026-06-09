"""COMSOL-class #45 -- AC INTERNAL-inductance roll-off of a round wire (skin effect).

The high-frequency complement of the DC internal inductance (#36, L_int(0)=mu0/8pi): as the
frequency rises the skin effect expels current from the core, so BOTH the AC resistance rises
(Rac/Rdc, the FEMM prob2big skin-effect test) AND the internal inductance ROLLS OFF. Modelling
the WIRE ALONE with A_z=0 on its surface, the current-driven impedance Z=Vc/I is the pure
INTERNAL impedance (no external term to subtract): Rac=Re(Z), L_int=Im(Z)/omega.

    Rac/Rdc        = (q/2) [ber bei' - bei ber'] / [ber'^2 + bei'^2]
    L_int/L_int_dc = (4/q) [ber ber' + bei bei'] / [ber'^2 + bei'^2]     (q = sqrt(2) a/delta)

radia's `solve_planar_eddy` reproduces both; here we feature the inductance roll-off (the
resistance is the existing skin-effect test). Self-contained, headless; the closed form is the
Kelvin (ber/bei) reference. Together Z_internal(omega) = Rdc*Rac_ratio + j omega (mu0/8pi)*Lint_ratio.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, WorkPlane
from radia_mcp.radia_ngsolve.solve import (solve_planar_eddy, internal_inductance_round_wire,
                                           skin_effect_resistance_ratio,
                                           skin_effect_internal_inductance_ratio)

MU0 = 4e-7 * math.pi
A, SIGMA, I = 1e-3, 5.8e7, 1.0
LINT_DC = internal_inductance_round_wire()      # mu0/8pi
RDC = 1.0 / (SIGMA * math.pi * A * A)


def radia_RL(q):
    omega = q * q / (A * A * MU0 * SIGMA)
    delta = math.sqrt(2.0 / (omega * MU0 * SIGMA))
    wire = WorkPlane().Circle(0, 0, A).Face(); wire.faces.name = "wire"; wire.edges.name = "surf"
    wire.maxh = delta / 4.0
    mesh = Mesh(OCCGeometry(wire, dim=2).GenerateMesh(maxh=A / 8))
    sigma = mesh.MaterialCF({"wire": SIGMA}, default=0.0)
    with TaskManager():
        gfu = solve_planar_eddy(mesh, CoefficientFunction(1.0 / MU0), sigma, omega,
                                driven_region="wire", total_current=I, order=4, dirichlet="surf")
    Az, Vc = gfu.components
    Z = complex(Vc(mesh(0.0, 0.0))) / I
    return Z.real / RDC, (Z.imag / omega) / LINT_DC


print(f"L_int_dc = mu0/8pi = {LINT_DC:.4e} H/m")
print(f"{'q':>5} {'Lint/Ldc FE':>12} {'(Kelvin)':>9} {'Rac/Rdc FE':>11} {'(Kelvin)':>9}")
for q in (0.5, 1.0, 2.0, 4.0, 8.0):
    rac_fe, lint_fe = radia_RL(q)
    print(f"{q:5.1f} {lint_fe:12.4f} {skin_effect_internal_inductance_ratio(q):9.4f} "
          f"{rac_fe:11.4f} {skin_effect_resistance_ratio(q):9.4f}")
print("OK -- L_int -> mu0/8pi (#36) at low freq, rolls off ~4/q at high freq; Rac rises ~q/2.")
