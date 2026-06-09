"""COMSOL-class #67 -- nonlinear magnetic-circuit operating point (saturating iron-frame electromagnet).

The saturating counterpart of the constant-mu reluctance circuits (#27 gap force, #39 inductance,
#42 PM load line). For a closed iron loop wound with N turns carrying I, Ampere's law gives the
operating flux density B from the iron B->H curve:

    NI = H_iron(B) * l_fe + (B/mu0) * gap .

With no gap every ampere-turn drops in the iron, so B sits on the BH KNEE; a gap adds the large
linear term B*gap/mu0 that de-saturates the iron fast. Here the iron follows H = 51 B + 2.5 B^15
[A/m] (mu_r ~ 15600 unsaturated, a hard saturation above ~1.5 T). The operating B ~ 1.56 T at
NI = 32 A-turns over a ~15.6 mm mean path is a stored regression reference for this electromagnet.
Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import magnetic_circuit_bh_operating_point, MU0

H_of_B = lambda B: 51.0 * B + 2.5 * B ** 15           # the iron BH curve [A/m], B in T
l_fe = 0.0156                                          # mean iron path of the closed frame [m]

print(f"iron BH: H = 51 B + 2.5 B^15  (mu_r(0) = {1.0/(MU0*51):.0f}, saturates above ~1.5 T)")
print("\nsaturation knee (closed frame, no gap):  NI[At]   B[T]   H(B)[A/m]   residual")
for NI in (8, 16, 32, 64, 128, 256):
    B = magnetic_circuit_bh_operating_point(NI, l_fe, H_of_B)
    res = H_of_B(B) * l_fe - NI
    star = "   <- operating point (1.56 T ref)" if NI == 32 else ""
    print(f"   {NI:6.0f}   {B:.4f}   {H_of_B(B):8.1f}   {res:+.1e}{star}")

# constant-mu limit reproduces the linear reluctance formula exactly
mu_r = 1000.0
lin = lambda B: B / (MU0 * mu_r)
for gap in (0.0, 1e-3):
    B = magnetic_circuit_bh_operating_point(100.0, l_fe, lin, gap=gap)
    B_exact = 100.0 / (l_fe / (MU0 * mu_r) + gap / MU0)
    print(f"\nconstant mu_r={mu_r:.0f}, NI=100, gap={gap*1e3:.0f} mm: B={B:.5f} T vs analytic "
          f"{B_exact:.5f} T (err {abs(B/B_exact-1):.1e})")

# a gap de-saturates the iron (air-gap reluctance dominates)
print("\ngap sweep (saturating iron, NI=200):  gap[mm]   B[T]")
for g in (0.0, 0.2e-3, 0.5e-3, 1.0e-3):
    print(f"   {g*1e3:5.1f}   {magnetic_circuit_bh_operating_point(200.0, l_fe, H_of_B, gap=g):.4f}")

print("\nOK -- closed-frame operating B sits on the BH knee (1.56 T at NI=32, l_fe=15.6 mm); "
      "constant-mu limit = analytic; a gap de-saturates.")
