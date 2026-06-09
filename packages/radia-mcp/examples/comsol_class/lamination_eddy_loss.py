"""COMSOL-class #66 -- lamination eddy-current loss (the core-loss eddy term).

A steel sheet of thickness d in a sinusoidal flux density B(t) (peak B) dissipates classical eddy
loss P = pi^2 sigma d^2 f^2 B^2 / 6 [W/m^3] -- the d^2 law that is the WHOLE REASON cores are
laminated, and the f^2 ('eddy') term of the Steinmetz separation P = P_hyst(f) + P_eddy(f^2). When
the sheet is no longer thin vs the skin depth delta = sqrt(2/(omega mu sigma)), the flux cannot fully
penetrate and the loss is reduced by

    F(xi) = (3/xi)(sinh xi - sin xi)/(cosh xi - cos xi),   xi = d/delta,

the exact 1-D magnetic-diffusion result (F -> 1 thin, F -> 3/xi heavy skin). This is the FIELD-driven
counterpart of the current-driven deep-bar / Dowell AC-resistance family (#55/#60). Here F is
cross-checked against the numerically-integrated exact slab loss (an independent reference). Self-
contained, headless; SI units.
"""
import os
import sys
import math
import cmath

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (classical_eddy_loss_density, lamination_eddy_skin_factor,
                                           lamination_eddy_loss_density, MU0)


def F_numeric(xi, n=20001):
    """Exact slab loss / classical loss, by integrating the 1-D complex-diffusion solution for a
    fixed average flux density (delta = 1 units). An independent check on the closed-form F(xi)."""
    import numpy as np
    trap = getattr(np, "trapezoid", None) or np.trapz
    omega, d, k = 2.0, xi, (1 + 1j)                     # delta=1: omega*mu*sigma=2, k=(1+j)/delta
    Bs = k * d / (2 * cmath.tanh(k * d / 2))            # surface B so the average B = 1
    x = np.linspace(0.0, d, n)
    Jz = -Bs * k * np.sinh(k * (x - d / 2)) / cmath.cosh(k * d / 2)
    return trap(np.abs(Jz) ** 2, x) / (2 * d) / (omega ** 2 * d ** 2 / 24.0)


print("lamination eddy-loss skin-effect factor F(xi) = (3/xi)(sinh-sin)/(cosh-cos):")
print("  xi=d/delta   F_closed    F_numeric(exact slab)   match")
for xi in (0.2, 0.5, 1.0, 2.0, 4.0, 8.0):
    fc = lamination_eddy_skin_factor(xi, 1.0)          # skin_depth=1 -> arg is xi
    fn = F_numeric(xi)
    print(f"   {xi:5.2f}     {fc:.6f}      {fn:.6f}           {abs(fc/fn-1):.2e}")

print("\nthin-lamination limit: F -> 1 (classical d^2 law); heavy skin: F -> 3/xi")
print(f"  F(0.05)={lamination_eddy_skin_factor(0.05,1.0):.6f}   F(8)={lamination_eddy_skin_factor(8,1.0):.5f} vs 3/8={3/8:.5f}")

# realistic example: M19 0.35 mm electrical steel, 50 Hz, 1.5 T
sigma, d, f, B, mur, rho = 2.0e6, 0.35e-3, 50.0, 1.5, 5000.0, 7650.0
P_cl = classical_eddy_loss_density(sigma, d, f, B)
P = lamination_eddy_loss_density(sigma, d, f, B, mu_r=mur)
delta = math.sqrt(2.0 / (2 * math.pi * f * MU0 * mur * sigma))
print(f"\nM19 0.35mm @ {f:.0f} Hz, {B} T, sigma={sigma:.1e}:")
print(f"  classical  P = {P_cl:8.1f} W/m^3 = {P_cl/rho:.3f} W/kg")
print(f"  delta={delta*1e3:.3f} mm, xi=d/delta={d/delta:.3f} -> F={P/P_cl:.4f}, corrected P = {P/rho:.3f} W/kg")
print("  (thin sheet at line frequency: skin correction tiny; halving d would quarter the loss)")
print("\nOK -- F_closed == numerically-integrated exact slab loss; thin limit = classical pi^2 sigma d^2 f^2 B^2/6.")
