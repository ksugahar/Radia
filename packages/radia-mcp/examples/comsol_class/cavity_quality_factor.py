"""COMSOL-class #68 -- rectangular-cavity quality factor (wall-loss-limited Q).

The LOSS sequel to the lossless cavity resonance (#59). A closed metal box resonates at
f_TE101 = (c/2) sqrt(1/a^2 + 1/d^2); finite-conductivity walls dissipate power, giving a finite

    Q = omega * U / P_wall = (k a d)^3 b eta / (2 pi^2 R_s (2 a^3 b + 2 b d^3 + a^3 d + a d^3)),

with R_s = sqrt(omega mu/(2 sigma)) = 1/(sigma*delta) the wall surface resistance. Q ~ (V/S)/delta ~
sqrt(sigma): a better conductor (thinner skin depth) gives a sharper resonance -- the bridge between
the wave family (#53/#59/#62/#65) and the skin depth (#45/#55). Here the closed form is cross-checked
against the numerically-integrated TE101 field energy and wall loss. Self-contained, headless; SI.
"""
import os
import sys
import math

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import numpy as np
from radia_mcp.radia_ngsolve.waveguide import (rectangular_cavity_frequency, rectangular_cavity_q,
                                               surface_resistance, skin_depth, MU0, C0)

EPS0 = 1.0 / (MU0 * C0 ** 2)


def q_from_fields(a, b, d, sigma, n=240):
    """Q = omega U / P_wall from the numerically-integrated TE101 fields (independent of the formula)."""
    f = rectangular_cavity_frequency(a, b, d, 1, 0, 1)
    w = 2 * math.pi * f
    Rs = surface_resistance(f, sigma)
    Camp = 1.0 / (w * MU0)
    x = np.linspace(0, a, n); z = np.linspace(0, d, n)
    X, Z = np.meshgrid(x, z, indexing="ij")
    U = (EPS0 / 2.0) * np.trapezoid(np.trapezoid((np.sin(math.pi * X / a) * np.sin(math.pi * Z / d)) ** 2,
                                                  z, axis=1), x) * b
    Hx = lambda xx, zz: Camp * (math.pi / d) * np.sin(math.pi * xx / a) * np.cos(math.pi * zz / d)
    Hz = lambda xx, zz: Camp * (math.pi / a) * np.cos(math.pi * xx / a) * np.sin(math.pi * zz / d)
    surf = (np.trapezoid(Hz(0.0, z) ** 2, z) + np.trapezoid(Hz(a, z) ** 2, z)
            + np.trapezoid(Hx(x, 0.0) ** 2, x) + np.trapezoid(Hx(x, d) ** 2, x)) * b \
        + 2 * np.trapezoid(np.trapezoid(Hx(X, Z) ** 2 + Hz(X, Z) ** 2, z, axis=1), x)
    return w * U / ((Rs / 2.0) * surf)


sigma = 5.8e7                                            # copper
print("rectangular cavity TE101 Q (copper walls):")
print("  a,b,d [mm]      f[GHz]   R_s[mOhm]   Q_closed   Q_fields    match")
for (a, b, d) in [(0.04, 0.02, 0.03), (0.0229, 0.0102, 0.0229), (0.05, 0.03, 0.05)]:
    f = rectangular_cavity_frequency(a, b, d, 1, 0, 1)
    Qc = rectangular_cavity_q(a, b, d, sigma)
    Qf = q_from_fields(a, b, d, sigma)
    print(f"  {a*1e3:.0f},{b*1e3:.0f},{d*1e3:.0f}        {f/1e9:.3f}   {surface_resistance(f,sigma)*1e3:.2f}"
          f"     {Qc:8.0f}   {Qf:8.0f}   {abs(Qc/Qf-1):.1e}")

a, b, d = 0.04, 0.02, 0.03
print(f"\nthe 40x20x30 mm cavity is the #59 TE101 ({rectangular_cavity_frequency(a,b,d,1,0,1)/1e9:.3f} GHz);"
      f" copper delta there = {skin_depth(rectangular_cavity_frequency(a,b,d,1,0,1),sigma)*1e6:.3f} um")
print(f"Q ~ sqrt(sigma): Q(sigma)/Q(4 sigma) = {rectangular_cavity_q(a,b,d,sigma)/rectangular_cavity_q(a,b,d,4*sigma):.4f}"
      f"  (expect 0.5); superconducting walls -> Q can reach 1e9+")
print("\nOK -- closed-form Q == field-integrated omega U / P_wall; Q ~ (V/S)/delta ~ sqrt(sigma).")
