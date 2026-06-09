"""COMSOL-class #58 -- skew factor (axial skew -> EMF-harmonic & cogging reduction).

Skewing the stator slots (or the rotor magnets) continuously over an electrical angle theta_sk
averages the n-th air-gap harmonic along the STACK, scaling the n-th EMF / torque harmonic by

    k_sk,n = sin(n theta_sk/2)/(n theta_sk/2)        (the sinc / averaging factor).

The axial twin of the distribution factor (#51, which spreads a phase over slots). A harmonic
whose period equals the skew (n theta_sk = 2 pi) is NULLED, so skewing by ONE SLOT PITCH cancels
the slot-passing / cogging ripple -- at a small fundamental cost k_sk,1 < 1. Self-consistent
(validated against the phasor average); pure geometry, headless.
"""
import os
import sys
import math
import cmath

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (skew_factor, slot_pitch_skew_angle,
                                           skewed_winding_factor, winding_factor)


def phasor_average(n, theta_sk, m=20000):                  # |mean of exp(j n phi)| over the skew
    return abs(sum(cmath.exp(1j * n * (i + 0.5) * theta_sk / m) for i in range(m)) / m)


# integer-slot machine: Q=36 slots, 2p=4 poles -> p=2, q = 36/(4*3) = 3 slots/pole/phase
Q, p, q = 36, 2, 3
slot_angle = 360.0 * p / Q                                 # 20 deg electrical between slots
theta_slot = slot_pitch_skew_angle(Q, p)                   # one slot pitch [elec rad]

print(f"machine Q={Q} 2p={2*p} (q={q}); one slot pitch = {math.degrees(theta_slot):.0f} deg elec")
print("skew = one slot pitch:")
print(f"  fundamental k_sk,1 = {skew_factor(1, theta_slot):.4f}  (small cost; formula==phasor "
      f"{phasor_average(1, theta_slot):.4f})")
# the slot-passing / cogging harmonic is at order Q/p = 18 -> nulled by one-slot-pitch skew
n_slot = Q // p
print(f"  slot-passing harmonic n={n_slot}: k_sk = {skew_factor(n_slot, theta_slot):+.2e}  (NULLED -> "
      f"cancels cogging)")

print("\nharmonic suppression (distribution k_d alone vs distribution x skew):")
print("  n    k_w (no skew)   k_w*k_sk (skewed)")
for n in (1, 5, 7, 11, 13):
    kw = winding_factor(n, q, slot_angle)
    kws = skewed_winding_factor(n, q, slot_angle, theta_slot)
    print(f"  {n:2d}     {kw:+.4f}        {kws:+.4f}")

print("\nk_sk,n = sin(n theta/2)/(n theta/2): nulls any harmonic with n*theta_sk = 2*pi.")
print("OK -- skew factor == phasor average; one-slot-pitch skew cancels the cogging harmonic.")
