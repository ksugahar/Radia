"""COMSOL-class #72 -- winding MMF space harmonics & the rotating-field synthesis.

A polyphase winding lays a stepped MMF around the air gap whose n-th space harmonic is
F_phi,n = (4/pi)(k_w,n/n)(N_ph I/2p) -- the (4/pi)/n square-wave envelope reduced by the winding
factor (#51). In a balanced m-phase winding the harmonics combine into ROTATING waves of amplitude
(m/2)F_phi,n: the fundamental rotates FORWARD, the triplen (3rd, 9th) CANCEL, and the surviving
n = 6k+-1 split into FORWARD (7,13,..) and BACKWARD (5,11,..) waves -- the source of the parasitic
asynchronous torques and rotor-surface / magnet-eddy losses. The MMF F_1 here is exactly the one
that drives the magnetizing inductance (#69). Validated by the direct rotating-field trig sum.
Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (winding_factor, single_phase_mmf_harmonic,
                                           rotating_mmf_amplitude, mmf_harmonic_direction)

# 36-slot / 4-pole, double-layer 5/6 pitch, q=3, 120 turns/phase, 10 A peak
q, p, slot_ang, pitch, N_ph, I = 3, 2, 360.0 * 2 / 36, 5.0 / 6.0, 120, 10.0
print("3-phase MMF space harmonics (36-slot/4-pole, 5/6 pitch):")
print("  n   k_w,n    F_phi,n [A]   rotating (m/2)F [A]   sense")
F1 = single_phase_mmf_harmonic(1, winding_factor(1, q, slot_ang, pitch), N_ph, I, p)
for n in (1, 3, 5, 7, 9, 11, 13):
    kw = winding_factor(n, q, slot_ang, pitch)
    Fp = single_phase_mmf_harmonic(n, kw, N_ph, I, p)
    Frot = rotating_mmf_amplitude(n, 3, Fp)
    d = mmf_harmonic_direction(n, 3)
    sense = {1: "forward", -1: "backward", 0: "cancel"}[d]
    print(f"  {n:2d}  {kw:+.4f}   {Fp:8.2f}     {Frot:8.2f}            {sense}")

# rotating-field theorem: the m phases (space & time both 120 deg apart) sum to a rotating wave
def mmf_sum(n, theta, t, m=3):
    return sum(math.cos(t - k * 2 * math.pi / m) * math.cos(n * (theta - k * 2 * math.pi / m))
               for k in range(m))

print("\nrotating-field theorem check (sum of 3 phases == (m/2)cos(n.theta - dir.t)):")
th, t = 0.4, 1.1
for n in (1, 3, 5, 7):
    d = mmf_harmonic_direction(n, 3)
    pred = rotating_mmf_amplitude(n, 3, 1.0) * (math.cos(n * th - d * t) if d else 0.0)
    print(f"  n={n:2d}: sum={mmf_sum(n, th, t):+.5f}  predicted={pred:+.5f}  err={abs(mmf_sum(n,th,t)-pred):.1e}")

print(f"\nfundamental: F_phi,1={F1:.1f} A, rotating (3/2)F_phi,1={rotating_mmf_amplitude(1,3,F1):.1f} A "
      f"(= the MMF behind L_m=(3/2)L_mu, #69)")
print("OK -- F_phi,n=(4/pi)(kw_n/n)(N I/2p); 3-phase: triplen cancel, 6k+-1 survive (fwd 7th/bwd 5th).")
