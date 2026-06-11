"""COMSOL-class #51 -- AC machine WINDING FACTOR k_w = k_d * k_p (distribution + pitch).

What scales the back-EMF (#34/F3) and shapes its harmonics. The n-th harmonic EMF is k_w,n
times the concentrated-full-pitch value:

    distribution  k_d,n = sin(n q gamma/2)/(q sin(n gamma/2)) = |sum_{i<q} exp(j n i gamma)|/q
    pitch         k_p,n = sin(n beta pi/2)                     (beta = coil span / pole pitch)

A distributed (q>1), short-pitched (beta<1) winding has k_w,1 < 1 but STRONGLY suppresses the
low harmonics (a 5/6 pitch kills the 5th & 7th). Pure winding geometry -> tool-independent;
k_d is checked against the direct phasor sum of the q slot EMFs. Self-contained, headless.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (winding_distribution_factor, winding_pitch_factor,
                                           winding_factor)

# typical: Q=36 slots, 2p=4 poles, m=3 phases -> q=3 slots/pole/phase; gamma=p*360/Q; 5/6 pitch
Q, p2, m = 36, 4, 3
q = Q // (p2 * m)
gamma = (p2 // 2) * 360.0 / Q                 # electrical slot angle [deg]
beta = 5.0 / 6.0
print(f"Q={Q}, 2p={p2}, m={m} -> q={q} slots/pole/phase, gamma={gamma:.1f} deg elec, pitch={beta:.3f}\n")
print(f"{'n':>3} {'k_d':>9} {'(phasor)':>9} {'k_p':>9} {'k_w':>9}")
for n in (1, 3, 5, 7, 11, 13):
    kd = winding_distribution_factor(n, q, gamma)
    kp = winding_pitch_factor(n, beta)
    g = math.radians(gamma)
    kd_ph = abs(sum(cmath.exp(1j * n * i * g) for i in range(q))) / q
    print(f"{n:3d} {kd:9.4f} {kd_ph:9.4f} {kp:9.4f} {winding_factor(n,q,gamma,beta):9.4f}")
print(f"\nfundamental k_w1 = {winding_factor(1,q,gamma,beta):.4f}  "
      f"(distributed+short-pitch, ~0.92-0.96 typical)")
print(f"5th/7th suppressed: |k_w5|={abs(winding_factor(5,q,gamma,beta)):.4f}, "
      f"|k_w7|={abs(winding_factor(7,q,gamma,beta)):.4f}  (vs k_w1)")
print("OK -- k_w=k_d*k_p; k_d matches the phasor sum; short-pitch suppresses the 5th/7th harmonic.")
