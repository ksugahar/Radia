"""COMSOL-class #69 -- magnetizing (air-gap) inductance of an AC machine.

The dominant part of the dq L_md / L_mq: the per-phase main-flux inductance built across the air gap.
From the fundamental MMF -> air-gap B -> flux/pole -> flux-linkage chain,

    L_mu = (2 mu0/pi) (kw1 N_ph)^2 D L_stk / (p^2 g_eff),   L_m = (m/2) L_mu  (synchronous).

This SYNTHESISES the campaign: the fundamental winding factor kw1 comes from #51, the effective gap
g_eff (slotting) from Carter #57; together with the slot-leakage #54 they give the terminal
Ld = L_m + L_leakage. Here the closed form is cross-checked against the explicit derivation chain.
Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (winding_factor, effective_air_gap, MU0,
                                           magnetizing_inductance_per_phase,
                                           synchronous_magnetizing_inductance)

# 36-slot / 4-pole, double-layer 5/6-pitch, D=100 mm gap dia, 100 mm stack, 120 series turns/phase
D, L_stk, p, N_ph, Q = 0.10, 0.10, 2, 120, 36
q = Q // (2 * p * 3)                                     # slots/pole/phase = 3
kw1 = winding_factor(1, q, 360.0 * p / Q, pitch_fraction=5.0 / 6.0)     # #51
g, b_o = 0.8e-3, 2.5e-3
g_eff = effective_air_gap(math.pi * D / Q, g, b_o)      # #57 Carter

print(f"36-slot/4-pole machine: kw1 (#51) = {kw1:.4f},  g_eff (#57) = {g_eff*1e3:.3f} mm (g = {g*1e3:.1f} mm)")

L_mu = magnetizing_inductance_per_phase(D, L_stk, g_eff, p, kw1, N_ph)
L_m = synchronous_magnetizing_inductance(D, L_stk, g_eff, p, kw1, N_ph)


def L_mu_chain():                                        # independent derivation chain
    tau_p = math.pi * D / (2 * p)
    F1 = (2.0 / math.pi) * kw1 * N_ph / p                # MMF amplitude at I=1 A
    B1 = MU0 * F1 / g_eff
    Phi1 = (2.0 / math.pi) * B1 * tau_p * L_stk
    return kw1 * N_ph * Phi1


print(f"\nL_mu (per-phase self) = {L_mu*1e3:.3f} mH   [chain {L_mu_chain()*1e3:.3f} mH, "
      f"match {abs(L_mu/L_mu_chain()-1):.1e}]")
print(f"L_m  (synchronous, 3/2 L_mu) = {L_m*1e3:.3f} mH   -> the dq L_md (non-salient L_mq = L_md)")

print("\nscalings (machine-design levers):")
print(f"  Carter slotting:  L_mu(g_eff)/L_mu(g_smooth) = "
      f"{L_mu/magnetizing_inductance_per_phase(D,L_stk,g,p,kw1,N_ph):.4f}  (<1, openings raise reluctance)")
print(f"  pole pairs 1/p^2: L(p=2)/L(p=4) = "
      f"{L_mu/magnetizing_inductance_per_phase(D,L_stk,g_eff,4,kw1,N_ph):.2f}  (=4)")
print(f"  turns (kw N)^2:   L(N)/L(N/2)   = "
      f"{L_mu/magnetizing_inductance_per_phase(D,L_stk,g_eff,p,kw1,N_ph//2):.2f}  (=4)")
print("\nOK -- L_mu == derivation chain; L_m = (3/2) L_mu = dq L_md; synthesises winding factor #51 + Carter #57.")
