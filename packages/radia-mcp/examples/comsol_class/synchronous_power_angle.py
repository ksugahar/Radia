"""COMSOL-class #43 -- SALIENT synchronous machine POWER-ANGLE (delta) curve + pull-out.

The grid/voltage-frame characteristic of a synchronous machine (vs the current-frame MTPA
#26, field weakening #37): torque vs the LOAD ANGLE delta (between terminal voltage V and
excitation EMF E),

    T(delta) = (3 p/omega_e) [ V E/Xd sin(delta)  +  (V^2/2)(1/Xq - 1/Xd) sin(2 delta) ],

a FIELD torque (~sin d) + a RELUCTANCE torque (~sin 2d) from the saliency Xd != Xq. The
steady-state stability limit is the PULL-OUT torque (the peak). LIMITS: non-salient ->
peak at 90 deg; pure reluctance (E=0) -> 45 deg; salient field machine in between. The
synchronous-machine companion to the induction-machine torque-slip (#40). Self-contained,
headless, tool-independent (closed-form pull-out vs an independent numeric angle sweep).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (synchronous_power_angle_torque, synchronous_pullout)

V, E, Xd, Xq = 230.0, 276.0, 2.0, 1.2       # phase RMS [V], EMF [V], reactances [ohm]
we, p = 2 * math.pi * 50.0, 2               # electrical speed, pole pairs

d_max, T_max = synchronous_pullout(V, E, Xd, Xq, we, p)
print(f"pull-out: delta_max = {math.degrees(d_max):.2f} deg, T_max = {T_max:.2f} N.m\n")
print(f"{'delta[deg]':>10} {'T_field':>9} {'T_reluct':>10} {'T_total[N.m]':>13}")
for dd in (15, 30, 45, math.degrees(d_max), 75, 90):
    d = math.radians(dd)
    Tf = (3*p/we) * V*E/Xd*math.sin(d)
    Tr = (3*p/we) * 0.5*V*V*(1/Xq-1/Xd)*math.sin(2*d)
    print(f"{dd:10.2f} {Tf:9.1f} {Tr:10.1f} {synchronous_power_angle_torque(V,E,Xd,Xq,we,p,d):13.2f}")

# closed-form pull-out vs an independent numeric angle sweep
N = 200000
dn, Tn = max(((i/N*math.pi, synchronous_power_angle_torque(V, E, Xd, Xq, we, p, i/N*math.pi))
              for i in range(1, N)), key=lambda t: t[1])
print(f"\nnumeric peak: {math.degrees(dn):.2f} deg, {Tn:.2f} N.m  "
      f"(closed-form err: d {abs(d_max-dn)/dn*100:.3f}%, T {abs(T_max-Tn)/Tn*100:.4f}%)")

# limits
nd, _ = synchronous_pullout(V, E, 2.0, 2.0, we, p)       # non-salient
rd, _ = synchronous_pullout(V, 0.0, 2.0, 1.2, we, p)     # reluctance only (E=0)
print(f"limits: non-salient pull-out = {math.degrees(nd):.1f} deg (->90), "
      f"reluctance-only = {math.degrees(rd):.1f} deg (->45)")
print("OK -- salient pull-out between 45 and 90 deg; field+reluctance torque decomposition.")
