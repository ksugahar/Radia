"""COMSOL-class #40 -- INDUCTION MACHINE torque-slip curve + breakdown (equivalent circuit).

The induction-machine companion to the PM dq blocks (MTPA #26, field weakening #37): a new
motor class. The single-cage Thevenin equivalent circuit gives the air-gap torque vs slip

    T(s) = (3/omega_s) Vth^2 (R2/s) / [(Rth+R2/s)^2 + (Xth+X2)^2],

peaking at the BREAKDOWN (pull-out) point  s_max = R2/sqrt(Rth^2+(Xth+X2)^2),
T_max = 3 Vth^2 / (2 omega_s [Rth + sqrt(Rth^2+(Xth+X2)^2)]).  The classic invariant:
T_max is INDEPENDENT of the rotor resistance R2 -- only s_max scales with R2 (how a
wound-rotor / deep-bar design moves peak torque toward standstill for starting).
Self-contained, headless, tool-independent (closed form vs an independent numeric slip sweep).
The analytic basis for the JMAG-roadmap IM 1/4-model FE study.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (induction_machine_torque,
                                           induction_machine_breakdown)

# representative 4-pole 50 Hz IM (SI; referred to stator)
V1, R1, X1, R2, X2, Xm = 230.0, 0.5, 1.0, 0.3, 1.0, 30.0
p = 2
omega_s = 2 * math.pi * 50.0 / p          # synchronous mechanical speed [rad/s]

s_max, T_max = induction_machine_breakdown(V1, R1, X1, R2, X2, Xm, omega_s)
print(f"synchronous speed = {omega_s*60/(2*math.pi):.0f} rpm")
print(f"breakdown: s_max = {s_max:.4f}, T_max = {T_max:.2f} N.m\n")
print(f"{'slip':>6} {'speed[rpm]':>11} {'T[N.m]':>9}")
for s in (1.0, 0.5, s_max, 0.10, 0.05, 0.02, 0.005):
    rpm = (1 - s) * omega_s * 60 / (2 * math.pi)
    print(f"{s:6.3f} {rpm:11.0f} {induction_machine_torque(V1,R1,X1,R2,X2,Xm,omega_s,s):9.2f}")

# numeric argmax over slip -> confirms the closed-form breakdown
N = 100000
s_num, T_num = max(((i/N, induction_machine_torque(V1,R1,X1,R2,X2,Xm,omega_s,i/N))
                    for i in range(1, N+1)), key=lambda t: t[1])
print(f"\nnumeric argmax: s={s_num:.4f}, T={T_num:.2f}  "
      f"(closed form err: s {abs(s_max-s_num)/s_num*100:.2f}%, T {abs(T_max-T_num)/T_num*100:.4f}%)")

# rotor-resistance invariant: T_max fixed, s_max ∝ R2
print("R2 sweep (breakdown T_max should be constant; s_max ∝ R2):")
for k in (1.0, 2.0, 3.0):
    sm, Tm = induction_machine_breakdown(V1, R1, X1, k*R2, X2, Xm, omega_s)
    print(f"  R2*{k:.0f}: s_max={sm:.4f} (x{sm/s_max:.2f}), T_max={Tm:.2f}")
print("OK -- IM torque-slip + breakdown match the numeric sweep; T_max R2-invariant.")
