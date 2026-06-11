"""COMSOL-class #49 -- PM machine 3-phase SHORT-CIRCUIT current + braking torque + demag.

The fault / protection regime (vs the operating-point #46): with the terminals shorted
(v_d=v_q=0) the dq voltage equations give the steady short-circuit currents

    id = -we^2 Lq lm /(R^2 + we^2 Ld Lq),   iq = -we R lm /(R^2 + we^2 Ld Lq).

HIGH-SPEED limit: id -> -lm/Ld = -Ich (the CHARACTERISTIC current, pure d-axis DEMAGNETISING),
iq -> 0, so |Isc| -> Ich. The braking torque (dq_torque of these currents) rises, peaks near
the critical speed, then falls to 0. Two design concerns: the demag current Ich (vs the magnet
knee) and the peak braking torque. Self-contained, headless, tool-independent (pure dq theory).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (short_circuit_dq_currents, characteristic_current,
                                           dq_torque)

lm, Ld, Lq, R, p = 0.1, 0.5e-3, 1.5e-3, 0.05, 4
Ich = characteristic_current(lm, Ld)
print(f"characteristic current Ich = lm/Ld = {Ich:.1f} A  (high-speed short-circuit / demag current)\n")
print(f"{'f_e[Hz]':>8} {'id[A]':>9} {'iq[A]':>8} {'|Isc|[A]':>9} {'T_brake[N.m]':>13}")
peakT, peakf = 0.0, 0.0
for fe in (10, 25, 50, 100, 200, 500, 2000):
    we = 2 * math.pi * fe
    id_, iq = short_circuit_dq_currents(R, Ld, Lq, lm, we)
    T = dq_torque(lm, Ld, Lq, id_, iq, p)
    if abs(T) > peakT:
        peakT, peakf = abs(T), fe
    print(f"{fe:8d} {id_:9.2f} {iq:8.2f} {math.hypot(id_,iq):9.2f} {T:13.3f}")

idh, iqh = short_circuit_dq_currents(R, Ld, Lq, lm, 2 * math.pi * 1e5)
print(f"\nhigh-speed: id -> {idh:.3f} A (= -Ich = {-Ich:.3f}), iq -> {iqh:.4f} (->0); "
      f"|Isc| -> Ich within {abs(math.hypot(idh,iqh)-Ich)/Ich*100:.3f}%")
print(f"braking torque peaks near {peakf} Hz (|T|~{peakT:.1f} N.m), -> 0 at high speed.")
print("OK -- steady short-circuit -> Ich (demag check); braking torque has a peak then decays.")
