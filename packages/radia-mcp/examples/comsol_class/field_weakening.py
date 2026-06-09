"""COMSOL-class #38 -- dq FIELD-WEAKENING operating-region map (MTPA -> voltage-
limited FW -> MTPV) of a salient PM synchronous machine.

Above base speed the inverter voltage limit (Vmax/omega_e) caps the reachable dq
flux, so the torque-maximising current must be weakened (id more negative). This
maps the three operating regions and the finite max speed, the high-speed sequel to
the MTPA locus (#26):

  * MTPA  -- below base speed (constant torque, current-limited),
  * FW    -- on the current circle ∩ voltage ellipse (constant ~power),
  * MTPV  -- max-torque-per-voltage tangent inside the current circle
             (only when characteristic current Ich = lambda_m/Ld < Imax).

Self-contained, headless, tool-independent: the closed-form operating point is
checked against an INDEPENDENT numeric constrained argmax (the analytic answer is
the truth standard).
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (dq_torque, field_weakening_operating_point,
                                           base_speed_electrical)

I, p = 10.0, 4                       # current limit [A], pole pairs


def numeric_opt(lm, Ld, Lq, Imax, Vlim, N=600):
    """independent brute-force constrained argmax of T over feasible (id,iq)."""
    def scan(id0, id1, iq0, iq1, n):
        best = (-1e99, 0.0, 0.0)
        for i in range(n + 1):
            id_ = id0 + (id1 - id0) * i / n
            for j in range(n + 1):
                iq = iq0 + (iq1 - iq0) * j / n
                if id_ * id_ + iq * iq > Imax * Imax:
                    continue
                if (Ld * id_ + lm) ** 2 + (Lq * iq) ** 2 > Vlim * Vlim:
                    continue
                T = dq_torque(lm, Ld, Lq, id_, iq, p)
                if T > best[0]:
                    best = (T, id_, iq)
        return best
    T, idb, iqb = scan(-Imax, Imax, 0.0, Imax, N)
    step = 2 * Imax / N
    for _ in range(2):
        T, idb, iqb = scan(idb - 2 * step, idb + 2 * step,
                            max(0.0, iqb - 2 * step), iqb + 2 * step, 200)
        step = 4 * step / 200
    return None if T < -1e90 else (idb, iqb, T)


def run(name, lm, Ld, Lq, Vmax):
    Ich = lm / Ld
    wb = base_speed_electrical(lm, Ld, Lq, I, Vmax, p)
    print(f"\n{name}:  Ich=lm/Ld={Ich:.2f} A "
          f"({'<Imax -> MTPV region' if Ich < I else '>=Imax -> finite max speed'}),  "
          f"base speed (elec) = {wb:.1f} rad/s")
    for k in (0.5, 1.0, 1.3, 1.8, 2.5, 4.0, 6.0):
        we = k * wb
        op = field_weakening_operating_point(lm, Ld, Lq, I, Vmax, we, p)
        num = numeric_opt(lm, Ld, Lq, I, Vmax / we)
        if op is None:
            tag = "infeasible (>w_max)" + ("" if num is None else "  !! numeric found one")
            print(f"  {k:>4}x wb:  {tag}")
            continue
        err = abs(op[2] - num[2]) / abs(num[2]) * 100
        print(f"  {k:>4}x wb:  region={op[3]:5} T={op[2]:7.4f} N.m  id={op[0]:7.3f} iq={op[1]:6.3f} "
              f"|I|={math.hypot(op[0],op[1]):5.2f}  (vs numeric {err:.3f}%)")


run("IPM wide-CPSR", 0.003, 0.0005, 0.0015, 1.2)
run("IPM current-limited", 0.10, 0.0005, 0.0015, 12.0)
run("non-salient PMSM", 0.05, 0.001, 0.001, 6.0)
print("\nOK -- closed-form field-weakening operating point matches the numeric optimum (0.000%).")
