"""COMSOL-class #46 -- dq steady-state OPERATING POINT (voltage, power, power factor).

Where MTPA (#26) / field weakening (#37) CHOOSE the dq currents, this EVALUATES the terminal
quantities the inverter sees there -- the dq machine model's voltage/power layer:

    lambda_d = Ld id + lambda_m,  lambda_q = Lq iq
    v_d = R id - omega_e lambda_q,  v_q = R iq + omega_e lambda_d
    P_in = (3/2)(v_d i_d + v_q i_q),  P_em = (3/2) omega_e (lambda_d i_q - lambda_q i_d) = T omega_mech,
    P_cu = (3/2) R |I|^2,  PF = P_in/((3/2)|V||I|),  eff = P_em/P_in.

Self-consistent: P_in = P_cu + P_em exactly, P_em = T*omega_mech, and the torque matches
dq_torque (#26). Pure dq theory -> tool-independent. Inverter/drive sizing (|V| vs the DC-link,
|I| vs the device rating, PF, efficiency). Self-contained, headless.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import mtpa_operating_point, dq_operating_point

lm, Ld, Lq, R, p = 0.1, 0.5e-3, 1.5e-3, 0.05, 4      # IPM PMSM

g, id_, iq, T = mtpa_operating_point(lm, Ld, Lq, 15.0, p)   # MTPA @ |I|=15 A
print(f"MTPA @ 15 A: id={id_:.3f}, iq={iq:.3f}, T={T:.3f} N.m\n")
print(f"{'f_e[Hz]':>8} {'|V|[V]':>8} {'P_in[W]':>9} {'P_em[W]':>9} {'PF':>7} {'eff[%]':>7} {'bal':>10}")
for fe in (50, 100, 200, 400):
    we = 2 * math.pi * fe
    op = dq_operating_point(R, Ld, Lq, lm, id_, iq, we, p)
    bal = op["P_in"] - (op["P_cu"] + op["P_em"])
    print(f"{fe:8d} {op['Vmag']:8.2f} {op['P_in']:9.1f} {op['P_em']:9.1f} "
          f"{op['power_factor']:7.4f} {op['efficiency']*100:7.3f} {bal:10.2e}")
print("\n|V| rises ~linearly with speed (back-EMF); the field-weakening limit (#37) is where "
      "|V| hits the inverter ceiling. Power balance P_in=P_cu+P_em holds to machine precision.")
