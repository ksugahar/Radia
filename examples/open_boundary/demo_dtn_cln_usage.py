# -*- coding: utf-8 -*-
r"""demo_dtn_cln_usage.py -- using radia.open_boundary (the Zs-DtN-CLN open boundary)
================================================================================
The PRODUCTION-API counterpart of the research demos
examples/kelvin_transformation/DtN_spectrum/act6_02_cln_dtn_cauer.py (inline math).
Here the verified operator ships as `radia.open_boundary`; this shows how to use it.

WHAT it gives you, for a SEPARABLE (spherical) MQS truncation:
  * the EXACT exterior eddy/diffusion DtN eigenvalue per multipole n,
  * its Cauer Ladder Network (CLN) realisation -- EXACT at n+1 stages, well-conditioned,
  * the companion auxiliary-ODE poles (Re<0) for a transient passive Robin boundary.

SCOPE (honest): this WINS over a CFS-PML only on its island -- compact / quasi-spherical
MQS problems wanting an exact, DC-well-conditioned boundary + a compact passive circuit
ROM.  Elongated/arbitrary geometry -> a box CFS-PML hugs better; genuine wave radiation
is outside radia's MQS scope.  NOT novel (Grote-Keller / Hagstrom-Warburton / Warburg-Cauer
/ Kameari CLN).  See docs/open_boundary/OPEN_BOUNDARY_MAP.md.

Every printed 'ok' is gated on an executed assertion (no overclaim).
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import radia.open_boundary as ob

print("=" * 78)
print(" radia.open_boundary : exact Zs-DtN-CLN open boundary -- usage")
print("=" * 78)

# A copper-like spherical truncation: R0 = 30 mm, mu*sigma (mu0 * 58 MS/m).
R0 = 0.03
MU_SIGMA = 4.0e-7 * np.pi * 5.8e7
omega = np.logspace(1, 5, 40)      # 10 .. 1e5 rad/s

print(f"\n[1] exact eddy/diffusion DtN per multipole (sphere R0={R0} m):")
for n in (1, 2, 3):
    g50 = ob.eddy_dtn(n, 1j * 5.0e4, R0=R0, mu_sigma=MU_SIGMA)
    print(f"    n={n}: G_n(i*5e4) = {g50:.4f}")

print("\n[2] Cauer ladder is EXACT at n+1 stages and reproduces the symbol:")
for n in (1, 2, 3, 4, 5, 6):
    stages = ob.cauer_ladder(n)
    Zc = np.array([ob.eval_ladder(stages, 1j * w, R0, MU_SIGMA) for w in omega])
    Zr = np.array([ob.eddy_dtn(n, 1j * w, R0, MU_SIGMA) for w in omega])
    nrmse = float(np.sqrt(np.mean(np.abs(Zc - Zr) ** 2)) / np.sqrt(np.mean(np.abs(Zr) ** 2)))
    print(f"    n={n}: stages={len(stages)} (= n+1)   ladder-vs-symbol NRMSE={nrmse:.1e}")
    assert len(stages) == n + 1 and nrmse < 1e-9

print("\n[3] transient Robin realisation: one auxiliary ODE per companion pole (Re<0):")
for n in (1, 2, 3):
    poles = ob.companion_poles(n)
    print(f"    n={n}: {len(poles)} poles, max Re = {poles.real.max():+.3f}  "
          f"=> passive / unconditionally stable (Grote-Keller form)")
    assert np.all(poles.real < 0)

print("\n[4] sqrt(s) diffusion-memory element as a finite passive ladder:")
g, p, e = ob.sqrt_s_passive_ladder(omega, 12)
print(f"    K=12: NRMSE={e:.2e}, all poles -p<0 (stable), passive (g>=0): "
      f"{bool(np.all(g >= 0) and np.all(p > 0))}")
assert e < 5e-2 and np.all(p > 0) and np.all(g >= 0)

print("\nALL CHECKS PASSED.")
