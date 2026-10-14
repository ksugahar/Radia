# -*- coding: utf-8 -*-
r"""demo_dtn_continued_fraction_usage.py -- using radia.open_boundary (exact DtN open boundary)
============================================================================================
The PRODUCTION-API counterpart of the retired research demos, now maintained
through `radia.open_boundary` and `validation_test/open_boundary/test_dtn_continued_fraction.py`.
Here the verified operator ships as `radia.open_boundary`; this shows how to use it.

For spherical MQS truncation this evaluates the exact diffusion DtN symbol
and its terminating continued fraction in q = R0*sqrt(s*mu_sigma).
A finite fraction in sqrt(s) is not an exact finite-dimensional physical-time
ODE. The companion polynomial roots below are q-plane roots. A separate
passive pole fit approximates diffusion memory only over the sampled band.
This demonstration does not establish time-stepping stability or superiority
to a PML on general geometry. See docs/open_boundary/OPEN_BOUNDARY_MAP.md.
Every printed 'ok' is gated on an executed assertion (no overclaim).
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import radia.open_boundary as ob

print("=" * 78)
print(" radia.open_boundary : exact DtN open boundary -- usage")
print("=" * 78)

# A copper-like spherical truncation: R0 = 30 mm, mu*sigma (mu0 * 58 MS/m).
R0 = 0.03
MU_SIGMA = 4.0e-7 * np.pi * 5.8e7
omega = np.logspace(1, 5, 40)      # 10 .. 1e5 rad/s

print(f"\n[1] exact eddy/diffusion DtN per multipole (sphere R0={R0} m):")
for n in (1, 2, 3):
    g50 = ob.eddy_dtn(n, 1j * 5.0e4, R0=R0, mu_sigma=MU_SIGMA)
    print(f"    n={n}: G_n(i*5e4) = {g50:.4f}")

print("\n[2] the continued fraction terminates after n+1 quotients and reproduces the symbol:")
for n in (1, 2, 3, 4, 5, 6):
    stages = ob.continued_fraction_stages(n)
    Zc = np.array([ob.eval_continued_fraction(stages, 1j * w, R0, MU_SIGMA) for w in omega])
    Zr = np.array([ob.eddy_dtn(n, 1j * w, R0, MU_SIGMA) for w in omega])
    nrmse = float(np.sqrt(np.mean(np.abs(Zc - Zr) ** 2)) / np.sqrt(np.mean(np.abs(Zr) ** 2)))
    print(f"    n={n}: quotients={len(stages)} (= n+1)   fraction-vs-symbol NRMSE={nrmse:.1e}")
    assert len(stages) == n + 1 and nrmse < 1e-9

print("\n[3] companion polynomial roots in q (not physical-time ODE poles):")
for n in (1, 2, 3):
    poles = ob.companion_poles(n)
    print(f"    n={n}: {len(poles)} poles, max Re = {poles.real.max():+.3f}  "
          f"=> left-half q-plane roots; no timestep stability claim")
    assert np.all(poles.real < 0)

print("\n[4] sqrt(s) diffusion memory as a passive pole-residue sum (one ODE per pole):")
g, p, e = ob.sqrt_s_passive_poles(omega, 12)
print(f"    K=12: NRMSE={e:.2e}, all poles -p<0 (stable), passive (g>=0): "
      f"{bool(np.all(g >= 0) and np.all(p > 0))}")
assert e < 5e-2 and np.all(p > 0) and np.all(g >= 0)

print("\n[5] Kelvin BUILDS the DtN (material-aware / non-separable companion):")
print("    a radial Kelvin-FEM reproduces the closed-form eddy DtN with NO DC floor")
band = 1j * np.logspace(-4, 2, 30)
for n in (1, 2, 3):
    G = np.array([ob.kelvin_fem_radial_dtn(n, s) for s in band])
    Gx = np.array([ob.eddy_dtn(n, s) for s in band])
    nrmse = float(np.sqrt(np.mean(np.abs(G - Gx) ** 2)) / np.sqrt(np.mean(np.abs(Gx) ** 2)))
    print(f"    n={n}: Kelvin-FEM build vs closed-form, band NRMSE={nrmse:.1e}")
    assert nrmse < 5e-2
print("    (the arbitrary-shape / iron-exterior path = ob.kelvin_dtn_matrix + ob.steklov_spectrum,")
print("     verified O_h-split on a cube -- needs NGSolve; see tests/open_boundary/test_kelvin_dtn.py.")
print("     NOTE: material-in-exterior Kelvin is CLASSICAL (Freeman-Lowther 1988/89).)")

print("\nALL CHECKS PASSED.")
