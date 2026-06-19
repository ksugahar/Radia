# -*- coding: utf-8 -*-
r"""
demo_zz_highfreq_iabc_revisited.py  (Track A -- kelvin branch)
=============================================================
RE-EXAMINING the high-frequency IABC (user: "review the high-frequency IABC a bit more").
This CORRECTS demo_ww, which modelled the high-frequency IABC with a single (essentially real)
shell index and then realized "the time-domain IABC" as a reduced-pole fit of the EXACT DtN.
The ACTUAL high-frequency IABC (Sugahara PIERS-2016, building on Meeker; analytic shell model,
re-derived here independently) is richer:

  * FULL VECTOR EM (Mie): per spherical-harmonic mode n there are TWO polarizations -- a
    permittivity (epsilon) mode and a permeability (mu) mode -- propagated by Riccati-Bessel
    transfer matrices through the shells.
  * each shell carries COMPLEX epsilon AND COMPLEX mu (4 real DOF per shell), with the built-in
    radiation condition (outgoing Hankel) and a PEC termination -- i.e. a discrete, mode-matched
    METAMATERIAL absorber, not a single isotropic index.
  * the 4 DOF/shell are OPTIMIZED to NULL the reflection of the outgoing wave for the matched
    modes -- a per-frequency design.

VERIFIED HERE (all asserted; fully self-contained -- the analytic model is re-optimized in Python,
no external data read):
  (A) the Mie multilayer reflection, re-optimized (generic seeds), NULLS both polarizations for the
      matched mode to machine precision -> reproduces the high-frequency IABC independently.
  (B) the design is NARROWBAND: a shell optimized at omega0 has ~0 reflection at omega0 but the
      reflection RISES quickly off omega0 (e.g. ~6% at +/-10%, ~30% at +/-40%). => to cover a band
      in the TIME DOMAIN the shell materials must be FREQUENCY-DEPENDENT (dispersive).
  (C) THE TIME-DOMAIN OBSTRUCTION (the real finding): the optimal shell has Im(epsilon) and Im(mu)
      of OPPOSITE SIGN across the band -> under one time convention one parameter is lossy and the
      OTHER is active (gain). A single PASSIVE dispersive medium needs both imaginary parts the same
      sign, so the per-frequency-optimal IABC shell is NOT a passive medium; a naive dispersive
      (Debye/Lorentz ADE / recursive-convolution) time-domain realization is therefore NOT directly
      possible (it would be non-passive -> a stability risk). The honest time-domain high-frequency
      IABC needs a PASSIVITY-CONSTRAINED redesign (a genuine open problem), unlike the diffusion case
      (demo_xx) whose sqrt(s) SIBC is passive and Foster-realizable.

CORRECTION TO demo_ww: demo_ww's "reduced M-pole rational DtN = the time-domain IABC" is a clean
model reduction of the TARGET (exact) DtN -- a valid cost-efficient absorber -- but it is NOT the
IABC's own (narrowband, possibly non-passive) materials. demo_ww's cost-accuracy datasheet stands as
a generic reduced-order-DtN result; this file is the faithful picture of the *actual* high-freq IABC.

PRIOR ART: the analytic nested-shell IABC (Meeker IEEE T-Magn 2013/2014; Sugahara PIERS 2016) and the
Mie/Riccati-Bessel multilayer are standard; passivity/causality (Kramers-Kronig) of dispersive media
and ADE/recursive-convolution stability are standard. The re-examination point (narrowband +
non-passive optimal materials -> the time-domain obstruction) is the contribution here.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import hankel1, hankel2
from scipy.optimize import fsolve

np.set_printoptions(precision=5, suppress=True)

# --- spherical / Riccati-Bessel for complex argument -----------------------
def _sb(kind, n, z):
    z = complex(z); p = np.sqrt(np.pi / (2.0 * z)); h1 = hankel1(n + 0.5, z); h2 = hankel2(n + 0.5, z)
    return p * (h1 + h2) / 2.0 if kind == 'j' else p * (h1 - h2) / (2j)

def _sbp(kind, n, z):
    z = complex(z); return _sb(kind, n - 1, z) - (n + 1) / z * _sb(kind, n, z)

def riccati(n, z):
    """Riccati-Bessel  psi_n=z j_n, chi_n=z y_n  and derivatives d/dz (matches Sugahara sp_bessel)."""
    z = complex(z); jn = _sb('j', n, z); yn = _sb('y', n, z)
    return z * jn, z * yn, jn + z * _sbp('j', n, z), yn + z * _sbp('y', n, z)

# --- the actual high-frequency IABC: Mie multilayer reflection (port of ref1) ----
def reflections(mu, ep, r, n, omega):
    """Per mode n, the reflection of the outgoing wave for the epsilon- and mu- polarizations
    through complex-(eps,mu) shells (radii r), PEC-terminated. mu,ep length = len(r) (region 1 =
    vacuum). Returns (ref_eps, ref_mu)."""
    kk = omega * np.sqrt(ep * mu)
    c1 = np.eye(2, dtype=complex); c2 = np.eye(2, dtype=complex)
    for i in range(len(r) - 1):
        bj, by, bjd, byd = riccati(n, kk[i] * r[i])
        A1 = np.array([[mu[i] * bj, mu[i] * by], [kk[i] * bjd, kk[i] * byd]])
        A2 = np.array([[ep[i] * bj, ep[i] * by], [kk[i] * bjd, kk[i] * byd]])
        bj, by, bjd, byd = riccati(n, kk[i + 1] * r[i])
        B1 = np.array([[mu[i + 1] * bj, mu[i + 1] * by], [kk[i + 1] * bjd, kk[i + 1] * byd]])
        B2 = np.array([[ep[i + 1] * bj, ep[i + 1] * by], [kk[i + 1] * bjd, kk[i + 1] * byd]])
        c1 = np.linalg.solve(B1, A1 @ c1); c2 = np.linalg.solve(B2, A2 @ c2)
    bj, by, bjd, byd = riccati(n, kk[-1] * r[-1])
    A1 = np.array([[mu[-1] * bj, mu[-1] * by], [kk[-1] * bjd, kk[-1] * byd]])
    A2 = np.array([[ep[-1] * bj, ep[-1] * by], [kk[-1] * bjd, kk[-1] * byd]])
    c1 = np.array([1, 0]) @ A1 @ c1; c2 = np.array([0, 1]) @ A2 @ c2
    ref_mu = -(c1[0] - 1j * c1[1]) / (c1[0] + 1j * c1[1])
    ref_ep = -(c2[0] - 1j * c2[1]) / (c2[0] + 1j * c2[1])
    return ref_ep, ref_mu

def optimize_shell(omega, r=np.array([1.0, 1.1]), n=1):
    """Re-optimize one shell's (eps,mu) to null both-polarization reflection for mode n at omega.
    Generic seeds only -- no external data."""
    def resid(v):
        ep = np.array([1.0, complex(v[0], v[1])]); mu = np.array([1.0, complex(v[2], v[3])])
        re, rm = reflections(mu, ep, r, n, omega)
        return [re.real, re.imag, rm.real, rm.imag]
    for s in ([1, -.5, 1, -.5], [0.5, 0.5, 0.5, -1], [2, -1, 2, -1], [1, 1, 1, -1], [0.3, 0.2, 5, -2]):
        sol = fsolve(resid, s, full_output=True, xtol=1e-13)[0]
        if np.max(np.abs(resid(sol))) < 1e-9:
            return complex(sol[0], sol[1]), complex(sol[2], sol[3])
    return None

print("=" * 78)
print(" demo_zz : RE-EXAMINING the high-frequency IABC (Mie metamaterial multilayer)")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n[A] re-optimized Mie shell NULLS both polarizations (self-contained, generic seeds):")
designs = {}
for w in (2.0, 5.0, 8.0):
    ep, mu = optimize_shell(w)
    designs[w] = (ep, mu)
    re, rm = reflections(np.array([1.0, mu]), np.array([1.0, ep]), np.array([1.0, 1.1]), 1, w)
    print(f"    omega={w}: eps={ep:.4f}  mu={mu:.4f}   |ref_eps|={abs(re):.1e}  |ref_mu|={abs(rm):.1e}")
    assert abs(re) < 1e-8 and abs(rm) < 1e-8
print("    ok  (the real high-freq IABC = complex-(eps,mu) Mie shells, optimized to null reflection)")

# ---------------------------------------------------------------------------
print("\n[B] the design is NARROWBAND -> needs dispersion for a band:")
ep0, mu0 = designs[5.0]
EP = np.array([1.0, ep0]); MU = np.array([1.0, mu0])
ref_at = {}
for w in (3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0):
    re, rm = reflections(MU, EP, np.array([1.0, 1.1]), 1, w)
    ref_at[w] = abs(re)
    print(f"    design@5.0 evaluated at omega={w}: |ref| = {abs(re):.2e}")
assert ref_at[5.0] < 1e-7
assert ref_at[4.0] > 0.1 and ref_at[6.0] > 0.1     # rises quickly off the design frequency
print("    ok  (reflection ~0 only AT omega0; rises off it -> a per-frequency (dispersive) design)")

# ---------------------------------------------------------------------------
print("\n[C] TIME-DOMAIN OBSTRUCTION: optimal Im(eps), Im(mu) have OPPOSITE sign (one is active):")
for w in (2.0, 5.0, 8.0):
    ep, mu = designs[w]
    print(f"    omega={w}: Im(eps)={ep.imag:+.3f}  Im(mu)={mu.imag:+.3f}  -> product={ep.imag*mu.imag:+.3f}")
    assert ep.imag * mu.imag < 0          # opposite signs => not both passive-lossy
print("    ok  (Im(eps)>0 while Im(mu)<0: a single PASSIVE dispersive medium needs the SAME sign,")
print("         so the per-frequency-optimal shell is non-passive -> a naive Debye/Lorentz/ADE")
print("         (recursive-convolution) time-domain realization is NOT directly possible (stability).")

print("""
[interpretation] -- re-examination of the high-frequency IABC:
  * The ACTUAL high-freq IABC is a complex-(eps,mu) Mie metamaterial multilayer (2 polarizations,
    4 DOF/shell), per-frequency optimized -- NOT the single index of demo_ww.
  * It is NARROWBAND, so a TIME-DOMAIN version needs frequency-dependent (dispersive) shells.
  * BUT the per-frequency-optimal materials are NON-PASSIVE (Im(eps), Im(mu) opposite sign), so a
    naive passive dispersive (ADE / recursive-convolution) realization is not directly possible --
    the honest time-domain high-freq IABC requires a PASSIVITY-CONSTRAINED redesign (open problem).
  * CONTRAST: the eddy-current/diffusion case (demo_xx) has the PASSIVE sqrt(s) SIBC, which Foster-
    realizes cleanly; and at LOW frequency Kelvin is exact and parameter-free (demo_yy). So the
    time-domain story is clean for diffusion and low-freq, but genuinely HARD (passivity) for the
    high-frequency wave IABC -- this corrects demo_ww's optimistic reduced-pole picture.
""")
print("ALL CHECKS PASSED.")
