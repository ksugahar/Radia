# -*- coding: utf-8 -*-
# DEMO (gg) (verified): the exterior-DtN spectrum in the COMPLEX PLANE -- Kelvin vs PML vs FEM-BEM
# on one comparison axis, and what "high frequency" does to it.
#
# The user's point: "comparing the high-frequency Kelvin spectrum, the PML spectrum and the FEM-BEM
# spectrum might be interesting" -- and "in that case the spectrum becomes a COMPLEX-PLANE spectrum,
# right?"  YES. The exterior Dirichlet-to-Neumann (Steklov-Poincare) operator on the truncation
# sphere r=R has, per spherical-harmonic degree n, the EXACT eigenvalue (Helmholtz, wavenumber k,
# outgoing/Sommerfeld radiation condition):
#
#       Lambda_n(kR) = kR * h_n^(1)'(kR) / h_n^(1)(kR)        (units 1/R; R=1 so the argument IS kR)
#
# h_n^(1) = spherical Hankel (outgoing).  This is COMPLEX at finite kR: |Im(Lambda_n)| = RADIATION
# (its sign is the e^{-iwt}/h^(1) convention; here Im>0).  Static limit kR->0 collapses it onto the
# REAL negative ladder Lambda_n -> -(n+1)/R (the Track-A datasheet).  A single Argand plot of
# {Lambda_n} therefore puts every open-boundary method on one axis:
#
#   * FEM-BEM      : reproduces the EXACT complex Hankel locus (BEM = the exterior DtN, up to surface
#                    discretisation) -> the reference / gold standard, complex, but DENSE.
#   * Kelvin (static core, the SA/magnetostatics tool) : has NO frequency in it -> its spectrum is
#                    PINNED to the REAL axis at -(n+1), frequency-INDEPENDENT.  It is the kR->0
#                    operator.  Its deviation from the exact complex value is the QUASI-STATIC ERROR:
#                    the missing radiation Im PLUS a real-part shift ~ kR^2/(2n).  It -> 0 as kR->0
#                    (quasi-static exactness) and as n->inf at fixed kR (rate kR^2/2n, verified).
#                    => Kelvin's PERFORMANCE DOMAIN is the quasi-static spectrum; it is the CHEAPEST
#                       (closed-form / sparse thin ball, real, zero wavelength to resolve) but is
#                       accurate only for kR<~1 (the static-apparatus / SA regime).
#   * PML          : a complex-stretched layer; it DOES carry the wave, so it reproduces the exact
#                    complex spectrum -- a thick/strong PML matches it at all n (verified) -- with a
#                    characteristic error KNEE at the propagating<->evanescent transition n~kR
#                    (verified): a thin PML is most stressed exactly where n~kR.  Cost = a resolved
#                    layer (DoF grows with kR to resolve the wavelength).
#
# Honest conclusion (no false "complementary crossover"): at ANY finite kR a layer/integral method
# that carries k^2 (PML, FEM-BEM, or a Helmholtz/extended-Kelvin) is more accurate than the static
# Kelvin operator -- static Kelvin is strictly the kR->0 tool.  So the kR AXIS is the two-paper
# boundary (act5_08_sommerfeld_frequency_sweep): kR<~1 quasi-static = the SA Kelvin paper (Kelvin exact AND cheapest); kR>~1
# radiating = the exact complex operator (BEM / PML / extended-Kelvin sugahara2025).  The genuine
# "high-frequency Kelvin" is the Helmholtz/extended-Kelvin inversion (sugahara2025), whose spectrum
# WOULD track the exact complex locus until the inverted far-field oscillation out-resolves the FE
# (a peel-off) -- that is the natural follow-up; this demo verifies the three clean closed-form
# objects (exact, static-Kelvin, PML) that frame it.
#
# VERIFICATIONS (all asserted from computed values; no overclaim):
#  (1) my complex-argument spherical-Bessel helpers == scipy.special.spherical_jn/yn on real args.
#  (2) the EXACT DtN static limit Lambda_n(kR->0) -> -(n+1)  (real ladder).
#  (3) Lambda_n is genuinely COMPLEX at finite kR; |Im| large for n<kR, ->0 for evanescent n>kR.
#  (4) the PML radial solver, with a strong thick layer, reproduces the EXACT complex Lambda_n for
#      propagating modes (correctness) -- and its r~ integration returns to R (a self-check).
#  (5a) Kelvin is the quasi-static operator: |Kelvin-exact| -> 0 as kR->0 (~kR^2 at fixed n), and the
#       evanescent error law |Kelvin-exact|*2n/kR^2 -> 1 for n>>kR.
#  (5b) the PML error KNEE: a thin PML's per-degree error PEAKS at n~kR (the propagating/evanescent
#       transition), confirming the textbook PML stress point.
#
# Pure numpy/scipy (no NGSolve): this is the SPECTRUM (the per-mode operator eigenvalue), closed-form
# on the sphere for all three -- the cleanest apples-to-apples comparison.  The FEM realisation of
# these operators is act1_05_assemble_dtn_matrix/w (Kelvin) and act4_02_kelvin_approximates_bem/r (BEM reference).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import jv, yv, hankel1, spherical_jn, spherical_yn

R = 1.0  # work in units of the truncation radius; the argument of the special functions IS kR


# ----------------------------------------------------------------------------- spherical Bessel
def _sph(kind, n, z):
    """Spherical Bessel f_n(z) for COMPLEX z via half-integer ordinary Bessel.
    kind in {'j','y','h1'}.  f_n(z) = sqrt(pi/2z) * F_{n+1/2}(z)."""
    z = np.asarray(z, dtype=complex)
    pref = np.sqrt(np.pi / (2.0 * z))
    if kind == 'j':
        return pref * jv(n + 0.5, z)
    if kind == 'y':
        return pref * yv(n + 0.5, z)
    if kind == 'h1':
        return pref * hankel1(n + 0.5, z)
    raise ValueError(kind)


def _sph_prime(kind, n, z):
    """f_n'(z) via the recurrence f_n' = f_{n-1} - (n+1)/z f_n (holds for j, y, h1)."""
    z = np.asarray(z, dtype=complex)
    return _sph(kind, n - 1, z) - (n + 1.0) / z * _sph(kind, n, z)


def dtn_exact(n, kR):
    """EXACT exterior-Helmholtz DtN eigenvalue on degree n: kR h_n^(1)'(kR)/h_n^(1)(kR) (units 1/R).
    The operator a (converged) FEM-BEM reproduces; the gold-standard complex locus."""
    z = complex(kR)
    return (z * _sph_prime('h1', n, z) / _sph('h1', n, z)).item()


def dtn_kelvin_static(n):
    """The Kelvin open boundary as used in magnetostatics IS the STATIC operator: the real ladder
    -(n+1)/R, frequency-independent, no radiation (Im == 0).  Exact only as kR->0."""
    return -(n + 1.0) / R


# ----------------------------------------------------------------------------- PML radial model
def dtn_pml(n, kR, d=1.5, sigma0=16.0, nstep=6000):
    """Finite complex-stretched radial PML DtN on degree n, terminated by a Dirichlet wall.
    PML region r in [R, R+d], smooth absorption sigma(r)=sigma0*((r-R)/d)^2 (vanishes at r=R so NO
    interface jump), complex stretch s(r)=1+i sigma(r)/k, k=kR/R.  Integrate the mode-n radial
    Helmholtz inward from the wall (u=0) to r=R; Lambda_PML=(du/dr)/u|_R (physical derivative,
    s(R)=1).  Returns (Lambda_PML, |rtilde(R)-R|) -- the 2nd value is a self-check (->0)."""
    k = kR / R
    r0, r1 = R, R + d
    h = (r1 - r0) / nstep

    def s(r):
        return 1.0 + 1j * (sigma0 * ((r - R) / d) ** 2) / k

    def deriv(r, st):
        u, q, rt = st                          # q = du/dr~ ; rt = r~ (stretched radius)
        sr = s(r)
        du = sr * q
        dq = sr * (-(2.0 / rt) * q - (k * k - n * (n + 1.0) / (rt * rt)) * u)
        return np.array([du, dq, sr], dtype=complex)

    rg = np.linspace(r0, r1, nstep + 1)
    sv = s(rg)
    rt_outer = R + np.sum(0.5 * (sv[1:] + sv[:-1]) * np.diff(rg))   # r~ at the wall (trapezoid)

    st = np.array([0.0 + 0j, 1.0 + 0j, rt_outer], dtype=complex)    # wall: u=0, q=1
    r = r1
    for _ in range(nstep):                                          # RK4 inward to r=R
        k1 = deriv(r, st)
        k2 = deriv(r - 0.5 * h, st - 0.5 * h * k1)
        k3 = deriv(r - 0.5 * h, st - 0.5 * h * k2)
        k4 = deriv(r - h, st - h * k3)
        st = st - (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        r -= h
    u, q, rt = st
    return (s(R) * q / u).item(), abs(rt - R)                       # s(R)=1


# ============================================================================= VERIFICATIONS
print("=" * 78)
print("exterior-DtN spectrum in the COMPLEX PLANE: Kelvin vs PML vs FEM-BEM")
print("  Lambda_n(kR) = kR h_n^(1)'(kR)/h_n^(1)(kR)  (R=1; complex = radiation)")
print("=" * 78)

# (1) spherical-Bessel helper vs scipy on real args -----------------------------------------
err_j = err_y = 0.0
for n in range(0, 8):
    for z in (0.3, 1.0, 2.5, 7.0):
        rj, ry = spherical_jn(n, z), spherical_yn(n, z)
        err_j = max(err_j, abs(_sph('j', n, z).item() - rj) / (abs(rj) + 1e-300))
        err_y = max(err_y, abs(_sph('y', n, z).item() - ry) / (abs(ry) + 1e-300))
print("\n(1) complex sph-Bessel helper vs scipy.special.spherical_jn/yn (real args, RELATIVE):")
print("    max rel|j_n - scipy| = %.2e ,  max rel|y_n - scipy| = %.2e" % (err_j, err_y))
assert err_j < 1e-11 and err_y < 1e-11, "spherical-Bessel helper mismatch"

# (2) exact DtN static limit -> real ladder -(n+1) -------------------------------------------
print("\n(2) EXACT DtN static limit  Lambda_n(kR=1e-6) -> -(n+1)  (the Track-A real ladder):")
print("     n   Lambda_n(1e-6)        -(n+1)    |err|")
emax = 0.0
for n in (1, 2, 3, 4, 5):
    lam = dtn_exact(n, 1e-6)
    e = abs(lam.real - (-(n + 1)))
    emax = max(emax, e)
    print("     %d   %18.10f   %6.1f    %.2e" % (n, lam.real, -(n + 1), e))
assert emax < 1e-3 and abs(dtn_exact(3, 1e-6).imag) < 1e-6, "static-limit ladder failed"
print("    => static spectrum is REAL and on the -(n+1) ladder (Kelvin's home).")

# (3) finite kR -> genuinely COMPLEX; |Im| large for propagating n < kR ----------------------
print("\n(3) at kR=4 the spectrum is COMPLEX (Im = radiation; |Im| large for n<kR, ->0 for n>kR):")
print("     n   Re Lambda_n     Im Lambda_n   regime")
kR = 4.0
for n in (1, 2, 3, 4, 6, 8, 12):
    lam = dtn_exact(n, kR)
    reg = "propagating (n<kR)" if n < kR else "evanescent (n>kR)"
    print("    %2d   %10.4f   %12.4f   %s" % (n, lam.real, lam.imag, reg))
assert abs(dtn_exact(1, kR).imag) > 0.3, "propagating mode should carry radiation Im"
assert abs(dtn_exact(12, kR).imag) < 1e-2, "evanescent mode Im should be tiny"
print("    => Im!=0 confirms the complex-plane spectrum; evanescent modes hug the real axis.")

# (4) PML model correctness: strong thick PML == EXACT for propagating modes -----------------
print("\n(4) PML radial solver correctness (strong thick layer reproduces EXACT, propagating mode):")
print("     n  kR   Lambda_PML (Re,Im)            Lambda_exact (Re,Im)         rel.err   |r~-R|")
emax_pml = 0.0
for (n, kR) in ((1, 8.0), (2, 8.0), (3, 10.0)):
    lp, rterr = dtn_pml(n, kR, d=2.5, sigma0=28.0, nstep=8000)
    le = dtn_exact(n, kR)
    rel = abs(lp - le) / abs(le)
    emax_pml = max(emax_pml, rel)
    print("    %2d %4.1f  (%9.4f,%9.4f)   (%9.4f,%9.4f)   %.2e  %.1e"
          % (n, kR, lp.real, lp.imag, le.real, le.imag, rel, rterr))
assert emax_pml < 5e-3, "PML solver does not reproduce exact for propagating modes"
print("    => the PML radial model is correct (and r~ integrates back to R, a self-check).")

# (5a) Kelvin = the quasi-static operator: kR->0 exactness + evanescent kR^2/2n law -----------
print("\n(5a) Kelvin is the QUASI-STATIC operator (real axis, frequency-independent):")
print("   (i) |Kelvin-exact| -> 0 as kR->0 at fixed n=1 (quasi-static exactness, ~kR^2):")
print("        kR      |Kelvin-exact|   /kR^2")
for kRv in (0.05, 0.1, 0.2, 0.4):
    le = dtn_exact(1, kRv)
    ek = abs(dtn_kelvin_static(1) - le)
    print("        %.2f    %.4e      %.3f" % (kRv, ek, ek / kRv**2))
assert abs(dtn_kelvin_static(1) - dtn_exact(1, 0.05)) < 5e-3, "Kelvin not exact as kR->0"
print("   (ii) evanescent error law  |Kelvin-exact|*2n/kR^2 -> 1  for n>>kR (kR=4):")
print("        n      |Kelvin-exact|   *2n/kR^2")
kR = 4.0
law = None
for n in (6, 9, 12, 16, 20, 30):
    le = dtn_exact(n, kR)
    ek = abs(dtn_kelvin_static(n) - le)
    law = ek * 2 * n / kR**2
    print("        %2d     %.4e      %.3f" % (n, ek, law))
assert abs(law - 1.0) < 0.05, "evanescent Kelvin error should follow kR^2/2n"
print("    => Kelvin's finite-kR error = kR^2/(2n) (evanescent) + radiation (propagating); the")
print("       static Kelvin operator is exact only quasi-statically (kR->0) -- the SA regime.")

# (5b) the PML error KNEE: a thin PML is most stressed at the transition n ~ kR ---------------
print("\n(5b) PML error KNEE -- thin PML (d=0.3, sigma0=20) per-degree error PEAKS at n~kR=4:")
print("     n   regime        |PML-exact|")
kR = 4.0
errs = {}
for n in (1, 2, 3, 4, 5, 7, 9, 12):
    lp, _ = dtn_pml(n, kR, d=0.3, sigma0=20.0, nstep=6000)
    errs[n] = abs(lp - dtn_exact(n, kR))
    reg = "propagating" if n < kR else "evanescent"
    print("    %2d   %-11s   %.4e" % (n, reg, errs[n]))
n_peak = max(errs, key=errs.get)
assert n_peak in (3, 4, 5), "PML error should peak near the transition n~kR"
assert errs[12] < errs[n_peak] and errs[1] < errs[n_peak], "knee should be a peak at n~kR"
print("    => max PML error is at n=%d (~kR=%g): the propagating/evanescent transition is the PML's"
      % (n_peak, kR))
print("       stress point (a thick PML, (4), removes it at extra cost). Opposite end from Kelvin.")

# --------------------------------------------------------------------- Argand snapshot (printed)
print("\nComplex-plane snapshot at kR=4 (what the Argand figure shows):")
print("     n   exact(BEM) Re/Im        Kelvin (real axis)   PML Re/Im (thin)")
kR = 4.0
for n in (1, 2, 3, 5, 8):
    le = dtn_exact(n, kR)
    lk = dtn_kelvin_static(n)
    lp, _ = dtn_pml(n, kR, d=0.3, sigma0=20.0, nstep=6000)
    print("    %2d   (%8.3f,%8.3f)     %8.3f + 0j      (%8.3f,%8.3f)"
          % (n, le.real, le.imag, lk, lp.real, lp.imag))

print("\n" + "=" * 78)
print("RESULT: the exterior-DtN spectrum is the single axis that orders the open-boundary methods.")
print("At finite kR it lives in the COMPLEX plane (Im = radiation). Static Kelvin is PINNED to the")
print("real axis = the kR->0 quasi-static operator (cheapest; error kR^2/2n + radiation, exact only")
print("for kR<~1 = the SA paper). PML carries the complex spectrum with an error knee at n~kR;")
print("FEM-BEM is the exact-but-dense locus. The kR axis is the honest two-paper boundary; the")
print("genuine high-frequency Kelvin is the extended/Helmholtz-Kelvin inversion (sugahara2025).")
print("=" * 78)
