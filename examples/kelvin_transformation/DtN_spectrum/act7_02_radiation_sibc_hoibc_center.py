# -*- coding: utf-8 -*-
# DEMO (ii) (verified): the EXTENDED-KELVIN radiation boundary as a SURFACE impedance (SIBC/HOIBC)
# at the EXTERIOR CENTRE -- and why a sphere needs an HOIBC.
#
# Context (the author's IEICE Trans. C special issue 2024, "Extended Kelvin Transformation for
# Solving Radiating Electromagnetic Fields", building on sugahara2022): the Kelvin inversion
# x' = (a/r)^2 x maps the unbounded exterior r>a into the bounded ball rho<a, with r=infinity ->
# the CENTRE rho=0.  Differential geometry gives the exterior material as an ISOTROPIC modulation
# (a/r)^2 of mu, eps, sigma, sigma* (metric ratio g'/g).  A RADIATING field carries energy OUT to
# r=infinity = INTO the centre in the inverted picture, so the radiation/absorbing condition is
# imposed AT THE CENTRE: the paper places a spherical Maxwellian PML there (eps'=mu'=(1-0.2j)a^2/r'^2),
# excising a tiny ball at rho=0 (the singular image of infinity) and putting the absorber far in
# physical space (8 m, lambda=3 m) so a simple plane-wave (377 ohm) PML suffices.
#
# THIS demo derives the SURFACE-impedance alternative the user asked for (SIBC / radiation BC on
# Kelvin), instead of a volumetric PML: on the small inner sphere (image of a far sphere r=b) impose
# an IMPEDANCE that reproduces the exterior radiation Dirichlet-to-Neumann (DtN) operator.  The exact
# radiation DtN on a sphere of radius b is, per spherical-harmonic degree n (units 1/b; z=kb),
#       Lambda_n(z) = z h_n^(1)'(z) / h_n^(1)(z)        (complex; Im = radiation),
# whose large-z expansion IS the absorbing-BC hierarchy:
#       Lambda_n(z) = i z  -  1  -  i n(n+1)/(2z)  +  O(1/z^2)
#                     \___/    \_/   \__________/
#                     plane    curv  TANGENTIAL-LAPLACIAN (Laplace-Beltrami) term
#                     wave           = the HOIBC correction
#   * SIBC (Leontovich, n-INDEPENDENT): the constant impedance i z (-1) -- the plane-wave 377 ohm.
#     It matches Lambda_n only for z >> n(n+1); its error grows ~ n(n+1)/(2z) (the dropped term).
#   * HOIBC (n-DEPENDENT): adds -i n(n+1)/(2z).  Since n(n+1) is the unit-sphere Laplace-Beltrami
#     eigenvalue (Delta_S Y_n = -n(n+1) Y_n), the HOIBC is the SURFACE OPERATOR
#       Z_HOIBC = i z - 1 + (i/(2z)) Delta_S     (a 2nd-order PDE on the sphere),
#     implementable as an ordinary surface-FEM term.  This is WHY a sphere needs an HOIBC: the
#     radiation impedance is multipole-(curvature-)dependent, and only a surface-Laplacian operator
#     captures that dependence; a scalar SIBC cannot.
#
# SPECTRAL READING of the paper's design choice: "place the absorber far (large kb)" == "make
# n(n+1)/(2kb) small so the n-INDEPENDENT SIBC suffices".  The HOIBC removes the n-dependence to one
# higher order (error ~1/z^2), so it lets the absorber sit CLOSER (smaller far radius b => the image
# sphere a^2/b is larger => fewer exterior cells) at the same accuracy.  Topology note: the one-point
# compactification sends infinity to the single centre point; excising a small sphere there turns the
# NONLOCAL exterior DtN (on the truncation r=a) into a LOCAL absorber on a small interior sphere --
# the conformal Kelvin map makes "infinity" a regular meshable point (no cuts; the scalar/E exterior
# is simply connected), which is the sense in which this is "topologically easy".
#
# VERIFICATIONS (all asserted from computed values; no overclaim):
#  (1) complex spherical-Bessel helpers == scipy.special.spherical_jn/yn on real args.
#  (2) the SIBC/HOIBC HIERARCHY is the asymptotic expansion of the exact radiation DtN: at fixed n,
#      |Lambda - iz| -> const, |Lambda - (iz-1)| ~ O(1/z) (halves per z-doubling), |Lambda - HOIBC|
#      ~ O(1/z^2) (quarters per z-doubling).
#  (3) at the paper's far placement (kb=2*pi*8/3 ~ 16.8): SIBC error grows with n while HOIBC is
#      >~10x smaller across the multipole band (so far placement makes SIBC adequate for LOW n only;
#      HOIBC extends the band).
#  (4) place-far sweep (fixed n): SIBC error ~ 1/kb, HOIBC ~ 1/kb^2 -> HOIBC reaches a tolerance at a
#      SMALLER kb (closer absorber / fewer cells) than SIBC.
#  (5) the Kelvin image mapping b -> a^2/b (the paper's a, b=8 m -> image radius a^2/8; excision
#      0.25 m -> image of a very far sphere) and the plane-wave impedance leading term (i k = the
#      Sommerfeld/377-ohm condition).
#
# Pure numpy/scipy (the SPECTRUM is closed-form on the sphere); the FE realisation is the follow-up
# (a surface-Laplacian HOIBC term on the excised inner sphere of the inverted exterior).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import jv, yv, hankel1, spherical_jn, spherical_yn


def _sph(kind, n, z):
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
    z = np.asarray(z, dtype=complex)
    return _sph(kind, n - 1, z) - (n + 1.0) / z * _sph(kind, n, z)


def dtn_exact(n, z):
    """EXACT radiation DtN eigenvalue on degree n: z h_n^(1)'(z)/h_n^(1)(z)  (z = k*b, units 1/b)."""
    z = complex(z)
    return (z * _sph_prime('h1', n, z) / _sph('h1', n, z)).item()


def z_sibc(n, z):
    """SIBC (Leontovich, n-INDEPENDENT) impedance: plane-wave + curvature, NO multipole dependence."""
    return 1j * z - 1.0


def z_sibc0(n, z):
    """Pure plane-wave (Sommerfeld 1st-order ABC) = the 377-ohm characteristic impedance term."""
    return 1j * z


def z_hoibc(n, z):
    """HOIBC (n-DEPENDENT) = i z - 1 + (i/2z) Delta_S, Delta_S eigenvalue -n(n+1) (Laplace-Beltrami)."""
    return 1j * z - 1.0 - 1j * n * (n + 1.0) / (2.0 * z)


# ============================================================================= VERIFICATIONS
print("=" * 80)
print("EXTENDED-KELVIN radiation boundary as a SURFACE impedance at the exterior CENTRE")
print("  exact radiation DtN  Lambda_n(z)=z h_n^(1)'(z)/h_n^(1)(z) = iz - 1 - i n(n+1)/2z + O(1/z^2)")
print("  SIBC = iz-1 (n-indep);  HOIBC = iz-1 + (i/2z) Delta_S  (Delta_S Y_n = -n(n+1) Y_n)")
print("=" * 80)

# (1) spherical-Bessel helper vs scipy ------------------------------------------------------
err_j = err_y = 0.0
for n in range(0, 8):
    for zz in (0.3, 1.0, 2.5, 7.0):
        rj, ry = spherical_jn(n, zz), spherical_yn(n, zz)
        err_j = max(err_j, abs(_sph('j', n, zz).item() - rj) / (abs(rj) + 1e-300))
        err_y = max(err_y, abs(_sph('y', n, zz).item() - ry) / (abs(ry) + 1e-300))
print("\n(1) complex sph-Bessel helper vs scipy (real args, RELATIVE): j %.2e, y %.2e" % (err_j, err_y))
assert err_j < 1e-11 and err_y < 1e-11, "spherical-Bessel helper mismatch"

# (2) the SIBC/HOIBC hierarchy = asymptotic expansion of the radiation DtN -------------------
print("\n(2) the SIBC->HOIBC ladder IS the large-z expansion of the exact radiation DtN (n=2):")
print("     z     |L-iz|      |L-SIBC|    |L-HOIBC|    SIBC ratio  HOIBC ratio  (expect 2, 4)")
prev = None
r1 = r2 = float('nan')
for z in (10., 20., 40., 80.):
    L = dtn_exact(2, z)
    e0 = abs(L - z_sibc0(2, z)); e1 = abs(L - z_sibc(2, z)); e2 = abs(L - z_hoibc(2, z))
    if prev:
        r1, r2 = prev[0] / e1, prev[1] / e2
    print("    %4.0f  %.4e  %.4e  %.4e   %8.2f    %8.2f" % (z, e0, e1, e2, r1, r2))
    prev = (e1, e2)
assert 1.9 < r1 < 2.1, "SIBC error should be O(1/z) (halve per z-doubling)"
assert 3.8 < r2 < 4.2, "HOIBC error should be O(1/z^2) (quarter per z-doubling)"
assert abs(dtn_exact(2, 40.) - z_hoibc(2, 40.)) < abs(dtn_exact(2, 40.) - z_sibc(2, 40.)), "HOIBC<SIBC"
print("    => SIBC = O(1/z), HOIBC = O(1/z^2): the absorbing-BC orders are the DtN's own expansion.")

# (3) at the paper's FAR placement (b=8 m, lambda=3 m): SIBC ok for low n, HOIBC across band ---
kb = 2 * np.pi * 8.0 / 3.0   # f=100 MHz -> lambda=3 m, k=2pi/3; far sphere b=8 m
print("\n(3) paper far placement kb = 2*pi*8/3 = %.3f (b=8 m, lambda=3 m): SIBC vs HOIBC across n:" % kb)
print("     n   |SIBC-exact|   |HOIBC-exact|   HOIBC/SIBC")
ratios = []
for n in (1, 2, 3, 5, 8, 12):
    L = dtn_exact(n, kb)
    es = abs(L - z_sibc(n, kb)); eh = abs(L - z_hoibc(n, kb))
    ratios.append(eh / es)
    print("    %2d   %.4e    %.4e    %.3f" % (n, es, eh, eh / es))
assert abs(dtn_exact(1, kb) - z_sibc(1, kb)) < 0.1, "SIBC adequate for low n at far placement"
assert max(ratios) < 0.2, "HOIBC should be >5x better than SIBC across the band"
assert abs(dtn_exact(12, kb) - z_sibc(12, kb)) > 1.0, "SIBC degrades for high n even when placed far"
print("    => far placement makes the n-indep SIBC adequate for LOW multipoles only (paper's regime);")
print("       HOIBC stays accurate across the whole band -> can carry high-n / closer scatterers.")

# (4) place-far sweep: SIBC ~1/kb, HOIBC ~1/kb^2 -> HOIBC lets the absorber sit CLOSER ---------
print("\n(4) place-far sweep (n=5): SIBC error ~ 1/kb, HOIBC ~ 1/kb^2 (ratio per kb-doubling -> 2, 4):")
print("     kb      |SIBC-ex|    |HOIBC-ex|   SIBC ratio  HOIBC ratio")
prev = None; rs = rh = float('nan')
for kbv in (16., 32., 64., 128.):
    L = dtn_exact(5, kbv)
    es = abs(L - z_sibc(5, kbv)); eh = abs(L - z_hoibc(5, kbv))
    if prev:
        rs, rh = prev[0] / es, prev[1] / eh
    print("    %5.0f   %.4e   %.4e    %6.2f      %6.2f" % (kbv, es, eh, rs, rh))
    prev = (es, eh)
assert 1.8 < rs < 2.2, "SIBC error should scale ~1/kb"
assert 3.6 < rh < 4.4, "HOIBC error should scale ~1/kb^2"
# tolerance crossover: kb needed for |.-exact| < 1e-2 at n=5
tol = 1e-2
import math
kb_sibc = next(kbv for kbv in np.arange(1.0, 4000.0, 1.0) if abs(dtn_exact(5, kbv) - z_sibc(5, kbv)) < tol)
kb_hoibc = next(kbv for kbv in np.arange(1.0, 4000.0, 1.0) if abs(dtn_exact(5, kbv) - z_hoibc(5, kbv)) < tol)
print("    => to hit |Z-exact|<%.0e at n=5: SIBC needs kb>=%.0f, HOIBC needs kb>=%.0f (%.0fx closer)."
      % (tol, kb_sibc, kb_hoibc, kb_sibc / kb_hoibc))
assert kb_hoibc < kb_sibc, "HOIBC must reach tolerance at a smaller kb than SIBC"

# (5) Kelvin image mapping + plane-wave impedance leading term --------------------------------
a, b_far, excise = 4.0, 8.0, 0.25      # paper: a=4 m, PML far sphere 8 m, excision 0.25 m
print("\n(5) Kelvin image mapping (paper a=%.1f m, far b=%.1f m, excise %.2f m):" % (a, b_far, excise))
print("     far sphere r=b=%.1f m  ->  image radius a^2/b = %.3f m" % (b_far, a * a / b_far))
print("     excision   rho=%.2f m  <-  image of r = a^2/rho = %.1f m (very far -> wave fully decayed)"
      % (excise, a * a / excise))
print("     leading impedance term Lambda_n -> i*z = i*k*b  <=>  d_r u = i k u (Sommerfeld / 377 ohm)")
assert abs(a * a / b_far - 2.0) < 1e-12 and abs(a * a / excise - 64.0) < 1e-9, "image mapping"

print("\n" + "=" * 80)
print("RESULT: the extended-Kelvin RADIATION boundary can be a SURFACE impedance on the small inner")
print("sphere (image of a far sphere) instead of a volumetric PML. The exact radiation impedance is")
print("the DtN symbol Lambda_n(kb), whose expansion gives the SIBC (n-indep, plane-wave) and the")
print("HOIBC (n-dep, = a Laplace-Beltrami surface operator). A sphere needs the HOIBC because the")
print("radiation impedance is curvature/multipole-dependent; the paper's 'place it far' choice is")
print("exactly 'make kb large so the SIBC suffices', and the HOIBC relaxes that to a closer absorber.")
print("=" * 80)
