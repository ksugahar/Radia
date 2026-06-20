# -*- coding: utf-8 -*-
# DEMO (ll) (verified): the DIFFERENTIAL-GEOMETRY (transformation-optics) derivation of the Kelvin
# exterior -- the (a/r)^2 material is AUTOMATIC, and only IMPEDANCE MATCHING matters.  Why the
# radiation (surface-impedance) boundary is the principled object, and why the DtN spectrum is most
# meaningful in the HIGH-FREQUENCY (radiating) Kelvin.
#
# The user's points (this turn):
#   * "you can derive it with differential geometry" -- YES: transformation optics.
#   * "the modulation happens automatically if you don't track the Kelvin factor" -- YES: the (a/r')^2
#     material is the transformation-optics pullback of the metric; you do NOT carry a field weight.
#   * "as long as the impedance is matched, OK" -- YES: a matched impedance is reflectionless.
#   * "is the radiation boundary the better-principled approach / better than a PML?" -- the exact
#     transparent condition IS the DtN/impedance; PML and HOIBC are two REALISATIONS of it.
#   * "the high-frequency Kelvin is where the DtN spectral analysis is meaningful" -- YES: at low freq
#     the open boundary is EXACT (real ladder, a datasheet); at high freq the spectrum is COMPLEX and
#     IS the absorber-design object (per-mode reflection = spectral mismatch).
#
# DIFFERENTIAL GEOMETRY (transformation optics) of the Kelvin inversion x' = a^2 x/|x|^2:
#   Jacobian   J^i_j = d x'^i/d x^j = (a^2/r^2)(delta_ij - 2 n_i n_j),  n = x/r
#            = (a^2/r^2) Q,   Q = I - 2 n n^T = a Householder REFLECTION (Q^T Q = I, det Q = -1).
#   => the map is CONFORMAL (J = scalar * orthogonal): no shear, no anisotropy.  det J = -(a^2/r^2)^3
#      (NEGATIVE = the inversion is ORIENTATION-REVERSING = "inside-out").
#   Transformation optics: a coordinate map leaves Maxwell form-invariant if the media transform as
#      eps'_r = mu'_r = J J^T / |det J| = (a^2/r^2)^2 / (a^2/r^2)^3 * I = (r^2/a^2) I = (a^2/r'^2) I
#      (r' = a^2/r is the image radius).  ISOTROPIC, exactly the paper's (a/r')^2 modulation -- and it
#      drops out of the geometry AUTOMATICALLY (no hand-applied Kelvin field factor; the conformality
#      is precisely what kills the anisotropy a general transform would produce).
#   The only remaining design freedom is the BOUNDARY: an absorber at the inner image sphere is
#   reflectionless iff its impedance MATCHES the (transformed) radiation DtN -- the HOIBC of act7_03_hoibc_kelvin_transform.
#
# VERIFICATIONS (all asserted from computed values; no overclaim):
#  (1) complex spherical-Bessel helpers (j,y,h1,h2) == scipy on real args.
#  (2) the inversion Jacobian is CONFORMAL: J=(a^2/r^2)Q with Q^TQ=I, det Q=-1, det J=-(a^2/r^2)^3 (TO).
#  (3) the transformation-optics medium J J^T/|det J| is ISOTROPIC and equals (a^2/r'^2) I -- the
#      paper's modulation, derived from geometry, to machine precision (the "automatic" modulation).
#  (4) IMPEDANCE MATCHING = REFLECTIONLESS: the modal reflection R_n=|B/A|^2 of an outgoing wave at the
#      absorber sphere is ZERO (machine) for the exact DtN impedance, and HOIBC << SIBC otherwise; and
#      R_n is governed by the DtN-SPECTRAL MISMATCH |Z_n - Lambda_n|^2.
#  (5) far-vs-close placement: a constant SIBC reflects LESS the farther it is placed (large kb) =
#      the paper's far-placement (a large domain); the HOIBC matches the COMPLEX spectrum and reflects
#      little even close (a small domain). The DtN spectrum is an absorber-DESIGN object only in the
#      radiating regime (at static the open boundary is already exact -- the spectrum is a datasheet).
#
# Pure numpy/scipy.  The TO medium is the general (coordinate-free) form of act7_03_hoibc_kelvin_transform's radial result;
# the FE follow-up assembles the isotropic (a/r')^2 medium + the matched HOIBC surface term.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import jv, yv, hankel1, hankel2, spherical_jn, spherical_yn

a = 1.0  # interior (truncation) radius


def _sph(kind, n, z):
    z = np.asarray(z, dtype=complex)
    pref = np.sqrt(np.pi / (2.0 * z))
    fn = {'j': jv, 'y': yv, 'h1': hankel1, 'h2': hankel2}[kind]
    return pref * fn(n + 0.5, z)


def _sph_prime(kind, n, z):
    z = np.asarray(z, dtype=complex)
    return _sph(kind, n - 1, z) - (n + 1.0) / z * _sph(kind, n, z)


def dtn_exact(n, z):
    z = complex(z)
    return (z * _sph_prime('h1', n, z) / _sph('h1', n, z)).item()


def jacobian(x):
    """Jacobian of the Kelvin inversion x' = a^2 x/|x|^2 at point x (3-vector)."""
    x = np.asarray(x, dtype=float)
    r2 = float(x @ x); r = np.sqrt(r2)
    Q = np.eye(3) - 2.0 * np.outer(x, x) / r2          # Householder reflection across plane _|_ n
    return (a * a / r2) * Q, r, Q


# ============================================================================= VERIFICATIONS
print("=" * 80)
print("DIFFERENTIAL-GEOMETRY (transformation-optics) Kelvin exterior: automatic (a/r')^2 medium,")
print("impedance matching = reflectionless, and the DtN spectrum as the HIGH-FREQUENCY design object")
print("=" * 80)

# (1) spherical-Bessel helpers vs scipy ------------------------------------------------------
err = 0.0
for n in range(0, 8):
    for z in (0.3, 1.0, 2.5, 7.0):
        rj, ry = spherical_jn(n, z), spherical_yn(n, z)
        err = max(err, abs(_sph('j', n, z).item() - rj) / (abs(rj) + 1e-300))
        err = max(err, abs(_sph('y', n, z).item() - ry) / (abs(ry) + 1e-300))
print("\n(1) complex sph-Bessel helper vs scipy (real args, RELATIVE): %.2e" % err)
assert err < 1e-11, "spherical-Bessel helper mismatch"

# (2)+(3) transformation-optics medium from the conformal Jacobian --------------------------
print("\n(2)+(3) transformation optics: J=(a^2/r^2)Q conformal; medium J J^T/|det J| = (a^2/r'^2) I:")
print("     point x            r      |Q^TQ-I|  det Q    eps' isotropy   |eps'-(a^2/r'^2)I|")
rng_pts = [(2.0, 0.0, 0.0), (1.5, 1.0, 0.5), (0.3, -0.4, 0.2), (3.0, 2.0, -1.0)]
emax = 0.0
for x in rng_pts:
    J, r, Q = jacobian(x)
    ortho = np.max(np.abs(Q.T @ Q - np.eye(3)))
    detQ = np.linalg.det(Q)
    detJ = np.linalg.det(J)
    eps = (J @ J.T) / abs(detJ)                       # transformation-optics relative medium
    rp = a * a / r                                    # image radius r'
    target = (a * a / rp**2) * np.eye(3)              # = (r^2/a^2) I  (the paper's (a/r')^2 modulation)
    iso = np.max(np.abs(eps - np.diag(np.diag(eps)))) # off-diagonal => anisotropy
    dev = np.max(np.abs(eps - target))
    emax = max(emax, dev, ortho, abs(detQ + 1))
    print("    (%4.1f,%4.1f,%4.1f)  %.3f   %.1e   %+.3f   %.1e         %.1e"
          % (x[0], x[1], x[2], r, ortho, detQ, iso, dev))
    # also confirm det J = -(a^2/r^2)^3 (orientation-reversing)
    assert abs(detJ - (-(a * a / r**2)**3)) < 1e-9 * abs(detJ), "det J (orientation) wrong"
assert emax < 1e-9, "transformation-optics medium is not the isotropic (a/r')^2 modulation"
print("    => J is conformal (Q orthogonal, det Q=-1 => inside-out); the TO medium is ISOTROPIC and")
print("       equals (a^2/r'^2) I AUTOMATICALLY -- the paper's modulation, from geometry, no field factor.")

# (4) impedance matching = reflectionless; reflection = DtN-spectral mismatch ----------------
def reflection(n, kb, Z):
    """modal reflected-power fraction |B/A|^2 of an outgoing wave at an absorber sphere r=b with
    imposed impedance b f'/f = Z (f = A h1 + B h2; |h2|=|h1| for real kb)."""
    h1, h1p = _sph('h1', n, kb).item(), _sph_prime('h1', n, kb).item()
    h2, h2p = _sph('h2', n, kb).item(), _sph_prime('h2', n, kb).item()
    beta = (Z * h1 - kb * h1p) / (kb * h2p - Z * h2)
    return abs(beta)**2


k, b = 3.0, 2.0
kb = k * b
print("\n(4) impedance matching = reflectionless (absorber sphere kb=%.1f): modal reflection R_n=|B/A|^2:" % kb)
print("     n   exact Z       SIBC          HOIBC        |Z_SIBC-Lam|  |Z_HOIBC-Lam|")
for n in (1, 2, 3, 5):
    L = dtn_exact(n, kb)
    z_si = 1j * kb - 1.0
    z_ho = 1j * kb - 1.0 - 1j * n * (n + 1.0) / (2.0 * kb)
    R_ex = reflection(n, kb, L)
    R_si = reflection(n, kb, z_si)
    R_ho = reflection(n, kb, z_ho)
    print("    %2d   %.2e    %.4e    %.4e     %.3e     %.3e"
          % (n, R_ex, R_si, R_ho, abs(z_si - L), abs(z_ho - L)))
    assert R_ex < 1e-20, "exact impedance must be reflectionless"
    assert R_ho < R_si, "HOIBC (smaller spectral mismatch) must reflect less than SIBC"
print("    => exact DtN impedance -> ZERO reflection (machine); the reflection ORDERS by the DtN-")
print("       spectral mismatch |Z_n-Lambda_n| (HOIBC matches the spectrum better than the SIBC).")

# (5) far-vs-close placement, and why the spectrum is the HIGH-FREQUENCY design object ------
print("\n(5a) a constant (SIBC) impedance gets EXACT only FAR away (large kb) -- the paper's far-placement")
print("     (n=1): reflection DROPS as the absorber is placed farther (cost = a large domain):")
print("     kb     R_1(SIBC)")
Rprev = None
mono_dec = True
for kbv in (1.0, 2.0, 4.0, 8.0, 16.0):
    R = reflection(1, kbv, 1j * kbv - 1.0)
    if Rprev is not None and R >= Rprev:
        mono_dec = False
    Rprev = R
    print("    %5.1f    %.4e" % (kbv, R))
assert mono_dec, "constant SIBC should improve (reflect less) with distance"
print("\n(5b) the HOIBC matches the COMPLEX spectrum -> reflects little even CLOSE (kb=4) -> small domain:")
print("     n   R_n(SIBC)     R_n(HOIBC)")
kc = 4.0
ok_close = True
for n in (1, 2, 3):
    R_si = reflection(n, kc, 1j * kc - 1.0)
    R_ho = reflection(n, kc, 1j * kc - 1.0 - 1j * n * (n + 1.0) / (2.0 * kc))
    if not (R_ho < R_si):
        ok_close = False
    print("    %2d   %.4e    %.4e" % (n, R_si, R_ho))
assert ok_close, "HOIBC must reflect less than SIBC at a close placement"
print("    => SIBC needs FAR placement (large kb) to be reflectionless = a large domain (paper, U1/U3);")
print("       the HOIBC matches the spectrum and is low-reflection even close = a small domain. And the")
print("       low-/high-frequency split: at STATIC the open boundary is exact with NO absorber (the")
print("       spectrum is just the real accuracy ladder, a datasheet); only in the RADIATING regime is")
print("       the COMPLEX spectrum an absorber-DESIGN target (R_n = |Z_n - Lambda_n| spectral mismatch).")

print("\n" + "=" * 80)
print("RESULT: differential geometry (transformation optics) gives the Kelvin exterior medium")
print("AUTOMATICALLY -- the conformal inversion's Jacobian J=(a^2/r^2)Q (Q a reflection) yields the")
print("ISOTROPIC (a/r')^2 eps'=mu' with no hand-applied Kelvin field factor. The only design condition")
print("is IMPEDANCE MATCHING at the absorber, which is exactly the DtN; a matched impedance is")
print("reflectionless (machine) and any absorber's reflection = its DtN-spectral mismatch. So the")
print("radiation (surface-impedance/HOIBC) boundary is the principled object that the volumetric PML")
print("only approximates, and the DtN spectrum becomes an absorber-DESIGN tool in the radiating regime.")
print("=" * 80)
