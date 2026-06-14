# -*- coding: utf-8 -*-
# DEMO (tt) (verified): evanescent / low-frequency -- does Kelvin beat a plain PML? An HONEST,
# mode-resolved answer (and a CORRECTION of a too-quick earlier claim).
#
# Intuition (and an earlier over-quick claim of mine): "for evanescent / low-frequency modes Kelvin+PML
# beats a plain PML, because Kelvin gives the evanescent field infinite room to decay." Tested fully on
# the radial exterior DtN, the truth is mode-dependent and narrower:
#   * a PLAIN PML does NOT absorb evanescent waves (the imaginary stretch is for propagating waves); for
#     evanescent modes only the natural DECAY over the layer helps. So a THIN plain PML reflects the
#     residual evanescent field -> large error. (This is the unfair comparison that made "Kelvin wins"
#     look true.) But a plain PML with ADEQUATE THICKNESS (reaching far enough that the evanescent field
#     decays) is excellent at the same DoF.
#   * the right Kelvin tool for evanescent is KELVIN-ONLY (the quasi-static compactification, NO PML):
#     adding a PML to Kelvin for evanescent HURTS (the complex stretch distorts the smooth decaying
#     field). Kelvin-only is EXACT + PARAMETER-FREE for DEEPLY evanescent modes (n >> ka), but FAILS for
#     transition/low-n modes (n ~ ka): there the field decays slowly / carries radiation and reaches the
#     excision wall -> needs absorption a bare compactification lacks.
#   * NET, per mode (ka=0.5): Kelvin-only WINS (and is parameter-free) for DEEPLY evanescent n/ka >~ 5;
#     a thickness-tuned plain PML WINS for the transition / low-n modes (n ~ ka) -- which are the
#     physically dominant ones (dipole, quadrupole). So "Kelvin beats plain PML for evanescent" is NOT a
#     blanket win; it holds only in the deeply-evanescent, parameter-free corner.
#   * (Spherical vs box PML, the related question: a spherical PML IS better than a box -- no corner
#     reflections, conforms to the spherical wavefronts, diagonal in the spherical-harmonic DtN -- but
#     that favours ANY spherical PML; Kelvin's spherical truncation inherits it, it is not a Kelvin-vs-
#     PML differentiator. Box-vs-sphere needs a genuine 3D test, out of scope of this 1D radial demo.)
#
# VERIFIED (a=1, ka=0.5; matched DoF M=120):
#  (1) a THIN plain PML (d=0.5) is bad for evanescent but a THICK one (d=4..6) is excellent -> the thin
#      comparison was unfair; with adequate thickness the plain PML is competitive.
#  (2) Kelvin-only (parameter-free) vs plain-PML-best-d, per mode: plain wins for n~ka (n=1,2), Kelvin
#      wins for deeply evanescent (n>=3, here) -> a CROSSOVER near n/ka ~ 5.
#  (3) adding a PML to Kelvin for a deeply-evanescent mode HURTS vs Kelvin-only (the stretch distorts the
#      smooth field) -> use Kelvin-ONLY for evanescent.
#
# Pure numpy/scipy.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import jv, yv, hankel1, spherical_jn


def _sph(kind, n, z):
    z = np.asarray(z, dtype=complex)
    pref = np.sqrt(np.pi / (2.0 * z))
    return pref * {'j': jv, 'y': yv, 'h1': hankel1}[kind](n + 0.5, z)


def _sph_prime(kind, n, z):
    z = np.asarray(z, dtype=complex)
    return _sph(kind, n - 1, z) - (n + 1.0) / z * _sph(kind, n, z)


def dtn_exact(n, ka):
    z = complex(ka); return (z * _sph_prime('h1', n, z) / _sph('h1', n, z)).item()


_GX = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)]); _GW = np.array([5 / 9, 8 / 9, 5 / 9])


def plain_pml(n, k, a, d, M, alpha, p=2):
    """plain PML layer [a, a+d], wall at a+d, DtN at a."""
    def s(r): return 1 + 1j * alpha * ((r - a) / d) ** p
    def rt(r): dr = r - a; return a + (dr + 1j * alpha / d**p * dr**(p + 1) / (p + 1))
    nod = np.linspace(a, a + d, M + 1); A = np.zeros((M + 1, M + 1), complex)
    for e in range(M):
        r0, r1 = nod[e], nod[e + 1]; h = r1 - r0; dp = np.array([-1 / h, 1 / h])
        K = np.zeros((2, 2), complex); C = np.zeros((2, 2), complex); Mm = np.zeros((2, 2), complex)
        for gp, gw in zip(_GX, _GW):
            r = .5 * (r0 + r1) + .5 * h * gp; w = .5 * h * gw
            ph = np.array([(r1 - r) / h, (r - r0) / h]); sr = s(r); rr = rt(r)
            K += np.outer(dp, dp) * (rr * rr / sr) * w
            C += np.outer(ph, ph) * sr * w; Mm += np.outer(ph, ph) * sr * rr * rr * w
        A[e:e + 2, e:e + 2] += K + n * (n + 1) * C - k * k * Mm
    u = np.zeros(M + 1, complex); u[0] = 1.0; u[M] = 0.0; idx = list(range(1, M))
    u[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -A[np.ix_(idx, [0])][:, 0])
    return -(A[0, :] @ u) / a


def plain_pml_best(n, k, a, M, alpha):
    return min(abs(plain_pml(n, k, a, d, M, alpha) - dtn_exact(n, k * a)) / abs(dtn_exact(n, k * a))
               for d in (0.5, 1.0, 2.0, 4.0, 6.0))


def kelvin_only(n, k, a, rho_e, M):
    """Kelvin quasi-static compactification, NO PML, hard wall at the excision rho_e (parameter-free)."""
    node = np.linspace(rho_e, a, M + 1); A = np.zeros((M + 1, M + 1), complex)
    for e in range(M):
        r0, r1 = node[e], node[e + 1]; h = r1 - r0; dp = np.array([-1 / h, 1 / h])
        K = np.outer(dp, dp) * (a * a * h); C = np.zeros((2, 2), complex); Mm = np.zeros((2, 2), complex)
        for gp, gw in zip(_GX, _GW):
            r = .5 * (r0 + r1) + .5 * h * gp; w = .5 * h * gw
            ph = np.array([(r1 - r) / h, (r - r0) / h])
            C += np.outer(ph, ph) * (a * a / r**2) * w; Mm += np.outer(ph, ph) * (a**6 / r**4) * w
        A[e:e + 2, e:e + 2] += K + n * (n + 1) * C - k * k * Mm
    u = np.zeros(M + 1, complex); u[M] = 1.0; u[0] = 0.0; idx = list(range(1, M))
    u[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -A[np.ix_(idx, [M])][:, 0])
    return -(A[M, :] @ u) / a


def kelvin_with_pml(n, k, a, rho_e, rho_pml, M, alpha, p=2):
    """Kelvin medium + a PML stretch in [rho_e, rho_pml] (to show the PML HURTS for evanescent)."""
    L = rho_pml - rho_e
    def s(r): return (1.0 + 0j) if r >= rho_pml else 1 + 1j * alpha * ((rho_pml - r) / L) ** p
    def rt(r):
        if r >= rho_pml: return complex(r)
        dr = rho_pml - r; return rho_pml - (dr + 1j * alpha / L**p * dr**(p + 1) / (p + 1))
    node = np.linspace(rho_e, a, M + 1); A = np.zeros((M + 1, M + 1), complex)
    for e in range(M):
        r0, r1 = node[e], node[e + 1]; h = r1 - r0; dp = np.array([-1 / h, 1 / h])
        K = np.zeros((2, 2), complex); C = np.zeros((2, 2), complex); Mm = np.zeros((2, 2), complex)
        for gp, gw in zip(_GX, _GW):
            r = .5 * (r0 + r1) + .5 * h * gp; w = .5 * h * gw
            ph = np.array([(r1 - r) / h, (r - r0) / h]); sr = s(r); rr = rt(r)
            K += np.outer(dp, dp) * (a * a / sr) * w
            C += np.outer(ph, ph) * sr * (n * (n + 1) * a * a / rr**2) * w
            Mm += np.outer(ph, ph) * sr * (a**6 / rr**4) * w
        A[e:e + 2, e:e + 2] += K + C - k * k * Mm
    u = np.zeros(M + 1, complex); u[M] = 1.0; u[0] = 0.0; idx = list(range(1, M))
    u[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -A[np.ix_(idx, [M])][:, 0])
    return -(A[M, :] @ u) / a


# ============================================================================= VERIFICATIONS
a, ka = 1.0, 0.5
k = ka / a
M = 120
print("=" * 82)
print("Evanescent / low-freq: does Kelvin beat a plain PML? An HONEST mode-resolved answer (a=1, ka=%.1f)" % ka)
print("=" * 82)

err = 0.0
for n in range(0, 4):
    rj = spherical_jn(n, 1.0); err = max(err, abs(_sph('j', n, 1.0).item() - rj) / (abs(rj) + 1e-300))
assert err < 1e-11

# (1) a thin plain PML is bad, a thick one is excellent -> the thin comparison was unfair ------
print("\n(1) plain PML thickness matters (evanescent n=1): thin reflects, thick (more decay room) is fine:")
L = dtn_exact(1, ka)
e_thin = abs(plain_pml(1, k, a, 0.5, M, 4.0) - L) / abs(L)
e_thick = abs(plain_pml(1, k, a, 6.0, M, 4.0) - L) / abs(L)
print("     plain PML d=0.5 (reaches r=1.5): %.3e   |   d=6 (reaches r=7): %.3e" % (e_thin, e_thick))
assert e_thick < e_thin / 50, "an adequately thick plain PML is far better than a thin one (decay room)"
print("    => the earlier 'thin plain PML loses to Kelvin' was an UNFAIR comparison; a thickness-tuned")
print("       plain PML is competitive. (A plain PML does not absorb evanescent -- only decay helps.)")

# (2) Kelvin-only (parameter-free) vs plain-PML-best-d, per mode: a CROSSOVER near n/ka ~ 5 ----
print("\n(2) Kelvin-ONLY (no PML, parameter-free) vs plain-PML(best d), per mode -- a CROSSOVER:")
print("     n   n/ka   Kelvin-only   plain best   winner")
win_low = win_high = None
for n in (1, 2, 3, 4, 6, 8):
    ek = abs(kelvin_only(n, k, a, 0.02, M) - dtn_exact(n, ka)) / abs(dtn_exact(n, ka))
    ep = plain_pml_best(n, k, a, M, 4.0)
    w = "Kelvin" if ek < ep else "plain"
    if n == 1: win_low = w
    if n == 8: win_high = w
    print("    %2d   %4.1f   %.2e     %.2e    %s" % (n, n / ka, ek, ep, w))
assert win_low == "plain", "for transition/low-n (n~ka) a tuned plain PML wins"
assert win_high == "Kelvin", "for DEEPLY evanescent (n>>ka) Kelvin-only wins AND is parameter-free"
print("    => plain PML wins for n~ka (the dominant dipole/quadrupole); Kelvin-only wins for DEEPLY")
print("       evanescent n>>ka (and needs no parameters) -> 'Kelvin beats PML for evanescent' is NARROW.")

# (3) adding a PML to Kelvin HURTS for evanescent -> use Kelvin-ONLY -------------------------
print("\n(3) for a deeply-evanescent mode (n=4), adding a PML to Kelvin HURTS vs Kelvin-only:")
n = 4; L = dtn_exact(n, ka)
e_only = abs(kelvin_only(n, k, a, 0.05, M) - L) / abs(L)
e_pml = abs(kelvin_with_pml(n, k, a, 0.2, 0.5, M, 4.0) - L) / abs(L)
print("     Kelvin-only: %.3e   |   Kelvin+PML: %.3e" % (e_only, e_pml))
assert e_only < e_pml, "the PML stretch distorts the smooth evanescent field -> use Kelvin-only"
print("    => the right Kelvin tool for evanescent is KELVIN-ONLY (the quasi-static compactification);")
print("       a PML there only distorts the smooth decaying field.")

print("\n" + "=" * 82)
print("ANSWER (honest, correcting a too-quick earlier claim): 'Kelvin+PML beats plain PML for evanescent'")
print("is NOT a blanket win. A thickness-tuned plain PML is competitive/better for the dominant low-n")
print("(n~ka) modes; Kelvin's genuine evanescent edge is NARROW -- DEEPLY evanescent modes (n>>ka),")
print("where KELVIN-ONLY (no PML) is both more accurate AND parameter-free (the exact quasi-static")
print("compactification). Adding a PML to Kelvin for evanescent hurts. (Spherical>box PML is true but")
print("favours any spherical PML, not Kelvin specifically.)")
print("=" * 82)
