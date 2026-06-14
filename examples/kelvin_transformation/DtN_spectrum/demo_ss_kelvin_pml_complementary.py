# -*- coding: utf-8 -*-
# DEMO (ss) (verified): Kelvin and PML are COMPLEMENTARY, not redundant -- and does a "fancier" PML
# help the combination? A genuine combined Kelvin+PML radial FE (inverted-exterior transformation-optics
# medium PLUS a complex-stretch PML in the near-centre = the image of the far field).
#
# The user's question: "or, if you use a more sophisticated PML (there are several PML types), does
# combining it with Kelvin give a benefit?" Tested directly:
#   * Kelvin's compactified centre (image of r=infinity) is where the inverted far-field oscillation
#     piles up (k_eff = k a^2/rho^2 -> infinity). If you merely EXCISE it and put a hard wall (no
#     absorber), the wall REFLECTS that oscillation -> the truncation DtN is catastrophically wrong.
#     So the PML is NOT optional in the Kelvin exterior: Kelvin REQUIRES an absorber there. They are
#     COMPLEMENTARY (this is the IEICE-2024 design: a PML at the centre of the inverted exterior).
#   * A modest STANDARD (polynomial) PML in that near-centre region makes the truncation DtN ~exact.
#     There is a SWEET SPOT in strength: too weak = no absorption; too strong = the steep near-centre
#     stretch out-resolves the mesh and re-reflects (the classic discretized-PML failure).
#   * A "fancier" profile (graded to MATCH the diverging k_eff, ~1/rho^2) does NOT help here -- it
#     over-concentrates the stretch where the mesh is coarsest and degrades faster than a simple
#     polynomial. So for the Kelvin centre the key is adequate strength + near-centre RESOLUTION, not a
#     sophisticated sigma profile; a simple PML is the robust choice. (CFS-PML's real shift targets
#     EVANESCENT waves, but the Kelvin near-centre field is oscillatory, not evanescent, so CFS adds
#     little here either.)
#
# Combined FE: image shell [rho_e, a]; the PML stretches INWARD from rho_pml toward the excision rho_e
# (so [rho_pml, a] stays physical, rho~=rho); hard wall R(rho_e)=0; truncation Dirichlet R(a)=1;
# DtN = -(A R)|_a / a (the inversion sign-flip). s(rho)=1+i*alpha*((rho_pml-rho)/L)^p in the PML region.
#
# VERIFIED (a=1, ka=4):
#  (1) EXCISION + HARD WALL (no PML) is catastrophic: |DtN-exact| ~ O(1..40) (cavity resonances) -> a
#      PML/absorber is REQUIRED in the Kelvin exterior (Kelvin and PML are complementary).
#  (2) a modest standard PML makes it ~exact (alpha~2 -> ~1e-4) across modes n; there is a SWEET SPOT
#      (alpha too large -> the under-resolved steep stretch re-reflects, error grows again).
#  (3) a k_eff-matched / over-graded profile does NOT beat the simple polynomial (it degrades sooner):
#      for the Kelvin centre, strength + resolution matter, not a fancy profile.
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


def kelvin_pml(n, k, a, rho_e, rho_pml, M, alpha, p=2, extra_grade=False):
    """combined Kelvin (transformation-optics medium) + PML (complex stretch in [rho_e, rho_pml]).
    extra_grade=True multiplies the stretch by (rho_pml/rho)^2 (a k_eff-matched, more-concentrated
    profile) -- used to show that the fancier profile does NOT help."""
    L = rho_pml - rho_e

    def s(r):
        if r >= rho_pml:
            return 1.0 + 0j
        g = ((rho_pml - r) / L) ** p
        if extra_grade:
            g = g * (rho_pml / max(r, 1e-9)) ** 2
        return 1 + 1j * alpha * g

    def rt(r):                                  # rho~: physical for r>=rho_pml, complex-stretched inward
        if r >= rho_pml:
            return complex(r)
        if extra_grade:                         # integrate s inward numerically (graded profile)
            tt = np.linspace(r, rho_pml, 32)
            return rho_pml - np.trapezoid(np.array([s(t) for t in tt]), tt)
        dr = rho_pml - r                        # analytic for the plain polynomial profile
        return rho_pml - (dr + 1j * alpha / L**p * dr**(p + 1) / (p + 1))

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
    u = np.zeros(M + 1, complex); u[M] = 1.0; u[0] = 0.0       # wall R(rho_e)=0 ; truncation R(a)=1
    idx = list(range(1, M))
    u[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -A[np.ix_(idx, [M])][:, 0])
    return -(A[M, :] @ u) / a


# ============================================================================= VERIFICATIONS
a, k = 1.0, 4.0
ka = k * a
M = 400
print("=" * 82)
print("Kelvin & PML are COMPLEMENTARY; does a fancier PML help the combination? (a=1, ka=%.0f)" % ka)
print("=" * 82)

err = 0.0
for n in range(0, 4):
    for z in (0.5, 2.0):
        rj = spherical_jn(n, z); err = max(err, abs(_sph('j', n, z).item() - rj) / (abs(rj) + 1e-300))
assert err < 1e-11

# (1) excision + hard wall (no PML) is catastrophic -> a PML is REQUIRED -----------------------
print("\n(1) EXCISION + hard wall, NO PML (n=1): the inverted far-field reflects -> DtN catastrophically wrong:")
print("     rho_e   physical wall a^2/rho_e   |DtN - exact|")
bad = 0.0
for rho_e in (0.5, 0.3, 0.2):
    d = kelvin_pml(1, k, a, rho_e, rho_e + 1e-9, M, 0.0)
    bad = max(bad, abs(d - dtn_exact(1, ka)))
    print("     %.2f    %5.1f                    %.3e" % (rho_e, a * a / rho_e, abs(d - dtn_exact(1, ka))))
assert bad > 1.0, "excision + hard wall must be catastrophic (PML required)"
print("    => Kelvin's compactified centre MUST be absorbed; excision alone fails. Kelvin NEEDS a PML.")

# (2) a modest standard PML -> ~exact; there is a SWEET SPOT in strength -----------------------
print("\n(2) standard polynomial PML in [0.3,0.6] (n=1): a SWEET SPOT in strength alpha:")
print("     alpha   |DtN - exact|")
e_best = 1e9; e_strong = 0.0
for alpha in (0.5, 1.0, 2.0, 4.0, 8.0, 20.0):
    e = abs(kelvin_pml(1, k, a, 0.3, 0.6, M, alpha) - dtn_exact(1, ka))
    e_best = min(e_best, e)
    if alpha == 20.0:
        e_strong = e
    print("     %4.1f    %.3e" % (alpha, e))
assert e_best < 1e-3, "a well-tuned PML makes the Kelvin truncation DtN ~exact"
assert e_strong > 10 * e_best, "too-strong PML re-reflects (under-resolved steep stretch) = the sweet spot"
print("    => a modest PML (alpha~2) -> ~exact; too weak (no absorption) or too strong (under-resolved")
print("       steep stretch re-reflects) are both worse -> the PML must be TUNED (complementary design).")

# multi-mode check at the sweet spot
print("\n    sweet-spot PML (alpha=2) across modes n:  |DtN-exact| =",
      ", ".join("%d:%.1e" % (n, abs(kelvin_pml(n, k, a, 0.3, 0.6, M, 2.0) - dtn_exact(n, ka))) for n in (1, 2, 3)))

# (3) a fancier (k_eff-matched, over-graded) profile does NOT help --------------------------
print("\n(3) does a FANCIER profile help? plain polynomial vs k_eff-matched (~1/rho^2) grading (n=1):")
print("     alpha   plain poly     k_eff-matched")
worse = False
for alpha in (1.0, 2.0, 4.0, 8.0):
    ep = abs(kelvin_pml(1, k, a, 0.3, 0.6, M, alpha, extra_grade=False) - dtn_exact(1, ka))
    eg = abs(kelvin_pml(1, k, a, 0.3, 0.6, M, alpha, extra_grade=True) - dtn_exact(1, ka))
    if alpha >= 4.0 and eg > ep:
        worse = True
    print("     %4.1f    %.3e     %.3e" % (alpha, ep, eg))
assert worse, "the k_eff-matched profile over-concentrates the stretch -> degrades sooner than plain"
print("    => the fancier (k_eff-matched) profile is NOT better: it over-stretches the under-resolved")
print("       centre and degrades sooner. For the Kelvin centre, STRENGTH + RESOLUTION matter, not a")
print("       sophisticated sigma profile -- a simple PML is the robust choice.")

print("\n" + "=" * 82)
print("ANSWER: Kelvin and PML are COMPLEMENTARY, NOT redundant -- excising Kelvin's centre and walling")
print("it off is catastrophic; the PML is REQUIRED to absorb the inverted far-field there (the IEICE-2024")
print("design). But a SOPHISTICATED PML adds little to the COMBINATION: a simple polynomial PML, tuned to")
print("the sweet spot, already makes the truncation DtN ~exact; fancier (k_eff-matched) or stronger")
print("profiles over-stress the under-resolved centre and re-reflect. So 'a better PML' helps a PML")
print("generally, but the Kelvin-specific need is just adequate strength + near-centre resolution.")
print("=" * 82)
