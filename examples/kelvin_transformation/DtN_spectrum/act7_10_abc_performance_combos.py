# -*- coding: utf-8 -*-
# DEMO (qq) (verified): evaluating ABC PERFORMANCE via the DtN -- radiation BC vs PML, and what the
# KELVIN placement does for each. Answers: "the radiation BC performs poorly (known); how do Kelvin+PML
# and Kelvin+radiation-BC compare; does effectively-farther placement help?"
#
# The DtN spectrum is the performance yardstick for ANY absorbing/open boundary: a boundary's quality =
# how well its truncation DtN matches the exact Lambda_n(ka) across the multipole band. Two ingredients:
#   ORDER of the termination (how many multipoles its impedance matches) and PLACEMENT (how far out).
#   * radiation BC (1st-order ABC / Leontovich SIBC, d_n u = ik u): a FIXED-ORDER impedance -- matches
#     only the leading DtN term; its per-mode error GROWS with n.
#   * PML (complex stretch + wall): an ALL-MULTIPOLE absorber -- its termination impedance matches
#     Lambda_n(kb) for EVERY n (it absorbs the outgoing wave regardless of mode).
#   * Kelvin sets the PLACEMENT: it maps the whole exterior into the ball, so an absorber at the image
#     centre sits effectively at b -> infinity, cheaply (a finite image shell).
#
# VERIFIED FINDINGS (a=1, ka=4; the honest, slightly non-obvious answer):
#  (1) TERMINATION ORDER at b: the radiation BC (SIBC) matches Lambda_n(kb) only for low n (error grows
#      0.12->3.05 over n=1..6) = FIXED-order; the PML matches all n (~6e-4 flat) = ALL-order.
#  (2) ABC PERFORMANCE (truncation DtN error vs Lambda_n(ka), per mode):
#        radiation BC AT the truncation (b=a) -> POOR (0.06..1.1, grows with n) = the known fact;
#        Kelvin+radiation-BC (placed far) -> much better than radBC@a, but a FIXED-order FLOOR remains;
#        Kelvin+HOIBC -> better (2nd order); Kelvin+PML -> ~exact (all-order) = the paper's method,
#        on par with a well-resolved plain PML.
#  (3) WHAT DISTANCE DOES (the key nuance): placing the absorber FARTHER (larger b, what Kelvin buys)
#        HELPS a radiation BC a lot (error ~1/kb: 0.15->0.009 over b=1.5..6, since the higher multipoles
#        have decayed so a fixed-order impedance suffices) but barely changes a PML (~flat 1e-4: an
#        all-order absorber already works at any distance).
#  => the user's reasoning, calibrated: (a) radiation BC alone is poor (fixed-order, at b=a); (b)
#     Kelvin+radiation-BC IS limited -- a fixed-order BC needs far placement, which Kelvin provides but
#     not for free (near-centre resolution, act7_03_hoibc_kelvin_transform), and a floor remains -> worse than Kelvin+PML;
#     (c) Kelvin+PML is the best (all-order, effectively at infinity), but its gain over a plain PML is
#     QUALITATIVE (infinity baked in -> no truncation-distance choice, AND it carries exterior
#     scatterers, the IEICE-2024 use case) rather than a large vacuum-accuracy gain (a PML's all-order
#     absorption already makes placement distance nearly irrelevant).
#
# Pure numpy/scipy (radial FEs / closed-form propagation; the assemblies are act7_07_threeway_kelvin_pml_bemfem's).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import jv, yv, hankel1, hankel2, spherical_jn, spherical_yn


def _sph(kind, n, z):
    z = np.asarray(z, dtype=complex)
    pref = np.sqrt(np.pi / (2.0 * z))
    return pref * {'j': jv, 'y': yv, 'h1': hankel1, 'h2': hankel2}[kind](n + 0.5, z)


def _sph_prime(kind, n, z):
    z = np.asarray(z, dtype=complex)
    return _sph(kind, n - 1, z) - (n + 1.0) / z * _sph(kind, n, z)


def dtn_exact(n, z):
    z = complex(z); return (z * _sph_prime('h1', n, z) / _sph('h1', n, z)).item()


def dtn_trunc(n, ka, kb, inner):
    """truncation DtN at r=a given an inner impedance 'inner' (=b f'/f) imposed at radius b."""
    h1b, h1pb = _sph('h1', n, kb).item(), _sph_prime('h1', n, kb).item()
    h2b, h2pb = _sph('h2', n, kb).item(), _sph_prime('h2', n, kb).item()
    beta = (inner * h1b - kb * h1pb) / (kb * h2pb - inner * h2b)
    h1a, h1pa = _sph('h1', n, ka).item(), _sph_prime('h1', n, ka).item()
    h2a, h2pa = _sph('h2', n, ka).item(), _sph_prime('h2', n, ka).item()
    return ka * (h1pa + beta * h2pa) / (h1a + beta * h2a)


_GX = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)]); _GW = np.array([5 / 9, 8 / 9, 5 / 9])


def pml_impedance(n, k, b, d, M, s0):
    """the impedance a PML layer [b, b+d] presents at its inner radius b (= the DtN there)."""
    def s(r): return 1 + 1j * s0 * (r - b)**2 / (k * d * d)
    def rt(r): return r + 1j * s0 * (r - b)**3 / (3 * k * d * d)
    nod = np.linspace(b, b + d, M + 1); A = np.zeros((M + 1, M + 1), complex)
    for e in range(M):
        r0, r1 = nod[e], nod[e + 1]; h = r1 - r0; dp = np.array([-1 / h, 1 / h])
        K = np.zeros((2, 2), complex); C = np.zeros((2, 2), complex); Mm = np.zeros((2, 2), complex)
        for gp, gw in zip(_GX, _GW):
            r = .5 * (r0 + r1) + .5 * h * gp; w = .5 * h * gw
            ph = np.array([(r1 - r) / h, (r - r0) / h]); sr = s(r); rr = rt(r)
            K += np.outer(dp, dp) * (rr * rr / sr) * w
            C += np.outer(ph, ph) * sr * w; Mm += np.outer(ph, ph) * sr * rr * rr * w
        A[e:e + 2, e:e + 2] += K + n * (n + 1) * C - k * k * Mm
    u = np.zeros(M + 1, complex); u[0] = 1.0; idx = list(range(1, M))
    u[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -A[np.ix_(idx, [0])][:, 0])
    return -(A[0, :] @ u) / b


def sibc(kb):
    return 1j * kb - 1.0


def hoibc(n, kb):
    return 1j * kb - 1.0 - 1j * n * (n + 1.0) / (2.0 * kb)


# ============================================================================= VERIFICATIONS
a, k = 1.0, 4.0
ka = k * a
b, d, M, s0 = 2.0, 1.0, 200, 15.0
kb = k * b
print("=" * 84)
print("ABC PERFORMANCE via the DtN: radiation BC vs PML, and Kelvin's placement (a=1, ka=%.0f, b=%.0f)" % (ka, b))
print("=" * 84)

err = 0.0
for n in range(0, 6):
    for z in (0.5, 2.0, 6.0):
        rj = spherical_jn(n, z)
        err = max(err, abs(_sph('j', n, z).item() - rj) / (abs(rj) + 1e-300))
assert err < 1e-11
print("\n(1) TERMINATION ORDER at b=%.0f (kb=%.0f): match to exact Lambda_n(kb) across n:" % (b, kb))
print("     n   |radiationBC - Lam(kb)|   |PML - Lam(kb)|")
e_si_n, e_pml_n = [], []
for n in (1, 2, 3, 4, 6):
    Lb = dtn_exact(n, kb)
    es = abs(sibc(kb) - Lb); ep = abs(pml_impedance(n, k, b, d, M, s0) - Lb)
    e_si_n.append(es); e_pml_n.append(ep)
    print("    %2d   %.3e                %.3e" % (n, es, ep))
assert e_si_n[-1] > 10 * e_si_n[0], "radiation BC is FIXED-order: error grows with n"
assert max(e_pml_n) < 5e-3, "PML is ALL-order: matches Lambda_n(kb) for every n"
print("    => radiation BC = FIXED-order (matches only low n); PML = ALL-order (matches every n).")

print("\n(2) ABC PERFORMANCE: truncation DtN error vs exact Lambda_n(ka=%.0f), per mode:" % ka)
print("     n   radBC@a    Kelvin+SIBC   Kelvin+HOIBC   Kelvin+PML    plain PML")
worst_rad = 0.0
for n in (1, 2, 3, 4, 6):
    L = dtn_exact(n, ka)
    e_rad = abs(sibc(ka) - L) / abs(L)                                   # radiation BC AT the truncation
    e_ksi = abs(dtn_trunc(n, ka, kb, sibc(kb)) - L) / abs(L)             # Kelvin + radiation BC at b
    e_kho = abs(dtn_trunc(n, ka, kb, hoibc(n, kb)) - L) / abs(L)         # Kelvin + HOIBC at b
    e_kpml = abs(dtn_trunc(n, ka, kb, pml_impedance(n, k, b, d, M, s0)) - L) / abs(L)  # Kelvin + PML at b
    e_ppml = abs(pml_impedance(n, k, a, d, M, s0) - L) / abs(L)          # plain PML at the truncation
    worst_rad = max(worst_rad, e_rad)
    print("    %2d   %.2e   %.2e      %.2e       %.2e      %.2e" % (n, e_rad, e_ksi, e_kho, e_kpml, e_ppml))
    assert e_ksi < e_rad, "Kelvin's far placement helps the radiation BC (vs radBC at the truncation)"
    assert e_kpml < e_kho < e_ksi, "ordering: Kelvin+PML < Kelvin+HOIBC < Kelvin+SIBC"
assert worst_rad > 0.3, "the bare 1st-order radiation BC performs poorly (the known fact)"
print("    => radBC@a POOR (known); Kelvin helps it but a fixed-order FLOOR remains; Kelvin+PML ~exact")
print("       (all-order), on par with a well-resolved plain PML.")

print("\n(3) WHAT DISTANCE DOES (n=3): farther placement helps a radiation BC, NOT a PML:")
print("     b    kb    Kelvin+radBC(SIBC)   Kelvin+PML")
n = 3; L = dtn_exact(n, ka); e_si_b, e_pml_b = [], []
for bb in (1.5, 2.0, 3.0, 4.0, 6.0):
    kbb = k * bb
    es = abs(dtn_trunc(n, ka, kbb, sibc(kbb)) - L) / abs(L)
    ep = abs(dtn_trunc(n, ka, kbb, pml_impedance(n, k, bb, d, M, s0)) - L) / abs(L)
    e_si_b.append(es); e_pml_b.append(ep)
    print("    %.1f  %4.1f   %.3e            %.3e" % (bb, kbb, es, ep))
assert e_si_b[0] > 5 * e_si_b[-1], "radiation BC error must DROP with distance (farther helps)"
assert max(e_pml_b) / min(e_pml_b) < 2.0, "PML error ~flat with distance (all-order; distance irrelevant)"
print("    => radiation BC error ~1/kb (farther helps -> Kelvin's placement genuinely improves it);")
print("       PML error ~flat (all-order absorbs at any distance -> Kelvin gives PML no accuracy gain).")

print("\n" + "=" * 84)
print("ANSWER: the DtN spectrum evaluates every ABC. A radiation BC is a FIXED-ORDER impedance (poor at")
print("the truncation = the known fact); a PML is an ALL-ORDER absorber. Kelvin sets the PLACEMENT.")
print("Kelvin+radiation-BC: placement helps (error ~1/kb) but a fixed-order FLOOR remains and far")
print("placement costs near-centre resolution -> NOT as good as Kelvin+PML. Kelvin+PML: all-order +")
print("effectively at infinity -> ~exact; its edge over a plain PML is QUALITATIVE (infinity baked in,")
print("no truncation distance, and it carries exterior scatterers) -- a PML's all-order absorption")
print("already makes distance nearly irrelevant for vacuum accuracy.")
print("=" * 84)
