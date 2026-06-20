# -*- coding: utf-8 -*-
# DEMO (rr) (verified): what does Kelvin+PML actually buy? -- NOT "far placement" (a PML is distance-
# insensitive, act7_10_abc_performance_combos); the genuine merit is carrying EXTERIOR SCATTERERS / MATERIAL.
#
# The user's sharp question: "if it's only about being able to place the PML farther, does Kelvin+PML
# have much merit?" Honest answer, two parts:
#   PART 1 (VACUUM): NO -- for a compact source in vacuum, Kelvin+PML and a plain PML right at the
#     truncation give the same DtN. A PML is an all-order absorber whose accuracy is nearly distance-
#     independent (act7_10_abc_performance_combos), so Kelvin's effective infinite placement buys almost nothing in vacuum;
#     a plain PML is simpler (no inverted medium, no centre singularity/excision). The skepticism is
#     correct for vacuum.
#   PART 2 (EXTERIOR SCATTERER): YES, decisively -- the real merit. When scatterers / material / a
#     ground plane sit OUTSIDE the truncation (the IEICE-2024 antenna-with-mountains use case), the
#     Kelvin map conformally compactifies the WHOLE exterior -- INCLUDING those objects -- into a
#     bounded mesh of FIXED cost (independent of how far the objects are). A plain PML is a VACUUM
#     absorber: it cannot contain a scatterer; to include the objects you must extend the INTERIOR mesh
#     to enclose them (cost grows with their distance). So Kelvin+PML reproduces the correct (shifted)
#     exterior DtN that a vacuum PML simply misses.
#
# VERIFIED (a=1, k0a=4; scalar Helmholtz; scatterer shell physical [c1,c2]=[1.5,2], m=sqrt(eps_r)=2):
#  (1) VACUUM (m=1): the Kelvin radial FE (+PML-quality inner) and a plain PML both reproduce the exact
#      vacuum DtN Lambda_n(k0 a) to ~1e-4 -> Kelvin+PML ~ plain PML (no vacuum merit beyond placement).
#  (2) SCATTERER (m=2): the exterior shell SHIFTS the exact truncation DtN by O(1) (|scat - vacuum| ~
#      1.4..2.5, comparable to the DtN magnitude); the Kelvin FE carrying the scatterer in the inverted
#      shell REPRODUCES the exact with-scatterer DtN to ~1e-4, while a plain (vacuum) PML gives the
#      vacuum DtN -- WRONG by the full O(1) shift (it cannot represent the scatterer at all).
#
# Pure numpy/scipy. The exact with-scatterer DtN is a 3-region radial transfer (vacuum / shell /
# outgoing); the Kelvin FE is act7_05_fe_kelvin_hoibc's with the scatterer's wavenumber mapped into the inverted shell.
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


def dtn_vacuum(n, ka):
    z = complex(ka); return (z * _sph_prime('h1', n, z) / _sph('h1', n, z)).item()


def dtn_scatterer(n, k0, a, c1, c2, m):
    """exact truncation DtN at r=a with an exterior shell [c1,c2] of wavenumber k1=m*k0 (vacuum else,
    outgoing beyond c2): 3-region radial matching (u, u' continuous)."""
    k1 = m * k0

    def uv(kind, k, r):
        return _sph(kind, n, k * r).item(), k * _sph_prime(kind, n, k * r).item()
    u3, up3 = uv('h1', k0, c2)                                   # outgoing region III at c2
    jc2, jpc2 = uv('j', k1, c2); yc2, ypc2 = uv('y', k1, c2)
    A2, B2 = np.linalg.solve(np.array([[jc2, yc2], [jpc2, ypc2]]), np.array([u3, up3]))
    jc1, jpc1 = uv('j', k1, c1); yc1, ypc1 = uv('y', k1, c1)
    u2, up2 = A2 * jc1 + B2 * yc1, A2 * jpc1 + B2 * ypc1
    Jc1, Jpc1 = uv('j', k0, c1); Yc1, Ypc1 = uv('y', k0, c1)
    A1, B1 = np.linalg.solve(np.array([[Jc1, Yc1], [Jpc1, Ypc1]]), np.array([u2, up2]))
    Ja, Jpa = uv('j', k0, a); Ya, Ypa = uv('y', k0, a)
    return a * (A1 * Jpa + B1 * Ypa) / (A1 * Ja + B1 * Ya)


_GX = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)]); _GW = np.array([5 / 9, 8 / 9, 5 / 9])


def fe_kelvin_scat(n, k0, a, b, M, inner, c1=None, c2=None, m=1.0):
    """Kelvin radial FE; an exterior shell [c1,c2] enters as a wavenumber boost m in its IMAGE region
    [a^2/c2, a^2/c1] of the inverted shell (scalar transformation optics)."""
    rb = a * a / b; rho = np.linspace(rb, a, M + 1); A = np.zeros((M + 1, M + 1), complex)
    rlo, rhi = (a * a / c2, a * a / c1) if c1 else (None, None)
    for e in range(M):
        r0, r1 = rho[e], rho[e + 1]; h = r1 - r0; dp = np.array([-1 / h, 1 / h])
        K = np.outer(dp, dp) * (a * a * h); C = np.zeros((2, 2), complex); Mm = np.zeros((2, 2), complex)
        for gp, gw in zip(_GX, _GW):
            r = .5 * (r0 + r1) + .5 * h * gp; w = .5 * h * gw
            ph = np.array([(r1 - r) / h, (r - r0) / h])
            mm = m if (rlo is not None and rlo <= r <= rhi) else 1.0
            C += np.outer(ph, ph) * (a * a / r**2) * w
            Mm += np.outer(ph, ph) * (mm * mm) * (a**6 / r**4) * w
        A[e:e + 2, e:e + 2] += K + n * (n + 1) * C - k0 * k0 * Mm
    A[0, 0] += -b * inner; u = np.zeros(M + 1, complex); u[M] = 1.0
    u[:M] = np.linalg.solve(A[:M, :M], -A[:M, M]); return -(A[M, :] @ u) / a


def fe_pml_vac(n, k, a, d, M, s0):
    """plain vacuum PML at the truncation [a, a+d]; DtN at a (it can only ever give the VACUUM DtN)."""
    def s(r): return 1 + 1j * s0 * (r - a)**2 / (k * d * d)
    def rt(r): return r + 1j * s0 * (r - a)**3 / (3 * k * d * d)
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
    u = np.zeros(M + 1, complex); u[0] = 1.0; idx = list(range(1, M))
    u[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -A[np.ix_(idx, [0])][:, 0])
    return -(A[0, :] @ u) / a


# ============================================================================= VERIFICATIONS
a, k0 = 1.0, 4.0
ka = k0 * a
b, M, d, s0 = 3.0, 400, 1.0, 15.0
kb = k0 * b
c1, c2, m = 1.5, 2.0, 2.0    # exterior scatterer shell, eps_r=4 -> m=2
print("=" * 82)
print("What does Kelvin+PML buy? NOT far placement -- it carries EXTERIOR SCATTERERS (a=1, k0a=%.0f)" % ka)
print("=" * 82)

err = 0.0
for n in range(0, 5):
    for z in (0.5, 2.0, 6.0):
        rj = spherical_jn(n, z); err = max(err, abs(_sph('j', n, z).item() - rj) / (abs(rj) + 1e-300))
assert err < 1e-11

# (1) VACUUM: Kelvin+PML ~ plain PML -> far placement buys ~nothing in vacuum -------------------
print("\n(1) VACUUM (no scatterer): Kelvin+PML vs plain PML, error vs exact vacuum DtN Lambda_n(k0 a):")
print("     n   Kelvin+PML err   plain PML err")
for n in (1, 2, 3, 4):
    L = dtn_vacuum(n, ka)
    ek = abs(fe_kelvin_scat(n, k0, a, b, M, dtn_vacuum(n, kb)) - L) / abs(L)   # no scatterer (m=1 default)
    ep = abs(fe_pml_vac(n, k0, a, d, M, s0) - L) / abs(L)
    print("    %2d   %.3e        %.3e" % (n, ek, ep))
    assert ek < 1e-3 and ep < 1e-3, "in vacuum both reproduce the exact DtN -> Kelvin gives no edge"
print("    => in VACUUM Kelvin+PML ~ plain PML (both ~exact): far placement buys ~nothing; a plain PML")
print("       is simpler (no inverted medium / centre excision). The skepticism is correct for vacuum.")

# (2) EXTERIOR SCATTERER: Kelvin carries it; a vacuum PML misses it entirely --------------------
print("\n(2) EXTERIOR SCATTERER shell [%.1f,%.1f], eps_r=%.0f (m=%.0f): the genuine merit:" % (c1, c2, m * m, m))
print("     n   exact w/scat        vacuum DtN          shift |scat-vac|   Kelvin+PML(carries scat) err   plainPML err vs truth")
for n in (1, 2, 3, 4):
    Lscat = dtn_scatterer(n, k0, a, c1, c2, m)
    Lvac = dtn_vacuum(n, ka)
    shift = abs(Lscat - Lvac)
    ek = abs(fe_kelvin_scat(n, k0, a, b, M, dtn_vacuum(n, kb), c1, c2, m) - Lscat) / abs(Lscat)
    ep = abs(fe_pml_vac(n, k0, a, d, M, s0) - Lscat) / abs(Lscat)             # vacuum PML vs the TRUE DtN
    print("    %2d   (%6.3f,%6.3f)   (%6.3f,%6.3f)   %.3f             %.2e                  %.2e"
          % (n, Lscat.real, Lscat.imag, Lvac.real, Lvac.imag, shift, ek, ep))
    assert shift > 1.0, "an exterior scatterer shifts the DtN by O(1)"
    assert ek < 1e-3, "Kelvin+PML carries the exterior scatterer (reproduces the true DtN)"
    assert ep > 0.2, "a vacuum PML misses the scatterer (large error vs the true DtN)"
print("    => the scatterer shifts the DtN by O(1); Kelvin+PML (scatterer mapped into the inverted")
print("       shell) reproduces it to ~1e-4; a plain VACUUM PML gives the vacuum DtN -> WRONG by the")
print("       full shift. A PML cannot contain a scatterer; Kelvin compactifies the exterior + objects.")

print("\n" + "=" * 82)
print("ANSWER: 'far placement' alone is NOT the merit of Kelvin+PML -- a PML is distance-insensitive,")
print("so in VACUUM Kelvin+PML ~ a (simpler) plain PML. The genuine merit is that Kelvin conformally")
print("compactifies the WHOLE EXTERIOR, INCLUDING scatterers / material / ground, into a bounded fixed-")
print("cost mesh: it reproduces the correct exterior DtN (verified ~1e-4 with an eps_r=4 shell that")
print("shifts the DtN by O(1)), which a vacuum PML fundamentally cannot represent. That -- the IEICE-2024")
print("antenna-with-surrounding-structures use case -- is why Kelvin+PML exists, not the placement.")
print("=" * 82)
