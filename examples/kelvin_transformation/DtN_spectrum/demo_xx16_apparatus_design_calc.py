# -*- coding: utf-8 -*-
r"""
demo_xx16_apparatus_design_calc.py  (Kelvin mesh-adequacy: an ENGINEERING apparatus)
====================================================================================
The bridge from the canonical demos (sphere / square / L-shape) to the paper's target
apparatus (static-apparatus / rotating machinery).  A reactor / transformer LEG winding
cross-section is run through the full Kelvin DESIGN CALCULATION: a common enclosing circle
Gamma, the optimal radius R, and the required element order p -- all in closed form from
the winding geometry, NO convergence study (the abstract's promise, applied).

The winding is a 2-D balanced current set {(z_k, I_k)}, sum I_k = 0 (a real leg: go +I /
return -I conductors).  Its exterior A_z multipole content about the centroid is the set
of CURRENT MOMENTS a_m = sum_k I_k z_k^m (the m-th moment); on Gamma (radius R) the
relative content is |a_m| R^{-m}, decaying geometrically with the winding compactness
a_app/R.  Two design numbers follow:
    R*    : optimal enclosing radius (min total-DoF proxy ~ (R/a_app)^2 * p*(R)^2)
    p*    : required element order = the highest m with |a_m| R^{-m} above the target eps

VERIFIED HERE (asserted; self-contained numpy):
  [1] the apparatus multipole content (current moments) == an FFT of A_z on Gamma (machine
      precision) -- and is COMPACT: a balanced winding has no monopole, and the spectrum
      decays geometrically with a_app/R (the criterion's source side, demo_xx12).
  [2] the design calc: R* ~ 3 a_app (the "air box 2-5x" optimum, demo_e) and p* from the
      content; a realistic leg needs only a MODEST order (no fine air mesh).
  [3] the criterion holds for the apparatus: reconstructing the exterior field with m<=p*
      modes matches the full field to <= eps (so p* is the design order).

NON-CLAIM: a linear iron core only RESCALES the moments a_m (it does not change the
geometric a_app/R decay), so the design order is unchanged; the FEM end-to-end confirmation
of the same chain is demo_xx14 (a point source) -- here the point is the DESIGN CALCULATION
on a realistic multi-conductor apparatus.

No overclaim: every printed 'ok' is gated on an executed numerical assertion.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

# ---- a reactor / transformer LEG winding cross-section (metres, balanced) ----
# go column (+I) and return column (-I), 3 turns each -- a compact leg window.
_TURNS_Y = np.array([-0.04, 0.0, 0.04])
_XCOL = 0.05
SOURCES = ([(-_XCOL + 0j + 1j * yy, +1.0) for yy in _TURNS_Y]
           + [(+_XCOL + 0j + 1j * yy, -1.0) for yy in _TURNS_Y])     # (z_k, I_k), sum I_k = 0


def centroid(sources):
    Z = np.array([z for z, _ in sources]); W = np.array([abs(I) for _, I in sources])
    return complex(np.sum(W * Z) / np.sum(W))


def current_moments(sources, c, nmax=24):
    """a_m = sum_k I_k (z_k - c)^m  (the exterior A_z multipole content, m=1..nmax)."""
    return np.array([sum(I * (z - c) ** m for z, I in sources) for m in range(1, nmax + 1)])


def Az_on_circle(sources, c, R, nth=4096):
    """A_z = -(1/2pi) sum I_k ln|z-z_k| sampled on |z-c|=R (mu0/2pi=1)."""
    th = 2 * np.pi * np.arange(nth) / nth
    z = c + R * np.exp(1j * th)
    A = np.zeros(nth)
    for z0, I in sources:
        A += -I * np.log(np.abs(z - z0))
    return A


print("=" * 78)
print(" demo_xx16 : Kelvin DESIGN CALCULATION for a reactor/transformer leg winding")
print("=" * 78)

c = centroid(SOURCES)
a_app = max(abs(z - c) for z, _ in SOURCES)              # minimal enclosing radius
print(f"\n  winding: {len(SOURCES)} conductors, balanced (sum I = {sum(I for _,I in SOURCES):+.1f});"
      f" centroid={c:.4f}, enclosing a_app={a_app:.4f} m")

print(f"\n[1] apparatus multipole content (current moments) == FFT of A_z on Gamma; COMPACT:")
R0 = 3 * a_app
a_mom = current_moments(SOURCES, c)
A = Az_on_circle(SOURCES, c, R0)
c_fft = np.fft.rfft(A) / (len(A) / 2)                    # real FFT -> cos(m th) coeffs
# A_z multipole on the circle: coeff of cos(m th) is Re(a_m)/m * R^{-m} (+ sin from Im)
fft_mag = np.abs(c_fft[1:13])
mom_mag = np.array([abs(a_mom[m - 1]) / m * R0 ** (-m) for m in range(1, 13)])
rel = np.abs(fft_mag / fft_mag[0] - mom_mag / mom_mag[0])
print(f"    |content_m|/|content_1| (m=1..6), current-moment law: "
      + " ".join(f"{x:.2e}" for x in (mom_mag / mom_mag[0])[:6]))
print(f"    same from an FFT of A_z on Gamma:                     "
      + " ".join(f"{x:.2e}" for x in (fft_mag / fft_mag[0])[:6]))
print(f"    (max dev current-moment vs FFT, m<=12: {rel.max():.1e})")
assert rel.max() < 1e-3, "the current-moment content must match the FFT of A_z on Gamma"
assert abs(sum(I for _, I in SOURCES)) < 1e-12, "a balanced winding has NO monopole"
# the SYMMETRIC leg (x->-x, I->-I) excites only ODD multipoles m=1,3,5,... (a selection
# rule, like the magnetized square's n=4k+1 in demo_e); decay over the NON-ZERO modes:
norm = mom_mag / mom_mag[0]
nz = np.where(norm > 1e-12)[0]                            # non-zero (odd) mode indices
nz_mags = norm[nz]
print(f"    excited modes m = {[int(m+1) for m in nz[:5]]} (ODD only -- the leg's x<->-x symmetry rule)")
tail_ratio = float(np.mean(nz_mags[2:5] / nz_mags[1:4]))  # ratio over consecutive odd modes (tail)
print(f"    geometric decay over odd modes (tail ratio) = {tail_ratio:.3f}  (< 1 -> compact)")
assert norm[2] < 0.1 and nz_mags[-1] < norm[0], "the apparatus exterior content must be compact (decaying)"
assert tail_ratio < 0.95, "the odd-mode content must decay geometrically (compact source)"

print(f"\n[2] the design calc -- optimal enclosing R* (min total-DoF) and required order p*:")


def p_star(R, eps):
    a = np.array([abs(a_mom[m - 1]) / m * R ** (-m) for m in range(1, 25)])
    a /= a[0]
    above = np.where(a > eps)[0]
    return int(above[-1] + 1) if len(above) else 1       # highest m above eps


eps = 1e-4
print(f"    target eps={eps:.0e};  sweep enclosing R/a_app, total-DoF proxy ~ (R/a)^2 * p*^2:")
print(f"    {'R/a_app':>8} {'p*':>4} {'DoF proxy':>10}")
best = None
for ra in (1.5, 2.0, 2.78, 3.5, 5.0, 8.0):
    R = ra * a_app; ps = p_star(R, eps); dof = (ra) ** 2 * ps ** 2
    print(f"    {ra:8.2f} {ps:4d} {dof:10.1f}")
    if best is None or dof < best[2]:
        best = (ra, ps, dof)
print(f"    => optimal R*/a_app = {best[0]} (p*={best[1]}); the 'air box ~3x' rule, for this leg")
assert 2.0 <= best[0] <= 5.0, "the DoF-optimal enclosing radius must land in the air-box 2-5x window"
assert best[1] <= 8, "a compact leg needs only a MODEST element order (no fine air mesh)"

print(f"\n[3] the criterion holds for the apparatus: m<=p* reconstruction matches full field:")
Rdes = best[0] * a_app; ps = best[1]
# exterior field on a check circle just outside Gamma: full vs truncated at m<=p*
Rchk = 1.1 * Rdes
th = 2 * np.pi * np.arange(2048) / 2048
zc = c + Rchk * np.exp(1j * th)
full = np.zeros(2048)
for z0, I in SOURCES:
    full += -I * np.log(np.abs(zc - z0))
full -= full.mean()                                      # drop the (gauge) constant
# truncated reconstruction from m<=p* moments: A_z ~ Re[ sum_{m<=p*} (a_m/m) (z-c)^{-m} ]
trunc = np.zeros(2048)
for m in range(1, ps + 1):
    trunc += np.real((a_mom[m - 1] / m) * (zc - c) ** (-m))
trunc -= trunc.mean()
err = np.linalg.norm(trunc - full) / np.linalg.norm(full)
print(f"    R_design/a_app={best[0]}, p*={ps}: exterior-field reconstruction err (m<=p*) = {err:.2e}  (<= eps={eps:.0e})")
assert err <= eps, "the m<=p* Kelvin closure must reproduce the apparatus exterior field to eps"

print("\n[verdict]")
print("  A reactor/transformer leg run through the Kelvin DESIGN CALCULATION: its exterior")
print("  content is the winding CURRENT MOMENTS (== FFT of A_z on Gamma), a balanced winding")
print(f"  has no monopole and a geometrically compact spectrum, so a common enclosing circle at")
print(f"  R* ~ {best[0]} a_app with a MODEST order p*={ps} closes the open boundary to {eps:.0e} -- the")
print("  required R and p straight from the winding geometry, NO convergence study (a linear iron")
print("  core only rescales the moments).  The mesh-adequacy criterion (demo_xx12/13), in practice.")
print("\nALL CHECKS PASSED.")
