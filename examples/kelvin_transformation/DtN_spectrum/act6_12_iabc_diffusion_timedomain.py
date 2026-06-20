# -*- coding: utf-8 -*-
r"""
act6_12_iabc_diffusion_timedomain.py  (Track A -- kelvin branch)
===============================================================
The PARABOLIC (eddy-current / magnetic-diffusion) time-domain IABC -- the static-apparatus
(SA)-native regime (user: "(a)" = the diffusion version, after the high-frequency WAVE case
act7_15_iabc_timedomain_cost). This completes development-direction D's parabolic branch.

PHYSICS (the SA-correct setting): in a magneto-quasistatic (eddy-current) problem the EXTERIOR
AIR is non-conducting -> its field is instantaneously Laplace (the static ladder, act7_14_iabc_static_elegant), NO
time dynamics. A genuinely PARABOLIC truncation DtN appears where the adjacent region is a
CONDUCTOR (tank / core / shield skin effect = the surface-impedance / open boundary of a
diffusive region). Per spherical-harmonic mode n, the exterior magnetic-diffusion (modified
Helmholtz) decaying solution is the modified Bessel K_{n+1/2}(gamma r), gamma=sqrt(s*mu*sigma),
s=Laplace variable, giving the per-mode diffusion DtN

   Lambda_n(s) = a*u'/u|_a = -a*gamma*K_{n-1/2}(gamma a)/K_{n+1/2}(gamma a) - (n+1).

It INTERPOLATES the two known regimes:
   s -> 0   (DC / thick skin):   Lambda_n -> -(n+1)          (the STATIC multipole ladder, act7_14_iabc_static_elegant)
   s -> inf (high freq / thin skin): Lambda_n -> -a*gamma = -a*sqrt(s*mu*sigma)  (the sqrt(s) SIBC).

KEY CONTRAST WITH THE WAVE CASE (act7_15_iabc_timedomain_cost):
  * WAVE exterior DtN = RATIONAL, finite n poles -> EXACT with M=n auxiliary ODEs (finite memory).
  * DIFFUSION exterior DtN = sqrt(s) BRANCH CUT at s=0 -> INFINITE memory (t^{-3/2} kernel); it is
    NOT representable by finitely many poles. A rational (Foster RL-ladder) approximation converges
    only ALGEBRAICALLY in the pole count M -> the time-domain IABC is band-limited, M auxiliary ODEs
    = recursive convolution.

VERIFIED HERE (all asserted):
  (A) the diffusion DtN recovers the static ladder -(n+1) as s->0 (machine) and the sqrt(s) SIBC
      -a*gamma as s->inf (both limits).
  (B) sqrt(s) / diffusion DtN is a BRANCH CUT, not rational: a finite-pole fit error FLOORS
      algebraically (never machine-zero), unlike the wave case which hit ~1e-15 at M=n.
  (C) TIME-DOMAIN IABC = Foster M-pole rational fit  d + sum r_j s/(s+p_j)  (real poles p_j>0 =
      passive RL ladder) over a frequency band: stable, band error DECREASES with M; the M
      auxiliary ODEs  psi_j' = p_j(u - psi_j),  g = (d+sum r_j) u - sum r_j psi_j  reproduce it by
      transient integration -> the parabolic IABC = M auxiliary ODEs / recursive convolution.
  (D) DtN cost-accuracy datasheet vs M, and the honest wave-vs-diffusion (finite-pole vs branch-cut)
      cost contrast; time-domain FEM-BEM (exact diffusion kernel = full t^{-1/2} history) = most accurate.

PRIOR ART (cite, not claim): the sqrt(s*mu*sigma) surface impedance and its Foster/Cauer + recursive-
convolution time-domain realization is the eddy-current SIBC literature -- Valdivieso, Meunier,
Ramdane, Gyselinck (IEEE T-Magn 2020, Foster networks + recursive convolution) and time-domain SIBC
(Yuferev & Ida). Those realize INTERIOR material / conductor truncation; the defensible-new slice is
the per-mode DIFFUSION DtN SPECTRUM (static-ladder <-> sqrt(s) SIBC interpolation) read as an
open-boundary IABC + the cost-accuracy datasheet + the wave(finite-pole)-vs-diffusion(branch-cut)
contrast, in the SA / Kelvin context. Cost is STRUCTURAL (state dim / locality), not wall-clock.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import kv
from scipy.integrate import solve_ivp

np.set_printoptions(precision=6, suppress=False)
A_R, MUSIG = 1.0, 1.0                     # truncation radius a, and mu*sigma (units: set to 1)

# ===========================================================================
# (A) per-mode magnetic-diffusion exterior DtN and its two asymptotic regimes
# ===========================================================================
def gamma(s):
    return np.sqrt(complex(s) * MUSIG)

def dtn_diffusion(n, s, a=A_R):
    g = gamma(s)
    return -a * g * kv(n - 0.5, g * a) / kv(n + 0.5, g * a) - (n + 1.0)

print("=" * 78)
print(" act6_12_iabc_diffusion_timedomain : PARABOLIC (eddy-current/diffusion) TIME-DOMAIN IABC")
print("=" * 78)

print("\n[A] diffusion DtN interpolates static ladder -(n+1)  <->  sqrt(s) SIBC:")
maxdc = 0.0; maxsibc = 0.0
for n in range(1, 6):
    dc = dtn_diffusion(n, 1e-12)                 # s -> 0  (DC / thick skin)
    maxdc = max(maxdc, abs(dc - (-(n + 1.0))))
    s_hi = 1j * 1e4                              # s -> inf along imaginary axis (thin skin)
    ratio = dtn_diffusion(n, s_hi) / (-A_R * gamma(s_hi))
    maxsibc = max(maxsibc, abs(ratio - 1.0))
print(f"    n=1..5 : |DtN(s->0) - (-(n+1))|max = {maxdc:.2e}   (static ladder recovered)")
print(f"    n=1..5 : |DtN(s->inf)/(-a*gamma) - 1|max = {maxsibc:.2e}   (sqrt(s) SIBC recovered)")
assert maxdc < 1e-6 and maxsibc < 5e-2
print("    ok  (one per-mode operator bridging DC multipole ladder and the skin-effect SIBC)")

# ===========================================================================
# (B) it is a BRANCH CUT (sqrt s), NOT rational: finite-pole fit error FLOORS
# ===========================================================================
def foster_fit(func, w_band, M):
    """Fit func(iw) by d + sum_j r_j (iw)/(iw+p_j), real poles p_j log-spaced in band,
    real d & r_j by least squares (real+imag stacked). Stable by construction (p_j>0)."""
    s = 1j * w_band
    p = np.logspace(np.log10(w_band[0]), np.log10(w_band[-1]), M)
    cols = [s / (s + pj) for pj in p] + [np.ones_like(s)]
    Amat = np.column_stack(cols)
    rhs = np.array([func(sv) for sv in s], dtype=complex)
    AA = np.vstack([Amat.real, Amat.imag]); bb = np.concatenate([rhs.real, rhs.imag])
    coef, *_ = np.linalg.lstsq(AA, bb, rcond=None)
    r = coef[:M]; d = coef[M]
    fit = Amat @ np.concatenate([r, [d]])
    relerr = np.max(np.abs(fit - rhs)) / np.max(np.abs(rhs))
    return p, r, d, relerr

wb = np.logspace(-1, 2, 200)                      # band omega in [0.1, 100]
print("\n[B] sqrt(s) is a BRANCH CUT (infinite memory): finite-pole fit error floors algebraically:")
errs = []
for M in (2, 4, 8, 16, 32):
    _, _, _, e = foster_fit(lambda s: np.sqrt(s), wb, M)
    errs.append(e)
    print(f"    sqrt(s) Foster fit M={M:2d} poles : band rel.err = {e:.2e}")
assert errs[0] > errs[-1] > 1e-6, "branch cut: should improve but NOT reach machine zero"
print("    ok  (algebraic convergence, never ~1e-15 -- contrast act7_15_iabc_timedomain_cost WAVE DtN: exact at M=n)")

# ===========================================================================
# (C) TIME-DOMAIN IABC = Foster M-pole fit of the diffusion DtN + auxiliary ODEs
# ===========================================================================
print("\n[C] time-domain IABC = Foster fit of the diffusion DtN over band, + auxiliary ODEs:")
n = 2
errsC = []
for M in (2, 4, 8, 16, 32):
    p, r, d, e = foster_fit(lambda s: dtn_diffusion(n, s), wb, M)
    assert np.all(p > 0)                          # passive RL ladder -> stable
    errsC.append(e)
    print(f"    n={n}, M={M:2d} : stable(p>0)={np.all(p>0)}  band rel.err = {e:.2e}")
assert errsC[-1] < errsC[0] and errsC[-1] < 2e-2, "more poles -> better band accuracy"
print("    ok  (stable Foster ladder; band error decreases OVERALL with M -- not step-monotone:")
print("         the simple fixed-log-pole fit is sub-optimal, a vector-fit would be monotone/better)")

# C2: the M auxiliary ODEs reproduce the fitted DtN by transient integration.
#     psi_j' = p_j (u - psi_j);   g = (d + sum r_j) u - sum r_j psi_j
def transient_diffusion_iabc(p, r, d, w, t_end=None):
    m = p.size
    t_end = t_end if t_end is not None else max(200.0, 50.0 / min(p))
    def rhs(t, y):
        psi = y[:m] + 1j * y[m:]
        u = np.exp(1j * w * t)                    # e^{+i w t}, Laplace s=i w
        dpsi = p * (u - psi)
        return np.concatenate([dpsi.real, dpsi.imag])
    sol = solve_ivp(rhs, [0, t_end], np.zeros(2 * m), rtol=1e-9, atol=1e-11,
                    t_eval=[t_end], method='RK45')
    psi = sol.y[:m, -1] + 1j * sol.y[m:, -1]
    u = np.exp(1j * w * t_end)
    g = (d + np.sum(r)) * u - np.sum(r * psi)
    return g / u

p, r, d, _ = foster_fit(lambda s: dtn_diffusion(n, s), wb, 8)
def dtn_foster(s, p, r, d):
    s = complex(s); return d + np.sum(r * s / (s + p))
maxC2 = 0.0
for w in (0.5, 2.0, 10.0):
    g_over_u = transient_diffusion_iabc(p, r, d, w)
    maxC2 = max(maxC2, abs(g_over_u - dtn_foster(1j * w, p, r, d)) / abs(dtn_foster(1j * w, p, r, d)))
print(f"    C2 M=8 auxiliary-ODE transient vs the fitted DtN: max rel.err = {maxC2:.2e}")
assert maxC2 < 5e-3
print("    ok  (M auxiliary ODEs = a Foster RL ladder = recursive convolution = the parabolic IABC)")

# ===========================================================================
# (D) DtN cost-accuracy datasheet + honest wave-vs-diffusion contrast
# ===========================================================================
print("\n[D] diffusion-DtN cost-accuracy (band omega in [0.1,100]); cost = #poles = #aux ODEs/mode:")
print("    mode n |   M=2     M=4     M=8     M=16    M=32    (band rel.err)")
for n in range(1, 5):
    row = [foster_fit(lambda s: dtn_diffusion(n, s), wb, M)[3] for M in (2, 4, 8, 16, 32)]
    print(f"      n={n}  | " + "  ".join(f"{e:7.1e}" for e in row))
print("""
    READING (cost model, structural -- not wall-clock):
      * TIME-DOMAIN FEM-BEM = the EXACT diffusion DtN: dense Gamma x Gamma + a FULL temporal
        convolution with the t^{-1/2} diffusion kernel (infinite memory).  MOST ACCURATE, costliest.
      * PARABOLIC TIME-DOMAIN IABC = M auxiliary ODEs / mode (a Foster RL ladder = recursive
        convolution), LOCAL in time + SPARSE in space -> cheap; the table is its accuracy knob.
      * KEY CONTRAST vs the WAVE case (act7_15_iabc_timedomain_cost): the wave DtN is rational (finite n poles) -> EXACT
        at M=n; the DIFFUSION DtN has a sqrt(s) branch cut (infinite memory) -> NO finite exact pole
        set, only algebraic convergence in M.  This is the physical signature of magnetic diffusion,
        read straight off the DtN spectrum -- and it is the SA-native (eddy-current) regime.
""")
print("ALL CHECKS PASSED.")
