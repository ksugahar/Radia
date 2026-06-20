# -*- coding: utf-8 -*-
r"""
act7_17_dispersive_iabc.py  (Track A -- kelvin branch)
======================================================
TIME-DOMAIN high-frequency IABC by JOINTLY optimizing a PASSIVE, CAUSAL,
DISPERSIVE shell material (URN / Debye-relaxation form) together with the
shell THICKNESS -- and the answer to "do we need the Hankel continued fraction?"

This refines act7_16_wideband_iabc (constant passive shell, ~8.7% band reflection, "no
dispersion") and act7_19_highfreq_iabc_revisited (the exact-null per-omega optimum is non-passive).
The user's question (msg): can we jointly optimize {URN fit + layer thickness +
complex eps,mu} to DERIVE a time-domain IABC, and is the Hankel-function
continued-fraction expansion required?

SETUP (generic published Mie / transfer-matrix method -- no internal data read).
A single MATCHED spherical shell (eps = mu = m, so the planar impedance Z=1 for
every omega) of inner radius a=1, thickness d, PEC-backed, in vacuum.  Per mode n
the reflection of the outgoing spherical wave is computed by Riccati-Bessel
transfer matrices (same kernel as act7_19_highfreq_iabc_revisited/act7_16_wideband_iabc).  We MINIMIZE the band-max
|reflection| by weighted least squares (least_squares over a band of kR).

THE MATERIAL is parametrized in the URN passive-relaxation basis:

      m(w) = einf  +  sigma/(j w)  +  sum_k  dEps_k / (1 + j w tau_k),
      einf >= 1 real,  sigma >= 0,  dEps_k >= 0,  tau_k >= 0.

Every term is PASSIVE by construction (Im(m) <= 0 in the e^{+jwt} convention) and
CAUSAL -- each maps to a local-in-time auxiliary ODE (an ADE / circuit element):
  * einf   : instantaneous  D = einf E
  * sigma  : Ohmic conduction  dD/dt = sigma E              (1 local term)
  * Debye  : tau dP/dt + P = dEps E                          (1 ODE per pole)

PROVEN HERE (all asserted):
  (1) PASSIVITY: the fitted m(w) has Im(m(w)) <= 0 at every band sample
      (so it is realizable as a passive medium).
  (2) CAUSAL BEATS ACAUSAL: a passive +1-Debye shell reaches a LOWER band-max
      reflection than the best ACAUSAL "constant complex" shell (m = a - j b,
      constant Im => Kramers-Kronig-violating, NOT a real time-domain material).
      => dispersion is not a tax for going to the time domain; it HELPS.
  (3) THE RELAXATION POLE IS WHAT HELPS: a passive einf + sigma/(jw) shell
      (a conductor/Drude term, NO relaxation pole) is WORSE than the constant
      shell; adding ONE Debye pole flips it to clearly better.  So the win is the
      relaxation SHAPE -- exactly the URN basis.
  (4) TIME-DOMAIN REALIZATION: the Debye term's auxiliary ODE
      tau dP/dt + P = dEps E, integrated from rest under a steady tone E=e^{jwt},
      reproduces the frequency-domain susceptibility dEps/(1+jw tau) (rel.err
      small) and is STABLE (rate 1/tau > 0).  => the optimized shell is a
      ready-to-run ADE/FETD material; no recursive convolution, no instability.
  (5) HANKEL CONTINUED FRACTION -- needed?  The EXACT lossless exterior DtN is a
      FINITE rational function whose poles are the zeros of the spherical Hankel
      h_n^{(1)} (the reverse-Bessel-polynomial / Bessel-filter set, act6_10_iabc_time_domain); its
      partial fraction == the Hankel continued fraction and is reproduced here to
      ~machine.  That is the EXACT lossless-radiation reference.  The LOSSY
      wideband IABC above does NOT expand it: it fits a low-order passive
      relaxation model to the per-mode reflection instead.  So the Hankel CF is
      the reference target, NOT a required ingredient of the time-domain IABC.

HONEST SCOPE.  Single matched shell, mode n=1, e^{+jwt} convention.  The optimizer
is multi-start least_squares with FIXED seeds (deterministic); asserts use robust
margins, not the exact optimum value (the landscape is non-convex).  Prior art:
ADE/Debye-Lorentz dispersive PML & FDTD materials (Gedney, Taflove-Hagness);
graded/Chebyshev multilayer absorbers (Collin, Orfanidi); Grote-Keller exact
local NRBC (act6_10_iabc_time_domain).  New angle (knowledge module, not asserted as theorem): the
matched-shell IABC realized in the URN passive-relaxation basis beats the
constant passive shell AND is time-domain trivial.  Ports only published analytic
method; no internal optimization tables embedded.

No overclaim: every printed 'ok' is gated on an executed numerical assertion.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import hankel1, hankel2
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp

np.seterr(all="ignore")
np.set_printoptions(precision=6, suppress=False)
N_FAIL = 0
def check(name, cond, detail=""):
    global N_FAIL
    tag = "ok " if cond else "FAIL"
    if not cond: N_FAIL += 1
    print(f"  [{tag}] {name}{('  -- ' + detail) if detail else ''}")

# ===========================================================================
# (0) Riccati-Bessel matched-shell reflection  (generic Mie transfer matrix)
# ===========================================================================
def _sb(kind, n, z):
    z = complex(z); p = np.sqrt(np.pi / (2 * z))
    h1 = hankel1(n + 0.5, z); h2 = hankel2(n + 0.5, z)
    return p * (h1 + h2) / 2 if kind == 'j' else p * (h1 - h2) / (2j)
def _sbp(kind, n, z):
    z = complex(z); return _sb(kind, n - 1, z) - (n + 1) / z * _sb(kind, n, z)
def spb(n, z):
    z = complex(z); jn = _sb('j', n, z); yn = _sb('y', n, z)
    return z * jn, z * yn, jn + z * _sbp('j', n, z), yn + z * _sbp('y', n, z)

def refl(m, d, n, w):
    """|reflection| of outgoing mode-n wave from a matched (eps=mu=m) PEC-backed
    shell: inner radius a=1, outer radius 1+d.  Outgoing radiation [1;-j]."""
    r = np.array([1.0, 1.0 + d])
    mu = np.array([1.0, m]); ep = mu
    kk = w * np.sqrt(ep * mu)
    c1 = np.eye(2, dtype=complex)
    for i in range(len(r) - 1):
        bj, by, bjd, byd = spb(n, kk[i] * r[i])
        A1 = np.array([[mu[i] * bj, mu[i] * by], [kk[i] * bjd, kk[i] * byd]])
        bj, by, bjd, byd = spb(n, kk[i + 1] * r[i])
        B1 = np.array([[mu[i + 1] * bj, mu[i + 1] * by], [kk[i + 1] * bjd, kk[i + 1] * byd]])
        c1 = np.linalg.solve(B1, A1 @ c1)
    bj, by, bjd, byd = spb(n, kk[-1] * r[-1])
    A1 = np.array([[mu[-1] * bj, mu[-1] * by], [kk[-1] * bjd, kk[-1] * byd]])
    c1 = np.array([1, 0]) @ A1 @ c1
    return -(c1[0] - 1j * c1[1]) / (c1[0] + 1j * c1[1])

N_MODE = 1
def band_profile(mfun, d, band):
    out = []
    for w in band:
        try:
            v = abs(refl(mfun(w), d, N_MODE, w))
        except Exception:
            v = 10.0
        out.append(v if np.isfinite(v) else 10.0)
    return np.array(out)

# ===========================================================================
# (1) URN passive-relaxation material  +  acausal-constant benchmark
# ===========================================================================
def passive_material(params, ndeb):
    """m(w) = einf + sigma/(jw) + sum dEps_k/(1+jw tau_k).  params >=0 => Im<=0."""
    einf, sigma = params[0], params[1]
    deb = [(params[2 + 2 * k], params[3 + 2 * k]) for k in range(ndeb)]
    def m(w):
        val = complex(einf, 0.0) + sigma / (1j * w)
        for dE, ta in deb:
            val = val + dE / (1 + 1j * w * ta)
        return val
    return m

def fit_passive(band, ndeb, tries, seed):
    lb = [0.5, 0.0] + [0.0, 0.01] * ndeb + [0.02]
    ub = [20.0, 50.0] + [50.0, 50.0] * ndeb + [3.0]
    rng = np.random.default_rng(seed)
    best, best_bm = None, 9e9
    for _ in range(tries):
        x0 = rng.uniform([0.6, 0.0] + [0.2, 0.1] * ndeb + [0.3],
                         [3.0, 5.0] + [3.0, 5.0] * ndeb + [1.0])
        try:
            sol = least_squares(lambda p: band_profile(passive_material(p, ndeb), p[-1], band),
                                x0, bounds=(lb, ub), xtol=1e-12, ftol=1e-12, max_nfev=4000).x
        except Exception:
            continue
        bm = band_profile(passive_material(sol, ndeb), sol[-1], band).max()
        if bm < best_bm:
            best_bm, best = bm, sol
    return best_bm, best

def fit_acausal(band, tries, seed):
    """Constant complex m = a - j b (b free).  Constant Im => violates Kramers-
    Kronig => NOT a realizable time-domain material; the unphysical benchmark."""
    lb = [0.01, 0.0, 0.02]; ub = [20.0, 50.0, 3.0]
    rng = np.random.default_rng(seed)
    best, best_bm = None, 9e9
    for _ in range(tries):
        x0 = rng.uniform([0.3, 0.1, 0.3], [3.0, 3.0, 1.0])
        try:
            sol = least_squares(lambda p: band_profile(lambda w: complex(p[0], -p[1]), p[2], band),
                                x0, bounds=(lb, ub), xtol=1e-12, ftol=1e-12, max_nfev=4000).x
        except Exception:
            continue
        bm = band_profile(lambda w: complex(sol[0], -sol[1]), sol[2], band).max()
        if bm < best_bm:
            best_bm, best = bm, sol
    return best_bm, best

print("=" * 74)
print(" act7_17_dispersive_iabc: time-domain IABC by passive dispersive (URN/Debye) shell")
print("=" * 74)

BAND = np.linspace(2.0, 8.0, 25)   # primary band kR in [2,8] (4:1) -- robust
print("\n[A] Joint WLS optimization over band kR in [2,8], matched shell, n=1:")
bm_ac, p_ac = fit_acausal(BAND, tries=40, seed=0)
bm_d0, p_d0 = fit_passive(BAND, ndeb=0, tries=40, seed=10)   # einf + sigma/(jw)
bm_d1, p_d1 = fit_passive(BAND, ndeb=1, tries=60, seed=20)   # + 1 Debye pole
print(f"    ACAUSAL constant (a - j b)          band-max|ref| = {bm_ac:.4f}   (unrealizable floor)")
print(f"    PASSIVE einf + sigma/(jw)  [0 Debye] band-max|ref| = {bm_d0:.4f}   d={p_d0[-1]:.3f}")
print(f"    PASSIVE einf + sigma + 1 Debye       band-max|ref| = {bm_d1:.4f}   "
      f"einf={p_d1[0]:.3f} sigma={p_d1[1]:.3f} dEps={p_d1[2]:.3f} tau={p_d1[3]:.3f} d={p_d1[-1]:.3f}")

# ---- (1) passivity of the fitted dispersive material -----------------------
print("\n[1] Passivity of the fitted +1-Debye material over the band:")
m1 = passive_material(p_d1, 1)
ims = np.array([m1(w).imag for w in BAND])
res = np.array([m1(w).real for w in BAND])
check("Im(m(w)) <= 0 at every band sample (passive)", bool(np.all(ims <= 1e-9)),
      f"max Im = {ims.max():.2e}, min Im = {ims.min():.2e}")
check("Re(m(w)) > 0 over band", bool(np.all(res > 0)), f"min Re = {res.min():.3f}")

# ---- (2) causal (passive dispersive) beats acausal constant ----------------
print("\n[2] Causal passive +1-Debye vs acausal constant floor:")
check("passive +1-Debye band-max < acausal-constant band-max",
      bm_d1 < bm_ac, f"{bm_d1:.4f} < {bm_ac:.4f}  (ratio {bm_ac / bm_d1:.1f}x)")

# ---- (3) the relaxation pole is what helps ---------------------------------
print("\n[3] The relaxation pole is the active ingredient:")
check("einf+sigma/(jw) alone (no pole) is WORSE than constant",
      bm_d0 > bm_ac, f"0-Debye {bm_d0:.4f} > constant {bm_ac:.4f}")
check("adding ONE Debye pole flips it to clearly better than constant",
      bm_d1 < 0.6 * bm_ac, f"1-Debye {bm_d1:.4f} < 0.6 * {bm_ac:.4f}")

# ===========================================================================
# (4) TIME-DOMAIN realization: the Debye auxiliary ODE reproduces the
#     susceptibility, integrated from rest under a steady tone (causal, stable).
# ===========================================================================
print("\n[4] Time-domain ADE for the Debye term  tau dP/dt + P = dEps * E :")
dEps, tau = p_d1[2], p_d1[3]
w_test = float(BAND[len(BAND) // 2])          # a band frequency
chi_freq = dEps / (1 + 1j * w_test * tau)     # target susceptibility chi(w)
# integrate the auxiliary polarization ODE from rest with complex tone E=e^{jwt}
def rhs(t, P):                                  # P complex, split re/im for solver
    Pc = P[0] + 1j * P[1]
    E = np.exp(1j * w_test * t)
    dP = (dEps * E - Pc) / tau
    return [dP.real, dP.imag]
T = 60.0 * tau                                  # many relaxation times -> steady
sol = solve_ivp(rhs, [0.0, T], [0.0, 0.0], rtol=1e-10, atol=1e-12, dense_output=True)
Pend = sol.y[0, -1] + 1j * sol.y[1, -1]
chi_time = Pend / np.exp(1j * w_test * T)       # steady complex amplitude P/E
rel = abs(chi_time - chi_freq) / abs(chi_freq)
print(f"    chi(w) target     = {chi_freq:.6f}")
print(f"    chi from ADE      = {chi_time:.6f}   (transient from rest, {T/tau:.0f} tau)")
check("ADE steady tone reproduces Debye susceptibility", rel < 1e-3, f"rel.err = {rel:.2e}")
check("relaxation rate 1/tau > 0 (stable / causal)", (1.0 / tau) > 0, f"1/tau = {1.0/tau:.4f}")

# ===========================================================================
# (5) HANKEL CONTINUED FRACTION = exact lossless DtN reference (finite n poles)
# ===========================================================================
print("\n[5] Exact lossless exterior DtN = Hankel continued fraction (reference):")
# Lambda_n(z) = z h_n^{(1)'}(z)/h_n^{(1)}(z); for n=1 closed form = (z^2 - 1 - i z)/(... )
# Verify the partial-fraction (pole at zero of h_1^{(1)}) against direct scipy.
def dtn_lossless(n, z):
    bj = lambda nn: _sb('j', nn, z); by = lambda nn: _sb('y', nn, z)
    h = bj(n) + 1j * by(n)
    hp = (bj(n - 1) + 1j * by(n - 1)) - (n + 1) / z * h
    return z * hp / h
def dtn_n1_poleform(z):
    # h_1^{(1)} zero at z0 = i (theta_1(x)=x+1 -> root -1 -> z0 = i*(-1)*? ); use exact: -i.
    # Lambda_1(z) = i z - 1 + z0/(z - z0) with z0 = -i (the single zero of h_1^{(1)}).
    z0 = -1j
    return 1j * z - 1.0 + z0 / (z - z0)
zs = np.array([1.3 + 0.0j, 2.7, 4.1, 5.9, 8.3])
err = max(abs(dtn_lossless(1, z) - dtn_n1_poleform(z)) for z in zs)
print(f"    Lambda_1 partial-fraction (pole at z0=-i) vs scipy DtN: max err = {err:.2e}")
check("lossless DtN = finite Hankel-pole partial fraction (n=1)", err < 1e-10, f"err = {err:.2e}")
print("    => Hankel CF is the EXACT lossless reference (finite n poles, act6_10_iabc_time_domain);")
print("       the lossy wideband IABC above fits a passive relaxation model instead,")
print("       so the Hankel continued-fraction expansion is NOT required for it.")

# ---- secondary report: wider band [1,12] (no strict assert; non-convex) ----
print("\n[+] Secondary report -- wider band kR in [1,12] (no strict assert):")
B2 = np.linspace(1.0, 12.0, 25)
bm_ac2, _ = fit_acausal(B2, tries=40, seed=1)
bm_d12, p12 = fit_passive(B2, ndeb=1, tries=60, seed=30)
m12 = passive_material(p12, 1); ims2 = np.array([m12(w).imag for w in B2])
print(f"    acausal constant  band-max|ref| = {bm_ac2:.4f}")
print(f"    passive +1 Debye  band-max|ref| = {bm_d12:.4f}  passive(Im<=0)={bool(np.all(ims2<=1e-9))}")

print("\n" + "=" * 74)
if N_FAIL == 0:
    print(" ALL CHECKS PASSED -- passive dispersive (URN/Debye) shell is a")
    print(" time-domain-realizable wideband IABC that beats the constant shell;")
    print(" the Hankel continued fraction is the lossless reference, not required.")
else:
    print(f" {N_FAIL} CHECK(S) FAILED")
print("=" * 74)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
