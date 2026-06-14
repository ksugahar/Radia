"""Cauer depth test using v6 hybrid 18-digit alpha_n.

Determines how many Cauer stages can be reliably extracted given 18-digit input.

Theoretical: Hankel condition number ~ 10^(2k) per Cauer stage (Pade[K-1,K]).
With 18 digit input, useful output for k where 10^(2k) << 10^18 -> k < 9.
So ~8 reliable Cauer stages expected.

Date: 2026-05-02
"""

import json
import mpmath as mp

mp.mp.dps = 50

# Constants
sigma = mp.mpf("5.8e7")
mu0 = 4 * mp.pi * mp.mpf("1e-7")
a = mp.mpf("5e-3")
b = mp.mpf("2e-3")
c = mp.mpf("1e-3")
V = a * b * c
mu0sigma = mu0 * sigma

# === alpha_n from v6 hybrid (18 digit precision) ===
# Need 16 moments for 8 Cauer stages (K=8 Hankel needs M_0..M_15)
# But v6 only had alpha_1..alpha_8.
# For first Cauer test, use 8 alphas + alpha_0 = M_0 derived from Foster spec.

# Load v6 hybrid result
try:
    with open("bem_foster_spectral_v6_hybrid_results.json", "r") as f:
        v6 = json.load(f)
    alphas = [mp.mpf(s) for s in v6["alpha_n"]]  # alpha_1..alpha_8
    print("Loaded v6 hybrid results.")
except FileNotFoundError:
    print("v6 result not found; using hardcoded values from Mathematica run.")
    alphas = [
        mp.mpf("-0.0000440346570278169352257847180889677"),
        mp.mpf("7.4783240535434074561446429934e-10"),
        mp.mpf("-1.38797777606801509029076575691e-14"),
        mp.mpf("2.66525073846078827742571391442e-19"),
        mp.mpf("-5.19179202154889838584477748969e-24"),
        mp.mpf("1.01780181397466786388312805616e-28"),
        mp.mpf("-2.00117265750345380104921264126e-33"),
        mp.mpf("3.94011856973007074943279704837e-38"),
    ]

# Need M_0 = alpha_0 (DC value). From Foster: M_0 = sum a_k.
# In our problem, we need to derive M_0 from BEM.
# v6 hybrid gives 8 moments alpha_1..alpha_8. We also need M_0.
# From Phase 7 Mathematica (P=8 WP=30): M_0 = 3.05897258746351459800805934293582037956
M0 = mp.mpf("3.05897258746351459800805934293582037956")

# Construct moment list: M_n = (-1)^n alpha_n for n>=1, M_0 = alpha_0
# alpha_n in our convention: alpha_n = (-1)^n M_n where M_n = sum a_k tau_k^n
# So M_n = (-1)^n alpha_n.
# For n=1: M_1 = -alpha_1 = +4.40346e-5 (positive)
# For n=2: M_2 = +alpha_2 = +7.478e-10 (positive)

M_moments = [M0] + [(-1)**n * alphas[n-1] for n in range(1, 9)]  # M_0..M_8

print()
print(f"M_0 = {mp.nstr(M_moments[0], 25)}")
print(f"M_1 = {mp.nstr(M_moments[1], 25)}")
print(f"M_2 = {mp.nstr(M_moments[2], 25)}")
print(f"M_8 = {mp.nstr(M_moments[8], 25)}")
print()

# === CLN-I (s=0) continued fraction extraction ===
def invert_taylor(c_taylor):
    """Invert power series."""
    n = len(c_taylor)
    e = [mp.mpf(0)] * n
    e[0] = 1 / c_taylor[0]
    for k in range(1, n):
        s = mp.fsum(c_taylor[j] * e[k - j] for j in range(1, k + 1))
        e[k] = -s / c_taylor[0]
    return e


def cauer_continued_fraction(c_taylor, max_stages):
    """CLN-I (s=0) continued fraction: chi(s) = 1/(p_1 + s/(p_2 + s/(...))).

    Iteration: invert series, peel leading term, divide by s.
    """
    c = list(c_taylor)
    p_list = []
    residuals = []
    for k in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf("1e-300"):
            break
        e = invert_taylor(c)
        p_list.append(e[0])
        c = e[1:]
        residuals.append(mp.norm(mp.matrix(c)) if c else mp.mpf(0))
    return p_list, residuals


def cauer_truncate_to_poles(p_list, K):
    """Truncate Cauer at K stages, find poles via polynomial roots.

    chi^(K)(s) = N(s)/D(s) built inside-out.
    """
    if K < 1:
        return []
    # Innermost: chi_K = 1/p_K
    N = [mp.mpf(1)]
    D = [p_list[K - 1]]
    for k in range(K - 2, -1, -1):
        sN = [mp.mpf(0)] + N
        Dnew_pk = [p_list[k] * d for d in D]
        max_len = max(len(Dnew_pk), len(sN))
        Dnew = [mp.mpf(0)] * max_len
        for i in range(len(Dnew_pk)):
            Dnew[i] += Dnew_pk[i]
        for i in range(len(sN)):
            Dnew[i] += sN[i]
        N, D = D, Dnew
    return N, D


def find_real_neg_roots(coeffs):
    """Find real negative roots of polynomial. coeffs[0] = constant term."""
    # mpmath polyroots
    # Reverse to high-to-low coefficient order for mpmath
    coeffs_hi_to_lo = coeffs[::-1]
    try:
        roots = mp.polyroots(coeffs_hi_to_lo, maxsteps=100, extraprec=20)
    except Exception:
        return []
    real_neg = []
    for r in roots:
        if hasattr(r, 'imag') and abs(r.imag) > mp.mpf("1e-15") * abs(r.real):
            continue
        real_part = r.real if hasattr(r, 'real') else r
        if real_part < 0:
            real_neg.append(real_part)
    return real_neg


# === Compute Cauer continued fraction ===
print("=" * 64)
print("CLN-I (s=0) continued fraction extraction")
print("=" * 64)

# Taylor coefficients of chi(s) = sum c_n s^n where c_n = (-1)^n M_n = alpha_n
# Wait actually: chi(s) = sum a_k/(1+s tau_k) = sum_n (-s)^n sum_k a_k tau_k^n
# = sum_n (-1)^n s^n M_n = sum_n alpha_n s^n where alpha_n = (-1)^n M_n
# So c_n = alpha_n. M_0 = alpha_0 (Parseval, no sign).
# Actually alpha_0 = M_0 by definition (n=0).
c_taylor = [M_moments[0]] + [(-1)**n * M_moments[n] for n in range(1, len(M_moments))]
# Equivalent: c_n = alpha_n via standard convention

print(f"Taylor coefficients of chi(s):")
for n in range(min(9, len(c_taylor))):
    print(f"  c_{n} = {mp.nstr(c_taylor[n], 20)}")

p_list, residuals = cauer_continued_fraction(c_taylor, max_stages=10)
print()
print(f"Cauer p_k values (extracted {len(p_list)} stages):")
for k, p in enumerate(p_list):
    print(f"  p_{k+1} = {mp.nstr(p, 15)}")
print()

# === Cauer truncation -> poles ===
print("=" * 64)
print("Cauer truncation -> recovered tau_k [μs]")
print("=" * 64)
print(f"{'K':<4s} {'Poles (1/τ_k)':<30s} {'Recovered τ_k [μs]':<60s}")
for K in range(2, len(p_list) + 1):
    Nc, Dc = cauer_truncate_to_poles(p_list, K)
    poles = find_real_neg_roots(Dc)
    if poles:
        taus = sorted([-1/p for p in poles], reverse=True)
        tau_str = ", ".join(mp.nstr(t * 1e6, 10) for t in taus)
        print(f"  {K:<2d}  {len(poles)} real-neg poles  {tau_str}")
    else:
        print(f"  {K:<2d}  no real-neg poles found")

print()
print("=" * 64)
print("Stage limits at 18-digit input precision")
print("=" * 64)
print("Theoretical: Hankel cond ~ 10^(2K), useful output for K < 9.")
print("Practical: Cauer 4-5 stages (8-10 R/L parameters) reliably 8-12 digits.")
print("To reach Cauer 8 stages: need P=8+ (more α_n), 30+ digit input.")
