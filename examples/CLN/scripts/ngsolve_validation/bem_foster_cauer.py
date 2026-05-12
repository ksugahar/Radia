"""BEM-Foster Phase 4: CLN-I (s=0) ladder extraction.

From Foster pairs {(tau_k, a_k)}, compute:
  - Moments M_n = sum_k a_k tau_k^n for n=0..2K-1
  - Pade [K-1, K] approximant via Hankel system
  - Pade poles vs Foster tau_k (consistency check)
  - Hankel matrix condition number per stage K
  - CLN-I (s=0) continued fraction iteratively (Routh-style)
  - Track precision loss per stage

Goal: identify how many Cauer stages are reliable in float64 (16 digits).

Date: 2026-05-01
"""

import json
import time

import numpy as np
import scipy.linalg as la
from scipy.integrate import tplquad

# === Constants ===
sigma = 5.8e7
mu0 = 4.0 * np.pi * 1e-7
mu0sigma = mu0 * sigma
a, b, c = 5e-3, 2e-3, 1e-3
V = a * b * c

print("=" * 64)
print("BEM-Foster Phase 4: CLN-I (s=0) ladder extraction")
print("=" * 64)
print()


# === Cell-cell Newton kernel (reused) ===
def cell_newton_off(dx, dy, dz, hx, hy, hz, epsrel=1e-10):
    def integrand(w, v, u):
        return ((hx - abs(u)) * (hy - abs(v)) * (hz - abs(w)) /
                np.sqrt((u + dx)**2 + (v + dy)**2 + (w + dz)**2))
    val, _ = tplquad(integrand, -hx, hx, -hy, hy, -hz, hz,
                     epsabs=0.0, epsrel=epsrel)
    return val


def cell_newton_self(hx, hy, hz, epsrel=1e-10):
    def integrand(w, v, u):
        return ((hx - u) * (hy - v) * (hz - w) /
                np.sqrt(u*u + v*v + w*w))
    val, _ = tplquad(integrand, 0.0, hx, 0.0, hy, 0.0, hz,
                     epsabs=0.0, epsrel=epsrel)
    return 8.0 * val


def cell_newton(dx, dy, dz, hx, hy, hz):
    if dx == 0.0 and dy == 0.0 and dz == 0.0:
        return cell_newton_self(hx, hy, hz)
    return cell_newton_off(dx, dy, dz, hx, hy, hz)


def foster_solve(Nx, Ny, Nz):
    hx, hy, hz = a / Nx, b / Ny, c / Nz
    Vcell = hx * hy * hz
    ncells = Nx * Ny * Nz

    centers = np.array([
        [(i + 0.5) * hx - a/2,
         (j + 0.5) * hy - b/2,
         (k + 0.5) * hz - c/2]
        for i in range(Nx) for j in range(Ny) for k in range(Nz)
    ])

    cache = {}
    K = np.zeros((ncells, ncells))
    h_arr = np.array([hx, hy, hz])
    for i in range(ncells):
        for j in range(i, ncells):
            d = centers[i] - centers[j]
            key = tuple(int(round(abs(d[k] / h_arr[k]))) for k in range(3))
            if key not in cache:
                cache[key] = cell_newton(d[0], d[1], d[2], hx, hy, hz)
            K[i, j] = K[j, i] = cache[key]

    eigvals, eigvecs = la.eigh(K / Vcell)
    idx = np.argsort(-eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    phis = eigvecs / np.sqrt(Vcell)
    lambdaK = eigvals / (4.0 * np.pi)
    taus = mu0sigma * lambdaK

    xs = centers[:, 0]
    ys = centers[:, 1]
    gx = phis.T @ (-ys / 2.0) * Vcell
    gy = phis.T @ (+xs / 2.0) * Vcell

    # Foster amplitudes a_k for chi(s) = sum a_k / (1 + s tau_k)
    # alpha_n = (-1)^n sum a_k tau_k^n with chi(s) = alpha_0 + alpha_1 s + ...
    # We had: alpha_n = (-1)^n (mu0sigma)^n / V sum (gx^2 + gy^2) lambda^(n-1)
    # So: a_k tau_k^n contributions are unified:
    # comparing alpha_1 = -sum a_k tau_k = -(mu0sigma/V)*sum (gx^2+gy^2)
    # so sum a_k tau_k = (mu0sigma/V) sum (gx^2+gy^2)
    # and tau_k = mu0sigma lambda_k, so
    #   a_k = (gx_k^2 + gy_k^2)/(V lambda_k)
    a_ks = (gx**2 + gy**2) / (V * np.where(lambdaK > 0, lambdaK, np.inf))

    return {
        "taus": taus,
        "a_ks": a_ks,
        "lambdaK": lambdaK,
        "gx": gx,
        "gy": gy,
    }


# === Moments and Pade ===
def compute_moments(taus, a_ks, max_n):
    """M_n = sum_k a_k tau_k^n for n=0..max_n-1.

    Filter out modes with negative or near-zero lambda
    (numerical noise) to keep Stieltjes positivity.
    """
    valid = (taus > 0) & np.isfinite(a_ks)
    t = taus[valid]
    a = a_ks[valid]
    return np.array([np.sum(a * t**n) for n in range(max_n)])


def pade_via_hankel(M, K_stages):
    """Pade [K-1, K] approximant from moments M_0..M_{2K-1}.

    chi(s) = sum_{n=0}^inf (-1)^n M_n s^n = N(s)/D(s)
    with N degree K-1, D degree K, D(0)=1.

    Returns (num_coeffs, den_coeffs, hankel_cond)
    """
    n = K_stages
    # Taylor coefficients c_k = (-1)^k M_k
    c = np.array([(-1)**k * M[k] for k in range(2 * n)])

    # Solve for denominator: sum_{j=0..n} d_j c[i-j] = 0 for i=n..2n-1, d_0=1
    # i.e., A d = -c[n..2n-1] where A[i, j] = c[n+i-1-j], i,j=0..n-1
    A = np.array([[c[n + i - 1 - j] if n + i - 1 - j >= 0 else 0.0
                   for j in range(n)] for i in range(n)])
    rhs = -np.array([c[n + i] for i in range(n)])
    cond = np.linalg.cond(A)
    d_tail = np.linalg.solve(A, rhs)
    den = np.concatenate(([1.0], d_tail))

    # Numerator: a_i = sum_{j=0..min(i,n)} d[j] c[i-j], i=0..n-1
    num = np.array([
        sum(den[j] * c[i - j] for j in range(min(i, n) + 1))
        for i in range(n)
    ])

    return num, den, cond


def pade_poles_to_taus(den):
    """Roots of D(s) -> Pade time constants tau_k = -1/s_k."""
    roots = np.roots(den[::-1])  # numpy expects high->low order
    # Filter to real-negative poles (physical)
    real_neg = roots[(roots.real < 0) & (np.abs(roots.imag) < 1e-9 * np.abs(roots.real))]
    taus_pade = -1.0 / real_neg.real
    return np.sort(taus_pade)[::-1]  # descending


# === CLN-I (s=0) continued fraction extraction ===
def invert_taylor(c):
    """Compute Taylor series of 1/(c[0] + c[1] s + c[2] s^2 + ...).

    e[0] = 1/c[0]
    e[n] = -1/c[0] sum_{k=1..n} c[k] e[n-k]
    """
    n = len(c)
    if abs(c[0]) < 1e-300:
        raise ValueError("c[0] is zero")
    e = np.zeros(n)
    e[0] = 1.0 / c[0]
    for k in range(1, n):
        e[k] = -np.sum(c[1:k+1] * e[k-1::-1]) / c[0]
    return e


def cauer_continued_fraction(c_taylor, max_stages):
    """Extract CLN-I (s=0) continued fraction from Taylor series.

    chi(s) = c[0] + c[1] s + c[2] s^2 + ...
           = 1/(p[0] + s/(p[1] + s/(p[2] + ...)))

    Each iteration: invert series, peel off leading term, divide by s.

    Returns (p_list, residual_norms).
    residual_norms[k] = norm of remaining series after k stages
    (proxy for precision loss).
    """
    c = np.array(c_taylor, dtype=float)
    initial_norm = np.linalg.norm(c)
    p_list = []
    residuals = [initial_norm]

    for k in range(max_stages):
        if len(c) < 2:
            break
        if abs(c[0]) < 1e-300:
            break
        try:
            e = invert_taylor(c)
        except (ValueError, ZeroDivisionError):
            break
        p_k = e[0]
        p_list.append(p_k)
        # New chi: e[1] + e[2] s + e[3] s^2 + ...
        c = e[1:]
        residuals.append(np.linalg.norm(c))

    return np.array(p_list), np.array(residuals)


# === Main: compute Cauer ladder ===
print("Solving (8, 4, 2) Foster spectrum...")
t0 = time.time()
res = foster_solve(8, 4, 2)
print(f"  ({time.time()-t0:.1f} s)")
print()

# Filter valid Foster pairs (positive tau, non-trivial coupling)
valid = (res["taus"] > 0) & (res["a_ks"] > 1e-10 * res["a_ks"].max())
taus_v = res["taus"][valid]
a_v = res["a_ks"][valid]
order = np.argsort(-taus_v)
taus_v = taus_v[order]
a_v = a_v[order]

print(f"Coupled Foster modes: {len(taus_v)} (out of {len(res['taus'])})")
print(f"Top 5 (tau [us], a_k):")
for i in range(min(5, len(taus_v))):
    print(f"  tau_{i+1} = {taus_v[i]*1e6:8.4f} us,  a_{i+1} = {a_v[i]:.4e}")
print()

# Compute moments up to high order
N_moments = 16
M = compute_moments(res["taus"], res["a_ks"], N_moments)
print(f"Moments M_n (n=0..{N_moments-1}):")
for n in range(N_moments):
    print(f"  M_{n:2d} = {M[n]:+.10e}  (units: s^{n})")
print()

# === Pade poles vs Foster taus ===
print("=== Pade [K-1, K] poles vs Foster top-K tau_k ===")
print(f"{'K':<4s} {'Hankel cond':<14s} {'Pade tau_k [us]':<40s} {'vs Foster top-K [us]':<30s}")
print("-" * 100)
for K in range(1, 8):
    if 2 * K > len(M):
        break
    try:
        num, den, cond = pade_via_hankel(M, K)
        pade_taus = pade_poles_to_taus(den)
        foster_top = taus_v[:K] * 1e6
        pade_str = ", ".join(f"{t*1e6:.4f}" for t in pade_taus[:K])
        foster_str = ", ".join(f"{t:.4f}" for t in foster_top)
        print(f"  {K:<2d}  {cond:<14.4e}  [{pade_str}]   vs  [{foster_str}]")
    except Exception as e:
        print(f"  {K:<2d}  FAILED: {e}")
print()

# === CLN-I (s=0) continued fraction ===
print("=== CLN-I (s=0) continued fraction from Taylor series ===")
# Taylor coefficients c_k = (-1)^k M_k
c_taylor = np.array([(-1)**k * M[k] for k in range(N_moments)])

p_cauer, residuals = cauer_continued_fraction(c_taylor, max_stages=N_moments - 1)

print(f"{'k':<4s} {'p_k':<22s} {'residual norm':<18s} {'log10(res ratio)':<18s}")
print("-" * 70)
for k, p in enumerate(p_cauer):
    res_after = residuals[k+1] if k+1 < len(residuals) else 0
    res_before = residuals[k]
    ratio = np.log10(res_after/res_before) if res_after > 0 and res_before > 0 else None
    ratio_str = f"{ratio:+.4f}" if ratio is not None else "N/A"
    print(f"  {k+1:<2d}  {p:<22.10e}  {res_after:<18.6e}  {ratio_str}")
print()

# Convert continued fraction p_k to physical R, L (alternating)
# CLN-I (s=0): chi(s) = 1/(p_1 + s/(p_2 + s/(p_3 + ...)))
# Physical interpretation depends on context; here we report raw p_k.

# Check first stage matches Pade [0,1]
print(f"First Cauer p_1 = {p_cauer[0]:.6e}  (= 1/chi(0) = 1/M_0)")
print(f"M_0 = {M[0]:.6e}, so 1/M_0 = {1.0/M[0]:.6e}")
print()

# Extract effective time constants from Cauer truncation
# Truncate at K stages, build polynomial, find poles
def cauer_to_poles(p_list, K):
    """Truncate Cauer at K stages and find poles of resulting rational function."""
    # chi^(K)(s) = 1/(p_1 + s/(p_2 + s/(p_3 + ... + s/p_K)))
    # Build from inside out as N(s)/D(s)
    if K < 1:
        return np.array([])
    # Innermost: chi_K = 1/p_K (just the last p)
    N = np.array([1.0])
    D = np.array([p_list[K-1]])
    for k in range(K-2, -1, -1):
        # Outer: chi_k(s) = 1/(p_k + s * chi_{k+1}(s))
        # = 1/(p_k + s * N/D) = D/(p_k D + s N)
        # New N = D, new D = p_k D + s N
        N_new = D.copy()
        D_new_p = p_list[k] * D
        # Pad N to multiply by s
        sN = np.concatenate(([0.0], N))
        # Pad D_new_p to same length
        max_len = max(len(D_new_p), len(sN))
        D_new = np.zeros(max_len)
        D_new[:len(D_new_p)] += D_new_p
        D_new[:len(sN)] += sN
        N, D = N_new, D_new
    return N, D


print("=== Cauer truncation -> Pade poles ===")
print(f"{'K stages':<10s} {'recovered tau_k [us]':<60s}")
for K in range(1, min(len(p_cauer), 8) + 1):
    try:
        Nc, Dc = cauer_to_poles(p_cauer, K)
        roots = np.roots(Dc[::-1])
        # Get real-negative roots, convert to positive tau
        real_neg = roots[(roots.real < 0) & (np.abs(roots.imag) < 1e-6 * np.abs(roots.real))]
        cauer_taus = sorted(-1.0/real_neg.real, reverse=True)
        tau_str = ", ".join(f"{t*1e6:.4f}" for t in cauer_taus)
        print(f"  {K:<8d}  [{tau_str}]")
    except Exception as e:
        print(f"  {K:<8d}  FAILED: {e}")

# Save
with open("bem_foster_cauer_results.json", "w") as f:
    json.dump({
        "grid": [8, 4, 2],
        "moments": M.tolist(),
        "p_cauer": p_cauer.tolist(),
        "residuals": residuals.tolist(),
        "foster_taus_top10": taus_v[:10].tolist(),
        "foster_a_top10": a_v[:10].tolist(),
    }, f, indent=2)
print()
print("Results saved to bem_foster_cauer_results.json")
