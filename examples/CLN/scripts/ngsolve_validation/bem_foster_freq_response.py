"""BEM-Foster Phase 5: frequency response comparison.

Compute chi(j omega) on a frequency sweep via:
  - Foster sum (full BEM eigendecomposition - "ground truth")
  - Cauer K-stage truncation (K = 2, 3, 4, 5, 6, 8)
  - Single-pole closed-form Pade[1,1] (closed alpha_1, alpha_2)

Plots: magnitude, phase, real/imag, relative error vs Foster.

Goal: identify the frequency range where each Cauer truncation is
physically usable. K=3 should be exact at DC, K=2 already represents
single-pole response.

Date: 2026-05-01
"""

import json
import time

import numpy as np
import scipy.linalg as la
from scipy.integrate import tplquad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# === Constants ===
sigma = 5.8e7
mu0 = 4.0 * np.pi * 1e-7
mu0sigma = mu0 * sigma
a, b, c = 5e-3, 2e-3, 1e-3
V = a * b * c

alpha1_closed = -mu0sigma * (a*a + b*b) / 48.0
alpha2_closed = 7.47832405353e-10
tau1_closed = alpha2_closed / (-alpha1_closed)
M0_closed = -alpha1_closed**2 / alpha2_closed  # DC value from Pade[1,1]

print("=" * 64)
print("BEM-Foster Phase 5: frequency response chi(j omega)")
print("=" * 64)
print()


# === Reused: cell newton, foster_solve ===
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

    a_ks = (gx**2 + gy**2) / (V * np.where(lambdaK > 0, lambdaK, np.inf))

    return {"taus": taus, "a_ks": a_ks, "lambdaK": lambdaK}


# === Cauer continued fraction (reused from Phase 4) ===
def invert_taylor(c):
    n = len(c)
    e = np.zeros(n)
    e[0] = 1.0 / c[0]
    for k in range(1, n):
        e[k] = -np.sum(c[1:k+1] * e[k-1::-1]) / c[0]
    return e


def cauer_continued_fraction(c_taylor, max_stages):
    c = np.array(c_taylor, dtype=float)
    p_list = []
    for k in range(max_stages):
        if len(c) < 2 or abs(c[0]) < 1e-300:
            break
        try:
            e = invert_taylor(c)
        except Exception:
            break
        p_list.append(e[0])
        c = e[1:]
    return np.array(p_list)


# === Frequency response evaluators ===
def chi_foster(omegas, taus, a_ks):
    """chi(j omega) = sum a_k / (1 + j omega tau_k) over coupled modes."""
    s = 1j * omegas
    # Vectorized over omega and modes
    chi = np.zeros(len(omegas), dtype=complex)
    for k in range(len(taus)):
        chi += a_ks[k] / (1.0 + s * taus[k])
    return chi


def chi_cauer(omegas, p_cauer, K):
    """chi(j omega) from K-stage Cauer truncation:
       chi = 1/(p_1 + s/(p_2 + s/(p_3 + ... + s/p_K)))
    Build from inside out.
    """
    s = 1j * omegas
    # Innermost: chi_inner = 1/p_K  (a constant, since chi_{K+1} = 0)
    # Actually if we truncate at K: chi^(K) = 1/(p_1 + s/(p_2 + ... + s/p_K))
    # Build from inside.
    inner = 1.0 / p_cauer[K-1]   # this is chi_K (just one term)
    # Wait, that's wrong. Let me redo.
    # chi^(K) = 1/(p_1 + s/(p_2 + s/(... + s/(p_{K-1} + s/p_K))))
    # Innermost level = p_{K-1} + s/p_K (treating p_K as terminating)
    # Actually for cleanest truncation: just stop at p_K with chi_{K+1}=0.
    #   chi^(K) such that the K-th level is 1/p_K (treating remainder as 0).
    # But that gives chi^(K) only when K is interpreted as "truncate after K levels".

    # Implementation: start with x = 1/p_K (innermost), then for k = K-2 down to 0:
    #    x = 1/(p_k + s*x_{prev})... no wait.
    # The continued fraction is:
    #   chi(s) = 1/(p_1 + s/(p_2 + s/(p_3 + s/(... + s/p_K))))
    # = 1/(p_1 + s * y_1) where y_1 = 1/(p_2 + s * y_2), etc.
    # Innermost y_{K-1} = 1/p_K.
    # Then y_{K-2} = 1/(p_{K-1} + s * y_{K-1})
    #     y_{K-3} = 1/(p_{K-2} + s * y_{K-2})
    #     ...
    #     y_0 = 1/(p_1 + s * y_1) = chi(s)
    # So loop: y = 1/p_K, then for k = K-1, K-2, ..., 1:
    #   y = 1/(p_k + s*y)
    if K < 1:
        return np.zeros_like(s)
    y = np.full_like(s, 1.0 / p_cauer[K-1])
    for k in range(K-2, -1, -1):
        y = 1.0 / (p_cauer[k] + s * y)
    return y


def chi_pade11(omegas, alpha1, alpha2):
    """Single-pole Pade [1,1]: chi(s) ~ M_0 / (1 + s*tau_1) where M_0=alpha1^2/alpha2."""
    M0 = alpha1**2 / alpha2  # = (sum a_k tau_k)^2 / sum a_k tau_k^2
    tau1 = alpha2 / (-alpha1)
    s = 1j * omegas
    return M0 / (1.0 + s * tau1)


# === Run ===
print("Solving (8, 4, 2) Foster spectrum...")
t0 = time.time()
res = foster_solve(8, 4, 2)
print(f"  ({time.time()-t0:.1f} s)")
print()

# Filter coupled modes
valid = (res["taus"] > 0) & (res["a_ks"] > 1e-10 * res["a_ks"].max())
taus_v = res["taus"][valid]
a_v = res["a_ks"][valid]
order = np.argsort(-taus_v)
taus_v = taus_v[order]
a_v = a_v[order]

print(f"Coupled modes: {len(taus_v)}")
print(f"tau range: {taus_v.max()*1e6:.4f} us  ...  {taus_v.min()*1e9:.4f} ns")
print()

# Moments
N_moments = 16
M = np.array([np.sum(a_v * taus_v**n) for n in range(N_moments)])
print(f"M_0 (sum a_k)        = {M[0]:.6e}")
print(f"M_1 (sum a_k tau_k)  = {M[1]:.6e}  (= -alpha_1)")
print(f"M_2                  = {M[2]:.6e}  (= alpha_2)")
print()

# Cauer params
c_taylor = np.array([(-1)**k * M[k] for k in range(N_moments)])
p_cauer = cauer_continued_fraction(c_taylor, max_stages=N_moments - 1)
print(f"Cauer extracted {len(p_cauer)} stages")
print()

# === Frequency sweep ===
# tau_1 = 17 us → omega_1 = 1/tau_1 ≈ 5.9e4 rad/s
# Span 5 decades around it: 1e2 to 1e8 rad/s
omegas = np.logspace(2, 8, 400)

chi_F = chi_foster(omegas, taus_v, a_v)
chi_P11 = chi_pade11(omegas, alpha1_closed, alpha2_closed)

Ks = [2, 3, 4, 5, 6, 8]
chi_Cs = {K: chi_cauer(omegas, p_cauer, K) for K in Ks}

# === Plot ===
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# Magnitude
ax = axs[0, 0]
ax.loglog(omegas, np.abs(chi_F), 'k-', label='Foster (full BEM, 16 modes)', lw=2.5)
ax.loglog(omegas, np.abs(chi_P11), 'r:', label=r'Pade[1,1] closed ($\tau_1$=16.98 us)', lw=2)
colors = plt.cm.viridis(np.linspace(0, 0.85, len(Ks)))
for K, color in zip(Ks, colors):
    ax.loglog(omegas, np.abs(chi_Cs[K]), '--', color=color, label=f'Cauer K={K}', lw=1.2)
# Vertical line at omega_1 = 1/tau_1
ax.axvline(1/tau1_closed, color='gray', ls=':', alpha=0.5)
ax.text(1/tau1_closed*1.1, 0.01, r'$\omega_1=1/\tau_1$', color='gray')
ax.set_xlabel(r'$\omega$ [rad/s]')
ax.set_ylabel(r'$|\chi(j\omega)|$')
ax.legend(loc='lower left', fontsize=9)
ax.grid(True, which='both', alpha=0.3)
ax.set_title(r'Magnitude of $\chi(j\omega)$')

# Phase
ax = axs[0, 1]
ax.semilogx(omegas, np.angle(chi_F)*180/np.pi, 'k-', label='Foster', lw=2.5)
ax.semilogx(omegas, np.angle(chi_P11)*180/np.pi, 'r:', label='Pade[1,1] closed', lw=2)
for K, color in zip(Ks, colors):
    ax.semilogx(omegas, np.angle(chi_Cs[K])*180/np.pi, '--', color=color, label=f'Cauer K={K}', lw=1.2)
ax.axvline(1/tau1_closed, color='gray', ls=':', alpha=0.5)
ax.set_xlabel(r'$\omega$ [rad/s]')
ax.set_ylabel('phase [deg]')
ax.legend(loc='lower left', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_title(r'Phase of $\chi(j\omega)$')

# Relative error vs Foster
ax = axs[1, 0]
for K, color in zip(Ks, colors):
    err = np.abs(chi_Cs[K] - chi_F) / np.abs(chi_F)
    ax.loglog(omegas, err, color=color, label=f'Cauer K={K}', lw=1.5)
err_P11 = np.abs(chi_P11 - chi_F) / np.abs(chi_F)
ax.loglog(omegas, err_P11, 'r:', label='Pade[1,1] closed', lw=2)
ax.axvline(1/tau1_closed, color='gray', ls=':', alpha=0.5)
ax.set_xlabel(r'$\omega$ [rad/s]')
ax.set_ylabel(r'$|\chi_K - \chi_{Foster}|/|\chi_{Foster}|$')
ax.legend(loc='lower left', fontsize=9)
ax.grid(True, which='both', alpha=0.3)
ax.set_title('Relative error vs full Foster')

# Imaginary part (proportional to loss)
ax = axs[1, 1]
ax.loglog(omegas, np.abs(np.imag(chi_F)), 'k-', label='Foster', lw=2.5)
ax.loglog(omegas, np.abs(np.imag(chi_P11)), 'r:', label='Pade[1,1] closed', lw=2)
for K, color in zip(Ks, colors):
    ax.loglog(omegas, np.abs(np.imag(chi_Cs[K])), '--', color=color, label=f'Cauer K={K}', lw=1.2)
ax.axvline(1/tau1_closed, color='gray', ls=':', alpha=0.5)
ax.set_xlabel(r'$\omega$ [rad/s]')
ax.set_ylabel(r'$|Im\,\chi(j\omega)|$  (loss component)')
ax.legend(loc='lower left', fontsize=9)
ax.grid(True, which='both', alpha=0.3)
ax.set_title(r'Imaginary part (proportional to power loss)')

plt.tight_layout()
plt.savefig('bem_foster_freq_response.png', dpi=120)
print("Plot saved to bem_foster_freq_response.png")

# === Print error at key frequencies ===
print()
print("=== Relative error at key frequencies ===")
key_omegas = [1e3, 1/tau1_closed, 1e6, 1e7]
print(f"{'omega':<12s} {'Pade[1,1] closed':<18s}", end="")
for K in Ks:
    print(f"{'Cauer K=' + str(K):<14s}", end="")
print()
for ow in key_omegas:
    idx = np.argmin(np.abs(omegas - ow))
    print(f"  {omegas[idx]:.2e}  ", end="")
    err_p = np.abs(chi_P11[idx] - chi_F[idx]) / np.abs(chi_F[idx])
    print(f"{err_p:<18.4e}", end="")
    for K in Ks:
        err = np.abs(chi_Cs[K][idx] - chi_F[idx]) / np.abs(chi_F[idx])
        print(f"{err:<14.4e}", end="")
    print()

# Save data
with open("bem_foster_freq_response.json", "w") as f:
    json.dump({
        "omegas": omegas.tolist(),
        "chi_foster_re": chi_F.real.tolist(),
        "chi_foster_im": chi_F.imag.tolist(),
        "Ks": Ks,
        "chi_cauer_re": {str(K): chi_Cs[K].real.tolist() for K in Ks},
        "chi_cauer_im": {str(K): chi_Cs[K].imag.tolist() for K in Ks},
        "chi_pade11_re": chi_P11.real.tolist(),
        "chi_pade11_im": chi_P11.imag.tolist(),
        "tau1_closed": tau1_closed,
        "p_cauer": p_cauer.tolist(),
    }, f, indent=2)
print("Data saved to bem_foster_freq_response.json")
