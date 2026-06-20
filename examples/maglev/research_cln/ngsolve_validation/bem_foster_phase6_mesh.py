"""BEM-Foster Phase 6: mesh refinement effect on Cauer & freq response.

Run BEM-Foster + Cauer + chi(j omega) at multiple grid sizes:
- (8, 4, 2)  =  64 cells (baseline)
- (12, 6, 3) = 216 cells (medium)
- (16, 8, 4) = 512 cells (fine)

Compare Cauer truncation accuracy and frequency response error
across mesh resolutions.

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


def chi_foster(omegas, taus, a_ks):
    s = 1j * omegas
    chi = np.zeros(len(omegas), dtype=complex)
    for k in range(len(taus)):
        chi += a_ks[k] / (1.0 + s * taus[k])
    return chi


def chi_cauer(omegas, p_cauer, K):
    s = 1j * omegas
    if K < 1:
        return np.zeros_like(s)
    y = np.full_like(s, 1.0 / p_cauer[K-1])
    for k in range(K-2, -1, -1):
        y = 1.0 / (p_cauer[k] + s * y)
    return y


# === Run on multiple grids ===
configs = [(8, 4, 2), (12, 6, 3), (16, 8, 4)]
data = {}

for cfg in configs:
    print(f"=== Grid {cfg} ===")
    t0 = time.time()
    res = foster_solve(*cfg)
    print(f"  Foster solve: {time.time()-t0:.1f} s")

    valid = (res["taus"] > 0) & (res["a_ks"] > 1e-10 * res["a_ks"].max())
    taus_v = res["taus"][valid]
    a_v = res["a_ks"][valid]
    order = np.argsort(-taus_v)
    taus_v = taus_v[order]
    a_v = a_v[order]

    N_moments = 16
    M = np.array([np.sum(a_v * taus_v**n) for n in range(N_moments)])
    c_taylor = np.array([(-1)**k * M[k] for k in range(N_moments)])
    p_cauer = cauer_continued_fraction(c_taylor, max_stages=N_moments - 1)

    # alpha_n = (-1)^n M_n; alpha_1 < 0, alpha_2 > 0; tau_1 = alpha_2/(-alpha_1) = M[2]/M[1]
    alpha1_BEM = -M[1]
    alpha2_BEM = M[2]
    tau1_BEM = M[2] / M[1]
    print(f"  alpha_1 = {alpha1_BEM:+.6e}  (err = {abs((alpha1_BEM - alpha1_closed)/alpha1_closed):.2e})")
    print(f"  alpha_2 = {alpha2_BEM:+.6e}  (err = {abs((alpha2_BEM - alpha2_closed)/alpha2_closed):.2e})")
    print(f"  tau_1   = {tau1_BEM*1e6:.6f} us  (err = {abs((tau1_BEM - tau1_closed)/tau1_closed):.2e})")
    print(f"  Coupled modes: {len(taus_v)}/{len(res['taus'])}")
    print()

    data[str(cfg)] = {
        "taus_v": taus_v.tolist(),
        "a_v": a_v.tolist(),
        "M": M.tolist(),
        "p_cauer": p_cauer.tolist(),
        "tau1_BEM": tau1_BEM,
        "alpha1_BEM": alpha1_BEM,
        "alpha2_BEM": alpha2_BEM,
        "ncoupled": int(len(taus_v)),
    }

# === Frequency response comparison ===
omegas = np.logspace(2, 8, 400)

fig, axs = plt.subplots(2, 3, figsize=(18, 10))
Ks_to_show = [2, 4, 6, 8]
colors_K = plt.cm.viridis(np.linspace(0, 0.85, len(Ks_to_show)))

for col, cfg in enumerate(configs):
    d = data[str(cfg)]
    taus_v = np.array(d["taus_v"])
    a_v = np.array(d["a_v"])
    p_cauer = np.array(d["p_cauer"])

    chi_F = chi_foster(omegas, taus_v, a_v)

    # Magnitude plot (top row)
    ax = axs[0, col]
    ax.loglog(omegas, np.abs(chi_F), 'k-', label=f'Foster (full, {d["ncoupled"]} modes)', lw=2.5)
    for K, color in zip(Ks_to_show, colors_K):
        if K <= len(p_cauer):
            chi_C = chi_cauer(omegas, p_cauer, K)
            ax.loglog(omegas, np.abs(chi_C), '--', color=color,
                      label=f'Cauer K={K}', lw=1.2)
    ax.axvline(1/tau1_closed, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel(r'$\omega$ [rad/s]')
    ax.set_ylabel(r'$|\chi(j\omega)|$')
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_title(f'Grid {cfg} - Magnitude')

    # Relative error plot (bottom row)
    ax = axs[1, col]
    for K, color in zip(Ks_to_show, colors_K):
        if K <= len(p_cauer):
            chi_C = chi_cauer(omegas, p_cauer, K)
            err = np.abs(chi_C - chi_F) / np.abs(chi_F)
            ax.loglog(omegas, err, color=color, label=f'Cauer K={K}', lw=1.5)
    ax.axvline(1/tau1_closed, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel(r'$\omega$ [rad/s]')
    ax.set_ylabel(r'rel. error vs Foster')
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_title(f'Grid {cfg} - Cauer error vs Foster')

plt.tight_layout()
plt.savefig('bem_foster_phase6_mesh.png', dpi=120)
print("Plot saved to bem_foster_phase6_mesh.png")

# === Summary table ===
print()
print("=== Mesh refinement summary ===")
print(f"{'Grid':<14s} {'N_cells':<8s} {'tau_1 [us]':<14s} {'tau_1 err':<12s} "
      f"{'alpha_2 err':<12s} {'N_coupled':<10s}")
for cfg in configs:
    d = data[str(cfg)]
    ncells = cfg[0]*cfg[1]*cfg[2]
    err_t = abs((d["tau1_BEM"] - tau1_closed)/tau1_closed)
    err_a = abs((d["alpha2_BEM"] - alpha2_closed)/alpha2_closed)
    print(f"  {str(cfg):<12s} {ncells:<8d} {d['tau1_BEM']*1e6:<14.8f} "
          f"{err_t:<12.4e} {err_a:<12.4e} {d['ncoupled']:<10d}")

with open("bem_foster_phase6_mesh_results.json", "w") as f:
    json.dump({
        "configs": [list(c) for c in configs],
        "data": data,
        "tau1_closed": tau1_closed,
        "alpha1_closed": alpha1_closed,
        "alpha2_closed": alpha2_closed,
    }, f, indent=2)
print("Data saved to bem_foster_phase6_mesh_results.json")
