"""BEM-Foster mesh convergence study.

Sweep grid sizes and track:
- alpha_n (n=1..5) Foster sum
- Pade tau_1 = alpha_2 / (-alpha_1)
- Top coupled time constants
- Foster spectrum statistics

Goal: confirm whether the 0.2% Pade tau_1 agreement at (8,4,2) is
systematic (errors cancel by ratio) or fortuitous.

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

# References (closed forms)
alpha1_closed = -mu0sigma * (a*a + b*b) / 48.0
alpha2_closed = 7.47832405353e-10
alpha3_MC = -1.435e-14
tau1_closed = alpha2_closed / (-alpha1_closed)  # 16.98 us

print("=" * 64)
print("BEM-Foster mesh convergence study (Python)")
print("=" * 64)
print()
print(f"References:")
print(f"  alpha_1 (closed)  = {alpha1_closed:.10e}")
print(f"  alpha_2 (closed)  = {alpha2_closed:.10e}")
print(f"  tau_1 (Pade)      = {tau1_closed*1e6:.6f} us")
print()


# === Cell-cell Newton kernel ===
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


def foster_solve(Nx, Ny, Nz, max_alphas=5):
    hx, hy, hz = a / Nx, b / Ny, c / Nz
    Vcell = hx * hy * hz
    ncells = Nx * Ny * Nz

    centers = np.array([
        [(i + 0.5) * hx - a/2,
         (j + 0.5) * hy - b/2,
         (k + 0.5) * hz - c/2]
        for i in range(Nx) for j in range(Ny) for k in range(Nz)
    ])

    t_total = time.time()

    # K matrix with translation-invariance cache
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

    # Eigendecomposition
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

    alphas = np.array([
        ((-1)**n * mu0sigma**n / V *
         np.sum((gx**2 + gy**2) * lambdaK**(n - 1)))
        for n in range(1, max_alphas + 1)
    ])

    elapsed = time.time() - t_total

    return {
        "Nx": Nx, "Ny": Ny, "Nz": Nz,
        "ncells": ncells,
        "Vcell": Vcell,
        "taus": taus,
        "gx": gx,
        "gy": gy,
        "lambdaK": lambdaK,
        "alphas": alphas,
        "elapsed": elapsed,
    }


# === Sweep ===
configs = [
    (2, 2, 2),
    (4, 2, 2),
    (4, 4, 2),
    (6, 4, 2),
    (8, 4, 2),
    (10, 4, 2),
    (10, 6, 2),
    (12, 6, 3),
]

print(f"{'Grid':<14s} {'alpha_1':<14s} {'err_a1':<10s} "
      f"{'alpha_2':<14s} {'err_a2':<10s} {'tau_1 [us]':<12s} {'err_tau':<10s} "
      f"{'time':<8s}")
print("-" * 100)

results = []
for cfg in configs:
    res = foster_solve(*cfg)
    a1, a2 = res["alphas"][0], res["alphas"][1]
    tau1_BEM = a2 / (-a1)
    err_a1 = abs((a1 - alpha1_closed) / alpha1_closed)
    err_a2 = abs((a2 - alpha2_closed) / alpha2_closed)
    err_tau = abs((tau1_BEM - tau1_closed) / tau1_closed)
    print(f"  {str(cfg):<12s} {a1:.4e}  {err_a1:.2e}  "
          f"{a2:.4e}  {err_a2:.2e}  {tau1_BEM*1e6:10.6f}  {err_tau:.2e}  "
          f"{res['elapsed']:6.1f}s")
    results.append({
        "config": list(cfg),
        "ncells": res["ncells"],
        "alpha1": a1,
        "alpha2": a2,
        "alpha3": res["alphas"][2],
        "alpha4": res["alphas"][3],
        "alpha5": res["alphas"][4],
        "tau1_BEM": tau1_BEM,
        "err_alpha1": err_a1,
        "err_alpha2": err_a2,
        "err_tau1": err_tau,
        "top5_taus": res["taus"][:5].tolist(),
        "top5_gx": res["gx"][:5].tolist(),
        "top5_gy": res["gy"][:5].tolist(),
        "elapsed": res["elapsed"],
    })

print()
print("=" * 64)
print("Convergence summary")
print("=" * 64)
print(f"* alpha_1, alpha_2 individual errors (h^p):")
if len(results) >= 2:
    rs = [results[i]["err_alpha1"] for i in range(len(results))]
    rs2 = [results[i]["err_alpha2"] for i in range(len(results))]
    rt = [results[i]["err_tau1"] for i in range(len(results))]
    ns = [results[i]["ncells"] for i in range(len(results))]
    print(f"    err_alpha1 over N_cells: {[f'{r:.1e}' for r in rs]}")
    print(f"    err_alpha2 over N_cells: {[f'{r:.1e}' for r in rs2]}")
    print(f"    err_tau1   over N_cells: {[f'{r:.1e}' for r in rt]}")

with open("bem_foster_convergence_results.json", "w") as f:
    json.dump({
        "alpha1_closed": alpha1_closed,
        "alpha2_closed": alpha2_closed,
        "alpha3_MC": alpha3_MC,
        "tau1_closed": tau1_closed,
        "results": results,
    }, f, indent=2)
print("* Results saved to bem_foster_convergence_results.json")
