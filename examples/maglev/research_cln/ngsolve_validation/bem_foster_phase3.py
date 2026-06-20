"""BEM-Foster Phase 3: full Foster eigendecomposition (Python).

Mirror of bem_foster_phase3.wls in float64.

Strategy:
- Build cell-cell Newton matrix K (Phase 2).
- Solve generalized eigenproblem K phi = mu V_cell phi via scipy.linalg.eigh.
  => Foster eigenvalues lambda_k = mu_k / (4 pi),
     time constants tau_k = mu0 sigma lambda_k.
- M-orthonormalize eigenvectors.
- Project A_ext = (-y/2, x/2, 0) onto modes.
- Foster sum: alpha_n = (-1)^n * mu0sigma^n / V * sum_k (gx^2 + gy^2) lambda^(n-1).
- Compare with closed-form alpha_1, alpha_2 and Born MC alpha_3.

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

# References
alpha1_closed = -mu0sigma * (a*a + b*b) / 48.0
alpha2_closed = 7.47832405353e-10
alpha3_MC = -1.435e-14
tau1_closed = 16.98281e-6

print("=" * 64)
print("BEM-Foster Phase 3: full Foster eigendecomposition (Python)")
print("=" * 64)
print()
print(f"Material : Cu, sigma = {sigma:.4g} S/m")
print(f"Geometry : 5 x 2 x 1 mm  (V = {V:.4e} m^3)")
print()
print("References:")
print(f"  alpha_1 (closed)    = {alpha1_closed:.10e}")
print(f"  alpha_2 (closed)    = {alpha2_closed:.10e}")
print(f"  alpha_3 (3-dig MC)  = {alpha3_MC:.4e}")
print(f"  tau_1   (closed)    = {tau1_closed:.6e} s = 16.98 us")
print()


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

    print(f"  Grid ({Nx}, {Ny}, {Nz}): {ncells} cells, V_cell = {Vcell:.4e}")

    # Build K with translation-invariance cache
    print("  Building K...")
    t0 = time.time()
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
    print(f"    K assembly: {time.time()-t0:.2f} s, {len(cache)} unique values")

    # Eigendecomposition: solve K c = mu V_cell c
    print(f"  Eigendecomposition ({ncells} modes)...")
    t0 = time.time()
    eigvals, eigvecs = la.eigh(K / Vcell)
    print(f"    done: {time.time()-t0:.2f} s")
    # eigh returns eigenvectors as columns, ascending eigenvalues
    # Reorder descending
    idx = np.argsort(-eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # M-orthonormalize: <phi_k, phi_k>_M = sum_i phi^2 V_cell = 1
    # eigh returns standard-orthonormal: phi^T phi = 1
    # M-orthonormal: phi^T (V_cell I) phi = V_cell -> divide by sqrt(V_cell)
    phis = eigvecs / np.sqrt(Vcell)  # columns are M-orthonormal eigenvectors

    # Newton eigenvalues lambda_k = mu_k / (4 pi)
    lambdaK = eigvals / (4.0 * np.pi)
    taus = mu0sigma * lambdaK

    # Project A_ext (B_0 = 1)
    xs = centers[:, 0]
    ys = centers[:, 1]
    # g_k^x = sum_i phi_{k,i} A_ext^x_i V_cell, with A_ext^x = -y/2
    gx = phis.T @ (-ys / 2.0) * Vcell   # shape (ncells,)
    gy = phis.T @ (+xs / 2.0) * Vcell

    # Foster sum
    alphas = []
    for n in range(1, 6):
        a_n = ((-1)**n * mu0sigma**n / V *
               np.sum((gx**2 + gy**2) * lambdaK**(n - 1)))
        alphas.append(a_n)

    return {
        "taus": taus,
        "alphas": np.array(alphas),
        "gx": gx,
        "gy": gy,
        "lambdaK": lambdaK,
        "centers": centers,
    }


# === Run on (8, 4, 2) ===
print("=== Solving (8, 4, 2) grid ===")
res = foster_solve(8, 4, 2)
print()

# === Foster spectrum: top 10 modes ===
print("=== Foster spectrum: top 10 modes by tau_k ===")
print(f"  k    tau_k [us]        g_k^x         g_k^y         a_k")
ord_idx = np.argsort(-res["taus"])
for k in range(min(10, len(res["taus"]))):
    ki = ord_idx[k]
    tk = res["taus"][ki]
    gxk = res["gx"][ki]
    gyk = res["gy"][ki]
    ak = (gxk**2 + gyk**2) / (V * res["lambdaK"][ki])
    print(f"  {k+1:2d}   {tk*1e6:14.6e}   {gxk:12.4e}  {gyk:12.4e}  {ak:12.4e}")
print()

# === alpha_n via Foster sum ===
print("=== alpha_n via Foster sum (B_0 = 1) ===")
refs = [alpha1_closed, alpha2_closed, alpha3_MC, None, None]
ref_labels = ["closed", "closed", "MC", "?", "?"]
for n in range(5):
    a_n = res["alphas"][n]
    ref = refs[n]
    print(f"  alpha_{n+1} = {a_n:.10e}", end="")
    if ref is not None:
        err = abs((a_n - ref) / ref)
        print(f"   ref ({ref_labels[n]}) = {ref:.6e}   rel.err = {err:.4e}")
    else:
        print()
print()

# === Summary ===
print("=" * 64)
print("Summary: BEM-Foster on (8,4,2) grid")
print("=" * 64)
print(f"* tau_1 (BEM dominant)  = {res['taus'][ord_idx[0]]*1e6:.6f} us")
print(f"* tau_1 (closed form)   = {tau1_closed*1e6:.6f} us")
print(f"* alpha_1 BEM/closed    = {res['alphas'][0]/alpha1_closed:.6f}")
print(f"* alpha_2 BEM/closed    = {res['alphas'][1]/alpha2_closed:.6f}")
print(f"* alpha_3 BEM/MC-ref    = {res['alphas'][2]/alpha3_MC:.6f}")
print()
print(f"* tau_k range: {res['taus'].max()*1e6:.4f} us  ...  "
      f"{res['taus'].min()*1e9:.4f} ns")
print(f"* Number of modes: {len(res['taus'])}")
print()

# Save
with open("bem_foster_phase3_results_py.json", "w") as f:
    json.dump({
        "grid": [8, 4, 2],
        "taus": res["taus"].tolist(),
        "gx": res["gx"].tolist(),
        "gy": res["gy"].tolist(),
        "alphas": res["alphas"].tolist(),
        "alpha1_closed": alpha1_closed,
        "alpha2_closed": alpha2_closed,
        "alpha3_MC": alpha3_MC,
    }, f, indent=2)
print("* Results saved to bem_foster_phase3_results_py.json")
