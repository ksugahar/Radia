"""BEM-Foster Phase 2: alpha_2 validation on Cartesian hex grid (Python).

Mirror of bem_foster_alpha2.wls in float64.

Strategy:
- Discretize cuboid Omega = [-a/2,a/2] x [-b/2,b/2] x [-c/2,c/2]
  into Nx x Ny x Nz uniform hex cells with piecewise-constant basis.
- Cell-cell Newton matrix K_{ij} = int_{Omega_i} int_{Omega_j} 1/|r-r'| dVdV'
  reduced to 3D via change of variables u = r - r':
    K_ij(d) = int_{|u_a|<h_a} prod_a(h_a - |u_a|) / |u + d| du,  d = r_i - r_j
- alpha_2^BEM = (mu0 sigma)^2 / (16 pi V) * sum_{ij} (x_i x_j + y_i y_j) K_{ij}
- Compare to closed form alpha_2 = 7.4783e-10 (Mathematica reference).

Date: 2026-05-01
"""

import json
import time

import numpy as np
from scipy.integrate import tplquad

# === Constants ===
sigma = 5.8e7  # S/m, Cu
mu0 = 4.0 * np.pi * 1e-7  # H/m
a, b, c = 5e-3, 2e-3, 1e-3  # m
V = a * b * c

alpha2_closed = 7.47832405353e-10  # closed form reference

print("=" * 64)
print("BEM-Foster Phase 2: alpha_2 validation (Python, float64)")
print("=" * 64)
print()
print(f"Material : Cu, sigma = {sigma:.4g} S/m")
print(f"Geometry : {a*1e3:.3g} x {b*1e3:.3g} x {c*1e3:.3g} mm  (V = {V:.4e} m^3)")
print(f"alpha_2 (closed form): {alpha2_closed:.10e}")
print()


def cell_newton_off(dx, dy, dz, hx, hy, hz, epsabs=0.0, epsrel=1e-10):
    """Off-diagonal cell-cell Newton integral (smooth integrand)."""
    def integrand(w, v, u):
        return ((hx - abs(u)) * (hy - abs(v)) * (hz - abs(w)) /
                np.sqrt((u + dx)**2 + (v + dy)**2 + (w + dz)**2))
    val, _ = tplquad(integrand,
                     -hx, hx, -hy, hy, -hz, hz,
                     epsabs=epsabs, epsrel=epsrel)
    return val


def cell_newton_self(hx, hy, hz, epsabs=0.0, epsrel=1e-10):
    """Self-cell Newton integral (1/r corner singularity, octant symmetry)."""
    def integrand(w, v, u):
        return ((hx - u) * (hy - v) * (hz - w) /
                np.sqrt(u*u + v*v + w*w))
    val, _ = tplquad(integrand,
                     0.0, hx, 0.0, hy, 0.0, hz,
                     epsabs=epsabs, epsrel=epsrel)
    return 8.0 * val


def cell_newton(dx, dy, dz, hx, hy, hz):
    if dx == 0.0 and dy == 0.0 and dz == 0.0:
        return cell_newton_self(hx, hy, hz)
    return cell_newton_off(dx, dy, dz, hx, hy, hz)


# === Sanity check: unit cube self-energy ===
print("=== Sanity check ===")
t0 = time.time()
K_unit = cell_newton_self(1.0, 1.0, 1.0)
print(f"  K_self(unit cube) = {K_unit:.12f}  ({time.time()-t0:.2f} s)")
print(f"  Expected         ~ 1.88231264338  (literature)")
print()


def build_and_compute(Nx, Ny, Nz):
    hx, hy, hz = a / Nx, b / Ny, c / Nz
    ncells = Nx * Ny * Nz

    centers = np.array([
        [(i + 0.5) * hx - a/2,
         (j + 0.5) * hy - b/2,
         (k + 0.5) * hz - c/2]
        for i in range(Nx) for j in range(Ny) for k in range(Nz)
    ])

    print(f"  Grid ({Nx}, {Ny}, {Nz}): {ncells} cells, "
          f"{ncells*(ncells+1)//2} pairs")
    t0 = time.time()

    # Translation-invariance cache: K depends on (|delta_i|, |delta_j|, |delta_k|)
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
    dt = time.time() - t0
    print(f"    K assembly: {dt:.2f} s, {len(cache)} unique kernel values")

    xs = centers[:, 0]
    ys = centers[:, 1]
    Dxx = xs @ K @ xs
    Dyy = ys @ K @ ys
    alpha2_BEM = (mu0 * sigma)**2 / (16.0 * np.pi * V) * (Dxx + Dyy)
    err = abs((alpha2_BEM - alpha2_closed) / alpha2_closed)
    return alpha2_BEM, err, dt


# === Convergence study ===
print("=== Convergence study ===")
print()
print(f"{'Grid':<18s} {'alpha_2^BEM':<22s} {'rel. error':<12s}")
print("-" * 58)

configs = [(2, 2, 2), (4, 2, 2), (4, 4, 2), (6, 4, 2), (8, 4, 2)]
results = []
for cfg in configs:
    a2, err, dt = build_and_compute(*cfg)
    print(f"  {str(cfg):<16s} {a2:.12e}  {err:.4e}  ({dt:.2f} s)")
    results.append({
        "config": list(cfg),
        "alpha2": a2,
        "rel_err": err,
        "time": dt,
    })

print()
print("=" * 64)
print("Summary")
print("=" * 64)
print(f"* alpha_2 (closed form): {alpha2_closed:.12e}")
print(f"* alpha_2 (BEM, finest): {results[-1]['alpha2']:.12e}")
print(f"* relative error       : {results[-1]['rel_err']:.4e}")

if len(results) >= 2:
    r1, r2 = results[-2], results[-1]
    h1 = (V / np.prod(r1['config']))**(1/3)
    h2 = (V / np.prod(r2['config']))**(1/3)
    p = np.log(r1['rel_err'] / r2['rel_err']) / np.log(h1 / h2)
    print(f"* Estimated convergence rate (err ~ h^p): p = {p:.3f}")

with open("bem_foster_alpha2_results_py.json", "w") as f:
    json.dump({
        "alpha2_closed": alpha2_closed,
        "results": results,
    }, f, indent=2)
print("* Results saved to bem_foster_alpha2_results_py.json")
