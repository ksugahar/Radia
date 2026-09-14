"""
Cylinder mixed Galerkin: rank-N bulk Krylov + 1-DOF surface envelope.

Phase 3 of the 2026-05-28 -> 2026-06-12 research sprint, corrected
post Phase 8b artifact discovery.

Tests whether adding bulk Krylov modes (phi_k satisfying recursive
-Laplacian phi_{k+1} = phi_k with phi_k(a) = 0) helps the wall-band
accuracy beyond the 1-DOF surface envelope baseline.

## Bulk Krylov basis

The Krylov modes are obtained symbolically (sympy) by recursive Poisson
inversion:

    phi_0(r) = (a^2 - r^2) / 4              (-Laplacian phi_0 = 1)
    phi_1(r) = ...                           (-Laplacian phi_1 = phi_0)
    ...

Compare this bulk-rank sweep with the surface-rank sweep in
02_senior_tower_truncation.py. Their relative effectiveness must be measured;
it is not established by the number or names of basis functions.

## Measurement status

The former table predicted rank-independent errors before the current
reference/projection corrections. Do not treat that prediction or the shared
WIP's replacement improvement factors as accepted evidence. `summary()` owns
the complex-relative-error measurement interface; numerical results remain HOLD
pending the numerical owner's independent validation on a compute host.
"""

from __future__ import annotations

import math
import cmath

import numpy as np
import sympy as sp
from scipy.integrate import quad

from radia.maglev.mixed_galerkin.references import (
    K_SIBC_cylinder,
    Y_DC_cylinder,
    Y_exact_cylinder,
)

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
A_NUM = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cylinder(A_NUM, SIGMA)
K_SIBC = K_SIBC_cylinder(A_NUM, SIGMA, MU)


# -----------------------------------------------------------------------
# Symbolic bulk Krylov basis: -Laplacian_2D phi_{k+1} = phi_k.
# -----------------------------------------------------------------------
N_MAX = 5
r_sym, a_sym = sp.symbols("r a", positive=True)
phi = [(a_sym**2 - r_sym**2) / 4]  # phi_0
for _ in range(N_MAX - 1):
    prev = phi[-1]
    rr = sp.symbols("rr", positive=True)
    inner = sp.integrate(rr * prev.subs(r_sym, rr), (rr, 0, r_sym))
    raw = sp.integrate(-inner / r_sym, r_sym)
    C = -raw.subs(r_sym, a_sym)
    phi.append(sp.expand(raw + C))


def _to_num(expr):
    return float(expr.subs(a_sym, A_NUM))


# Symbolic bulk-bulk matrices.
K0_bb_sym = sp.zeros(N_MAX, N_MAX)
K1_bb_sym = sp.zeros(N_MAX, N_MAX)
b_b_sym = sp.zeros(N_MAX, 1)
for i in range(N_MAX):
    dphi_i = sp.diff(phi[i], r_sym)
    b_b_sym[i] = sp.integrate(phi[i] * r_sym, (r_sym, 0, a_sym)) * 2 * sp.pi
    for j in range(N_MAX):
        dphi_j = sp.diff(phi[j], r_sym)
        K0_bb_sym[i, j] = sp.integrate(dphi_i * dphi_j * r_sym, (r_sym, 0, a_sym)) * 2 * sp.pi
        K1_bb_sym[i, j] = sp.integrate(phi[i] * phi[j] * r_sym, (r_sym, 0, a_sym)) * 2 * sp.pi

K0_BB_N = np.array([[_to_num(K0_bb_sym[i, j]) for j in range(N_MAX)] for i in range(N_MAX)])
K1_BB_N = np.array([[_to_num(K1_bb_sym[i, j]) for j in range(N_MAX)] for i in range(N_MAX)])
B_B_N = np.array([_to_num(b_b_sym[i]) for i in range(N_MAX)])

phi_funcs = [sp.lambdify(r_sym, p.subs(a_sym, A_NUM), "numpy") for p in phi]
dphi_funcs = [sp.lambdify(r_sym, sp.diff(p, r_sym).subs(a_sym, A_NUM), "numpy") for p in phi]


def _integrate_complex(f, lo, hi, t=None, limit=400):
    points = None
    if t is not None and abs(t) > 1.0:
        skin = 1.0 / abs(t)
        if skin < hi - lo:
            points = list(np.geomspace(skin / 100, min(20 * skin, (hi - lo) * 0.99), 30))
            points = [p for p in points if lo < p < hi]
    r_real, _ = quad(lambda x: f(x).real, lo, hi, limit=limit, points=points)
    r_imag, _ = quad(lambda x: f(x).imag, lo, hi, limit=limit, points=points)
    return complex(r_real, r_imag)


def Y_mixed_galerkin(s, N_bulk: int):
    t = cmath.sqrt(s * MS)

    K0_ss = _integrate_complex(
        lambda u: t**2 * cmath.exp(-2 * u * t) * 2 * math.pi * (A_NUM - u), 0, A_NUM, t=t)
    K1_ss = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1)**2 * 2 * math.pi * (A_NUM - u), 0, A_NUM, t=t)
    b_s = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1) * 2 * math.pi * (A_NUM - u), 0, A_NUM, t=t)

    K0_bs = np.zeros(N_bulk, dtype=complex)
    K1_bs = np.zeros(N_bulk, dtype=complex)
    for i in range(N_bulk):
        K0_bs[i] = _integrate_complex(
            lambda u, ii=i: dphi_funcs[ii](A_NUM - u) * t * cmath.exp(-u * t) * 2 * math.pi * (A_NUM - u),
            0, A_NUM, t=t)
        K1_bs[i] = _integrate_complex(
            lambda u, ii=i: phi_funcs[ii](A_NUM - u) * (cmath.exp(-u * t) - 1) * 2 * math.pi * (A_NUM - u),
            0, A_NUM, t=t)

    dim = N_bulk + 1
    K_mat = np.zeros((dim, dim), dtype=complex)
    b_vec = np.zeros(dim, dtype=complex)
    K_mat[:N_bulk, :N_bulk] = K0_BB_N[:N_bulk, :N_bulk] + s * MS * K1_BB_N[:N_bulk, :N_bulk]
    b_vec[:N_bulk] = B_B_N[:N_bulk]
    K_mat[N_bulk, N_bulk] = K0_ss + s * MS * K1_ss
    b_vec[N_bulk] = b_s
    cross = K0_bs + s * MS * K1_bs
    K_mat[:N_bulk, N_bulk] = cross
    K_mat[N_bulk, :N_bulk] = cross

    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / (math.pi * A_NUM**2)
    return Y_DC * (1 + v_avg)


SWEEP = np.logspace(0, 8, 81)
WALL_BAND_HZ = (1e4, 1e6)
METRIC = "abs(Y_exact - Y_mixed) / abs(Y_exact)"


def summary() -> dict:
    """Measure the existing solver; do not imply numerical acceptance."""
    wall = (SWEEP > WALL_BAND_HZ[0]) & (SWEEP < WALL_BAND_HZ[1])
    if not np.any(wall):
        raise ValueError("Sweep has no points inside the wall band")
    tower = {}
    for rank in (1, 2, 3, 4):
        errors = []
        for frequency in SWEEP:
            s = 2j * math.pi * frequency
            exact = Y_exact_cylinder(s, A_NUM, SIGMA, MU)
            mixed = Y_mixed_galerkin(s, rank)
            if not np.isfinite(exact) or not np.isfinite(mixed) or abs(exact) == 0:
                raise ValueError("Non-finite admittance or zero reference in sweep")
            errors.append(abs(exact - mixed) / abs(exact))
        errors = np.asarray(errors)
        if not np.all(np.isfinite(errors)):
            raise ValueError("Non-finite relative errors in sweep")
        tower[str(rank)] = {
            "n_unknowns": rank + 1,
            "max_error_pct": float(errors.max() * 100),
            "max_error_at_hz": float(SWEEP[errors.argmax()]),
            "wall_band_max_error_pct": float(errors[wall].max() * 100),
        }
    return {
        "case": "cylinder_bulk_tower", "validation_status": "HOLD",
        "description": "rank-N CLN bulk + one surface DOF; numerical acceptance pending",
        "metric": METRIC,
        "geometry": {"a_m": A_NUM, "sigma_S_per_m": SIGMA, "mu_H_per_m": MU},
        "sweep": {"f_lo_hz": float(SWEEP[0]), "f_hi_hz": float(SWEEP[-1]),
                  "n_points": int(SWEEP.size), "wall_band_hz": list(WALL_BAND_HZ),
                  "wall_band_endpoints": "excluded"},
        "by_n_bulk": tower,
    }


def main():
    print("=== Cylinder mixed Galerkin: rank-N bulk Krylov + 1-DOF surface ===")
    print(f"a = {A_NUM*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print()
    print("Bulk Krylov mode boundary values phi_k(a) (should all be 0):")
    for k, p in enumerate(phi):
        print(f"  phi_{k}(a) = {sp.expand(p.subs(r_sym, a_sym))}")
    print()
    print("Sweep:  rank-N bulk + 1-DOF surface (envelope intrinsic, no Senior tower)")
    measured = summary()
    print("Numerical acceptance: HOLD")
    print(f"{'N_bulk':>7}  {'basis':>5}  {'max anywhere':>15}  {'wall band max':>15}")
    for rank, row in measured["by_n_bulk"].items():
        print(f"  {rank}    {row['n_unknowns']:5d}  {row['max_error_pct']:13.5f}%  "
              f"{row['wall_band_max_error_pct']:13.5f}%")


if __name__ == "__main__":
    main()
