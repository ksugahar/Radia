"""1-D nonlinear magnetic diffusion: do deliberate surface modes still pay off?

sigma dB(H)/dt = d2H/dx2 on a half-slab x in [0, L]; H(0,t) = H0 sin(wt), dH/dx(L) = 0.
B(H) = mu0 H + Js tanh(mu0 mur H / Js).  Backward Euler + Newton on a fine grid gives the
periodic-steady reference.  Over the last period each H(., t) is projected onto
  (a) bulk only: Legendre polynomials,
  (b) bulk + analytic surface modes from linear theory, Re/Im exp(-(1+i) x/delta_lin),
  (c) bulk + surface POD learned from other amplitudes (bulk component removed first).
The metric is the relative space-time L2 projection error, i.e. the best any Galerkin
ROM on that basis could do.  Diagnostic evidence only.
"""
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import scipy

import numpy as np
from numpy.polynomial import legendre as leg
from scipy.linalg import solve_banded

MU0 = 4e-7 * np.pi
L, SIGMA, MUR, JS, FREQ = 0.010, 5.0e6, 1000.0, 1.6, 50.0
OMEGA = 2 * np.pi * FREQ
DELTA_LIN = np.sqrt(2.0 / (OMEGA * MU0 * MUR * SIGMA))
N = 2000                       # grid cells
STEPS_PER_PERIOD, PERIODS = 800, 4
x = np.linspace(0.0, L, N + 1)
h = x[1] - x[0]
w = np.full(N + 1, h); w[0] = w[-1] = 0.5 * h      # trapezoid weights


def b_of_h(H):
    return MU0 * H + JS * np.tanh(MU0 * MUR * H / JS)


def db_dh(H):
    return MU0 + MU0 * MUR / np.cosh(MU0 * MUR * H / JS) ** 2


def simulate(H0):
    dt = 1.0 / (FREQ * STEPS_PER_PERIOD)
    H = np.zeros(N + 1)
    snaps = []
    for step in range(1, STEPS_PER_PERIOD * PERIODS + 1):
        t = step * dt
        Hn = H.copy()
        Bn = b_of_h(Hn)
        H[0] = H0 * np.sin(OMEGA * t)
        for _ in range(50):
            # interior + Neumann (ghost node) rows 1..N; node 0 is Dirichlet.
            lap = np.zeros(N + 1)
            lap[1:N] = (H[0:N - 1] - 2 * H[1:N] + H[2:N + 1]) / h**2
            lap[N] = 2 * (H[N - 1] - H[N]) / h**2
            r = SIGMA * (b_of_h(H[1:]) - Bn[1:]) / dt - lap[1:]
            diag = SIGMA * db_dh(H[1:]) / dt + 2.0 / h**2
            upper = np.full(N, -1.0 / h**2); lower = np.full(N, -1.0 / h**2)
            lower[-1] = -2.0 / h**2                      # Neumann row couples to N-1 twice
            ab = np.zeros((3, N))
            ab[0, 1:] = upper[:-1]
            ab[1] = diag
            ab[2, :-1] = lower[1:]
            dH = solve_banded((1, 1), ab, -r)
            H[1:] += dH
            if np.max(np.abs(dH)) <= 1e-10 * max(1.0, abs(H0)):
                break
        else:
            raise RuntimeError(f"Newton did not converge at H0={H0}, step {step}")
        if step > STEPS_PER_PERIOD * (PERIODS - 1):
            snaps.append(H.copy())
    return np.array(snaps).T                             # (N+1, steps)


def orthonormal(columns):
    Q, _ = np.linalg.qr(np.sqrt(w)[:, None] * columns)
    return Q                                             # orthonormal in the weighted norm


def projection_error(snaps, basis_cols):
    Q = orthonormal(basis_cols)
    S = np.sqrt(w)[:, None] * snaps
    residual = S - Q @ (Q.T @ S)
    return float(np.linalg.norm(residual) / np.linalg.norm(S))


def bulk(n):
    t = 2 * x / L - 1
    return np.column_stack([leg.legval(t, [0] * k + [1]) for k in range(n)])


def analytic_surface():
    z = np.exp(-(1 + 1j) * x / DELTA_LIN)
    return np.column_stack([z.real, z.imag])


def surface_pod(training, n_bulk, m, extra=None):
    B = orthonormal(bulk(n_bulk) if extra is None else np.hstack([bulk(n_bulk), extra]))
    # Normalize each amplitude so small-amplitude (linear) shapes are not swamped.
    S = np.hstack([np.sqrt(w)[:, None] * s / np.linalg.norm(np.sqrt(w)[:, None] * s)
                   for s in training])
    R = S - B @ (B.T @ S)                                # remove the bulk content
    U, sv, _ = np.linalg.svd(R, full_matrices=False)
    return U[:, :m] / np.sqrt(w)[:, None], sv[: m + 3] / sv[0]


def main():
    started = time.perf_counter()
    train_amps = [30.0, 300.0, 3000.0, 30000.0]
    test_amps = [10.0, 100.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0]
    sims = {a: simulate(a) for a in sorted(set(train_amps + test_amps))}
    peak_b = {a: float(np.max(np.abs(b_of_h(s[0])))) for a, s in sims.items()}
    n_bulk = 5
    pod2, sv = surface_pod([sims[a] for a in train_amps], n_bulk, 2)
    pod4, _ = surface_pod([sims[a] for a in train_amps], n_bulk, 4)
    lin = analytic_surface()
    pod_after_lin2, _ = surface_pod([sims[a] for a in train_amps], n_bulk, 2, extra=lin)
    rows = []
    for a in test_amps:
        s = sims[a]
        row = {"H0_A_per_m": a, "surface_peak_B_T": peak_b[a],
               "in_training_set": a in train_amps,
               "bulk_only_n5": projection_error(s, bulk(5)),
               "bulk_only_n7": projection_error(s, bulk(7)),
               "bulk_only_n9": projection_error(s, bulk(9)),
               "bulk5_plus_analytic_surface2": projection_error(s, np.hstack([bulk(5), analytic_surface()])),
               "bulk5_plus_surface_pod2": projection_error(s, np.hstack([bulk(5), pod2])),
               "bulk5_plus_surface_pod4": projection_error(s, np.hstack([bulk(5), pod4])),
               "bulk5_plus_analytic2_plus_pod2": projection_error(
                   s, np.hstack([bulk(5), lin, pod_after_lin2]))}
        rows.append(row)
        print(f"H0={a:>8.0f} A/m  Bpk={peak_b[a]:.2f} T  "
              f"bulk5 {row['bulk_only_n5']:.1e}  bulk7 {row['bulk_only_n7']:.1e}  "
              f"bulk5+lin2 {row['bulk5_plus_analytic_surface2']:.1e}  "
              f"bulk5+pod2 {row['bulk5_plus_surface_pod2']:.1e}  bulk5+pod4 {row['bulk5_plus_surface_pod4']:.1e}  "
              f"bulk5+lin2+pod2 {row['bulk5_plus_analytic2_plus_pod2']:.1e}"
              + ("  (train)" if a in train_amps else ""), flush=True)
    out = Path(__file__).with_name("nonlinear_surface_modes_1d.json")
    out.write_text(json.dumps({
        "schema": "radia.validation.nonlinear-surface-modes-1d.v1",
        "host": platform.node(), "python": sys.version.split()[0],
        "numpy": np.__version__, "scipy": scipy.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "acceptance": "HOLD: projection (best-approximation) error, not a solved ROM",
        "parameters": {"L_m": L, "sigma_S_per_m": SIGMA, "mur": MUR, "Js_T": JS, "freq_Hz": FREQ,
                       "delta_linear_m": DELTA_LIN, "cells": N, "steps_per_period": STEPS_PER_PERIOD,
                       "periods": PERIODS, "bh": "B = mu0 H + Js tanh(mu0 mur H / Js)"},
        "training_amplitudes_A_per_m": train_amps, "surface_pod_singular_value_ratios": sv.tolist(),
        "metric": "relative space-time L2 projection error over the last period",
        "rows": rows, "runtime_s": time.perf_counter() - started}, indent=2), encoding="utf-8")
    print(f"delta_lin={DELTA_LIN*1e3:.3f} mm  runtime {time.perf_counter() - started:.1f}s")


if __name__ == "__main__":
    main()
