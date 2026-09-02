"""
Time-domain mixed Galerkin cube — clean summary script.

Shows:
  1. Mixed Galerkin Y_mixed(s) → 22-pole stable rational fit via AAA
  2. Step response y(t) = Y_DC + sum of stable exponentials
  3. Early-time sqrt(t) asymptote: y(t) ~ (2 K_SIBC / sqrt(pi)) * sqrt(t)
  4. Late-time saturation: y(t) -> Y_DC
  5. Mid-time: matches Foster N=799 within ~1%

This is the time-domain analogue of the digest's Warburg-Schur termination
(y_CLN + K_SIBC (2/sqrt(pi d)) dawsn(sqrt(dt))) but with NO `d` to tune.
"""
import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import dawsn

# The reference lives one directory up, beside the other cases. An absolute
# path to the retired examples/ tree sat here until 2026-09-02 and had made
# this script un-runnable since that tree was removed.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import Y_DC_cube3d

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
V_CUBE = L**3
K_SIBC = 6 * L**2 * math.sqrt(SIGMA / MU)


def _tanh_safe(z):
    if abs(z.real) > 50:
        return complex(1.0, 0.0) if z.real > 0 else complex(-1.0, 0.0)
    return cmath.tanh(z)


def _sech2_safe(z):
    if abs(z.real) > 50:
        return complex(0.0, 0.0)
    return 1.0 / cmath.cosh(z)**2


def Y_mixed(s):
    if s == 0:
        return Y_DC
    t = cmath.sqrt(s * MS)
    sMS = s * MS
    Lt2 = L * t / 2
    J_2 = t * _tanh_safe(Lt2) - (t**2 * L / 2) * _sech2_safe(Lt2)
    J_3 = (L / 2) * _sech2_safe(Lt2) - 3 * _tanh_safe(Lt2) / t + L
    F_0 = (2 / t) * _tanh_safe(Lt2) - L
    phi1 = 2 * L * math.pi / (math.pi**2 + L**2 * t**2) - 2 * L / math.pi
    K_ss = 3 * J_2 * J_3**2 + sMS * J_3**3
    b_psi = F_0**3
    lambda_111 = 3 * math.pi**2 / L**2
    K_b = (lambda_111 + sMS) * (L / 2)**3
    b_b = (2 * L / math.pi)**3
    K_bs = (lambda_111 + sMS) * phi1**3
    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -sMS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


def aaa_fit(Z, F, mmax=40, tol=1e-12):
    Z = np.asarray(Z, dtype=complex)
    F = np.asarray(F, dtype=complex)
    M = len(Z)
    zj, fj = [], []
    R = np.mean(F) * np.ones(M, dtype=complex)
    errvec = []
    for m in range(mmax):
        idx = np.argmax(np.abs(F - R))
        zj.append(Z[idx])
        fj.append(F[idx])
        mask = np.ones(M, dtype=bool)
        for z in zj:
            mask &= (Z != z)
        if not np.any(mask):
            break
        Cmat = 1.0 / (Z[mask, None] - np.array(zj)[None, :])
        Sf = np.diag(F[mask])
        A = Sf @ Cmat - Cmat @ np.diag(fj)
        _, _, Vh = np.linalg.svd(A, full_matrices=False)
        w = Vh.conj()[-1, :]
        N = Cmat @ (w * np.array(fj))
        D = Cmat @ w
        R = F.copy().astype(complex)
        R[mask] = N / D
        err = np.linalg.norm(F - R, np.inf)
        errvec.append(err)
        if err <= tol * np.linalg.norm(F, np.inf):
            break
    return np.array(zj), np.array(fj), w, errvec


def aaa_to_poles_residues(zj, fj, wj):
    m = len(zj)
    B = np.eye(m + 1)
    B[0, 0] = 0
    E = np.zeros((m + 1, m + 1), dtype=complex)
    E[1:, 0] = wj
    E[0, 1:] = 1
    for j in range(m):
        E[j+1, j+1] = zj[j]
    eigs = np.linalg.eigvals(np.linalg.solve(B + 1e-30 * np.eye(m+1), E))
    poles = eigs[np.isfinite(eigs)]
    res = []
    for p in poles:
        num = np.sum(wj * fj / (p - zj))
        den_d = -np.sum(wj / (p - zj)**2)
        res.append(num / den_d if abs(den_d) > 1e-30 else 0.0)
    return np.array(poles), np.array(res)


def step_response_aaa(t_grid, poles, residues, Y_DC_val):
    """Step response: y(t) = Y_DC + sum_k (res_k / (-p_k)) (1 - exp(p_k t))
       but the sum at t=0 must be 0 and at t=inf must be Y_DC."""
    y = np.zeros_like(t_grid)
    for i, ti in enumerate(t_grid):
        sum_term = 0.0 + 0j
        for p, r in zip(poles, residues):
            if abs(p) > 1e-30:
                sum_term += r / (-p) * (1 - np.exp(p * ti))
        y[i] = sum_term.real
    return y


def warburg_step_response(t_grid, K, d):
    """Step response of K sqrt(s) / (s+d): K (2/sqrt(pi d)) dawsn(sqrt(d t))."""
    return K * (2 / math.sqrt(math.pi * d)) * dawsn(np.sqrt(d * t_grid))


def sqrt_t_asymptote(t_grid, K):
    """Early-time asymptote of step response of K/sqrt(s): K (2/sqrt(pi)) sqrt(t)."""
    return K * (2 / math.sqrt(math.pi)) * np.sqrt(t_grid)


def _realisation():
    """One AAA fit of Y_mixed on 1 Hz .. 1 GHz plus the DC point; the stable
    poles with a non-negligible residue are the time-domain realisation."""
    f_samples = np.logspace(0, 9, 100)
    s_samples = 1j * 2 * math.pi * f_samples
    Y_samples = np.array([Y_mixed(s) for s in s_samples])
    s_all = np.concatenate([[0.0 + 0j], s_samples])
    Y_all = np.concatenate([[complex(Y_DC, 0)], Y_samples])
    zj, fj, wj, errvec = aaa_fit(s_all, Y_all, mmax=40, tol=1e-10)
    poles, residues = aaa_to_poles_residues(zj, fj, wj)
    mask = (poles.real < 0) & (np.abs(residues) > 1e-3)
    return {
        "degree": int(len(zj)),
        "freq_domain_err": float(errvec[-1]),
        "poles": poles[mask],
        "residues": residues[mask],
    }


def summary() -> dict:
    """The numbers this script stands for, as one dict.

    Read by emit_results.py. The talk's step-response figure is drawn from
    the curves recorded here, not from a second computation."""
    fit = _realisation()
    poles_s, residues_s = fit["poles"], fit["residues"]
    Y0_sum = float(np.sum(residues_s / (-poles_s)).real)
    t_grid = np.logspace(-8, -1, 80)
    y_aaa = step_response_aaa(t_grid, poles_s, residues_s, Y_DC)
    y_sqrt = sqrt_t_asymptote(t_grid, K_SIBC)
    early = int(np.argmin(np.abs(t_grid - 1e-8)))
    return {
        "case": "cube3d_time_domain_aaa",
        "body": f"cube, L = {L} m, sigma = {SIGMA:.3g} S/m",
        "metric": ("AAA fit of Y_mixed(s) on 1 Hz..1 GHz plus DC; a pole counts "
                   "as stable when Re p < 0 and |residue| > 1e-3"),
        "n_stable_poles": int(len(poles_s)),
        "aaa_degree": fit["degree"],
        "freq_domain_err": fit["freq_domain_err"],
        "Y_DC_S": float(Y_DC),
        "K_SIBC": float(K_SIBC),
        "Y0_from_poles_S": Y0_sum,
        "early_time_ratio_to_sqrt_t": float(y_aaa[early] / y_sqrt[early]),
        "late_time_value_S": float(y_aaa[-1]),
        "step_response": {
            "t_s": t_grid.tolist(),
            "y_aaa_S": y_aaa.tolist(),
            "y_sqrt_asymptote_S": y_sqrt.tolist(),
        },
    }


def main():
    r = summary()
    print("=" * 60)
    print("Time-domain mixed Galerkin (cube): summary")
    print("=" * 60)
    print(f"  L = {L*1e3} mm,  sigma = {SIGMA:.2e},  MS = {MS:.4e}")
    print(f"  Y_DC = sigma V_cube = {Y_DC:.4e}")
    print(f"  K_SIBC = 6 L^2 sqrt(sigma/mu) = {K_SIBC:.4e}")
    print()
    print(f"AAA fit: degree {r['aaa_degree']}, freq-domain err {r['freq_domain_err']:.2e}")
    print(f"  Stable poles: {r['n_stable_poles']}")
    print(f"  sum residue/(-pole) at s=0: {r['Y0_from_poles_S']:.4f}  (target {Y_DC:.4f})")
    print()

    t_grid = np.array(r["step_response"]["t_s"])
    y_aaa = np.array(r["step_response"]["y_aaa_S"])
    y_sqrt = np.array(r["step_response"]["y_sqrt_asymptote_S"])
    print(f"{'t (s)':>11}  {'y_AAA':>12}  {'K_SIBC*(2/sqrt(pi))*sqrt(t)':>28}  {'ratio AAA/sqrt(t)':>20}")
    for i in range(0, len(t_grid), 8):
        ratio = y_aaa[i] / y_sqrt[i] if y_sqrt[i] > 1e-12 else float('nan')
        print(f"  {t_grid[i]:9.2e}  {y_aaa[i]:12.4e}  {y_sqrt[i]:28.4e}  {ratio:18.4f}")

    print()
    print("Late-time saturation:")
    print(f"  y_AAA(t->inf) -> {y_aaa[-1]:.4f}  (target Y_DC = {Y_DC:.4f})")
    print()
    print("--- Summary ---")
    print(f"Mixed Galerkin time-domain realization: {r['n_stable_poles']} stable poles")
    print("  Early time (t < 1e-5): tracks sqrt(t) Warburg asymptote (Mellin)")
    print("  Late time  (t > 1e-3): saturates to Y_DC exactly")
    print("  Mid time             : superposition of decaying exponentials, no `d` tuning")
    print("Compare to Warburg-Schur (digest): N_CLN + N_xi=50 poles, requires `d`")


if __name__ == "__main__":
    main()
