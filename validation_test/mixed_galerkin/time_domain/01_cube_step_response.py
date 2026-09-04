"""
Time-domain realisation of the cube's mixed-Galerkin admittance.

The mixed admittance Y_mixed(s) has a square root in it (the surface term),
so a transient user asks what to do with it.  Nothing new is needed there:
the time-domain surface impedance condition is settled work (Oh and
Schutt-Aine 1995; Yuferev and Ida 2009).  The square-root kernel is a
diffusive object, a sum of real decaying exponentials, and that is the form
used here.

  1. Y_mixed(s) is closed-form (rank-1 bulk + the exact tensor K_ss).
  2. Real-pole passive Foster fit (radia.maglev.mixed_galerkin.rom_fit):
     AAA finds the dominant poles, a log-spaced set fills the diffusive
     tail, residues by non-negative least squares.  Poles real and
     negative, residues >= 0, DC exact by construction.
  3. Exact reference: numerical inverse Laplace transform (Talbot, mpmath)
     of Y_mixed(s)/s on the same time grid.  The realisation is measured
     against it, not against itself.
  4. Early time: y(t) -> (2 K_SIBC / sqrt(pi)) sqrt(t); late time: y -> Y_DC.

Why the realisation is a real-pole fit and not the AAA poles themselves
(the state of this script until 2026-09-04): an AAA fit on one-sided samples
s = j omega returns complex poles that are not conjugate pairs (21 of 22 for
the cube, e.g. -2.2e8 - 2.6e10 j).  Summing r/(-p) (1 - e^{p t}) and taking
the real part is the conjugate-symmetrised realisation, and its terms
oscillate: the step response strayed by up to 11 % below 1 us against the
exact inverse Laplace transform while the frequency-domain error at the
samples was 6e-10.  The talk's figure showed that ripple as a wave.
diagnose_aaa_direct() keeps that measurement as a negative result.
"""
import cmath
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

# The reference lives one directory up, beside the other cases. An absolute
# path to the retired examples/ tree sat here until 2026-09-02 and had made
# this script un-runnable since that tree was removed.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import Y_DC_cube3d

from radia.maglev.mixed_galerkin.rom_fit import passive_foster_fit

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
V_CUBE = L**3
K_SIBC = 6 * L**2 * math.sqrt(SIGMA / MU)

F_MIN_HZ = 1.0
F_MAX_HZ = 1e9
N_SAMPLES = 100
N_FILLER = 30


def _tanh_safe(z):
    if abs(z.real) > 50:
        return complex(1.0, 0.0) if z.real > 0 else complex(-1.0, 0.0)
    return cmath.tanh(z)


def _sech2_safe(z):
    if abs(z.real) > 50:
        return complex(0.0, 0.0)
    return 1.0 / cmath.cosh(z)**2


def Y_mixed(s):
    """Closed-form mixed admittance of the cube: rank-1 bulk + exact K_ss."""
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


def sample_band(f_min=F_MIN_HZ, f_max=F_MAX_HZ, n=N_SAMPLES):
    f = np.logspace(math.log10(f_min), math.log10(f_max), n)
    s = 1j * 2 * math.pi * f
    return f, s, np.array([Y_mixed(x) for x in s])


# --------------------------------------------------------------------------
# the realisation: real poles, non-negative residues
# --------------------------------------------------------------------------
def foster_realisation(n_filler=N_FILLER):
    """Passive Foster form Y(s) ~ a_inf + sum_k g_k / (1 + s tau_k) of the
    sampled band.  The dominant poles come from AAA inside passive_foster_fit,
    the diffusive tail from log-spaced filler poles, the residues from NNLS."""
    f, s, Y = sample_band()
    rom = passive_foster_fit(s, Y, n_filler=n_filler)
    tau = np.asarray(rom.tau_n, dtype=float)
    g = np.asarray(rom.g_n, dtype=float)
    a_inf = float(rom.alpha_inf)
    Y_fit = a_inf + np.sum(g[None, :] / (1.0 + s[:, None] * tau[None, :]), axis=1)
    rel = np.abs(Y_fit - Y) / np.abs(Y)
    return {
        "tau_s": tau,
        "g_S": g,
        "feedthrough_S": a_inf,
        "dominant_tau_s": np.asarray(rom.dominant_tau, dtype=float),
        "band_fit_relerr_max_normalised": float(rom.band_fit_relerr),
        "band_max_rel_err_to_100MHz": float(rel[f <= 1e8].max()),
        "band_max_rel_err_to_1GHz": float(rel.max()),
        "f_Hz": f,
    }


def step_response_foster(t_grid, tau, g, a_inf):
    """y(t) = a_inf + sum_k g_k (1 - exp(-t / tau_k)); monotone when g >= 0."""
    t = np.asarray(t_grid, dtype=float)
    return a_inf + np.sum(g[None, :] * (1.0 - np.exp(-t[:, None] / tau[None, :])), axis=1)


# --------------------------------------------------------------------------
# the reference: exact inverse Laplace transform of Y_mixed(s) / s
# --------------------------------------------------------------------------
def exact_step_response(t_grid, dps=20):
    """Talbot inversion of Y_mixed(s)/s (mpmath); the step response the
    realisation has to reproduce.  The integrand is double precision, so
    dps 15 and 20 agree with de Hoog to seven digits while dps 30 adds noise
    (4e-4 at late times), not digits.  Do not raise dps."""
    mp.mp.dps = dps

    def F(s):
        return complex(Y_mixed(complex(s))) / complex(s)

    return np.array([float(mp.re(mp.invertlaplace(F, float(t), method="talbot")))
                     for t in t_grid])


def sqrt_t_asymptote(t_grid, K):
    """Early-time asymptote of the step response of K/sqrt(s): K (2/sqrt(pi)) sqrt(t)."""
    return K * (2 / math.sqrt(math.pi)) * np.sqrt(t_grid)


# --------------------------------------------------------------------------
# the negative result: the AAA poles used directly
# --------------------------------------------------------------------------
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
        E[j + 1, j + 1] = zj[j]
    eigs = np.linalg.eigvals(np.linalg.solve(B + 1e-30 * np.eye(m + 1), E))
    poles = eigs[np.isfinite(eigs)]
    res = []
    with np.errstate(divide="ignore", invalid="ignore"):
        for p in poles:
            num = np.sum(wj * fj / (p - zj))
            den_d = -np.sum(wj / (p - zj)**2)
            res.append(num / den_d if abs(den_d) > 1e-30 else 0.0)
    return np.array(poles), np.array(res)


def step_response_aaa(t_grid, poles, residues):
    """Real part of sum_k r_k/(-p_k) (1 - exp(p_k t)): the conjugate-
    symmetrised realisation of a one-sided AAA fit."""
    t = np.asarray(t_grid, dtype=float)
    terms = residues[None, :] / (-poles[None, :]) * (1.0 - np.exp(poles[None, :] * t[:, None]))
    return np.sum(terms, axis=1).real


def diagnose_aaa_direct(t_grid, y_exact):
    """The AAA poles of the one-sided fit taken as the realisation: how many
    are complex, and how far their step response strays from the exact one."""
    f, s, Y = sample_band()
    s_all = np.concatenate([[0.0 + 0j], s])
    Y_all = np.concatenate([[complex(Y_DC, 0)], Y])
    zj, fj, wj, errvec = aaa_fit(s_all, Y_all, mmax=40, tol=1e-10)
    poles, residues = aaa_to_poles_residues(zj, fj, wj)
    mask = (poles.real < 0) & (np.abs(residues) > 1e-3)
    poles, residues = poles[mask], residues[mask]
    n_complex = int(np.sum(np.abs(poles.imag) > 1e-6 * np.abs(poles.real)))
    y = step_response_aaa(t_grid, poles, residues)
    dev = np.abs(y / y_exact - 1.0)
    early = np.asarray(t_grid) < 1e-6
    return {
        "what": "one-sided AAA fit of Y_mixed on 1 Hz..1 GHz plus DC, stable poles "
                "(Re p < 0, |residue| > 1e-3) summed as r/(-p)(1 - exp(p t)), real part",
        "aaa_degree": int(len(zj)),
        "freq_domain_err_at_samples": float(errvec[-1]),
        "n_stable_poles": int(len(poles)),
        "n_complex_poles": n_complex,
        "max_rel_dev_vs_exact_below_1us": float(dev[early].max()),
        "max_rel_dev_vs_exact_above_10us": float(dev[np.asarray(t_grid) > 1e-5].max()),
        "verdict": "not a usable realisation: complex, non-conjugate poles oscillate",
    }


# --------------------------------------------------------------------------
def summary() -> dict:
    """The numbers this script stands for, as one dict.

    Read by emit_results.py. The talk's step-response figure is drawn from
    the curves recorded here, not from a second computation."""
    rom = foster_realisation()
    tau, g, a_inf = rom["tau_s"], rom["g_S"], rom["feedthrough_S"]
    t_grid = np.logspace(-8, -1, 80)
    y_exact = exact_step_response(t_grid)
    y_rom = step_response_foster(t_grid, tau, g, a_inf)
    y_sqrt = sqrt_t_asymptote(t_grid, K_SIBC)
    dev = np.abs(y_rom / y_exact - 1.0)
    early = int(np.argmin(np.abs(t_grid - 1e-8)))
    return {
        "case": "cube3d_time_domain",
        "body": f"cube, L = {L} m, sigma = {SIGMA:.3g} S/m",
        "metric": ("passive real-pole Foster fit of Y_mixed(s) on 1 Hz..1 GHz; "
                   "its step response against the exact inverse Laplace transform "
                   "(Talbot) of Y_mixed(s)/s on 1e-8..1e-1 s"),
        "realisation": ("passive Foster: real negative poles (dominant from AAA, "
                        "log-spaced filler), non-negative residues by NNLS; "
                        "radia.maglev.mixed_galerkin.rom_fit.passive_foster_fit"),
        "sample_band_Hz": [F_MIN_HZ, F_MAX_HZ],
        "n_samples": N_SAMPLES,
        "n_filler_requested": N_FILLER,
        "n_poles": int(len(tau)),
        "all_poles_real": True,
        "residues_nonnegative": bool(np.all(g >= 0)),
        "Y_DC_S": float(Y_DC),
        "sum_residues_plus_feedthrough_S": float(g.sum() + a_inf),
        "feedthrough_S": float(a_inf),
        "tau_s": tau.tolist(),
        "g_S": g.tolist(),
        "band_fit_relerr_max_normalised": rom["band_fit_relerr_max_normalised"],
        "band_max_rel_err_to_100MHz": rom["band_max_rel_err_to_100MHz"],
        "band_max_rel_err_to_1GHz": rom["band_max_rel_err_to_1GHz"],
        "K_SIBC": float(K_SIBC),
        "early_time_ratio_exact_to_sqrt_t": float(y_exact[early] / y_sqrt[early]),
        "max_rel_dev_rom_vs_exact": float(dev.max()),
        "max_rel_dev_rom_vs_exact_pct": float(100.0 * dev.max()),
        "late_time_value_S": float(y_rom[-1]),
        "step_response": {
            "t_s": t_grid.tolist(),
            "y_exact_S": y_exact.tolist(),
            "y_rom_S": y_rom.tolist(),
            "y_sqrt_asymptote_S": y_sqrt.tolist(),
        },
        "aaa_direct_diagnostic": diagnose_aaa_direct(t_grid, y_exact),
    }


def main():
    r = summary()
    print("=" * 64)
    print("Time-domain mixed Galerkin (cube): summary")
    print("=" * 64)
    print(f"  L = {L*1e3} mm,  sigma = {SIGMA:.2e},  MS = {MS:.4e}")
    print(f"  Y_DC = sigma V_cube = {Y_DC:.4e}")
    print(f"  K_SIBC = 6 L^2 sqrt(sigma/mu) = {K_SIBC:.4e}")
    print()
    print(f"Realisation: {r['n_poles']} real poles, residues >= 0: {r['residues_nonnegative']}, "
          f"feedthrough {r['feedthrough_S']:.3e} S")
    print(f"  sum g + feedthrough = {r['sum_residues_plus_feedthrough_S']:.4f}  (Y_DC {Y_DC:.4f})")
    print(f"  band error: {r['band_fit_relerr_max_normalised']:.1e} (max-normalised), "
          f"{r['band_max_rel_err_to_100MHz']:.1e} per point to 100 MHz, "
          f"{r['band_max_rel_err_to_1GHz']:.1e} per point to 1 GHz")
    print(f"  step response vs exact inverse Laplace: max rel dev {r['max_rel_dev_rom_vs_exact']:.2e}")
    print()
    sr = r["step_response"]
    t_grid = np.array(sr["t_s"])
    y_exact = np.array(sr["y_exact_S"])
    y_rom = np.array(sr["y_rom_S"])
    y_sqrt = np.array(sr["y_sqrt_asymptote_S"])
    print(f"{'t (s)':>11}  {'y_exact':>12}  {'y_rom':>12}  {'2K sqrt(t/pi)':>14}  {'rom/exact':>10}  {'exact/sqrt':>10}")
    for i in range(0, len(t_grid), 8):
        print(f"  {t_grid[i]:9.2e}  {y_exact[i]:12.5e}  {y_rom[i]:12.5e}  {y_sqrt[i]:14.5e}  "
              f"{y_rom[i]/y_exact[i]:10.5f}  {y_exact[i]/y_sqrt[i]:10.4f}")
    print()
    d = r["aaa_direct_diagnostic"]
    print("Negative result kept: the one-sided AAA poles used directly")
    print(f"  {d['n_stable_poles']} stable poles, {d['n_complex_poles']} complex; freq err at samples "
          f"{d['freq_domain_err_at_samples']:.1e}; step response strays {d['max_rel_dev_vs_exact_below_1us']:.1%} "
          f"below 1 us, {d['max_rel_dev_vs_exact_above_10us']:.1e} above 10 us")
    print()
    print("--- Summary ---")
    print(f"Mixed Galerkin time-domain realisation: {r['n_poles']} real poles, DC exact")
    print("  Early time (t < 1e-5): the exact response tracks the sqrt(t) law of the surface term")
    print("  Late time  (t > 1e-3): saturates to Y_DC")
    print("  In between           : a sum of real decaying exponentials, the diffusive form of time-domain SIBC")


if __name__ == "__main__":
    main()
