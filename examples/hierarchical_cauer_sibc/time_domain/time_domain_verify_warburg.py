"""
Time-domain verification of the Warburg block: K_SIBC * sqrt(s) / (s + d).

This script isolates the verification of the diffusive Foster approximation
of the irrational Warburg block, comparing its step response (computed via
finite-dimensional LTI state-space) against the closed-form analytical
inverse Laplace transform.

Analytical impulse response of f(s) = sqrt(s)/(s+d):
    h(t) = L^{-1}{ sqrt(s) / (s + d) } = 1 / sqrt(pi t) - sqrt(d) * exp(d t) * erfc(sqrt(d t))

Step response (integral of h):
    H(t) = integral_0^t h(tau) dtau
         = 2 sqrt(t / pi) + (1/sqrt(d)) * [exp(d t) * erfc(sqrt(d t)) - 1]

We verify: the finite-dim diffusive Foster approximation gives a step response
that matches H(t) closely across many decades of time.

Usage:
    python time_domain_verify_warburg.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
from scipy import special


def setup_times_new_roman():
    try:
        for font in ['Times New Roman', 'Times']:
            try:
                font_manager.findfont(font, fallback_to_default=False)
                plt.rcParams['font.family'] = font
                plt.rcParams['font.serif'] = [font]
                plt.rcParams['mathtext.fontset'] = 'stix'
                return
            except Exception:
                continue
    except Exception:
        pass


def warburg_step_analytical(t, d):
    """
    Closed-form step response of f(s) = sqrt(s)/(s+d).

    Derivation: f(s) = sqrt(s)/(s+d). Partial fraction in sqrt(s):
        f(s) = (1/2) [1/(sqrt(s) - i sqrt(d)) + 1/(sqrt(s) + i sqrt(d))]
    Using L^{-1}{1/(sqrt(s) + a)} = 1/sqrt(pi t) - a e^{a^2 t} erfc(a sqrt(t))
    and combining the complex-conjugate pair gives the real impulse response
        h(t) = 1/sqrt(pi t) - sqrt(d) e^{-d t} erfi(sqrt(d t))
    The step response is H(t) = ∫_0^t h(tau) dtau, which evaluates to
        H(t) = e^{-d t} erfi(sqrt(d t)) / sqrt(d)
             = (2 / sqrt(pi d)) * dawsn(sqrt(d t))
    using the Dawson function dawsn(x) = e^{-x^2} ∫_0^x e^{u^2} du =
    (sqrt(pi)/2) e^{-x^2} erfi(x).

    Limits:
        t → 0+: dawsn(x) ~ x   =>  H(t) ~ 2 sqrt(t/pi)
        t → inf: dawsn(x) ~ 1/(2x)  =>  H(t) ~ 1/sqrt(pi d t)
        H(0) = 0, H(inf) = 0     (matches final-value theorem on f)
    """
    t = np.asarray(t)
    x = np.sqrt(d * t)
    H = (2.0 / np.sqrt(np.pi * d)) * special.dawsn(x)
    return H


def warburg_step_diffusive(t, d, xi, c):
    """
    Step response of the diffusive Foster approximation:
        f(s) ~= sum_k c_k / (s + xi_k)
    Step response (integral of impulse response):
        H_approx(t) = sum_k c_k * (1 - exp(-xi_k t)) / xi_k

    xi, c may be real or complex; output is real (sum of conjugate pairs).
    """
    t = np.asarray(t)
    H = np.zeros_like(t)
    for k in range(len(xi)):
        H += np.real(c[k] / xi[k] * (1 - np.exp(-xi[k] * t)))
    return H


if __name__ == "__main__":
    # Load Warburg diffusive fit
    fit_path = Path(__file__).parent / "warburg_diffusive_fit.npz"
    if not fit_path.exists():
        print(f"ERROR: run diffusive_quadrature.py first to produce {fit_path}")
        exit(1)
    data = np.load(fit_path)
    d = float(data['d'])
    xi = data['xi']
    c = data['c']
    N_xi = len(xi)
    freq_err = float(data['max_rel_err'])
    print(f"Loaded Warburg fit:")
    print(f"  d = {d:.3e} rad/s ({d/(2*np.pi):.2e} Hz)")
    print(f"  N_xi = {N_xi}")
    print(f"  freq-domain max rel err = {freq_err:.2e}")
    print()

    # Time grid: cover many decades
    # Characteristic time t_d = 1/d ~ 2.3e-6 s for d ~ 4.4e5
    # Cover 1e-9 to 1e-3 (6 decades)
    t = np.logspace(-9, -2, 200)

    # Compute step responses
    H_analytical = warburg_step_analytical(t, d)
    H_diffusive = warburg_step_diffusive(t, d, xi, c)

    # Error metrics
    rel_err = np.abs(H_diffusive - H_analytical) / np.maximum(
        np.abs(H_analytical), np.max(np.abs(H_analytical)) * 1e-12
    )
    max_rel_err = np.max(rel_err)
    median_err = np.median(rel_err)
    print(f"Time-domain step response:")
    print(f"  max rel err = {max_rel_err:.2e}")
    print(f"  median rel err = {median_err:.2e}")
    print()

    # Time-axis characteristic markers
    t_d = 1.0 / d
    print(f"  Wall-band time t_d = 1/d = {t_d:.2e} s")
    print(f"  H_analytical(t_d) = {warburg_step_analytical(t_d, d):.3e}")
    print()

    # Asymptotic checks: t -> 0+: H ~ 2 sqrt(t/pi); peak near t ~ 1/d; t -> inf: H -> 0
    H_peak_idx = np.argmax(H_analytical)
    print(f"  H_analytical peak: {H_analytical[H_peak_idx]:.3e} at t = {t[H_peak_idx]:.2e}")
    print(f"  H_analytical(t_max={t[-1]:.2e}) = {H_analytical[-1]:.3e}")
    print(f"  H_diffusive (t_max) = {H_diffusive[-1]:.3e}")
    print()

    # Plot
    setup_times_new_roman()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 3.0))

    ax1.semilogx(t, H_analytical, 'k-', linewidth=1.4,
                 label='analytical (Dawson)')
    ax1.semilogx(t, H_diffusive, 'r--', linewidth=1.0,
                 label=f'diffusive Foster ($N_\\xi={N_xi}$)')
    ax1.axvline(t_d, color='gray', linestyle=':', linewidth=0.6,
                label='$t = 1/d$')
    ax1.set_xlabel('time $t$ (s)')
    ax1.set_ylabel('step response $H(t)$')
    ax1.legend(fontsize=8, frameon=False, loc='upper right')
    ax1.grid(True, which='both', alpha=0.3)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    ax2.loglog(t, rel_err, 'b-', linewidth=1.0)
    ax2.axvline(t_d, color='gray', linestyle=':', linewidth=0.6)
    ax2.axhline(1e-2, color='gray', linestyle='--', linewidth=0.6,
                label='1% engineering')
    ax2.set_xlabel('time $t$ (s)')
    ax2.set_ylabel('relative error on $H(t)$')
    ax2.legend(fontsize=8, frameon=False, loc='upper left')
    ax2.grid(True, which='both', alpha=0.3)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "time_domain_warburg_step.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / 'time_domain_warburg_step.png'}")

    # Save data
    np.savez(out_dir / "warburg_step_data.npz",
             t=t, H_analytical=H_analytical, H_diffusive=H_diffusive,
             rel_err=rel_err, max_rel_err=max_rel_err,
             N_xi=N_xi, d=d)
    print(f"Saved: {out_dir / 'warburg_step_data.npz'}")
