"""
Realistic driving inputs for hierarchical Cauer transient analysis.

Three engineering-relevant input scenarios:
  (a) Frequency-swept sinusoid (chirp) covering the wall-band crossover d
  (b) PWM with dead time (gate-driver realistic with blanking)
  (c) Resonant transient: damped sinusoid (RLC ringing)

All three drive the cylinder's hierarchical Cauer LTI; we compare port-level
response across them and against analytical reference where available.

Usage:
    python realistic_inputs.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
from scipy import special
from scipy.integrate import solve_ivp


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


def cyl_modal(N, a, sigma, mu):
    j0 = special.jn_zeros(0, N)
    j1 = special.j1(j0)
    xi = j0**2 / (a**2 * mu * sigma)
    c = (4 * np.pi / mu) * np.ones(N)
    A = 2.0 / (j0 * j1)
    return j0, j1, xi, c, A


def integrate_cylinder(H_ext_func, t_end, N_modes, a, sigma, mu, n_eval=4000):
    """Integrate Foster bulk LTI; return (t, z, m, i_port)."""
    j0, j1, xi, c, A = cyl_modal(N_modes, a, sigma, mu)
    t_eval = np.linspace(0, t_end, n_eval)

    def rhs(t, z):
        return -xi * z + c * H_ext_func(t)

    sol = solve_ivp(rhs, [0, t_end], np.zeros(N_modes),
                    t_eval=t_eval, method='RK45',
                    rtol=1e-8, atol=1e-12)
    z = sol.y
    H_arr = np.array([H_ext_func(t) for t in t_eval])
    # Modal amplitudes
    m = np.zeros_like(z)
    for k in range(N_modes):
        m[k] = A[k] * (H_arr - xi[k] * z[k] / c[k])
    i_port = np.sum(z, axis=0)
    return t_eval, z, m, i_port, H_arr


# ----------------------------------------------------------------------
# Input signal builders
# ----------------------------------------------------------------------
def chirp_signal(t, f0, f1, T, amplitude=1.0):
    """Linear chirp from f0 to f1 over duration T."""
    k = (f1 - f0) / T
    phase = 2 * np.pi * (f0 * t + 0.5 * k * t**2)
    return amplitude * np.sin(phase)


def pwm_dead_time(t, f_sw, duty=0.5, dead_time=0.05, amplitude=1.0):
    """
    PWM with dead time: rectangular wave at f_sw, but a fraction
    `dead_time` of the period is held at zero around each switching edge.

    During dead time, the output is 0 (both gate drivers off).
    """
    T = 1.0 / f_sw
    phase = (t * f_sw) % 1.0
    # On phase: [dead_time/2, duty - dead_time/2]
    # Off (dead) zones: [0, dead_time/2] and [duty - dead_time/2, duty + dead_time/2]
    # ...
    # Simpler: dead bands at [phase < dt/2] and [|phase - duty| < dt/2] and [phase > 1 - dt/2]
    dt = dead_time / 2
    in_high = (phase > dt) & (phase < duty - dt)
    in_low = (phase > duty + dt) & (phase < 1.0 - dt)
    return np.where(in_high, amplitude,
                    np.where(in_low, -amplitude, 0.0))


def damped_ringing(t, f_ring, decay_const, amplitude=1.0):
    """Damped sinusoidal ringing: A * exp(-t/tau) * sin(2 pi f_ring t)."""
    return amplitude * np.exp(-t / decay_const) * np.sin(2 * np.pi * f_ring * t)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    d = 2 * np.pi * 7e4
    t_d = 1.0 / d
    N_modes = 30
    print(f"Cylinder a = {a*1e3:.1f} mm Cu, wall-band t_d = {t_d*1e6:.2f} us")
    print(f"Foster bulk modes: N = {N_modes}")
    print()

    # ------------------------------------------------------------------
    # Scenario (a): Chirp from 100 Hz to 1 MHz
    # ------------------------------------------------------------------
    print("=" * 50)
    print("(a) Chirp 100 Hz -> 1 MHz over 500 us")
    print("=" * 50)
    f0, f1, T_chirp = 100, 1e6, 500e-6
    H_chirp = lambda t: chirp_signal(t, f0, f1, T_chirp)
    t_c, z_c, m_c, ip_c, Hi_c = integrate_cylinder(
        H_chirp, T_chirp, N_modes, a, sigma_cu, mu_0)
    print(f"  Integrated {len(t_c)} samples")
    print(f"  Port current max = {np.max(np.abs(ip_c)):.3e}")
    print()

    # ------------------------------------------------------------------
    # Scenario (b): PWM with dead time
    # ------------------------------------------------------------------
    print("=" * 50)
    print("(b) PWM 30 kHz with 5% dead time")
    print("=" * 50)
    f_sw = 30e3
    dead = 0.05
    n_periods = 4
    t_end_pwm = n_periods / f_sw
    H_pwm = lambda t: pwm_dead_time(t, f_sw, duty=0.5, dead_time=dead)
    t_p, z_p, m_p, ip_p, Hi_p = integrate_cylinder(
        H_pwm, t_end_pwm, N_modes, a, sigma_cu, mu_0)
    print(f"  Integrated {len(t_p)} samples")
    print(f"  Port current max = {np.max(np.abs(ip_p)):.3e}")
    print()

    # ------------------------------------------------------------------
    # Scenario (c): Damped sinusoidal ringing (RLC transient)
    # ------------------------------------------------------------------
    print("=" * 50)
    print("(c) Damped ringing at 50 kHz, tau = 100 us")
    print("=" * 50)
    f_ring = 50e3
    tau_ring = 100e-6
    t_end_ring = 5 * tau_ring
    H_ring = lambda t: damped_ringing(t, f_ring, tau_ring)
    t_r, z_r, m_r, ip_r, Hi_r = integrate_cylinder(
        H_ring, t_end_ring, N_modes, a, sigma_cu, mu_0)
    print(f"  Integrated {len(t_r)} samples")
    print(f"  Port current max = {np.max(np.abs(ip_r)):.3e}")
    print()

    # ------------------------------------------------------------------
    # Plot all three
    # ------------------------------------------------------------------
    setup_times_new_roman()
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 7.5))

    panels = [
        ('(a) Chirp 100 Hz -> 1 MHz', t_c * 1e6, Hi_c, ip_c, t_c[-1] * 1e6),
        ('(b) PWM 30 kHz + 5\\% dead time', t_p * 1e6, Hi_p, ip_p, t_p[-1] * 1e6),
        ('(c) Damped ringing 50 kHz', t_r * 1e6, Hi_r, ip_r, t_r[-1] * 1e6),
    ]
    for row, (title, t_us, H, ip, t_max) in enumerate(panels):
        ax_input = axes[row, 0]
        ax_resp = axes[row, 1]
        ax_input.plot(t_us, H, 'k-', linewidth=0.7)
        ax_input.set_xlabel('time ($\\mu$s)')
        ax_input.set_ylabel(r'$H_{\mathrm{ext}}(t)$')
        ax_input.set_title(title + ' -- input', fontsize=9)
        ax_input.grid(True, alpha=0.3)
        ax_input.set_xlim(0, t_max)
        for spine in ax_input.spines.values():
            spine.set_linewidth(0.8)
        ax_input.tick_params(direction='in', which='both')

        ax_resp.plot(t_us, ip, 'r-', linewidth=0.7)
        ax_resp.set_xlabel('time ($\\mu$s)')
        ax_resp.set_ylabel(r'$i_{\mathrm{port}}(t)$ (rel.)')
        ax_resp.set_title(title + ' -- response', fontsize=9)
        ax_resp.grid(True, alpha=0.3)
        ax_resp.set_xlim(0, t_max)
        for spine in ax_resp.spines.values():
            spine.set_linewidth(0.8)
        ax_resp.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    out_path = out_dir / "realistic_inputs.png"
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {out_path}")

    np.savez(out_dir / "realistic_inputs_data.npz",
             chirp_t=t_c, chirp_H=Hi_c, chirp_i=ip_c,
             pwm_t=t_p, pwm_H=Hi_p, pwm_i=ip_p,
             ring_t=t_r, ring_H=Hi_r, ring_i=ip_r,
             N_modes=N_modes, a=a, sigma=sigma_cu, mu=mu_0)
    print(f"Saved: realistic_inputs_data.npz")
