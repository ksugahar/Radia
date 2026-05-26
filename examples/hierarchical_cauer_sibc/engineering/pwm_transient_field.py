"""
PWM switching transient with simultaneous field reconstruction.

Driving the hierarchical Cauer LTI (Foster expansion of cylinder admittance)
with a realistic engineering input (square wave 30 kHz, PWM-like) and
computing BOTH the port-level current i(t) AND the spatial field
H_z(r, t), J_phi(r, t) by tracking modal amplitudes z_k(t).

Mathematical setup:
  Input: H_ext(t) (the external field; for SPICE analogy this is the
                  port voltage driving the Foster admittance Y)
  State: z_k(t) = bulk Foster branch current for mode k
                  dz_k/dt = -xi_k z_k + c_k H_ext(t),  z_k(0) = 0
  with xi_k = j_{0,k}^2/(a^2 mu sigma), c_k = 4 pi / mu.

Port output:
  i_port(t) = sum_k z_k(t)   (total induced eddy current admittance response)

Modal amplitudes for field reconstruction:
  m_k(t) = A_k * (H_ext(t) - xi_k z_k(t) / c_k)
  where A_k = 2 / (j_{0,k} J_1(j_{0,k}))   (Fourier-Bessel coefficient of 1)

Field at any point r in (0, a):
  H_z(r, t) = H_ext(t) - sum_k m_k(t) * J_0(j_{0,k} r/a)
  J_phi(r, t) = -d H_z/dr = sum_k m_k(t) * (j_{0,k}/a) * J_1(j_{0,k} r/a)

Reference: convolve the analytical step response with the input via
linearity (square wave = sum of step pulses).

Usage:
    python pwm_transient_field.py
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


def cylinder_modal_data(N_modes, a, sigma, mu):
    """Foster modal data for the cylinder."""
    j0_zeros = special.jn_zeros(0, N_modes)
    j1_vals = special.j1(j0_zeros)
    xi = j0_zeros**2 / (a**2 * mu * sigma)
    c = (4 * np.pi / mu) * np.ones(N_modes)
    A = 2.0 / (j0_zeros * j1_vals)
    return j0_zeros, j1_vals, xi, c, A


def square_wave(t, f_sw, amplitude=1.0):
    """Square wave H_ext(t) at frequency f_sw, amplitude H_0."""
    phase = (t * f_sw) % 1.0
    return amplitude * np.where(phase < 0.5, 1.0, -1.0)


def lti_rhs(t, z, xi, c, H_ext_func):
    """RHS of state-space ODE: dz_k/dt = -xi_k z_k + c_k H_ext(t)."""
    H = H_ext_func(t)
    return -xi * z + c * H


def step_response_modal(t_arr, xi, A):
    """Modal amplitudes for unit step at t=0:
       m_k(t) = A_k * exp(-xi_k * t) for t > 0, 0 for t <= 0.
    """
    m = np.zeros((len(t_arr), len(xi)))
    for ti, t in enumerate(t_arr):
        if t > 0:
            m[ti] = A * np.exp(-xi * t)
    return m


def step_response_field(r, t, N_modes, a, sigma, mu):
    """Step response H_z(r, t)/H_0 (analytical reference)."""
    j0, j1v, xi, c, A = cylinder_modal_data(N_modes, a, sigma, mu)
    if t <= 0:
        return np.zeros_like(r)
    arg = j0[None, :] * r[:, None] / a
    J0_arr = special.j0(arg)
    decay = np.exp(-xi * t)
    return 1.0 - np.sum(A[None, :] * J0_arr * decay[None, :], axis=1)


def square_wave_field_by_superposition(r, t, edges, polarities, H_0,
                                        N_modes, a, sigma, mu):
    """
    Field response to square wave at time t, by superposition of step responses
    at each switching edge.

    edges[i]: time of i-th edge
    polarities[i]: amplitude change at edge i (+/- 2*H_0 typically)

    H_z(r, t) = sum_{i: edges[i] < t} polarities[i] * step_response(r, t - edges[i])
    """
    H = np.zeros_like(r)
    for edge, pol in zip(edges, polarities):
        if edge < t:
            H += pol * step_response_field(r, t - edge, N_modes, a, sigma, mu)
    return H


if __name__ == "__main__":
    # Cu cylinder
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    d = 2 * np.pi * 7e4
    t_d = 1.0 / d

    # Modal setup
    N_bulk = 30
    j0, j1v, xi, c, A = cylinder_modal_data(N_bulk, a, sigma_cu, mu_0)
    print(f"Cylinder a = {a*1e3:.1f} mm Cu")
    print(f"Bulk Foster modes: N = {N_bulk}")
    print(f"  xi_1  = {xi[0]:.3e} rad/s (lowest, slowest decay)")
    print(f"  xi_N  = {xi[-1]:.3e} rad/s (highest, fastest decay)")
    print(f"  wall-band 1/d = {t_d:.3e} s")
    print()

    # PWM input parameters
    H_0 = 1.0  # amplitude
    f_sw = 30e3  # 30 kHz switching frequency (typical IGBT)
    T_sw = 1.0 / f_sw  # 33.3 us period
    n_periods = 5
    t_end = n_periods * T_sw
    print(f"PWM input: square wave, f_sw = {f_sw/1e3:.1f} kHz, "
          f"T_sw = {T_sw*1e6:.1f} us")
    print(f"Simulation: {n_periods} periods, t_end = {t_end*1e6:.1f} us")
    print()

    # Define input
    def H_ext(t):
        return square_wave(t, f_sw, H_0)

    # Compute switching edges and polarity changes.
    # System starts from rest (H_ext = 0 for t < 0, z(0) = 0); first edge at
    # t = 0+ injects the initial H_ext value, then subsequent edges flip the
    # square-wave polarity.
    edges = []
    polarities = []
    H_prev = 0.0  # rest before t = 0
    n_edges = int(2 * t_end * f_sw) + 2
    for i in range(n_edges):
        edge = i * T_sw / 2
        if edge > t_end:
            break
        H_after = H_ext(edge + 1e-12)
        polarities.append(H_after - H_prev)
        edges.append(edge)
        H_prev = H_after
    edges = np.array(edges)
    polarities = np.array(polarities)
    print(f"Switching edges: {len(edges)}")
    print(f"  polarities: {polarities}")
    print()

    # ----------------------------------------------------------------------
    # Integrate LTI state-space: dz_k/dt = -xi_k z_k + c_k H_ext(t)
    # ----------------------------------------------------------------------
    print("Integrating LTI state-space (RK45, rtol=1e-8)...")
    t_eval = np.linspace(0, t_end, 4000)
    sol = solve_ivp(
        lti_rhs, [0, t_end], np.zeros(N_bulk),
        t_eval=t_eval, args=(xi, c, H_ext),
        method='RK45', rtol=1e-8, atol=1e-12
    )
    z = sol.y  # shape (N_bulk, n_t)
    print(f"  Integration done.  status = {sol.status}, n_steps = {sol.t.size}")
    print()

    # Port-level eddy current response
    i_port = np.sum(z, axis=0)
    # Note: i_port has units of (4pi/mu)/(unit xi) -- it's effectively the
    # net induced current density integrated over the cylinder cross-section.

    # Modal amplitudes m_k(t) = A_k (H_ext - xi_k z_k / c_k)
    H_input_arr = H_ext(t_eval)
    m = np.zeros_like(z)
    for k in range(N_bulk):
        m[k] = A[k] * (H_input_arr - xi[k] * z[k] / c[k])

    # ----------------------------------------------------------------------
    # Verify against analytical superposition
    # ----------------------------------------------------------------------
    # Port current can also be reconstructed analytically as sum of step responses
    # i.e., i_port(t) = sum over edges (polarity * Y_step(t-edge))
    # where Y_step(t) = sum_k (c_k/xi_k) (1 - exp(-xi_k t)) for t > 0
    Y_step = lambda tt: (np.sum((c/xi) * (1 - np.exp(-xi * tt)))
                          if tt > 0 else 0.0)
    i_port_analytical = np.array([
        sum(pol * Y_step(t - edge) for edge, pol in zip(edges, polarities)
            if edge < t)
        for t in t_eval
    ])

    port_err = np.max(np.abs(i_port - i_port_analytical))
    port_max = np.max(np.abs(i_port_analytical))
    print(f"Port-current verification:")
    print(f"  max |LTI - analytical|   = {port_err:.2e}")
    print(f"  max |analytical|         = {port_max:.2e}")
    print(f"  relative error           = {port_err/port_max:.2e}")
    print()

    # ----------------------------------------------------------------------
    # Field reconstruction at chosen snapshot times
    # ----------------------------------------------------------------------
    r = np.linspace(0, a, 150)
    # Snapshots inside half-periods (away from switching edges to avoid
    # boundary ambiguity between LTI and superposition reference)
    t_snap = np.array([0.25, 0.75, 1.25, 1.75, 2.75, 4.25]) * T_sw

    H_lti = np.zeros((len(t_snap), len(r)))
    H_ref = np.zeros((len(t_snap), len(r)))
    for si, t_s in enumerate(t_snap):
        # From LTI: find closest t_eval index, evaluate m_k at that time
        idx = np.argmin(np.abs(t_eval - t_s))
        m_at_t = m[:, idx]
        # Reconstruct H_z(r) = H_ext(t_s) - sum_k m_k J_0(j_{0,k} r/a)
        arg = j0[None, :] * r[:, None] / a
        J0_arr = special.j0(arg)
        H_lti[si] = H_ext(t_s) - np.sum(m_at_t[None, :] * J0_arr, axis=1)
        # Reference: square-wave superposition
        H_ref[si] = square_wave_field_by_superposition(
            r, t_s, edges, polarities, H_0, 200, a, sigma_cu, mu_0)

    # ----------------------------------------------------------------------
    # Plot
    # ----------------------------------------------------------------------
    setup_times_new_roman()
    fig = plt.figure(figsize=(8.5, 7.5))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.30)

    # (a) Input H_ext + port current i_port (top-left, wide)
    ax_input = fig.add_subplot(gs[0, :])
    ax_input2 = ax_input.twinx()
    ax_input.plot(t_eval * 1e6, H_input_arr, 'k-', linewidth=0.8,
                  label=r'$H_{\mathrm{ext}}(t)$ (PWM)')
    ax_input2.plot(t_eval * 1e6, i_port, 'r-', linewidth=1.0,
                   label=r'$i_{\mathrm{port}}(t)$ (LTI)')
    ax_input2.plot(t_eval * 1e6, i_port_analytical, 'b--', linewidth=0.5,
                   label='analytical superposition')
    for t_s in t_snap:
        ax_input.axvline(t_s * 1e6, color='gray', linewidth=0.4,
                         linestyle=':', alpha=0.6)
    ax_input.set_xlabel('time $t$ ($\\mu$s)')
    ax_input.set_ylabel(r'$H_{\mathrm{ext}}$ (A/m)', color='k')
    ax_input2.set_ylabel(r'$i_{\mathrm{port}}$ (rel.)', color='r')
    ax_input.tick_params(direction='in', axis='y', labelcolor='k')
    ax_input2.tick_params(direction='in', axis='y', labelcolor='r')
    ax_input.legend(loc='upper left', fontsize=8, frameon=False)
    ax_input2.legend(loc='upper right', fontsize=8, frameon=False)
    ax_input.set_title(f'(a) PWM input ({f_sw/1e3:.0f} kHz) and port response',
                       fontsize=9)
    for spine in ax_input.spines.values():
        spine.set_linewidth(0.8)

    # (b)-(d) Field snapshots at 3 selected times
    # use first 3 of t_snap for snapshots
    snap_pos = [(1, 0), (1, 1), (2, 0), (2, 1)]
    for si, (row, col) in enumerate(snap_pos):
        if si >= 4:
            break
        ax = fig.add_subplot(gs[row, col])
        idx = si if si < len(t_snap) else len(t_snap) - 1
        ax.plot(r * 1e3, H_lti[idx], 'r--', linewidth=1.2,
                label=f'LTI reconstruction ($N$={N_bulk})')
        ax.plot(r * 1e3, H_ref[idx], 'k-', linewidth=0.8, alpha=0.8,
                label='analytical superposition')
        ax.axhline(H_ext(t_snap[idx]), color='gray', linestyle=':',
                   linewidth=0.6, alpha=0.6,
                   label=f'$H_{{\\mathrm{{ext}}}}={H_ext(t_snap[idx]):+.0f}$')
        ax.set_xlabel('radial position $r$ (mm)')
        ax.set_ylabel(r'$H_z(r, t)/H_0$')
        ax.set_title(f'$t = {t_snap[idx]*1e6:.1f}$ $\\mu$s '
                     f'($t/T_{{sw}}={t_snap[idx]/T_sw:.2f}$)',
                     fontsize=9)
        ax.legend(fontsize=7, frameon=False, loc='center left')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-1.1, 1.1)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    out_path = out_dir / "pwm_transient_field.png"
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {out_path}")
    print()

    # Verification at each snapshot
    print("Field reconstruction error at snapshots:")
    for si, t_s in enumerate(t_snap):
        err = np.max(np.abs(H_lti[si] - H_ref[si]))
        print(f"  t = {t_s*1e6:5.1f} us (t/T_sw = {t_s/T_sw:.2f}): "
              f"max |LTI - ref| = {err:.2e}")
    print()

    # Skin penetration evolution (depth where H drops below 1 - 1/e of surface value)
    # for a steady-state PWM cycle (last cycle)
    np.savez(out_dir / "pwm_transient_field_data.npz",
             t=t_eval, H_input=H_input_arr, i_port=i_port,
             i_port_analytical=i_port_analytical,
             z=z, m=m,
             r=r, t_snap=t_snap, H_lti=H_lti, H_ref=H_ref,
             a=a, sigma=sigma_cu, mu=mu_0, f_sw=f_sw, N_bulk=N_bulk)
    print(f"Saved: {out_dir / 'pwm_transient_field_data.npz'}")
