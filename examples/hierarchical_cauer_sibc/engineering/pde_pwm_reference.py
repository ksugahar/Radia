"""
Independent reference: 1D radial diffusion PDE for cylinder PWM transient.

Method of lines:
  dH_z/dt = D (1/r) d/dr (r dH_z/dr),  D = 1/(mu sigma)
  BC: H_z(a, t) = H_ext(t),  d/dr H_z(0, t) = 0 (symmetry)
  IC: H_z(r, 0) = 0

Discretize on radial grid (N_r points), use scipy.integrate.solve_ivp to
step in time.  This gives an INDEPENDENT FEM-class reference (no Bessel
modal decomposition assumed) for cross-verifying the LTI hierarchical
Cauer of pwm_transient_field.py.

Usage:
    python pde_pwm_reference.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
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


def square_wave(t, f_sw, A=1.0):
    """Square wave starting at +A at t = 0."""
    if np.isscalar(t):
        return A if (t * f_sw) % 1.0 < 0.5 else -A
    return A * np.where((t * f_sw) % 1.0 < 0.5, 1.0, -1.0)


def build_pde_rhs(N_r, a, D):
    """
    Build the right-hand side function for the 1D radial diffusion PDE.

    Grid: r_i = i * dr, i = 0..N_r-1, dr = a/(N_r - 1).
    Interior: dH_i/dt = D * (1/r_i) * d/dr [r dH/dr]_i
                      = D * (r_{i+1/2}(H_{i+1}-H_i) - r_{i-1/2}(H_i - H_{i-1})) / (r_i dr^2)
    Boundary i=0 (r=0): use symmetry, d/dr H = 0, so H_1 = H_-1 (ghost),
                       => dH_0/dt = D * (2/dr^2) * (H_1 - H_0) (limit form)
    Boundary i=N_r-1 (r=a): Dirichlet, H = H_ext(t) (handled in time-stepper)
    """
    dr = a / (N_r - 1)
    r = np.linspace(0, a, N_r)
    # Precompute coefficients
    r_int = r[1:-1]
    r_plus = (r[:-2] + r[1:-1]) / 2 + dr/2   # r_{i+1/2}
    r_minus = (r[:-2] + r[1:-1]) / 2 - dr/2  # r_{i-1/2}
    # Actually let me redo: for i = 1..N_r-2 (interior):
    # r_{i+1/2} = (r[i] + r[i+1])/2 = r[i] + dr/2
    # r_{i-1/2} = (r[i] + r[i-1])/2 = r[i] - dr/2
    r_plus = r_int + dr/2
    r_minus = r_int - dr/2

    def rhs(t, H_interior, H_ext_func):
        """
        H_interior: array of length N_r - 1 (i = 0, 1, ..., N_r - 2).
        H_{N_r-1} = H_ext(t) handled separately.
        """
        H = np.empty(N_r)
        H[:-1] = H_interior
        H[-1] = H_ext_func(t)

        dH = np.zeros(N_r - 1)
        # i = 0 (axis): symmetry boundary
        dH[0] = D * 2 * (H[1] - H[0]) / dr**2
        # Interior: i = 1, ..., N_r - 2
        for i in range(1, N_r - 1):
            r_p = r[i] + dr/2
            r_m = r[i] - dr/2
            dH[i] = D * (r_p * (H[i+1] - H[i]) - r_m * (H[i] - H[i-1])) / (r[i] * dr**2)
        return dH

    return rhs, r


def cylinder_pde_pwm(f_sw, H_0, t_end, N_r, a, sigma, mu):
    """
    Solve the cylinder PWM diffusion PDE.
    Returns: t_eval (n_t,), r (N_r,), H_field (n_t, N_r).
    """
    D = 1.0 / (mu * sigma)
    H_ext_func = lambda t: square_wave(t, f_sw, H_0)
    rhs, r = build_pde_rhs(N_r, a, D)

    # Initial condition: H = 0 inside, H_ext(0+) at boundary
    H0_interior = np.zeros(N_r - 1)

    print(f"  Time integration (PDE: {N_r-1} radial DOFs)...")
    import time
    t0 = time.time()
    sol = solve_ivp(rhs, [0, t_end], H0_interior,
                    args=(H_ext_func,),
                    method='LSODA', rtol=1e-7, atol=1e-10,
                    t_eval=np.linspace(0, t_end, 400))
    print(f"  Done in {time.time()-t0:.1f}s, {sol.t.size} steps")

    # Assemble full field (add boundary value back)
    H_field = np.zeros((len(sol.t), N_r))
    H_field[:, :-1] = sol.y.T
    H_field[:, -1] = np.array([H_ext_func(ti) for ti in sol.t])

    return sol.t, r, H_field


def lti_field(t_arr, r, f_sw, H_0, N_modes, a, sigma, mu):
    """LTI hierarchical Cauer field for comparison."""
    from scipy import special
    j0 = special.jn_zeros(0, N_modes)
    j1v = special.j1(j0)
    xi = j0**2 / (a**2 * mu * sigma)
    c = (4 * np.pi / mu) * np.ones(N_modes)
    A = 2.0 / (j0 * j1v)

    H_ext_func = lambda t: square_wave(t, f_sw, H_0)

    def rhs(t, z):
        return -xi * z + c * H_ext_func(t)

    sol = solve_ivp(rhs, [0, t_arr[-1]], np.zeros(N_modes),
                    t_eval=t_arr, method='RK45', rtol=1e-8, atol=1e-12)
    z = sol.y
    H_input_arr = np.array([H_ext_func(ti) for ti in t_arr])

    # Modal amplitudes
    m = np.zeros_like(z)
    for k in range(N_modes):
        m[k] = A[k] * (H_input_arr - xi[k] * z[k] / c[k])

    # Field reconstruction
    arg = j0[None, :] * r[:, None] / a
    J0_arr = special.j0(arg)
    H_field = np.zeros((len(t_arr), len(r)))
    for ti in range(len(t_arr)):
        H_field[ti] = H_input_arr[ti] - np.sum(m[:, ti][None, :] * J0_arr,
                                                axis=1)
    return H_field


if __name__ == "__main__":
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    f_sw = 30e3
    T_sw = 1.0 / f_sw
    n_periods = 3
    t_end = n_periods * T_sw

    print(f"Cylinder a = {a*1e3:.1f} mm Cu")
    print(f"PWM: {f_sw/1e3:.0f} kHz, {n_periods} periods, t_end = {t_end*1e6:.1f} us")
    print()

    # Solve PDE
    N_r = 100
    print("Solving 1D radial diffusion PDE (FD method of lines, LSODA)...")
    t_pde, r_pde, H_pde = cylinder_pde_pwm(f_sw, 1.0, t_end, N_r,
                                             a, sigma_cu, mu_0)
    print()

    # Solve LTI
    print("Solving LTI hierarchical Cauer (30 Foster modes)...")
    H_lti = lti_field(t_pde, r_pde, f_sw, 1.0, 30, a, sigma_cu, mu_0)
    print()

    # Compare at selected snapshots
    t_snaps_idx = [
        np.argmin(np.abs(t_pde - 0.25 * T_sw)),
        np.argmin(np.abs(t_pde - 0.75 * T_sw)),
        np.argmin(np.abs(t_pde - 1.25 * T_sw)),
        np.argmin(np.abs(t_pde - 2.75 * T_sw)),
    ]
    print("Cross-verification at snapshots:")
    for ti in t_snaps_idx:
        err = np.max(np.abs(H_pde[ti] - H_lti[ti]))
        print(f"  t = {t_pde[ti]*1e6:6.2f} us "
              f"(t/T_sw = {t_pde[ti]/T_sw:.2f}): "
              f"max |PDE - LTI| = {err:.2e}")
    print()

    # Plot
    setup_times_new_roman()
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.0))
    for sj, ti in enumerate(t_snaps_idx):
        ax = axes[sj // 2, sj % 2]
        ax.plot(r_pde * 1e3, H_pde[ti], 'k-', linewidth=1.4,
                label='PDE (FD, LSODA)')
        ax.plot(r_pde * 1e3, H_lti[ti], 'r--', linewidth=1.0,
                label='LTI (30 Foster + RK45)')
        ax.set_xlabel('radial position $r$ (mm)')
        ax.set_ylabel(r'$H_z(r, t) / H_0$')
        ax.set_title(f'$t$ = {t_pde[ti]*1e6:.1f} $\\mu$s '
                     f'($t/T_{{sw}}$ = {t_pde[ti]/T_sw:.2f})',
                     fontsize=9)
        ax.legend(fontsize=8, frameon=False, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-1.1, 1.1)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "pde_pwm_reference.png", dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / 'pde_pwm_reference.png'}")

    np.savez(out_dir / "pde_pwm_reference_data.npz",
             t=t_pde, r=r_pde, H_PDE=H_pde, H_LTI=H_lti,
             a=a, sigma=sigma_cu, mu=mu_0, f_sw=f_sw, N_r=N_r)
    print(f"Saved: pde_pwm_reference_data.npz")
