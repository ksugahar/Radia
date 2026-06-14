"""export.py -- Export conductor Y(s) Foster ladder as MATLAB / Simulink LTI.

Foster representation:
    Y(s) = sum_n g_n / (1 + s tau_n) + K_SIBC/sqrt(s) + c_1/s

State-space realization with one state per Foster pole + n_warburg
diffusive-quadrature states + one integrator state for c_1/s.  The
resulting (A, B, C, D) can be loaded into MATLAB with `ss(A,B,C,D)`
and dropped into Simulink as a "State-Space" or "LTI System" block.
"""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np


def diffusive_quadrature(K_SIBC, n_aux=30, f_lo=1e-3, f_hi=1e10):
    """Approximate K_SIBC/sqrt(s) by n_aux first-order RC rungs.

    Identity: 1/sqrt(s) = (1/pi) * integral_0^inf 1/(sqrt(xi) (s + xi)) dxi.

    Logarithmic quadrature: xi_k log-spaced over [omega_lo, omega_hi]:
        K_SIBC/sqrt(s) ~ sum_k r_k / (s + xi_k),
        r_k = K_SIBC * sqrt(xi_k) * d(ln xi) / pi.

    Returns
    -------
    xi_arr : ndarray (n_aux,)
        Pole frequencies (rad/s).
    r_arr : ndarray (n_aux,)
        Per-pole residues.
    """
    omega_lo = 2 * math.pi * f_lo
    omega_hi = 2 * math.pi * f_hi
    xi_arr = np.logspace(math.log10(omega_lo), math.log10(omega_hi), n_aux)
    d_log = math.log(xi_arr[1] / xi_arr[0])
    r_arr = K_SIBC * np.sqrt(xi_arr) * d_log / math.pi
    return xi_arr, r_arr


def build_state_space(g_n, tau_n, V, sigma, K_SIBC, c1, n_warburg_rungs=30):
    """Build (A, B, C, D) state-space for the conductor admittance Y(s).

    Y(s) = sum_n g_n / (1 + s tau_n) + K_SIBC/sqrt(s) + c_1/s

    Returns
    -------
    A : ndarray (n_states, n_states)
    B : ndarray (n_states, 1)
    C : ndarray (1, n_states)
    D : ndarray (1, 1)
    n_foster, n_warburg, n_integrator : int
    """
    n_foster = len(g_n)
    xi_w, r_w = diffusive_quadrature(K_SIBC, n_aux=n_warburg_rungs)
    n_warburg = len(xi_w)
    n_int = 1

    n_states = n_foster + n_warburg + n_int
    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, 1))
    C = np.zeros((1, n_states))
    D = np.zeros((1, 1))

    # Foster poles: residue g_n at pole 1/tau_n
    for n in range(n_foster):
        A[n, n] = -1.0 / tau_n[n]
        B[n, 0] = 1.0
        C[0, n] = g_n[n] / tau_n[n]

    # Warburg rungs: K_SIBC/sqrt(s) via diffusive quadrature
    for k in range(n_warburg):
        idx = n_foster + k
        A[idx, idx] = -xi_w[k]
        B[idx, 0] = 1.0
        C[0, idx] = r_w[k]

    # Integrator for c_1/s
    idx_int = n_foster + n_warburg
    A[idx_int, idx_int] = 0.0
    B[idx_int, 0] = 1.0
    C[0, idx_int] = c1

    return A, B, C, D, n_foster, n_warburg, n_int


def save_mat(out_path, *, A, B, C, D, V, S, sigma, mu, K_SIBC, c_1,
              tau_foster, g_foster, n_foster, n_warburg, n_integrator,
              description="Conductor Y(s) Foster ladder LTI state-space."):
    """Save state-space + metadata to MATLAB-loadable .mat.

    Also writes a sibling .m helper script for Simulink loading.
    """
    from scipy.io import savemat
    out_path = Path(out_path)
    savemat(str(out_path), {
        "A": A, "B": B, "C": C, "D": D,
        "V_conductor": V, "S_conductor": S,
        "sigma": sigma, "mu": mu,
        "K_SIBC": K_SIBC, "c_1": c_1,
        "tau_foster": np.asarray(tau_foster),
        "g_foster": np.asarray(g_foster),
        "n_foster": n_foster, "n_warburg": n_warburg,
        "n_integrator": n_integrator,
        "description": description,
    })

    m_script = out_path.with_suffix(".m")
    with open(m_script, "w", encoding="utf-8") as fp:
        fp.write(f"""% {out_path.name} loader for Simulink
load('{out_path.name}');

% State-space LTI
sys = ss(A, B, C, D);
disp(sys);
fprintf('  states: %d (%d Foster + %d Warburg + %d integrator)\\n', ...
        size(A, 1), n_foster, n_warburg, n_integrator);
fprintf('  conductor: V = %.2f mm^3, S = %.2f mm^2, sigma = %.2e\\n', ...
        V_conductor*1e9, S_conductor*1e6, sigma);

figure;
bodeplot(sys);
title('Conductor Y(s) frequency response');

% Time-domain step response example
T = 0:1e-6:5e-3;
u_step = 0.1 * (T > 1e-3);
du_step = gradient(u_step, T(2)-T(1));
[y_out, t_out] = lsim(sys, du_step, T);
figure; plot(t_out, y_out);
xlabel('time (s)'); ylabel('Y(t) admittance');
title('Step response');
""")
    return str(out_path), str(m_script)
