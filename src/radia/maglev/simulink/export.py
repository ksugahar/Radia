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


def _append_matrix_residue(R, pole, poles, Brows, Ccols, tol_rel=1e-12):
    """Realize a symmetric matrix residue term R / (s - pole) as rank-1 states.

    Eigen-decomposes R = sum_j kappa_j v_j v_j^T and appends one first-order
    state per significant eigenpair:  dx_j/dt = pole x_j + v_j^T u,
    y += kappa_j v_j x_j, whose transfer contribution is
    kappa_j v_j v_j^T / (s - pole); summed over j this is exactly
    R / (s - pole).  Indefinite R (e.g. the edge term C1, which is negative
    for convex bodies) is handled -- eigenvalues may be negative.
    """
    R = np.asarray(R, dtype=float)
    if R.size == 0:
        return
    w, Vv = np.linalg.eigh(0.5 * (R + R.T))
    wmax = float(np.max(np.abs(w))) if w.size else 0.0
    if wmax <= 0.0:
        return
    for j in range(len(w)):
        if abs(w[j]) <= tol_rel * wmax:
            continue
        v = Vv[:, j]
        poles.append(pole)
        Brows.append(v.copy())
        Ccols.append(w[j] * v)


def build_state_space_mimo(G_n, tau_n, V, sigma, K_mat, C1_mat,
                           n_warburg_rungs=30):
    """MIMO (P-port) state-space (A, B, C, D) for the matrix admittance Y(s).

        Y(s)_{pq} = sum_n G_n[n]_{pq} / (1 + s tau_n)
                    + K_mat_{pq} / sqrt(s) + C1_mat_{pq} / s

    with G_n (n_eigen, P, P) from alpha.bulk_foster_matrix_via_eigen, K_mat
    from alpha.K_SIBC_matrix, and C1_mat from cad_edges.edge_moment_matrix.

    Realization: each Foster pole contributes ONE state (its residue
    G_n[n]/tau_n is rank 1, an outer product b b^T); the Warburg tail uses
    P x n_warburg_rungs diffusive-quadrature states; the edge term uses up to
    P integrator states.  The resulting (A, B, C, D) is a passive MIMO LTI
    that drops into MATLAB `ss(A,B,C,D)` / Simulink.  Its transfer-function
    matrix reproduces Y_matrix_mixed(s) exactly (up to the diffusive-
    quadrature approximation of 1/sqrt(s)).  P=1 reduces to the same transfer
    function as build_state_space.

    Returns
    -------
    A, B, C, D : ndarray
        Shapes (n_states, n_states), (n_states, P), (P, n_states), (P, P).
    n_foster, n_warburg, n_integrator : int
        State counts of each block.
    """
    G_n = np.asarray(G_n, dtype=float)
    tau_n = np.asarray(tau_n, dtype=float)
    K_mat = np.asarray(K_mat, dtype=float)
    C1_mat = np.asarray(C1_mat, dtype=float)
    P = K_mat.shape[0]

    poles, Brows, Ccols = [], [], []

    # Foster poles: residue G_n[n]/tau_n at pole -1/tau_n (rank 1 -> 1 state)
    for n in range(G_n.shape[0]):
        _append_matrix_residue(G_n[n] / tau_n[n], -1.0 / tau_n[n],
                               poles, Brows, Ccols)
    n_foster = len(poles)

    # Warburg tail K_mat/sqrt(s) via diffusive quadrature (unit-K weights)
    xi_w, w_w = diffusive_quadrature(1.0, n_aux=n_warburg_rungs)
    for k in range(len(xi_w)):
        _append_matrix_residue(w_w[k] * K_mat, -xi_w[k], poles, Brows, Ccols)
    n_warburg = len(poles) - n_foster

    # Edge term C1_mat/s at pole 0
    _append_matrix_residue(C1_mat, 0.0, poles, Brows, Ccols)
    n_int = len(poles) - n_foster - n_warburg

    n_states = len(poles)
    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, P))
    C = np.zeros((P, n_states))
    D = np.zeros((P, P))
    for i in range(n_states):
        A[i, i] = poles[i]
        B[i, :] = Brows[i]
        C[:, i] = Ccols[i]

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
