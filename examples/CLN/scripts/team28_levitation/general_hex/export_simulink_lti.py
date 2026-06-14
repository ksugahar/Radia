"""export_simulink_lti.py -- Export the conductor alpha(s) Foster ladder
as a Simulink/MATLAB-loadable LTI state-space model (.mat) + helper .m script.

Foster representation of alpha(s) / V:
    alpha(s) / V = 1 - Y(s) / (sigma V)
                 = 1 - sum_n g_n / (sigma V (1 + s tau_n))
                   - K_SIBC / (sigma V sqrt(s))            (Mellin c_0 tail)
                   - c_1 / (sigma V s)                      (Mellin c_1 tail)

State-space form:
    For each Foster pole n: state x_n with x_n_dot = -x_n/tau_n + u
                            contribution to y: -(g_n / sigma V tau_n) * x_n
                            (the negative sign + 1/tau_n factor reproduce the pole's residue)

    Warburg K_SIBC/sqrt(s) via diffusive quadrature -> M auxiliary states
        (logarithmically-spaced first-order RC rungs, see CLAUDE memory
         "SPICE_subcircuit_template_Foster_of_Cauers" and Nagamine 2026.)

    Pure integrator c_1/s -> 1 additional state.

The resulting (A, B, C, D) is what Simulink's "LTI System" block consumes
via the From Workspace pattern.

Run:
    python export_simulink_lti.py --vol plate_100x50x5.vol --out ecb_plate_lti.mat
Then in MATLAB:
    load ecb_plate_lti.mat
    sys = ss(A, B, C, D);    % continuous-time state-space LTI
    % Simulink: place LTI System block, set sys reference
    bode(sys);               % verify frequency response
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cln_sibc_general_hex import (
    bulk_foster_via_eigen,
    K_SIBC_total,
    c1_polyhedral,
    measure_total_area_and_edges,
)


SIGMA_AL = 3.5e7
MU_0 = 4 * math.pi * 1e-7


def diffusive_quadrature(K_SIBC, n_aux=30, f_lo=1e-3, f_hi=1e10):
    """Approximate K_SIBC/sqrt(s) by sum of n_aux first-order RC rungs.

    Identity: 1/sqrt(s) = (1/pi) * integral_0^inf 1/(sqrt(xi) (s + xi)) dxi.

    Logarithmic quadrature: xi_k log-spaced over [omega_lo, omega_hi],
    with d(ln xi) = constant -> dxi = xi_k * d(ln xi). Then:
        1/sqrt(s) ~ (1/pi) sum_k (1/sqrt(xi_k)) * (1/(s + xi_k)) * xi_k * d(ln xi)
                   = (1/pi) sum_k sqrt(xi_k) * d(ln xi) / (s + xi_k)

    So:  K_SIBC/sqrt(s) ~ sum_k r_k / (s + xi_k),
         r_k = K_SIBC * sqrt(xi_k) * d(ln xi) / pi.

    The state-space realization with A_kk = -xi_k, B_k = 1, C_k = r_k
    gives transfer r_k / (s + xi_k) per rung.

    Returns: (xi_arr, r_arr) where xi_arr are pole frequencies (rad/s) and
    r_arr are the per-pole residues.
    """
    omega_lo = 2 * math.pi * f_lo
    omega_hi = 2 * math.pi * f_hi
    xi_arr = np.logspace(math.log10(omega_lo), math.log10(omega_hi), n_aux)
    d_log = math.log(xi_arr[1] / xi_arr[0])
    r_arr = K_SIBC * np.sqrt(xi_arr) * d_log / math.pi
    return xi_arr, r_arr


def build_state_space(g_n, tau_n, V, sigma, K_SIBC, c1, n_warburg_rungs=20):
    """Build (A, B, C, D) state-space for the conductor admittance Y(s).

    Y(s) = sum_n g_n / (1 + s tau_n) + K_SIBC/sqrt(s) + c_1/s

    Each Foster pole gives 1 state: x_dot = -x/tau_n + u, y += g_n/tau_n * x.
    Each Warburg rung gives 1 state: x_dot = -x/tau_k + u, y += K_k/tau_k * x.
    The 1/s integrator gives 1 state: x_dot = u, y += c_1 * x.

    The TOTAL output y(t) is the eddy current admittance response.
    Input u(t) is the time-derivative of B_ext scalar drive.
    """
    n_foster = len(g_n)
    xi_w, r_w = diffusive_quadrature(K_SIBC, n_aux=n_warburg_rungs)
    n_warburg = len(xi_w)
    n_int = 1   # integrator for c_1/s

    n_states = n_foster + n_warburg + n_int

    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, 1))
    C = np.zeros((1, n_states))
    D = np.zeros((1, 1))

    # Foster poles: residue g_n at pole 1/tau_n.  Each pole gives
    #   transfer g_n / (1 + s tau_n) = (g_n/tau_n) / (s + 1/tau_n)
    # State-space: A_nn = -1/tau_n, B_n = 1, C_n = g_n/tau_n.
    for n in range(n_foster):
        A[n, n] = -1.0 / tau_n[n]
        B[n, 0] = 1.0
        C[0, n] = g_n[n] / tau_n[n]

    # Warburg rungs (diffusive quadrature of K_SIBC/sqrt(s)):
    #   Each pole transfer is r_k / (s + xi_k).
    # State-space: A_kk = -xi_k, B_k = 1, C_k = r_k.
    for k in range(n_warburg):
        idx = n_foster + k
        A[idx, idx] = -xi_w[k]
        B[idx, 0] = 1.0
        C[0, idx] = r_w[k]

    # Integrator for c_1/s: A = 0, B = 1, C = c1.  Transfer = c1/s.
    idx_int = n_foster + n_warburg
    A[idx_int, idx_int] = 0.0
    B[idx_int, 0] = 1.0
    C[0, idx_int] = c1

    return A, B, C, D, n_foster, n_warburg, n_int


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vol", default="plate_100x50x5.vol")
    parser.add_argument("--n-eigen", type=int, default=300)
    parser.add_argument("--n-warburg", type=int, default=20)
    parser.add_argument("--out", default="ecb_plate_lti.mat")
    parser.add_argument("--sigma", type=float, default=SIGMA_AL)
    args = parser.parse_args()

    print("=== Simulink LTI export ===")
    print(f"  vol = {args.vol}")
    print(f"  Foster eigenmodes: {args.n_eigen}")
    print(f"  Warburg rungs: {args.n_warburg}")
    print()

    from ngsolve import Mesh, TaskManager
    mesh = Mesh(args.vol)
    print(f"  ne = {mesh.ne}")
    with TaskManager():
        lam, tau, g_n, V = bulk_foster_via_eigen(mesh, args.sigma, MU_0,
                                                   n_eigen=args.n_eigen)
        S_total, edges = measure_total_area_and_edges(mesh)
    K_SIBC = K_SIBC_total(S_total, args.sigma, MU_0)
    c1 = c1_polyhedral(edges, MU_0)

    A, B, C, D, n_f, n_w, n_i = build_state_space(
        g_n, tau, V, args.sigma, K_SIBC, c1, n_warburg_rungs=args.n_warburg)

    print(f"  Volume V = {V*1e9:.2f} mm^3, S = {S_total*1e6:.2f} mm^2")
    print(f"  K_SIBC = {K_SIBC:.4e}, c_1 = {c1:.4e}")
    print(f"  State-space: {A.shape[0]} states ({n_f} Foster + {n_w} Warburg + {n_i} integrator)")
    print()

    # Verify by Bode at sample frequencies
    print(f"  {'f (Hz)':>10}  {'|Y(s) ss|':>13}  {'|Y(s) ref|':>13}  {'rel err':>10}")
    import scipy.linalg as la
    for f in [1, 10, 100, 1e3, 1e4]:
        s = 1j * 2 * math.pi * f
        # State-space Y(s) = C (sI - A)^{-1} B + D
        ss_resp = C @ la.solve(s*np.eye(A.shape[0]) - A, B) + D
        Y_ss = complex(ss_resp[0, 0])
        # Reference direct Foster + Mellin tail
        import cmath
        Y_ref = np.sum(g_n / (1.0 + s*tau)) + K_SIBC/cmath.sqrt(s) + c1/s
        err = abs(Y_ss - Y_ref) / abs(Y_ref) * 100
        print(f"  {f:10.2e}  {abs(Y_ss):13.4e}  {abs(Y_ref):13.4e}  {err:9.4f}%")

    # Export to .mat for MATLAB / Simulink
    try:
        from scipy.io import savemat
        savemat(args.out, {
            "A": A, "B": B, "C": C, "D": D,
            "V_conductor": V, "S_conductor": S_total,
            "sigma": args.sigma, "mu": MU_0,
            "K_SIBC": K_SIBC, "c_1": c1,
            "tau_foster": np.asarray(tau),
            "g_foster": np.asarray(g_n),
            "n_foster": n_f, "n_warburg": n_w, "n_integrator": n_i,
            "description": ("Conductor Y(s) Foster ladder LTI state-space. "
                            "Input u = dB_ext/dt scalar; Output y = induced "
                            "admittance current density (scalar A_z approx). "
                            "alpha(s)/V = 1 - Y(s)/(sigma V)."),
        })
        print(f"\n  Saved {args.out}")
    except ImportError:
        print("  scipy unavailable; cannot export .mat")
        return

    # Also write a MATLAB helper script
    m_script = Path(args.out).with_suffix(".m")
    with open(m_script, "w", encoding="utf-8") as fp:
        fp.write(f"""% {args.out} loader for Simulink
load('{Path(args.out).name}');

% State-space LTI: input u = dB_ext/dt scalar, output y = Y(s) admittance.
sys = ss(A, B, C, D);
disp(sys);
fprintf('  states: %d (%d Foster + %d Warburg + %d integrator)\\n', ...
        size(A, 1), n_foster, n_warburg, n_integrator);
fprintf('  conductor: V = %.2f mm^3, S = %.2f mm^2, sigma = %.2e\\n', ...
        V_conductor*1e9, S_conductor*1e6, sigma);

% Verify with Bode
figure;
bodeplot(sys);
title('Conductor Y(s) frequency response');

% Time-domain: drive with dB_ext/dt and integrate
% Example: step from 0 to B_ext = 0.1 T over 1 ms
T = 0:1e-6:5e-3;
u_step = 0.1 * (T > 1e-3);   % B_ext(t)
du_step = gradient(u_step, T(2)-T(1));   % dB/dt
[y_out, t_out] = lsim(sys, du_step, T);
figure; plot(t_out, y_out);
xlabel('time (s)'); ylabel('Y(t) admittance');
title('Step response');
""")
    print(f"  Saved {m_script}")


if __name__ == "__main__":
    main()
