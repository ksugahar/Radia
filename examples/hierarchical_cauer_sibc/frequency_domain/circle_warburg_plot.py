"""Generate circle_warburg.pdf for the IGTE 2026 digest (Fig. 1).

Plots |Y(jω)| log-log for a 2D Cu cylinder (a=5 mm), comparing:
- exact Bessel admittance Y_cyl(s) = π a² σ · 2I_1(γa)/[γa I_0(γa)]
- rank-10 CLN L-terminated (Foster sum, decays as 1/f at high f)
- rank-10 CLN R-terminated (Foster sum + plateau closure)
- Schur-augmented CLN basis 6 (additive form, Paper 1 Theorem 1):
  Y_R(s) = Y_Pade(s) + K_SIBC √s / (s + d)
  with d the bulk-surface crossover frequency.
- K_SIBC^cyl / √ω asymptote (slope -1/2)

TrueType fonts (pdf.fonttype=42) for IEEE/Elsevier journal submission.
"""
from __future__ import annotations
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
# Explicit Times New Roman registration (Windows path) so the font is
# actually applied, not silently fallen back to DejaVu Sans.
import matplotlib.font_manager as fm
for _p in (r"C:\Windows\Fonts\times.ttf", r"C:\Windows\Fonts\timesbd.ttf",
           r"C:\Windows\Fonts\timesi.ttf", r"C:\Windows\Fonts\timesbi.ttf"):
    try:
        fm.fontManager.addfont(_p)
    except Exception:
        pass
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
matplotlib.rcParams["mathtext.fontset"] = "stix"     # serif-style math, pairs with TNR
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import iv, jn_zeros
from sympy import symbols, factorial, Rational
import sympy as sp

# Cu cylinder parameters
SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
A = 5e-3
Y_DC = math.pi * A**2 * SIGMA
K_SIBC = 2 * math.pi * A * math.sqrt(SIGMA / MU)
T_char = MU * SIGMA * A**2


def Y_exact(s):
    gamma = np.sqrt(s * MU * SIGMA)
    gA = gamma * A
    result = np.empty_like(gA, dtype=complex)
    small = np.abs(gA) < 50
    if np.any(small):
        result[small] = (Y_DC * 2.0 * iv(1, gA[small])
                         / (gA[small] * iv(0, gA[small])))
    large = ~small
    if np.any(large):
        result[large] = Y_DC * 2.0 / gA[large]
    return result


def foster_spectrum(M):
    x = jn_zeros(0, M)
    tau = MU * SIGMA * A**2 / x**2
    g_sq = Y_DC * 4.0 / x**2
    return x, tau, g_sq


def Y_L_term(s, N):
    """L-terminated Cauer-I ladder, basis N = 2N+1 elements.
    Convention (Sugahara 2018, slab CF analysis): the ladder starts with
    series L_c, then N rungs of (R + L parallel pairs), terminated by
    the final L of the Nth rung.  Total: 1 initial L + N R + N L = 2N+1
    elements, with N+1 inductors = N+1 poles in Y(s).  At high s, the
    chain ends in an inductor so Y -> 0 like 1/s."""
    _, tau, g_sq = foster_spectrum(N + 1)   # N+1 modes
    return np.sum(g_sq[None, :] / (1.0 + s[:, None] * tau[None, :]), axis=1)


def Y_R_term(s, N):
    """R-terminated Cauer-II ladder, basis N = 2N elements.
    Same topology as L-term but truncated one step earlier, ending in
    the R of the Nth rung (the final L is omitted).  Total: 1 initial L
    + N R + (N-1) L = 2N elements, with N inductors = N poles, plus the
    terminating R that produces a plateau at high s.

    The plateau height absorbs the "missing tail mass":
        c_inf = Y_DC - sum_{k<=N} g_k^2

    So R-term has 1 LESS element than L-term (2N vs 2N+1) per the
    Sugahara 2018 convention.  L-term branches at f_{N+1} (higher);
    R-term branches at f_N (lower)."""
    _, tau, g_sq = foster_spectrum(N)        # N modes
    Y_L_N = np.sum(g_sq[None, :] / (1.0 + s[:, None] * tau[None, :]), axis=1)
    c_inf = Y_DC - np.sum(g_sq)
    return Y_L_N + c_inf


# Pade[N-1/N] in s (at s=0) of the cylinder
def taylor_coefs(order):
    z = symbols('z')
    I0 = sum(Rational(1, factorial(k)**2) * (z/2)**(2*k) for k in range(order + 2))
    I1 = sum(Rational(1, factorial(k)*factorial(k+1)) * (z/2)**(2*k+1) for k in range(order + 2))
    Y_norm = 2 * I1 / (z * I0)
    Y_taylor = sp.series(Y_norm, z, 0, 2*(order+1) + 1).removeO()
    Y_expanded = sp.expand(Y_taylor)
    return [Y_expanded.coeff(z, 2*k) for k in range(order + 1)]


def pade(coefs, L, M):
    A_mat = np.zeros((M, M))
    b_vec = np.zeros(M)
    for i in range(M):
        row_k = L + 1 + i
        b_vec[i] = -float(coefs[row_k])
        for j in range(M):
            col_idx = row_k - 1 - j
            if 0 <= col_idx <= L + M:
                A_mat[i, j] = float(coefs[col_idx])
    q_sol = np.linalg.solve(A_mat, b_vec)
    q_coefs = np.concatenate(([1.0], q_sol))
    p_coefs = np.zeros(L + 1)
    for k in range(L + 1):
        p = float(coefs[k])
        for j in range(1, min(k, M) + 1):
            p += float(coefs[k - j]) * q_coefs[j]
        p_coefs[k] = p
    return p_coefs, q_coefs


def Y_Pade(s, N):
    coefs_z2 = taylor_coefs(2 * N + 1)
    coefs_s = [float(c) * T_char**k for k, c in enumerate(coefs_z2)]
    p, q = pade(coefs_s, N - 1, N)
    s_arr = np.atleast_1d(s).astype(complex)
    P = sum(p[k] * s_arr**k for k in range(len(p)))
    Q = sum(q[k] * s_arr**k for k in range(len(q)))
    return Y_DC * P / Q


def Y_Schur_aug(s, N_Pade, d):
    """Schur-augmented CLN (Paper 1 Theorem 1 additive form, Pade base):

        Y_R(s) = Y_Pade(s) + K_SIBC * sqrt(s) / (s + d).

    The single tuning parameter d > 0 sets the bulk-surface crossover.
    Preserves DC exactly (sqrt(s)/(s+d) vanishes at s=0) and recovers
    the SIBC asymptote (sqrt(s)/(s+d) -> 1/sqrt(s) as s -> infty).
    """
    YP = Y_Pade(s, N_Pade)
    return YP + K_SIBC * np.sqrt(s) / (s + d)


def Y_CLN_modal(s, M0_order, M_modes=400):
    """Canonical Cauer Ladder Network in modal (Foster-eigenmode) space.

    Build K_0=I, M=diag(tau_n), b=g_n in the first M_modes Foster modes,
    then do a single-point Krylov reduction with M0_order vectors at s=0:
        V_k = (K_0^{-1} M)^k b = (tau^k) * g,  k = 0..M0_order-1
    K_0-orthonormalise (in <u, v> = sum u_n v_n) to obtain Q (M0 x M_modes).
    Reduced operators:
        K_0_r = Q Q^T = I_{M0 x M0}
        M_r   = Q diag(tau) Q^T
        b_r   = Q g
    Return Y_CLN(s) = b_r^T (K_0_r + s M_r)^{-1} b_r.
    M0_order = 3 reproduces the canonical CLN3 (Kameari 2018, Padé[2/3]
    at s=0).
    """
    _, tau, g_sq = foster_spectrum(M_modes)
    g = np.sqrt(g_sq)
    # Build raw single-point Krylov basis at s=0
    V = np.zeros((M0_order, M_modes), dtype=complex)
    for k in range(M0_order):
        V[k] = (tau ** k) * g
    # MGS in K_0 = I inner product
    Q = np.zeros_like(V)
    for i in range(M0_order):
        u = V[i].copy()
        for j in range(i):
            u = u - np.vdot(Q[j], u) * Q[j]
        Q[i] = u / np.sqrt(np.vdot(u, u).real)
    K0_r = np.eye(M0_order, dtype=complex)
    M_r  = (Q * tau[None, :]) @ Q.conj().T
    b_r  = Q @ g
    Y = np.empty(s.shape, dtype=complex)
    for k, sv in enumerate(s):
        x = np.linalg.solve(K0_r + sv * M_r, b_r)
        Y[k] = b_r @ x
    return Y


def Y_CLN3_Schur(s, d, M_modes=400):
    """Canonical Schur-augmented CLN: basis 4 = CLN3 bulk + sqrt(s) Schur DOF.

    This is the Paper 1 Theorem 1 single-DOF Schur augmentation applied to
    the canonical 3-Krylov CLN at s=0.  Total basis size 4 (3 bulk + 1
    surface).  Wall-band peak relative error on the Cu cylinder is ~17%
    with optimally-tuned d/(2pi) ~ 1.58e5 Hz.
    """
    return Y_CLN_modal(s, M0_order=3, M_modes=M_modes) + \
           K_SIBC * np.sqrt(s) / (s + d)


def main():
    f = np.logspace(0, 8, 240)
    omega = 2 * math.pi * f
    s = 1j * omega

    Y_ex = Y_exact(s)
    Y_L = Y_L_term(s, 10)
    Y_R = Y_R_term(s, 10)
    # Canonical Schur-augmented CLN3 (Paper 1 Theorem 1):
    # basis 4 = 3 single-point Krylov vectors at s=0 + 1 sqrt(s) Schur DOF.
    # d/(2pi) = 1.58e5 Hz tuned to minimise wall-band peak relative error.
    d_canon = 2 * math.pi * 1.58e5
    Y_Schur_basis4 = Y_CLN3_Schur(s, d=d_canon)

    _, tau10, _ = foster_spectrum(10)
    f_N = 1.0 / (2 * math.pi * tau10[-1])

    # Sanity: max rel-err for Schur-augmented CLN3 (basis 4)
    err_Schur = np.abs(Y_Schur_basis4 - Y_ex) / np.abs(Y_ex) * 100
    print(f"Schur-augmented CLN3 (basis 4) d/(2pi)={d_canon/(2*math.pi):.2e} Hz: "
          f"max rel-err over 1Hz-1e8Hz = {err_Schur.max():.3f}% "
          f"@ {f[err_Schur.argmax()]:.2e} Hz")

    fig, ax = plt.subplots(figsize=(3.45, 2.6))

    ax.loglog(f, np.abs(Y_ex), color="black", linestyle="-", lw=1.8,
              label=r"$Y_{\mathrm{exact}}$ (Bessel)", zorder=5)
    # Schur-augmented CLN3 as dotted line so the wall-band deviation
    # against the exact Bessel curve is visible through the gaps.
    ax.loglog(f, np.abs(Y_Schur_basis4), color="red", linestyle=":",
              lw=2.0, label="Schur-augmented CLN3 (basis 4)", zorder=10)
    ax.loglog(f, np.abs(Y_L), color="royalblue", linestyle=(0, (5, 3)),
              lw=1.3, label="CLN L-term (basis 10)", zorder=4)
    ax.loglog(f, np.abs(Y_R), color="darkorange", linestyle="-.", lw=1.3,
              label="CLN R-term (basis 10)", zorder=4)

    ax.axvline(f_N, color="gray", linestyle=":", lw=0.7, alpha=0.6)
    ax.text(f_N * 1.15, 2.0, "$f_N$", color="gray", fontsize=8,
            ha="left", va="bottom")

    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Y(j\omega)|$ (S/m)")
    ax.set_xlim(1e0, 1e8)
    ax.grid(True, which="both", linestyle=":", alpha=0.35)
    ax.legend(loc="lower left", fontsize=6.5, frameon=False,
              handlelength=2.0, labelspacing=0.25)

    fig.tight_layout(pad=0.25)
    out = Path(__file__).parent / "circle_warburg"
    fig.savefig(str(out) + ".pdf")
    fig.savefig(str(out) + ".png", dpi=220)
    print(f"wrote {out}.pdf, {out}.png")

    print(f"K_SIBC^cyl = {K_SIBC:.4e}")
    print(f"f_N (rank 10) = {f_N:.3e} Hz")


if __name__ == "__main__":
    main()
